#!/usr/bin/env bash
# Source this to set MPI environment and helpers.
# Modes:
#   auto     : detect cluster vs laptop; safe defaults
#   emulate  : force laptop-friendly Open MPI settings (loopback/TCP/oversubscribe)
#   cluster  : disable emulation; favor scheduler/HCAs
#
# Usage:
#   source mpi_env.sh [auto|emulate|cluster]
# After sourcing:
#   mpi_exec 4 ./your_mpi_binary
#
# It also exports CTest/CMake-friendly variables:
#   MPIEXEC_EXECUTABLE, MPIEXEC_NUMPROC_FLAG, MPIEXEC_PREFLAGS, MPIEXEC_POSTFLAGS

# --- must be sourced ---
if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
  echo "Please 'source' this file:  source ${BASH_SOURCE[0]} [mode]"
  exit 1
fi

# --- helpers ---
_have() { command -v "$1" >/dev/null 2>&1; }
_lower() { printf '%s' "$1" | tr '[:upper:]' '[:lower:]'; }

MPI_MODE="$(_lower "${1:-${MPI_MODE:-auto}}")"
MPI_TEST_NP="${MPI_TEST_NP:-2}"      # default ranks for mpi_exec when NP omitted
MPI_BIND="${MPI_BIND:-none}"         # bind policy for Open MPI when emulating
MPI_MAP_BY="${MPI_MAP_BY:-slot}"     # map policy for Open MPI when emulating

# Detect vendor
MPI_VENDOR="unknown"
if _have mpirun; then
  if mpirun --version 2>&1 | grep -qi 'Open MPI'; then
    MPI_VENDOR="openmpi"
  elif mpirun --version 2>&1 | grep -qiE 'MPICH|Intel\(R\) MPI|Hydra'; then
    MPI_VENDOR="mpich"
  fi
fi

# Choose launcher: prefer SLURM when in a job
MPI_LAUNCHER=""
MPI_LAUNCHER_OPTS=""
if [[ -n "${SLURM_JOB_ID:-}" ]] && _have srun; then
  MPI_LAUNCHER="srun"
  MPI_LAUNCHER_OPTS="--mpi=pmix"
elif _have mpirun; then
  MPI_LAUNCHER="mpirun"
elif _have mpiexec; then
  MPI_LAUNCHER="mpiexec"
fi

# Clean slate for emulation knobs
_unset_emulation() {
  unset OMPI_MCA_btl OMPI_MCA_btl_tcp_if_include OMPI_MCA_oob_tcp_if_include OMPI_MCA_pml
  unset UCX_TLS
}

_apply_emulation() {
  # Favor stability for laptop dev, esp. Open MPI.
  if [[ "$MPI_VENDOR" == "openmpi" ]]; then
    export OMPI_MCA_btl="self,tcp"
    export OMPI_MCA_btl_tcp_if_include="lo"
    export OMPI_MCA_oob_tcp_if_include="lo"
    export UCX_TLS="tcp"
    : "${OMPI_MCA_pml:=ob1}"
    # Launcher-side niceties
    if [[ "$MPI_LAUNCHER" == "mpirun" ]]; then
      MPI_LAUNCHER_OPTS+=" --oversubscribe --bind-to ${MPI_BIND} --map-by ${MPI_MAP_BY}"
    fi
  fi
}

case "$MPI_MODE" in
  auto)
    _unset_emulation
    # On clusters (SLURM detected) keep things clean; on laptops add OpenMPI emulation
    if [[ -z "${SLURM_JOB_ID:-}" ]]; then
      [[ "$MPI_VENDOR" == "openmpi" ]] && _apply_emulation
    fi
    ;;
  emulate)
    _unset_emulation
    _apply_emulation
    ;;
  cluster)
    _unset_emulation
    ;;
  *)
    echo "Unknown mode '$MPI_MODE' (use auto|emulate|cluster)"
    return 1
    ;;
esac

# Export CMake/CTest friendly variables, derived from the launcher choice
case "$MPI_LAUNCHER" in
  srun)
    export MPIEXEC_EXECUTABLE="srun"
    export MPIEXEC_NUMPROC_FLAG="-n"
    # shellcheck disable=SC2089
    export MPIEXEC_PREFLAGS="${MPI_LAUNCHER_OPTS}"
    ;;
  mpirun|mpiexec)
    export MPIEXEC_EXECUTABLE="${MPI_LAUNCHER}"
    # OpenMPI/mpiexec commonly use -np; mpiexec from MPICH also accepts -n, but -np is fine.
    export MPIEXEC_NUMPROC_FLAG="-np"
    export MPIEXEC_PREFLAGS="${MPI_LAUNCHER_OPTS}"
    ;;
  *)
    unset MPIEXEC_EXECUTABLE MPIEXEC_NUMPROC_FLAG MPIEXEC_PREFLAGS
    ;;
esac
# Optional: postflags placeholder (rarely needed, but supported by CTest/CMake)
: "${MPIEXEC_POSTFLAGS:=}"
export MPIEXEC_POSTFLAGS

# Convenience runner
mpi_exec() {
  local np="$1"
  shift || true
  if [[ -z "${np:-}" || "$np" =~ ^[^0-9]+$ ]]; then
    # NP omitted; treat first arg as command
    np="$MPI_TEST_NP"
  fi
  if [[ -z "${MPIEXEC_EXECUTABLE:-}" ]]; then
    echo "No MPI launcher found (srun/mpirun/mpiexec)."
    return 1
  fi
  if [[ "$MPIEXEC_EXECUTABLE" == "srun" ]]; then
    echo "+ srun -n $np ${MPIEXEC_PREFLAGS} $*"
    srun -n "$np" ${MPIEXEC_PREFLAGS} "$@"
  else
    echo "+ ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} $np ${MPIEXEC_PREFLAGS} $*"
    "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "$np" ${MPIEXEC_PREFLAGS} "$@"
  fi
}

# Friendly summary (one-liners, safe to pipe into logs)
echo "MPI mode=${MPI_MODE} vendor=${MPI_VENDOR} launcher=${MPI_LAUNCHER:-none} npflag=${MPIEXEC_NUMPROC_FLAG:-n/a}"
[[ -n "${MPIEXEC_PREFLAGS:-}" ]] && echo "MPI preflags=${MPIEXEC_PREFLAGS}"
[[ -n "${UCX_TLS:-}" ]] && echo "UCX_TLS=${UCX_TLS}"
