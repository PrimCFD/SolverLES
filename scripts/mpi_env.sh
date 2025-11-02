#!/usr/bin/env bash
set -euo pipefail
# Source this to set MPI environment + helpers that work on a laptop and scale to HPC.
# Modes:
#   auto     : detect SLURM/cluster vs. laptop; choose safe defaults
#   emulate  : force laptop-friendly Open MPI settings (loopback/TCP/oversubscribe)
#   cluster  : disable emulation; favor scheduler + high-speed fabrics
#
# Usage:
#   source scripts/mpi_env.sh [auto|emulate|cluster]
# After sourcing:
#   mpi_exec 4 ./your_mpi_binary args...

# --- must be sourced ---
if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
  echo "Please 'source' this file:  source ${BASH_SOURCE[0]} [mode]"
  exit 1
fi

# --- helpers ---
_have() { command -v "$1" >/dev/null 2>&1; }
_lower() { printf '%s' "$1" | tr '[:upper:]' '[:lower:]'; }

MPI_MODE="$(_lower "${1:-${MPI_MODE:-auto}}")"
MPI_TEST_NP="${MPI_TEST_NP:-8}"      # default np for mpi_exec if omitted
# Laptop-friendly binding (Open MPI) — overridden by scheduler on clusters
MPI_BIND="${MPI_BIND:-core}"
MPI_MAP_BY="${MPI_MAP_BY:-core}"
# Refuse to oversubscribe when enabled (ranks × PE must fit cores)
MPI_STRICT_PE="${MPI_STRICT_PE:-0}"

# Detect vendor/launcher
MPI_VENDOR="unknown"
if _have mpirun; then
  if mpirun --version 2>&1 | grep -qi 'Open MPI'; then
    MPI_VENDOR="openmpi"
  elif mpirun --version 2>&1 | grep -qiE 'MPICH|Intel\\(R\\) MPI|Hydra'; then
    MPI_VENDOR="mpich"
  fi
fi

MPI_LAUNCHER=""
MPI_LAUNCHER_OPTS=""
if [[ -n "${SLURM_JOB_ID:-}" ]] && _have srun; then
  # Prefer srun inside SLURM allocations; let PMIx/PMI wireup MPI libs
  MPI_LAUNCHER="srun"
  MPI_LAUNCHER_OPTS="--mpi=pmix"   # PMIx is the common fast path. :contentReference[oaicite:0]{index=0}
elif _have mpirun; then
  MPI_LAUNCHER="mpirun"
elif _have mpiexec; then
  MPI_LAUNCHER="mpiexec"
fi


# Open MPI: allow running as root inside CI/containers if needed. :contentReference[oaicite:1]{index=1}
if [[ "$MPI_VENDOR" == "openmpi" && "${EUID:-$(id -u)}" -eq 0 ]]; then
  export OMPI_ALLOW_RUN_AS_ROOT=1
  export OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
  MPI_LAUNCHER_OPTS+=" --allow-run-as-root"
fi

# --- make CMake & PETSc use the same MPI wrappers we detected ---
: "${EXTRA_CMAKE_ARGS:=}"

if command -v mpicc >/dev/null 2>&1;  then
  export MPI_C_COMPILER="$(command -v mpicc)"
fi
if command -v mpicxx >/dev/null 2>&1; 
then
  export MPI_CXX_COMPILER="$(command -v mpicxx)"
fi
# Fortran is optional; PETSc benefits if present
if command -v mpifort >/dev/null 2>&1; then
  export MPI_Fortran_COMPILER="$(command -v mpifort)"
elif command -v mpif90 >/dev/null 2>&1; then
  export MPI_Fortran_COMPILER="$(command -v mpif90)"
fi

# Hand these to CMake’s FindMPI (append only when set)
if [[ -n "${MPI_C_COMPILER:-}" ]];      then EXTRA_CMAKE_ARGS+=" -DMPI_C_COMPILER=${MPI_C_COMPILER}"; fi
if [[ -n "${MPI_CXX_COMPILER:-}" ]];    then EXTRA_CMAKE_ARGS+=" -DMPI_CXX_COMPILER=${MPI_CXX_COMPILER}"; fi
if [[ -n "${MPI_Fortran_COMPILER:-}" ]]; then EXTRA_CMAKE_ARGS+=" -DMPI_Fortran_COMPILER=${MPI_Fortran_COMPILER}"; fi
export EXTRA_CMAKE_ARGS

# Also useful for PETSc’s configure: force its CC/CXX/FC to the same wrappers
export PETSC_CC="${MPI_C_COMPILER:-cc}"
export PETSC_CXX="${MPI_CXX_COMPILER:-c++}"
if [[ -n "${MPI_Fortran_COMPILER:-}" ]]; then export PETSC_FC="${MPI_Fortran_COMPILER}"; fi

# Clean slate for emulation knobs
_unset_emulation() {
  unset OMPI_MCA_btl OMPI_MCA_btl_tcp_if_include OMPI_MCA_oob_tcp_if_include OMPI_MCA_pml UCX_TLS
}

# How many usable logical cores do we have on this node?
_available_cores() {
  if [[ -n "${SLURM_CPUS_ON_NODE:-}" ]]; then
    printf "%s" "${SLURM_CPUS_ON_NODE}"
  elif command -v nproc >/dev/null 2>&1; then
    nproc
  else
    # Fallback: count "processor" lines in /proc/cpuinfo
    grep -c '^processor' /proc/cpuinfo 2>/dev/null || echo 1
  fi
}

# Compute PE (threads per MPI rank) from OMP/SLURM, default 1.
_compute_pe() {
  local pe=""
  # Prefer the scheduler value if present
  if [[ -n "${SLURM_CPUS_PER_TASK:-}" ]]; then
    pe="${SLURM_CPUS_PER_TASK}"
    # If user didn't set OMP, mirror scheduler to OMP for library behavior
    [[ -z "${OMP_NUM_THREADS:-}" ]] && export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
  fi
  # Fallback to OMP_NUM_THREADS
  if [[ -z "${pe}" && -n "${OMP_NUM_THREADS:-}" ]]; then
    pe="${OMP_NUM_THREADS}"
  fi
  [[ -z "${pe}" ]] && pe=1
  printf "%s" "${pe}"
}

# Pick the right shared-memory BTL for this Open MPI
# OMPI 3.x/4.x -> "vader"; OMPI 5.x -> "sm" (vader remains an alias)
_detect_ompi_shm_btl() {
  local ver btl="vader"
  if command -v ompi_info >/dev/null 2>&1; then
    ver="$(ompi_info --version 2>/dev/null | awk -F'[ .]' '/Open MPI:/ {print $3}')"
    # If major version >= 5, prefer "sm"
    if [[ -n "$ver" && "$ver" -ge 5 ]]; then
      btl="sm"
    fi
    # If explicitly present, "vader" works across versions (kept as alias in 5.x)
    # We keep "vader" as safe default for 3.x/4.x.
  fi
  printf "%s" "$btl"
}

_apply_emulation_openmpi() {
  # Favor local-node perf: shared-memory first, then TCP loopback.
  local shm_btl; shm_btl="$(_detect_ompi_shm_btl)"
  local pe; pe="$(_compute_pe)"
  export OMPI_MCA_btl="self,${shm_btl},tcp"
  export OMPI_MCA_btl_tcp_if_include="lo"
  export OMPI_MCA_oob_tcp_if_include="lo"
  : "${OMPI_MCA_pml:=ob1}"                           # simple PML
  : "${UCX_TLS:=sm,self,tcp}"                        # prefer UCX shared-memory; TCP as fallback
  if [[ "$MPI_LAUNCHER" == "mpirun" || "$MPI_LAUNCHER" == "mpiexec" ]]; then
    # Hybrid done right on laptops: reserve PE cores per rank and bind to cores.
    # (--map-by socket:PE=... and --bind-to core are the OMPI-documented knobs) 
    MPI_LAUNCHER_OPTS+=" --oversubscribe --map-by socket:PE=${pe} --bind-to core"
    # (optional) debug bindings: MPI_LAUNCHER_OPTS+=" --report-bindings"
  fi
}

_apply_cluster_defaults() {
  # Let vendor defaults + scheduler decide. Still provide reasonable hints:
  if [[ "$MPI_VENDOR" == "openmpi" ]]; then
    # Prefer UCX PML on clusters when available (InfiniBand/RoCE etc.)
    : "${OMPI_MCA_pml:=ucx}"
    export OMPI_MCA_pml
  elif [[ "$MPI_VENDOR" == "mpich" ]]; then
    # MPICH/Intel MPI pinning knobs (users can override).
    : "${I_MPI_PIN:=1}"               # Intel MPI pin on
    : "${I_MPI_PIN_DOMAIN:=auto}"     # domain auto (socket/core)
    export I_MPI_PIN I_MPI_PIN_DOMAIN  # :contentReference[oaicite:4]{index=4}
    # MPICH rank reorder is useful under SLURM topologies; off by default. :contentReference[oaicite:5]{index=5}
    : "${MPICH_RANK_REORDER_METHOD:=0}"
    export MPICH_RANK_REORDER_METHOD
  fi

  # Enforce PE-aware placement depending on launcher:
  local pe; pe="$(_compute_pe)"
  if [[ "$MPI_LAUNCHER" == "srun" ]]; then
    # With SLURM, request/bind cores for hybrid tasks.
    # If the allocation doesn't already define cpus-per-task, supply it from OMP.
    if [[ -z "${SLURM_CPUS_PER_TASK:-}" && -n "${pe}" && "${pe}" -gt 0 ]]; then
      MPI_LAUNCHER_OPTS+=" --cpus-per-task=${pe}"
      # Keep OMP in sync for threaded libraries
      [[ -z "${OMP_NUM_THREADS:-}" ]] && export OMP_NUM_THREADS="${pe}"
    fi
    MPI_LAUNCHER_OPTS+=" --cpu-bind=cores"
  elif [[ "$MPI_VENDOR" == "openmpi" && ( "$MPI_LAUNCHER" == "mpirun" || "$MPI_LAUNCHER" == "mpiexec" ) ]]; then
    # On clusters when launching with mpirun, still do PE-aware mapping/binding (no oversubscribe).
    MPI_LAUNCHER_OPTS+=" --map-by socket:PE=${pe} --bind-to core"
  fi
}

case "$MPI_MODE" in
  auto)
    _unset_emulation
    if [[ -n "${SLURM_JOB_ID:-}" ]]; then
      _apply_cluster_defaults
    else
      [[ "$MPI_VENDOR" == "openmpi" ]] && _apply_emulation_openmpi
    fi
    ;;
  emulate)
    _unset_emulation
    [[ "$MPI_VENDOR" == "openmpi" ]] && _apply_emulation_openmpi
    ;;
  cluster)
    _unset_emulation
    _apply_cluster_defaults
    ;;
  *)
    echo "Unknown mode '$MPI_MODE' (use auto|emulate|cluster)"
    return 1
    ;;
esac

# Export CTest/CMake-friendly variables from the chosen launcher. :contentReference[oaicite:6]{index=6}
case "$MPI_LAUNCHER" in
  srun)
    export MPIEXEC_EXECUTABLE="srun"
    export MPIEXEC_NUMPROC_FLAG="-n"
    export MPIEXEC_PREFLAGS="${MPI_LAUNCHER_OPTS}"
    ;;
  mpirun|mpiexec)
    export MPIEXEC_EXECUTABLE="${MPI_LAUNCHER}"
    export MPIEXEC_NUMPROC_FLAG="-np"      # works for Open MPI and MPICH hydra. :contentReference[oaicite:7]{index=7}
    export MPIEXEC_PREFLAGS="${MPI_LAUNCHER_OPTS}"
    ;;
  *)
    unset MPIEXEC_EXECUTABLE MPIEXEC_NUMPROC_FLAG MPIEXEC_PREFLAGS
    ;;
esac
: "${MPIEXEC_POSTFLAGS:=}"
export MPIEXEC_POSTFLAGS

# Convenience runner
mpi_exec() {
  local np="$1"; shift || true
  if [[ -z "${np:-}" || "$np" =~ ^[^0-9]+$ ]]; then np="$MPI_TEST_NP"; fi
  # Enforce MPI_STRICT_PE (no oversubscription by cores)
  if [[ "${MPI_STRICT_PE}" == "1" ]]; then
    local pe cores need
    pe="$(_compute_pe)"
    cores="$(_available_cores)"
    need=$(( np * pe ))
    if (( need > cores )); then
      echo "❌ MPI_STRICT_PE=1: Refusing to run, ranks(${np}) × PE(${pe}) = ${need} exceeds available cores (${cores}). Set MPI_STRICT_PE=0 to override." >&2
      return 2
    fi
  fi
  if [[ -z "${MPIEXEC_EXECUTABLE:-}" ]]; then
    echo "No MPI launcher found (srun/mpirun/mpiexec)."; return 1; fi
  if [[ "$MPIEXEC_EXECUTABLE" == "srun" ]]; then
    echo "+ srun -n $np ${MPIEXEC_PREFLAGS} $*"
    srun -n "$np" ${MPIEXEC_PREFLAGS} "$@"
  else
    echo "+ ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} $np ${MPIEXEC_PREFLAGS} $*"
    "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "$np" ${MPIEXEC_PREFLAGS} "$@"
  fi
}

# Friendly summary for logs
echo "MPI mode=${MPI_MODE} vendor=${MPI_VENDOR} launcher=${MPI_LAUNCHER:-none} npflag=${MPIEXEC_NUMPROC_FLAG:-n/a}"
[[ -n "${MPIEXEC_PREFLAGS:-}" ]] && echo "MPI preflags=${MPIEXEC_PREFLAGS}"
[[ -n "${UCX_TLS:-}" ]] && echo "UCX_TLS=${UCX_TLS}"
echo "MPI strict PE:    $([[ \"${MPI_STRICT_PE}\" == \"1\" ]] && echo on || echo off)"
