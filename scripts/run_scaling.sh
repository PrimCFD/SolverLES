#!/usr/bin/env bash
# =============================================================================
# run_scaling.sh — sweep grids × ranks × threads and dump CSVs for scaling plots
# =============================================================================
# Location: scripts/run_scaling.sh
# Output:   build-scaling/build_info.csv, build-scaling/results_scaling.csv
#
# Usage examples:
#   scripts/run_scaling.sh
#   # Choose a near-cubic DMDA factorization (px,py,pz) for (N, R)
choose_dmda() {
  local N="$1"; local R="$2"
  local best_px=1 best_py=1 best_pz="$R" best_cost=1e18
  local ds=() d
  for ((d=1; d<=R; ++d)); do
    if (( R % d == 0 )); then ds+=("$d"); fi
  done
  for px in "${ds[@]}"; do
    if (( N % px != 0 )); then continue; fi
    for py in "${ds[@]}"; do
      if (( (R / px) % py != 0 )); then continue; fi
      local pz=$(( R / (px*py) ))
      if (( pz < 1 )); then continue; fi
      if (( N % py != 0 || N % pz != 0 )); then continue; fi
      local cost
      cost=$(awk -v a="$px" -v b="$py" -v c="$pz" 'BEGIN{m=(a+b+c)/3.0; print ((a-m)^2+(b-m)^2+(c-m)^2)}')
      awk -v c1="$cost" -v c2="$best_cost" 'BEGIN{exit !(c1 < c2)}' && {
        best_cost="$cost"; best_px="$px"; best_py="$py"; best_pz="$pz";
      }
    done
  done
  echo "$best_px $best_py $best_pz"
}

# SCALING_GRIDS="128 256" SCALING_RANKS="1 4 8" SCALING_THREADS="1 2" scripts/run_scaling.sh
#   PETSC_OPTIONS_EXTRA="-pc_type mg -pc_mg_levels 6 -pc_mg_galerkin pmat -mg_levels_ksp_type chebyshev -mg_levels_pc_type jacobi -ksp_type pipecg" scripts/run_scaling.sh
#
# Notes:
# - Grids are cubic (N^3). We pass DMDA sizes via -da_grid_{x,y,z} (PETSc).
# - We always collect KSP convergence info (-ksp_converged_reason) and a brief -ksp_view.
# - We try to parse mean_us/stddev_us from bench output; if not present, we fall back to wall time.
# - We save per-run stdout/stderr logs under build-scaling/logs/ for later inspection.
#
set -euo pipefail

# --------------------------- repo & environment detection --------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
# Build and stage everything in build-scaling, not build-perf
BUILD_DIR="${REPO_ROOT}/build-scaling"
SCALING_DIR="${REPO_ROOT}/build-scaling"
LOG_DIR="${SCALING_DIR}/logs"
mkdir -p "${SCALING_DIR}" "${LOG_DIR}"

# Try to match MPI env helper if present
export MPI_STRICT_PE="${MPI_STRICT_PE:-1}"
MPI_ENV="${SCRIPT_DIR}/mpi_env.sh"
if [[ -f "${MPI_ENV}" ]]; then
  # shellcheck disable=SC1090
  source "${MPI_ENV}"
fi

# Defaults (can be overridden by env)
: "${SCALING_GRIDS:=64 128 256}"
: "${SCALING_RANKS:=1 2 4}"
: "${SCALING_THREADS:=1 2}"
: "${SCALING_REPS:=3}"
# Build before running (1=yes/0=no)
: "${SCALING_BUILD_FIRST:=1}"
# Auto layout -da_processors_{x,y,z} close to cubic for each -np (1=yes/0=no)
: "${SCALING_SET_DA_PROCS:=1}"

# PETSc options: caller can extend via PETSC_OPTIONS_EXTRA
: "${PETSC_OPTIONS_EXTRA:=}"
# Always include useful diagnostics; caller can still override tolerances/types via EXTRA
BASE_PETSC_OPTS="-ksp_converged_reason -ksp_view"

# MPI launcher and flags
MPI_LAUNCHER="${MPI_LAUNCHER:-mpirun}"
MPI_PREFLAGS="${MPI_PREFLAGS:---oversubscribe --map-by core --bind-to core}"
MPI_NPFLAG="${MPI_NPFLAG:--np}"

# Executable (built by run_perf_tests.sh); allow override
BENCH_EXE="${BENCH_EXE:-${BUILD_DIR}/tests/performance/physics/fluids/bench_poisson_fluids}"

# --------------------------- presentation -----------------------------------
echo "📁 Repo root:        ${REPO_ROOT}"
echo "🏗️  Build dir:        ${BUILD_DIR}"
echo "📊 Scaling dir:       ${SCALING_DIR}"
echo "🧰 Launcher:          ${MPI_LAUNCHER}"
echo "🚀 MPI pref.:         ${MPI_PREFLAGS}"
echo "🧵 OMP threads set:   ${SCALING_THREADS}"
echo "🧮 Ranks to test:     ${SCALING_RANKS}"
echo "🧊 Grids (N^3):       ${SCALING_GRIDS}"
echo "🧪 Repetitions:       ${SCALING_REPS}"
echo "🐶 PETSc base opts:   ${BASE_PETSC_OPTS}"
echo "➕ PETSc extra opts:  ${PETSC_OPTIONS_EXTRA:-<none>}"
echo

# --------------------------- (optional) build -------------------------------
if [[ "${SCALING_BUILD_FIRST}" == "1" ]]; then
  # Build only perf tests (enables bench_poisson_fluids target)
  BUILD_DIR="${BUILD_DIR}" \
  CMAKE_BUILD_TYPE="${CMAKE_BUILD_TYPE:-Release}" \
  ENABLE_MPI="${ENABLE_MPI:-On}" \
  ENABLE_CUDA="${ENABLE_CUDA:-OFF}" \
  USE_CUDA_UM="${USE_CUDA_UM:-OFF}" \
  BUILD_TESTS=ON \
  TEST_TOGGLES="-DENABLE_TESTS_UNIT=OFF -DENABLE_TESTS_PERF=ON -DENABLE_TESTS_REGRESSION=OFF" \
  EXTRA_CMAKE_ARGS="${EXTRA_CMAKE_ARGS:+$EXTRA_CMAKE_ARGS }${TEST_TOGGLES}${EXTRA_CMAKE_ARGS_USER:+ ${EXTRA_CMAKE_ARGS_USER}}" \
  "${SCRIPT_DIR}/build.sh"
fi
[[ -x "${BENCH_EXE}" ]] || { echo "❌ Bench not found at ${BENCH_EXE}"; exit 1; }

# Git info (best-effort)
GIT_HASH="$(git -C "${REPO_ROOT}" rev-parse --short=12 HEAD 2>/dev/null || echo unknown)"
GIT_TAG="$(git -C "${REPO_ROOT}" describe --tags --dirty 2>/dev/null || echo unknown)"
HOST="$(hostname)"
DATE_ISO="$(date -Iseconds)"

# --------------------------- CSV headers ------------------------------------
BUILD_CSV="${SCALING_DIR}/build_info.csv"
RESULTS_CSV="${SCALING_DIR}/results_scaling.csv"

if [[ ! -s "${BUILD_CSV}" ]]; then
  echo "date_iso,host,git_hash,git_tag,bench_exe,mpi_launcher,mpi_preflags" > "${BUILD_CSV}"
fi
if [[ ! -s "${RESULTS_CSV}" ]]; then
  echo "date_iso,host,git_hash,grid_n,ranks,omp_threads,rep,mean_us,stddev_us,wall_s,ksp_iters,converged_reason,log_stdout,log_stderr" > "${RESULTS_CSV}"
fi

# Append one build row (idempotent if you run multiple times in same tree it's fine)
echo "${DATE_ISO},${HOST},${GIT_HASH},${GIT_TAG},${BENCH_EXE},${MPI_LAUNCHER},\"${MPI_PREFLAGS}\"" >> "${BUILD_CSV}"

# --------------------------- helpers ----------------------------------------
# choose (px,py,pz) s.t. px*py*pz=R and each divides N; prefer cubic
choose_dmda() {
  local N="$1" R="$2" ds=() d best_px=1 best_py=1 best_pz="$R" best_cost=1e18
  for ((d=1; d<=R; ++d)); do (( R % d==0 )) && ds+=("$d"); done
  for px in "${ds[@]}"; do
    (( N % px )) && continue
    for py in "${ds[@]}"; do
      (( (R/px) % py )) && continue
      local pz=$(( R/(px*py) ))
      (( pz<1 || N%py || N%pz )) && continue
      local cost; cost=$(awk -v a="$px" -v b="$py" -v c="$pz" 'BEGIN{m=(a+b+c)/3.0;print((a-m)^2+(b-m)^2+(c-m)^2)}')
      awk -v c1="$cost" -v c2="$best_cost" 'BEGIN{exit !(c1<c2)}' && { best_cost="$cost"; best_px="$px"; best_py="$py"; best_pz="$pz"; }
    done
  done
  echo "$best_px $best_py $best_pz"
}

parse_field() {
  # usage: parse_field mean_us <file1> [file2]
  local key="$1"; shift
  # accept integers, decimals, scientific notation
  grep -Eo "${key}=[+-]?[0-9]*\.?[0-9]+([eE][+-]?[0-9]+)?" "$@" \
    | sed -E "s/.*${key}=//" | tail -n1 || true
}

parse_iters_reason() {
  # Anchored to PETSc KSP summary lines; ignores smoother lines like "maximum iterations=3".
  # usage: parse_iters_reason <stdout_log> [stderr_log]
  local files=("$@")
  local line iters="" reason=""

  # Preferred: KSPConvergedReasonView
  line="$(grep -E 'Linear( [a-zA-Z0-9_+-]+)? solve (converged|did not converge) due to[[:space:]]+[A-Z_]+[[:space:]]*[;, ]?[[:space:]]*iterations[[:space:]]+[0-9]+' "${files[@]}" 2>/dev/null | tail -n1 || true)"
  if [[ -n "$line" ]]; then
    reason="$(sed -E 's/.*due to[[:space:]]+([A-Z_]+).*/\1/' <<<"$line")"
    iters="$(sed -E 's/.*iterations[[:space:]]+([0-9]+).*/\1/' <<<"$line")"
    echo "${iters},${reason}"
    return
  fi

  # Fallback: compact -ksp_view single-line summary
  # "KSP cg CONVERGED_RTOL iterations 37"  or  "... its=37 ..."
  line="$(grep -E 'KSP[[:space:]]+[a-z0-9_+-]+[[:space:]]+((CONVERGED|DIVERGED)_[A-Z_]+).*?(its|iterations)[[:space:]=:]*[0-9]+' "${files[@]}" 2>/dev/null | tail -n1 || true)"
  if [[ -n "$line" ]]; then
    reason="$(sed -E 's/.*KSP[[:space:]]+[a-z0-9_+-]+[[:space:]]+((CONVERGED|DIVERGED)_[A-Z_]+).*/\1/I' <<<"$line")"
    iters="$(sed -E 's/.*(its|iterations)[[:space:]=:]*([0-9]+).*/\2/I' <<<"$line")"
    echo "${iters},${reason}"
    return
  fi

  # No match
  echo ","
}

# --------------------------- sweep ------------------------------------------

run_one() {
  local N="$1" R="$2" T="$3" REP="$4"

  local ts="$(date +%Y%m%d-%H%M%S)"
  local tag="N${N}_r${R}_t${T}_rep${REP}_${ts}"
  local out_log="${LOG_DIR}/stdout_${tag}.log"
  local err_log="${LOG_DIR}/stderr_${tag}.log"

  echo "▶️  N=${N}^3 | ranks=${R} | OMP=${T} | rep=${REP}"
  # Manage OpenMP threads per run
  export OMP_NUM_THREADS="${T}"
  export OMP_PLACES="${OMP_PLACES:-cores}"
  export OMP_PROC_BIND="${OMP_PROC_BIND:-close}"
  export OMP_DYNAMIC="FALSE"

  # --- near-cubic DMDA proc layout (px*py*pz = R) ---
  LAYOUT_OPTS=""
  if [[ "${SCALING_SET_DA_PROCS:-1}" == "1" ]]; then
    read -r PX PY PZ < <(choose_dmda "$N" "$R")
    LAYOUT_OPTS="-da_processors_x ${PX} -da_processors_y ${PY} -da_processors_z ${PZ}"
  fi

  # --- per-run binding policy: relax when using OpenMP >1 thread ---
  local PREF="${MPI_PREFLAGS}"
  if [[ "${T}" -gt 1 ]]; then
    # remove any existing '--bind-to ...' then append our choice
    PREF="$(sed -E 's/--bind-to[[:space:]]+[a-zA-Z0-9_:-]+//g' <<<"${PREF}") --bind-to none"
  fi

  # Time the run in bash (seconds with fractional)
  local start end wall_s
  start="$(date +%s.%N)"
  set +e
  ${MPI_LAUNCHER} ${PREF} ${MPI_NPFLAG} "${R}" \
    "${BENCH_EXE}" \
      ${BASE_PETSC_OPTS} ${PETSC_OPTIONS_EXTRA} \
      ${LAYOUT_OPTS} \
      -nx "${N}" -ny "${N}" -nz "${N}" \
      -da_grid_x "${N}" -da_grid_y "${N}" -da_grid_z "${N}" \
      >"${out_log}" 2>"${err_log}"
  local rc=$?
  set -e
  end="$(date +%s.%N)"
  wall_s="$(python3 - <<PY
import decimal,sys
a=decimal.Decimal('${start}'); b=decimal.Decimal('${end}')
print((b-a).normalize())
PY
)"
  if [[ ${rc} -ne 0 ]]; then
    echo "❗ Run failed (rc=${rc}). See ${out_log} / ${err_log}"
  fi

  local mean_us stddev_us iters reason
  mean_us="$(parse_field 'mean_us' "${out_log}" "${err_log}")"
  stddev_us="$(parse_field 'stddev_us' "${out_log}" "${err_log}")"
  IFS=, read -r iters reason <<<"$(parse_iters_reason "${out_log}" "${err_log}")"

  echo "   ⏱️ mean_us=${mean_us:-NA} stddev_us=${stddev_us:-NA} wall_s=${wall_s} iters=${iters:-NA} reason=${reason:-NA}"

  echo "${DATE_ISO},${HOST},${GIT_HASH},${N},${R},${T},${REP},${mean_us},${stddev_us},${wall_s},${iters},${reason},${out_log},${err_log}" >> "${RESULTS_CSV}"
}


# --------------------------- main loop --------------------------------------
for N in ${SCALING_GRIDS}; do
  for R in ${SCALING_RANKS}; do
    for T in ${SCALING_THREADS}; do
      for REP in $(seq 1 "${SCALING_REPS}"); do
        run_one "${N}" "${R}" "${T}" "${REP}"
      done
    done
  done
done

echo
echo "✅ Scaling sweep complete."
echo "📄 Build CSV:    ${BUILD_CSV}"
echo "📄 Results CSV:  ${RESULTS_CSV}"
echo "🪵 Logs:         ${LOG_DIR}/stdout_*.log, ${LOG_DIR}/stderr_*.log"
