#!/bin/bash
# push_ci.sh — push + start self-hosted runner in tmux + live-stream all job logs
# (minimal changes from your original)
set -euo pipefail

# ---------- Config ----------
SESSION="${PUSHCI_SESSION:-gh-runner}"                    # tmux session name
RUNNER_DIR="${RUNNER_DIR:-/home/actions/actions-runner}"  # where ./config.sh was run (owned by 'actions')
RUNNER_USER="${RUNNER_USER:-actions}"                     # runner account
REMOTE_DEFAULT="${PUSHCI_REMOTE:-origin}"                 # default git remote
WORKFLOW="${PUSHCI_WORKFLOW:-*}"                          # workflow to watch: filename, workflow name, or "*" for any
POLL_TOTAL="${PUSHCI_POLL_TOTAL:-300}"                    # seconds to wait for a suitable run to appear
POLL_STEP="${PUSHCI_POLL_STEP:-3}"                        # seconds between polls
STREAM_POLL="${PUSHCI_STREAM_POLL:-2}"                    # NEW: seconds between log polls
TMUX_PATH="/usr/bin/tmux"                                 # path allowed in sudoers
# ----------------------------

die() { printf '%s\n' "Error: $*" >&2; exit 1; }
need() { command -v "$1" >/dev/null 2>&1 || die "Missing dependency: $1"; }
msg()  { printf '\033[1;36m==> %s\033[0m\n' "$*"; }

need git
need gh
[ -x "$TMUX_PATH" ] || die "tmux not found at $TMUX_PATH"
git rev-parse --is-inside-work-tree >/dev/null 2>&1 || die "Run this from inside a git repo."
[ -d "$RUNNER_DIR" ] || die "Runner dir '$RUNNER_DIR' not found."
[ -x "$RUNNER_DIR/run.sh" ] || die "run.sh missing in '$RUNNER_DIR' (did you run ./config.sh as $RUNNER_USER?)"
gh auth status >/dev/null 2>&1 || die "GitHub CLI not authenticated. Run: gh auth login"

# sudo helpers
SUDO_N="sudo -n -u $RUNNER_USER"
SUDO_I="sudo    -u $RUNNER_USER"
if $SUDO_N "$TMUX_PATH" -V >/dev/null 2>&1; then
  TMUX="$SUDO_N $TMUX_PATH"
else
  echo "(info) Non-interactive sudo not available; you may be prompted for your password."
  TMUX="$SUDO_I $TMUX_PATH"
fi

BRANCH="$(git rev-parse --abbrev-ref HEAD)"
HEAD_BEFORE="$(git rev-parse HEAD || true)"
UPSTREAM="$(git rev-parse --abbrev-ref --symbolic-full-name @\{u\} 2>/dev/null || true)"

# Ensure cleanup if we started a tmux session
STARTED=0
cleanup() {
  if [ "$STARTED" -eq 1 ]; then
    msg "Stopping runner session '$SESSION'…"
    $TMUX kill-session -t "$SESSION" >/dev/null 2>&1 || true
  fi
}
trap cleanup EXIT INT TERM

# Start runner in tmux if needed (same as your original)
if ! $TMUX has-session -t "$SESSION" 2>/dev/null; then
  msg "Starting runner session '$SESSION' as $RUNNER_USER…"
  $TMUX new -ds "$SESSION" -c "$RUNNER_DIR" "./run.sh"
  STARTED=1
else
  msg "Runner already active (session '$SESSION')."
fi

# Push (same as your original)
if [ -z "${UPSTREAM:-}" ]; then
  msg "No upstream; pushing HEAD to ${REMOTE_DEFAULT} and setting upstream…"
  git push -u "$REMOTE_DEFAULT" HEAD
else
  msg "Pushing to $UPSTREAM…"
  git push
fi

HEAD_SHA="$(git rev-parse HEAD)"
[ -n "${HEAD_BEFORE:-}" ] && [ "$HEAD_SHA" = "$HEAD_BEFORE" ] && \
  msg "Note: commit SHA unchanged; will watch the newest run on '$BRANCH'."

# Build the workflow filter (empty for "*")
WF_FLAG=()
[ "$WORKFLOW" != "*" ] && WF_FLAG=(--workflow "$WORKFLOW")

msg "Waiting for workflow run (workflow: ${WORKFLOW}, branch: $BRANCH, sha: $HEAD_SHA)…"
RUN_ID=""
deadline=$(($(date +%s) + POLL_TOTAL))
FIELDS="databaseId,workflowName,displayTitle,headSha,headBranch,status,conclusion,createdAt"

while [ -z "$RUN_ID" ] && [ "$(date +%s)" -lt "$deadline" ]; do
  # Selection strategy (same as your original):
  # 1) newest non-completed (queued|in_progress) run for this SHA+branch
  # 2) else newest completed run that isn't cancelled
  RUN_ID="$(gh run list "${WF_FLAG[@]}" --limit 40 --json $FIELDS \
    --jq '
      ( map(select(.headSha=="'"$HEAD_SHA"'" and .headBranch=="'"$BRANCH"'"))
        | sort_by(.createdAt) | reverse ) as $r
      | ( $r | map(select(.status!="completed")) | .[0] //
          ( $r | map(select(.status=="completed" and .conclusion!="cancelled")) | .[0] )
        ) // empty
      | .databaseId
    ' 2>/dev/null || true)"
  [ -n "$RUN_ID" ] || sleep "$POLL_STEP"
done

if [ -z "$RUN_ID" ]; then
  msg "No suitable run found yet. Tip: increase POLL_TOTAL or set PUSHCI_WORKFLOW=linux.yml"
  exit 1
fi

INFO="$(gh run view "$RUN_ID" --json workflowName,displayTitle,status,conclusion \
       --jq '.workflowName + " · " + .displayTitle + " (" + .status + (if .status=="completed" then "/" + .conclusion else "" end) + ")"' 2>/dev/null || true)"
[ -n "$INFO" ] && msg "Streaming logs for: $INFO (id: $RUN_ID)…" || msg "Streaming logs for run $RUN_ID…"

# ----------------------------
# Live streaming of all job logs with job-name prefixes
# ----------------------------
stream_all_logs() {
  local rid="$1"
  local tmpdir
  tmpdir="$(mktemp -d -t pushci-logs-XXXXXX)"
  msg "Log cache: $tmpdir"
  declare -A sizes   # bytes printed per job id
  declare -A names   # job id -> job name

  # little helper to color headings
  color_hdr() { printf '\033[1;35m[%s]\033[0m\n' "$1"; }

  local run_status=""
  local run_conclusion=""

  while :; do
    # List all jobs seen so far
    while IFS=$'\t' read -r jid jname jstatus || [ -n "${jid:-}" ]; do
      [ -n "${jid:-}" ] || continue
      names["$jid"]="$jname"

      # Fetch full current log snapshot for this job (may be empty if not started yet)
      if gh run view "$rid" --job "$jid" --log > "$tmpdir/$jid.new" 2>/dev/null; then
        local prev="${sizes[$jid]:-0}"
        local size
        size="$(wc -c < "$tmpdir/$jid.new" | tr -d ' ')"
        # If job restarted and log shrank, print from beginning again
        if [ "$size" -lt "$prev" ]; then
          prev=0
          color_hdr "$jname (restarted)"
        fi
        if [ "$size" -gt "$prev" ]; then
          # Print only new bytes, prefixing each line with the job name
          tail -c +$((prev+1)) "$tmpdir/$jid.new" | sed -u "s/^/[$jname] /"
          sizes["$jid"]="$size"
        fi
        mv -f "$tmpdir/$jid.new" "$tmpdir/$jid.cur" 2>/dev/null || true
      fi
    done < <(gh run view "$rid" --json jobs \
              --jq '.jobs[]? | "\(.databaseId)\t\(.name)\t\(.status)"' 2>/dev/null || true)

    run_status="$(gh run view "$rid" --json status --jq .status 2>/dev/null || echo "")"
    if [ "$run_status" = "completed" ]; then
      run_conclusion="$(gh run view "$rid" --json conclusion --jq .conclusion 2>/dev/null || echo "")"
      break
    fi
    sleep "$STREAM_POLL"
  done

  # Final flush: if anything new appeared between last poll and completion
  while IFS=$'\t' read -r jid jname _ || [ -n "${jid:-}" ]; do
    [ -n "${jid:-}" ] || continue
    names["$jid"]="$jname"
    if gh run view "$rid" --job "$jid" --log > "$tmpdir/$jid.final" 2>/dev/null; then
      local prev="${sizes[$jid]:-0}"
      local size
      size="$(wc -c < "$tmpdir/$jid.final" | tr -d ' ')"
      if [ "$size" -gt "$prev" ]; then
        tail -c +$((prev+1)) "$tmpdir/$jid.final" | sed -u "s/^/[$jname] /"
      fi
    fi
  done < <(gh run view "$rid" --json jobs \
            --jq '.jobs[]? | "\(.databaseId)\t\(.name)\t\(.status)"' 2>/dev/null || true)

  # Summary + exit with the run's conclusion
  run_conclusion="$(gh run view "$rid" --json conclusion --jq .conclusion 2>/dev/null || echo "")"
  case "${run_conclusion:-}" in
    success|"") msg "Run $rid finished: ${run_conclusion:-success}"; return 0 ;;
    *)          msg "Run $rid finished: $run_conclusion"; return 1 ;;
  esac
}

# Replace the previous 'gh run watch --exit-status' with the merged streaming
STATUS=0
if ! stream_all_logs "$RUN_ID"; then
  STATUS=1
fi

exit "$STATUS"
