#!/bin/bash
# push_ci.sh — push + start self-hosted runner in tmux + status-only CI watcher
# Quiet stream: only job stage transitions and final summary.
set -euo pipefail

# ---------- Config ----------
SESSION="${PUSHCI_SESSION:-gh-runner}"                    # tmux session name
RUNNER_DIR="${RUNNER_DIR:-/home/actions/actions-runner}"  # where ./config.sh was run (owned by 'actions')
RUNNER_USER="${RUNNER_USER:-actions}"                     # runner account
REMOTE_DEFAULT="${PUSHCI_REMOTE:-origin}"                 # default git remote
WORKFLOW="${PUSHCI_WORKFLOW:-*}"                          # workflow to watch: filename, workflow name, or "*" for any
POLL_TOTAL="${PUSHCI_POLL_TOTAL:-300}"                    # seconds to wait for a suitable run to appear
POLL_STEP="${PUSHCI_POLL_STEP:-3}"                        # seconds between polls
STREAM_POLL="${PUSHCI_STREAM_POLL:-3}"                    # seconds between status polls
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

# Start runner in tmux if needed
if ! $TMUX has-session -t "$SESSION" 2>/dev/null; then
  msg "Starting runner session '$SESSION' as $RUNNER_USER…"
  $TMUX new -ds "$SESSION" -c "$RUNNER_DIR" "./run.sh"
  STARTED=1
else
  msg "Runner already active (session '$SESSION')."
fi

# Push
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
[ -n "$INFO" ] && msg "Watching: $INFO (id: $RUN_ID)…" || msg "Watching run $RUN_ID…"

# -------- status-only watcher (quiet tmux output) --------
watch_status() {
  local rid="$1" interval="$2"
  declare -A last
  echo "Status-only stream (every ${interval}s):"
  while :; do
    while IFS=$'\t' read -r jid name status concl || [ -n "$jid" ]; do
      [ -n "$jid" ] || continue
      cur="$status/${concl:-}"
      if [[ "${last[$jid]:-}" != "$cur" ]]; then
        if [[ "$status" == "completed" && -n "$concl" ]]; then
          printf '[%s] %s → %s\n' "$name" "$status" "$concl"
        else
          printf '[%s] %s\n' "$name" "$status"
        fi
        last[$jid]="$cur"
      fi
    done < <(gh run view "$rid" --json jobs \
              --jq '.jobs[]? | "\(.databaseId)\t\(.name)\t\(.status)\t\(.conclusion//"")"' 2>/dev/null || true)

    run_status="$(gh run view "$rid" --json status --jq .status 2>/dev/null || echo "")"
    [[ "$run_status" == "completed" ]] && break
    sleep "$interval"
  done

  echo
  echo "----- Workflow summary -----"
  gh run view "$rid" --summary || true
  echo "----------------------------"

  conclusion="$(gh run view "$rid" --json conclusion --jq .conclusion 2>/dev/null || echo "")"
  [[ "$conclusion" == "success" || -z "$conclusion" ]]
}
# --------------------------------------------------------

if watch_status "$RUN_ID" "$STREAM_POLL"; then
  exit 0
else
  exit 1
fi
