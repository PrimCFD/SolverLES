#!/bin/bash
set -eu

# ---------- Config ----------
SESSION="${PUSHCI_SESSION:-gh-runner}"                    # tmux session name
RUNNER_DIR="${RUNNER_DIR:-/home/actions/actions-runner}"  # where ./config.sh was run (owned by 'actions')
RUNNER_USER="${RUNNER_USER:-actions}"                     # runner account
REMOTE_DEFAULT="${PUSHCI_REMOTE:-origin}"                 # default git remote
WORKFLOW="${PUSHCI_WORKFLOW:-*}"                          # workflow to watch: filename, workflow name, or "*" for any
POLL_TOTAL="${PUSHCI_POLL_TOTAL:-300}"                    # seconds to wait for a suitable run to appear
POLL_STEP="${PUSHCI_POLL_STEP:-3}"                        # seconds between polls
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

# Start runner in tmux if needed
STARTED=0
if ! $TMUX has-session -t "$SESSION" 2>/dev/null; then
  msg "Starting runner session '$SESSION' as $RUNNER_USER…"
  $TMUX new -ds "$SESSION" -c "$RUNNER_DIR" "./run.sh"
  STARTED=1
else
  msg "Runner already active (session '$SESSION')."
fi

# Push (set upstream if missing)
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
  # Grab the latest ~40 runs (optionally filtered by workflow)
  # Selection strategy:
  #   1) Prefer newest NON-completed run for this SHA+branch (queued|in_progress)
  #   2) Else prefer newest COMPLETED run whose conclusion != "cancelled"
  #   3) Ignore cancelled runs (keep polling) to avoid latching on to concurrency-cancelled runs
  RUN_ID="$(gh run list "${WF_FLAG[@]}" --limit 40 --json $FIELDS \
    --jq '
      # all runs for this SHA+branch, newest first
      ( map(select(.headSha=="'"$HEAD_SHA"'" and .headBranch=="'"$BRANCH"'"))
        | sort_by(.createdAt) | reverse ) as $r
      | ( $r | map(select(.status!="completed")) | .[0] //    # queued/in_progress
          ( $r | map(select(.status=="completed" and .conclusion!="cancelled")) | .[0] ) # finished but not cancelled
        ) // empty
      | .databaseId
    ' 2>/dev/null || true)"

  [ -n "$RUN_ID" ] || sleep "$POLL_STEP"
done

if [ -z "$RUN_ID" ]; then
  msg "No suitable run found yet. Tip: increase POLL_TOTAL or set PUSHCI_WORKFLOW=linux.yml"
  [ "$STARTED" -eq 1 ] && $TMUX kill-session -t "$SESSION" || true
  exit 1
fi

# Nice log header about what we're going to watch
INFO="$(gh run view "$RUN_ID" --json workflowName,displayTitle,status,conclusion \
       --jq '.workflowName + " · " + .displayTitle + " (" + .status + (if .status=="completed" then "/" + .conclusion else "" end) + ")"' 2>/dev/null || true)"
[ -n "$INFO" ] && msg "Streaming logs for: $INFO (id: $RUN_ID)…" || msg "Streaming logs for run $RUN_ID…"

# Stream logs; exit with run status
gh run watch "$RUN_ID" --exit-status || STATUS=$?; STATUS=${STATUS:-0}

# Stop runner if we started it
if [ "$STARTED" -eq 1 ]; then
  msg "Stopping runner session '$SESSION'…"
  $TMUX kill-session -t "$SESSION" || true
fi

exit "$STATUS"