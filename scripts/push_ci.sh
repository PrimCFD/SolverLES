#!/bin/bash
set -eu

# ---------- Config ----------
SESSION="${PUSHCI_SESSION:-gh-runner}"                   # tmux session name for the runner
RUNNER_DIR="${RUNNER_DIR:-/home/actions/actions-runner}" # where ./config.sh was run (owned by 'actions')
RUNNER_USER="${RUNNER_USER:-actions}"                    # runner account
REMOTE_DEFAULT="${PUSHCI_REMOTE:-origin}"                # default git remote
POLL_TOTAL="${PUSHCI_POLL_TOTAL:-240}"                   # seconds to wait for run to appear
POLL_STEP="${PUSHCI_POLL_STEP:-3}"                       # seconds between polls
TMUX_PATH="/usr/bin/tmux"                                # path allowed in sudoers
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

# sudo helpers (prefer non-interactive; fall back to interactive if sudoers not set)
SUDO_N="sudo -n -u $RUNNER_USER"
SUDO_I="sudo    -u $RUNNER_USER"
if $SUDO_N "$TMUX_PATH" -V >/dev/null 2>&1; then
  TMUX="$SUDO_N $TMUX_PATH"
else
  echo "(info) Non-interactive sudo not available; you may be prompted for your password."
  TMUX="$SUDO_I $TMUX_PATH"
fi

# Determine branch/upstream & snapshot HEAD before push
BRANCH="$(git rev-parse --abbrev-ref HEAD)"
HEAD_BEFORE="$(git rev-parse HEAD || true)"
UPSTREAM="$(git rev-parse --abbrev-ref --symbolic-full-name @\{u\} 2>/dev/null || true)"

# Start runner if needed (tmux runs as 'actions'; -c sets working dir)
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

# Find run for this commit (prefer exact SHA)
HEAD_SHA="$(git rev-parse HEAD)"
if [ -n "${HEAD_BEFORE:-}" ] && [ "$HEAD_SHA" = "$HEAD_BEFORE" ]; then
  msg "Note: commit SHA unchanged; will watch most recent run on '$BRANCH'."
fi

msg "Waiting for workflow run (branch: $BRANCH, sha: $HEAD_SHA)…"
RUN_ID=""
deadline=$(($(date +%s) + POLL_TOTAL))
while [ -z "$RUN_ID" ] && [ "$(date +%s)" -lt "$deadline" ]; do
  RUN_ID="$(gh run list --limit 20 --json databaseId,headSha,headBranch \
            --jq 'map(select(.headSha=="'"$HEAD_SHA"'" and .headBranch=="'"$BRANCH"'")) | .[0].databaseId' 2>/dev/null || true)"
  if [ -z "$RUN_ID" ]; then
    RUN_ID="$(gh run list --branch "$BRANCH" --limit 1 --json databaseId --jq '.[0].databaseId' 2>/dev/null || true)"
  fi
  [ -n "$RUN_ID" ] || sleep "$POLL_STEP"
done

if [ -z "$RUN_ID" ]; then
  msg "No run found yet. Try: gh run list --branch \"$BRANCH\""
  if [ "$STARTED" -eq 1 ]; then $TMUX kill-session -t "$SESSION" || true; fi
  exit 0
fi

# Stream logs; exit with run status
msg "Streaming logs for run $RUN_ID…"
if gh run watch "$RUN_ID" --exit-status; then
  STATUS=0
else
  STATUS=$?
fi

# Stop runner if we started it
if [ "$STARTED" -eq 1 ]; then
  msg "Stopping runner session '$SESSION'…"
  $TMUX kill-session -t "$SESSION" || true
fi

exit "$STATUS"