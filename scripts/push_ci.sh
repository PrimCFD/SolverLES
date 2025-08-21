#!/bin/bash
# push_ci.sh — push + start self-hosted runner in tmux; no log streaming
# Prints the run URL (and opens it if possible). No tmux/stdout streaming.
set -euo pipefail

# ---------- Config ----------
SESSION="${PUSHCI_SESSION:-gh-runner}"                    # tmux session name
RUNNER_DIR="${RUNNER_DIR:-/home/actions/actions-runner}"  # where ./config.sh was run (owned by 'actions')
RUNNER_USER="${RUNNER_USER:-actions}"                     # runner account
REMOTE_DEFAULT="${PUSHCI_REMOTE:-origin}"                 # default git remote
WORKFLOW="${PUSHCI_WORKFLOW:-*}"                          # workflow to watch: filename, workflow name, or "*" for any
POLL_TOTAL="${PUSHCI_POLL_TOTAL:-300}"                    # seconds to wait for a suitable run to appear
POLL_STEP="${PUSHCI_POLL_STEP:-3}"                        # seconds between polls
OPEN_WEB="${PUSHCI_OPEN_WEB:-1}"                          # open run page in browser (1=yes, 0=no)
TMUX_PATH="/usr/bin/tmux"
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

SUDO_N="sudo -n -u $RUNNER_USER"
SUDO_I="sudo    -u $RUNNER_USER"
if $SUDO_N "$TMUX_PATH" -V >/dev/null 2>&1; then TMUX="$SUDO_N $TMUX_PATH"; else TMUX="$SUDO_I $TMUX_PATH"; fi

BRANCH="$(git rev-parse --abbrev-ref HEAD)"
UPSTREAM="$(git rev-parse --abbrev-ref --symbolic-full-name @\{u\} 2>/dev/null || true)"

# Start runner in tmux if needed
if ! $TMUX has-session -t "$SESSION" 2>/dev/null; then
  msg "Starting runner session '$SESSION' as $RUNNER_USER…"
  $TMUX new -ds "$SESSION" -c "$RUNNER_DIR" "./run.sh"
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
# Build the workflow filter (empty for "*")
WF_FLAG=()
[ "$WORKFLOW" != "*" ] && WF_FLAG=(--workflow "$WORKFLOW")

# Find the run id
msg "Locating the newest run for $BRANCH @ $HEAD_SHA…"
RUN_ID=""
deadline=$(($(date +%s) + POLL_TOTAL))
FIELDS="databaseId,workflowName,displayTitle,headSha,headBranch,status,conclusion,createdAt,url"

while [ -z "$RUN_ID" ] && [ "$(date +%s)" -lt "$deadline" ]; do
  RUN_INFO="$(gh run list "${WF_FLAG[@]}" --limit 40 --json $FIELDS \
    --jq '
      ( map(select(.headSha=="'"$HEAD_SHA"'" and .headBranch=="'"$BRANCH"'"))
        | sort_by(.createdAt) | reverse ) as $r
      | ( $r | map(select(.status!="completed")) | .[0] //
          ( $r | map(select(.status=="completed" and .conclusion!="cancelled")) | .[0] )
        ) // empty
    ' 2>/dev/null || true)"
  RUN_ID="$(jq -r '.databaseId // empty' <<<"$RUN_INFO" 2>/dev/null || true)"
  [ -n "$RUN_ID" ] || sleep "$POLL_STEP"
done

if [ -z "$RUN_ID" ]; then
  msg "No suitable run found. You can open the Actions tab to monitor."
  exit 0
fi

RUN_URL="$(jq -r '.url' <<<"$RUN_INFO" 2>/dev/null || gh run view "$RUN_ID" --json url --jq .url)"
msg "Run URL: $RUN_URL"
if [ "$OPEN_WEB" = "1" ]; then
  gh run view "$RUN_ID" --web >/dev/null 2>&1 || true
fi

exit 0
