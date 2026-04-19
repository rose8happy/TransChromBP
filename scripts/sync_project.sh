#!/bin/bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
BRANCH="master"

REMOTE_6000_USER="zhoujiazhen"
REMOTE_6000_HOST="127.0.0.1"
REMOTE_6000_PORT="6000"
REMOTE_6000_ROOT="/data1/zhoujiazhen/bylw_atac"
OFFICIAL_6000_DIR="$REMOTE_6000_ROOT/chromBPNet"
RUNTIME_6000_DIR="$REMOTE_6000_ROOT/TransChromBP"

REMOTE_6002_USER="zhengwei"
REMOTE_6002_HOST="127.0.0.1"
REMOTE_6002_PORT="6002"
REMOTE_6002_KEY="/home/zhengwei/.ssh/codex_6002_ed25519"
REMOTE_6002_ROOT="/home/zhengwei/bylw_atac"
RUNTIME_6002_DIR="$REMOTE_6002_ROOT/TransChromBP"

SSH_RSYNC_6000="ssh -p $REMOTE_6000_PORT"
SSH_RSYNC_6002="ssh -i $REMOTE_6002_KEY -p $REMOTE_6002_PORT"

ACTION="${1:-help}"
if [[ $# -gt 0 ]]; then
    shift
fi

DRY_RUN=0
for arg in "$@"; do
    case "$arg" in
        --dry-run)
            DRY_RUN=1
            ;;
        *)
            echo "Unknown option: $arg" >&2
            exit 1
            ;;
    esac
done

RSYNC_ARGS=("-avz" "--progress")
if [[ "$DRY_RUN" -eq 1 ]]; then
    RSYNC_ARGS+=("--dry-run")
fi

PUBLISH_EXCLUDES=(
    "--exclude=.agents/"
    "--exclude=.claude/"
    "--exclude=.codex"
    "--exclude=.codex_remote_edit/"
    "--exclude=.git/"
    "--exclude=.gitignore"
    "--exclude=.idea/"
    "--exclude=.pytest_cache/"
    "--exclude=.venv/"
    "--exclude=.venv-report/"
    "--exclude=__pycache__/"
    "--exclude=chrombpnet_data/"
    "--exclude=data/"
    "--exclude=logs/"
    "--exclude=outputs/"
    "--exclude=temp/"
    "--exclude=tmp/"
    "--exclude=tmp_remote_edit/"
    "--exclude=trained_models/"
    "--exclude=wandb/"
    "--exclude=references/local-only/"
    "--exclude=*.bam"
    "--exclude=*.bed"
    "--exclude=*.bed.gz"
    "--exclude=*.bigWig"
    "--exclude=*.bw"
    "--exclude=*.fa"
    "--exclude=*.fai"
    "--exclude=*.h5"
    "--exclude=*.pdf"
    "--exclude=*.pt"
    "--exclude=*.pyc"
)

REPORT_EXCLUDES=(
    "--exclude=*.aux"
    "--exclude=*.fdb_latexmk"
    "--exclude=*.fls"
    "--exclude=*.log"
    "--exclude=*.out"
    "--exclude=*.pdf"
    "--exclude=*.toc"
    "--exclude=*.xdv"
)

ssh_6000() {
    ssh -p "$REMOTE_6000_PORT" "$REMOTE_6000_USER@$REMOTE_6000_HOST" "$@"
}

ssh_6002() {
    ssh -i "$REMOTE_6002_KEY" -p "$REMOTE_6002_PORT" "$REMOTE_6002_USER@$REMOTE_6002_HOST" "$@"
}

print_help() {
    cat <<EOF
Usage: $0 <command> [--dry-run]

Repository topology is documented in docs/env/repository_governance.md.

Commands:
  help                  Show this help text
  status                Show local git status and origin
  status-all            Show local status plus the canonical 6000/6002 roles
  pull                  Git pull from origin/$BRANCH
  push                  Git commit (optional) and push to origin/$BRANCH
  publish-runtime-6000  Publish the local canonical repo to 6000 runtime ($RUNTIME_6000_DIR)
  publish-runtime-6002  Publish the local canonical repo to 6002 runtime ($RUNTIME_6002_DIR)
  pull-results-6000     Pull logs and report sources from 6000 runtime
  pull-results-6002     Pull logs and report sources from 6002 runtime

Deprecated commands:
  deploy                replaced by publish-runtime-6000 / publish-runtime-6002
  download_results      replaced by pull-results-6000 / pull-results-6002
EOF
}

print_local_status() {
    echo "== Local canonical trunk =="
    git -C "$REPO_ROOT" status --short --branch
    echo "origin: $(git -C "$REPO_ROOT" remote get-url origin)"
}

print_remote_status() {
    local label=$1
    local role=$2
    local cmd=$3
    local path=$4

    echo "== $label =="
    echo "role: $role"
    if "$cmd" "test -d '$path' && printf 'path: %s\n' '$path' || printf 'missing: %s\n' '$path'"; then
        :
    else
        echo "unreachable"
    fi
}

publish_runtime() {
    local label=$1
    local ssh_cmd=$2
    local target=$3

    echo "== Publish $label =="
    echo "source: $REPO_ROOT"
    echo "target: $target"
    if [[ "$DRY_RUN" -eq 1 ]]; then
        echo "mode: dry-run"
    fi

    rsync "${RSYNC_ARGS[@]}" \
        --delete \
        "${PUBLISH_EXCLUDES[@]}" \
        -e "$ssh_cmd" \
        "$REPO_ROOT/" "$target/"
}

pull_results() {
    local label=$1
    local ssh_cmd=$2
    local remote_root=$3
    local machine_dir=$4
    local logs_dest="$REPO_ROOT/logs/$machine_dir"
    local reports_dest="$REPO_ROOT/tmp_remote_edit/results/$machine_dir/reports"

    echo "== Pull results $label =="
    echo "remote root: $remote_root"
    echo "logs -> $logs_dest"
    echo "reports -> $reports_dest"
    if [[ "$DRY_RUN" -eq 1 ]]; then
        echo "mode: dry-run"
    fi

    mkdir -p "$logs_dest" "$reports_dest"

    rsync "${RSYNC_ARGS[@]}" -e "$ssh_cmd" \
        "$remote_root/logs/" "$logs_dest/"

    rsync "${RSYNC_ARGS[@]}" \
        "${REPORT_EXCLUDES[@]}" \
        -e "$ssh_cmd" \
        "$remote_root/reports/" "$reports_dest/"
}

deprecated_command() {
    local legacy=$1
    local replacement=$2
    echo "Command '$legacy' is deprecated. Use $replacement instead." >&2
    exit 2
}

case "$ACTION" in
    help)
        print_help
        ;;
    status)
        print_local_status
        ;;
    status-all)
        print_local_status
        print_remote_status "6000 official lookup" "official ChromBPNet lookup/repro" ssh_6000 "$OFFICIAL_6000_DIR"
        print_remote_status "6000 runtime" "TransChromBP runtime workspace" ssh_6000 "$RUNTIME_6000_DIR"
        print_remote_status "6002 runtime" "TransChromBP runtime mirror" ssh_6002 "$RUNTIME_6002_DIR"
        ;;
    pull)
        git -C "$REPO_ROOT" pull origin "$BRANCH"
        ;;
    push)
        if [[ -n "$(git -C "$REPO_ROOT" status -s)" ]]; then
            echo "You have uncommitted changes:"
            git -C "$REPO_ROOT" status -s
            read -r -p "Do you want to commit them now? (y/n) " reply
            if [[ "$reply" =~ ^[Yy]$ ]]; then
                read -r -p "Enter commit message: " msg
                git -C "$REPO_ROOT" add .
                git -C "$REPO_ROOT" commit -m "$msg"
            fi
        fi
        git -C "$REPO_ROOT" push origin "$BRANCH"
        ;;
    publish-runtime-6000)
        publish_runtime "6000 runtime" "$SSH_RSYNC_6000" "$REMOTE_6000_USER@$REMOTE_6000_HOST:$RUNTIME_6000_DIR"
        ;;
    publish-runtime-6002)
        publish_runtime "6002 runtime" "$SSH_RSYNC_6002" "$REMOTE_6002_USER@$REMOTE_6002_HOST:$RUNTIME_6002_DIR"
        ;;
    pull-results-6000)
        pull_results "6000 runtime" "$SSH_RSYNC_6000" \
            "$REMOTE_6000_USER@$REMOTE_6000_HOST:$RUNTIME_6000_DIR" "6000"
        ;;
    pull-results-6002)
        pull_results "6002 runtime" "$SSH_RSYNC_6002" \
            "$REMOTE_6002_USER@$REMOTE_6002_HOST:$RUNTIME_6002_DIR" "6002"
        ;;
    deploy)
        deprecated_command "deploy" "publish-runtime-6000 or publish-runtime-6002"
        ;;
    download_results)
        deprecated_command "download_results" "pull-results-6000 or pull-results-6002"
        ;;
    *)
        echo "Unknown command: $ACTION" >&2
        print_help >&2
        exit 1
        ;;
esac
