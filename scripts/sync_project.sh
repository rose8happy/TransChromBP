#!/bin/bash
# Sync script for TransChromBP project
# Usage: ./scripts/sync_project.sh [pull|push|deploy|download_results|status]

ACTION=${1:-status}
BRANCH="master"  # or main, adjust as needed

# Remote Configuration
REMOTE_USER="zhoujiazhen"
REMOTE_HOST="127.0.0.1"
REMOTE_PORT="6000"
REMOTE_DIR="/data1/zhoujiazhen/bylw_atac/chromBPNet"
SSH_CMD="ssh -p $REMOTE_PORT"

echo "=== TransChromBP Project Sync ($ACTION) ==="
echo "Current directory: $PWD"

if [ "$ACTION" == "pull" ]; then
    echo "⬇️  [Git] Pulling latest changes from origin/$BRANCH..."
    git pull origin $BRANCH

elif [ "$ACTION" == "push" ]; then
    echo "⬆️  [Git] Pushing local changes to origin/$BRANCH..."
    # Check for uncommitted changes
    if [[ -n $(git status -s) ]]; then
        echo "⚠️  You have uncommitted changes:"
        git status -s
        read -p "Do you want to commit them now? (y/n) " -n 1 -r
        echo
        if [[ $REPLY =~ ^[Yy]$ ]]; then
            read -p "Enter commit message: " msg
            git add .
            git commit -m "$msg"
        fi
    fi
    git push origin $BRANCH

elif [ "$ACTION" == "deploy" ]; then
    echo "🚀 [Rsync] Deploying local code to remote server (Fast Sync)..."
    echo "Target: $REMOTE_USER@$REMOTE_HOST:$REMOTE_DIR"
    
    # Exclude heavy/generated files to speed up sync and avoid overwriting server data
    rsync -avz --progress \
        --exclude='.git/' \
        --exclude='.gitignore' \
        --exclude='logs/' \
        --exclude='outputs/' \
        --exclude='data/' \
        --exclude='chrombpnet_data/' \
        --exclude='trained_models/' \
        --exclude='*.h5' \
        --exclude='*.bam' \
        --exclude='*.bigWig' \
        --exclude='*.pdf' \
        --exclude='__pycache__/' \
        --exclude='*.pyc' \
        --exclude='wandb/' \
        -e "$SSH_CMD" \
        ./ $REMOTE_USER@$REMOTE_HOST:$REMOTE_DIR/
        
    echo "✅ Deploy complete. You can now run the code on the server."

elif [ "$ACTION" == "download_results" ]; then
    echo "XB [Rsync] Downloading training results (logs & reports) from server..."
    
    # Sync logs
    echo "📄 Syncing logs..."
    mkdir -p logs
    rsync -avz --progress -e "$SSH_CMD" \
        $REMOTE_USER@$REMOTE_HOST:$REMOTE_DIR/logs/ ./logs/
        
    # Sync reports
    echo "QC Syncing reports..."
    mkdir -p reports
    rsync -avz --progress -e "$SSH_CMD" \
        $REMOTE_USER@$REMOTE_HOST:$REMOTE_DIR/reports/ ./reports/
        
    echo "✅ Download complete. Check ./logs and ./reports directories."

elif [ "$ACTION" == "status" ]; then
    echo "ℹ️  Git Status:"
    git status
    echo ""
    echo "Remote 'origin': $(git remote get-url origin)"
    echo "Deploy Target: $REMOTE_USER@$REMOTE_HOST:$REMOTE_DIR"

else
    echo "Usage: $0 [pull|push|deploy|download_results|status]"
    echo ""
    echo "Commands:"
    echo "  pull             : Git pull from remote repository"
    echo "  push             : Git commit (optional) and push to remote repository"
    echo "  deploy           : Rsync local code to server (Fast, no commit needed)"
    echo "  download_results : Rsync logs and reports from server to local"
    echo "  status           : Show current git status"
    exit 1
fi

echo "✅ Done."
