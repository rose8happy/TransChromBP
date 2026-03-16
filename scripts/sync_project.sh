#!/bin/bash
# Sync script for TransChromBP project
# Usage: ./scripts/sync_project.sh [pull|push|status]

ACTION=${1:-status}
BRANCH="master"  # or main, adjust as needed

# Function to check if we are on local or remote
is_remote() {
    if [[ "$(hostname)" == *"6000"* ]] || [[ "$PWD" == *"/data1/"* ]]; then
        return 0
    else
        return 1
    fi
}

echo "=== TransChromBP Project Sync ($ACTION) ==="
echo "Current directory: $PWD"

if [ "$ACTION" == "pull" ]; then
    echo "⬇️  Pulling latest changes from origin/$BRANCH..."
    git pull origin $BRANCH
elif [ "$ACTION" == "push" ]; then
    echo "⬆️  Pushing local changes to origin/$BRANCH..."
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
elif [ "$ACTION" == "status" ]; then
    echo "ℹ️  Git Status:"
    git status
    echo ""
    echo "Remote 'origin': $(git remote get-url origin)"
else
    echo "Usage: $0 [pull|push|status]"
    exit 1
fi

echo "✅ Done."
