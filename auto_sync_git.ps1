# Auto-sync script for Git repository
# This script monitors file changes and automatically commits and pushes to GitHub

$repoPath = "d:\work"
$checkInterval = 30  # seconds

function Test-GitChanges {
    cd $repoPath
    $status = git status --porcelain
    return $status.Length -gt 0
}

function Invoke-GitCommitAndPush {
    cd $repoPath
    $timestamp = Get-Date -Format "yyyy-MM-dd HH:mm:ss"
    
    # Add all changes
    git add .
    
    # Commit with timestamp
    git commit -m "Auto-sync: $timestamp"
    
    # Push to remote (using token in remote URL)
    git push origin master
}

Write-Host "Starting auto-sync monitoring for $repoPath..." -ForegroundColor Green
Write-Host "Press Ctrl+C to stop" -ForegroundColor Yellow

while ($true) {
    if (Test-GitChanges) {
        Write-Host "Changes detected at $(Get-Date -Format 'HH:mm:ss')" -ForegroundColor Cyan
        Invoke-GitCommitAndPush
        Write-Host "Sync completed" -ForegroundColor Green
    }
    
    Start-Sleep -Seconds $checkInterval
}
