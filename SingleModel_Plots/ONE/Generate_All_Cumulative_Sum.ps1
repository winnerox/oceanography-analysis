# Batch generate cumulative sum plots for all datasets and states

# Dataset list
$datasets = @('EN4', 'IAP', 'Ishii')

# State list
$states = @('Average', 'StdRef')

# MATLAB script path
$matlabScript = "cd('d:\work\SingleModel_Plots\ONE');"

# Loop to generate all plots
foreach ($dataset in $datasets) {
    foreach ($state in $states) {
        Write-Host "Generating cumulative sum plots for ${dataset} ${state}..." -ForegroundColor Green
        
        # Build MATLAB command
        $matlabCommand = "${matlabScript} Plot_Cumulative_Sum('${dataset}', '${state}'); exit;"
        
        # Run MATLAB
        & matlab -r $matlabCommand
        
        Write-Host "Completed: ${dataset} ${state}" -ForegroundColor Cyan
        Write-Host "----------------------------------------"
    }
}

Write-Host "All plots generated successfully!" -ForegroundColor Yellow
