# PowerShell script to run the first order analysis
Write-Host "Starting first order analysis..."

# Set working directory
Set-Location "d:\work\Task_Convergence"

# Create log file
$logFile = "FirstOrder_Analysis_Log.txt"
"$(Get-Date): Starting analysis" | Out-File -FilePath $logFile -Append

# Run MATLAB script
Write-Host "Running Plot_FirstOrder_Spatial.m..."
"$(Get-Date): Running Plot_FirstOrder_Spatial.m" | Out-File -FilePath $logFile -Append

try {
    # Run MATLAB with log output
    & "C:\Program Files\MATLAB\R2024a\bin\matlab.exe" -nodisplay -r "try, run('Plot_FirstOrder_Spatial.m'), catch ME, disp(['Error: ' ME.message]), end, exit;" | Out-File -FilePath $logFile -Append
    
    Write-Host "MATLAB script executed. Checking output files..."
    "$(Get-Date): MATLAB script executed" | Out-File -FilePath $logFile -Append
    
    # Check if output files were created
    if (Test-Path "FirstOrder_Spatial_Settings.mat") {
        Write-Host "✓ FirstOrder_Spatial_Settings.mat created successfully"
        "$(Get-Date): FirstOrder_Spatial_Settings.mat created successfully" | Out-File -FilePath $logFile -Append
    } else {
        Write-Host "✗ FirstOrder_Spatial_Settings.mat not found"
        "$(Get-Date): FirstOrder_Spatial_Settings.mat not found" | Out-File -FilePath $logFile -Append
    }
    
    if (Test-Path "FirstOrder_Analysis.txt") {
        Write-Host "✓ FirstOrder_Analysis.txt created successfully"
        "$(Get-Date): FirstOrder_Analysis.txt created successfully" | Out-File -FilePath $logFile -Append
        # Display the content
        Write-Host "\nContent of FirstOrder_Analysis.txt:"
        Get-Content "FirstOrder_Analysis.txt"
    } else {
        Write-Host "✗ FirstOrder_Analysis.txt not found"
        "$(Get-Date): FirstOrder_Analysis.txt not found" | Out-File -FilePath $logFile -Append
    }
    
} catch {
    Write-Host "Error running MATLAB: $($_.Exception.Message)"
    "$(Get-Date): Error running MATLAB: $($_.Exception.Message)" | Out-File -FilePath $logFile -Append
}

Write-Host "\nAnalysis completed. Check $logFile for details."
"$(Get-Date): Analysis completed" | Out-File -FilePath $logFile -Append
