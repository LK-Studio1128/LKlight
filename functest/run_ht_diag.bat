@echo off
echo === start %date% %time% === > c:\lkwork\ht_diag_run.log 2>&1
powershell -NoProfile -ExecutionPolicy Bypass -File c:\lkwork\win_drive_gpu.ps1 -Glow 1000 -Steps 1000 -Method dna -Out c:\lkwork\ht_diag >> c:\lkwork\ht_diag_run.log 2>&1
echo === diag done rc=%errorlevel% %date% %time% === >> c:\lkwork\ht_diag_run.log 2>&1
