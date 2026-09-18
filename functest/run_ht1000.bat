@echo off
echo === HT1000 start %date% %time% === > c:\lkwork\ht_1000_run.log 2>&1
powershell -NoProfile -ExecutionPolicy Bypass -File c:\lkwork\win_drive.ps1 -Glow 1000 -Steps 1000 -Method dna -Out c:\lkwork\ht_1000 >> c:\lkwork\ht_1000_run.log 2>&1
echo === HT1000 done rc=%errorlevel% %date% %time% === >> c:\lkwork\ht_1000_run.log 2>&1
