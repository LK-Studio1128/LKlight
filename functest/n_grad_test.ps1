$ErrorActionPreference = "Continue"
foreach ($n in 64,128,256,384,512,768,1000) {
    $d = "c:\lkwork\ht_n$n"
    Remove-Item -Recurse -Force $d -ErrorAction SilentlyContinue
    New-Item -ItemType Directory $d | Out-Null
    Get-Content c:\lkwork\initial_positions_0.dat -TotalCount $n | Set-Content "$d\initial_positions_0.dat"
    $s = @{seed=324324;anm_seed=324324;swarms=1;glowworms=$n;starting_points_seed=324324;use_anm=$false;noh=$false;membrane=$false;noxt=$false;now=$false;receptor_pdb="1azp_receptor.pdb";ligand_pdb="1azp_ligand.pdb";receptor_restraints=$null;ligand_restraints=$null;restraints=$null;ftdock_file=$null;verbose_parser=$false;anm_rec=0;anm_lig=0}
    ($s | ConvertTo-Json -Compress) | Set-Content "$d\setup.json"
    Copy-Item c:\lkwork\lightdock_1azp_receptor.pdb $d
    Copy-Item c:\lkwork\lightdock_1azp_ligand.pdb $d
    $log = "$d\err.log"
    $p = Start-Process -FilePath "c:\lkwork\v120\target\release\LKlight.exe" -ArgumentList @("run","$d\setup.json","$d\initial_positions_0.dat","10","dna") -WorkingDirectory $d -RedirectStandardError $log -RedirectStandardOutput "$d\out.log" -NoNewWindow -Wait -PassThru
    $err = if (Test-Path $log) { Get-Content $log -Raw } else { "" }
    $hasErr = $err -match "batch_score error"
    $active = $err -match "BATCH scoring ACTIVE"
    $wall = $p.ExitCode
    Write-Host ("N={0} batch_err={1} batch_active={2} rc={3}" -f $n, $hasErr, $active, $wall)
}
Write-Host "NGRAD_DONE"
