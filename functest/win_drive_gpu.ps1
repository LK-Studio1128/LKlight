# win_drive_gpu.ps1 — GPU-only GSO driver (diagnostic focus: capture cuda_batch_score stage errors).
# Usage: powershell -ExecutionPolicy Bypass -File win_drive_gpu.ps1 -Glow 1000 -Steps 1000 -Method dna -Out c:\lkwork\ht_diag
# Reuses pre-staged lightdock_ PDBs + initial_positions_0.dat from -PDBdir.

param(
    [int]$Glow = 1000,
    [int]$Steps = 1000,
    [string]$Method = "dna",
    [string]$Out = "c:\lkwork\ht_diag",
    [string]$PDBdir = "c:\lkwork"
)

$ErrorActionPreference = "Continue"
$bin = Join-Path $PDBdir "v120\target\release\LKlight.exe"
if (-not (Test-Path $bin)) { Write-Error "missing $bin"; exit 3 }

Write-Host "=== GPU-only GSO driver (diagnostic) ==="
Write-Host "Glow=$Glow Steps=$Steps Method=$Method Out=$Out"

# -- 1) Work tree --
Remove-Item -Recurse -Force $Out -ErrorAction SilentlyContinue
New-Item -ItemType Directory -Path $Out | Out-Null

# -- 2) setup.json --
$setup = @{
  seed = 324324; anm_seed = 324324; ftdock_file = $null; noh = $false
  anm_rec = 0; anm_lig = 0; swarms = 1
  starting_points_seed = 324324; verbose_parser = $false; noxt = $false; now = $false
  restraints = $null; use_anm = $false; glowworms = $Glow; membrane = $false
  receptor_pdb = "1azp_receptor.pdb"; ligand_pdb = "1azp_ligand.pdb"
  receptor_restraints = $null; ligand_restraints = $null
} | ConvertTo-Json -Compress
Set-Content -Path (Join-Path $Out "setup.json") -Value $setup

# -- 3) initial positions: pre-staged reuse (5090 has no python) --
$initPath = Join-Path $Out "initial_positions_0.dat"
$stagedInit = Join-Path $PDBdir "initial_positions_0.dat"
if (Test-Path $stagedInit) {
    Copy-Item $stagedInit $initPath -Force
    Write-Host "reused pre-staged $stagedInit"
} else {
    Write-Error "no pre-staged initial_positions_0.dat at $PDBdir (generate on Mac, upload first)"; exit 4
}

# -- 4) lightdock_-prefixed PDBs into work dir --
Copy-Item (Join-Path $PDBdir "lightdock_1azp_receptor.pdb") $Out -Force
Copy-Item (Join-Path $PDBdir "lightdock_1azp_ligand.pdb")   $Out -Force

function Best-Energy($dir) {
    $best = $null
    Get-ChildItem -Path (Join-Path $dir "swarm_0\gso_*.out") -ErrorAction SilentlyContinue | ForEach-Object {
        foreach ($ln in (Get-Content $_ -ErrorAction SilentlyContinue)) {
            if ($ln -match "^\s*#" -or $ln -notmatch "\S") { continue }
            $tok = ($ln -split "\s+")[-1]
            try { $v = [double]$tok } catch { continue }
            if ($null -eq $best -or $v -lt $best) { $best = $v }
        }
    }
    return $best
}

Write-Host ">>> gpu  ($Method, ${Glow}glow x ${Steps}steps) ..."
$sw = [System.Diagnostics.Stopwatch]::StartNew()
$p  = Start-Process -FilePath $bin -ArgumentList @("run",(Join-Path $Out "setup.json"),(Join-Path $Out "initial_positions_0.dat"),"$Steps",$Method) `
       -WorkingDirectory $Out -RedirectStandardOutput (Join-Path $Out "stdout.log") -RedirectStandardError (Join-Path $Out "stderr.log") -NoNewWindow -Wait -PassThru
$sw.Stop()
$best = Best-Energy $Out
Write-Host ("    rc={0}  wall={1}s  best={2}" -f $p.ExitCode, [math]::Round($sw.Elapsed.TotalSeconds,2), $best)
Write-Host "GPU_ONLY_DONE dir=$Out method=$Method glow=$Glow steps=$Steps rc=$($p.ExitCode)"
