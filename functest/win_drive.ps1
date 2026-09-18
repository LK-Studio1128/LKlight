# win_drive.ps1 — three-engine high-throughput GSO driver for LKlight family on Windows.
# Usage: powershell -ExecutionPolicy Bypass -File win_drive.ps1 -Glow 1000 -Steps 1000 -Method dna -Out c:\lkwork\ht_work
# Requires: 5090 (or compatible) with c:\lkwork\{exact120,grid120,v120}\target\release\LKlight.exe built.
# Requires: Python on PATH for initial_positions generation (install Python 3.9+ or pre-generate).

param(
    [int]$Glow = 1000,
    [int]$Steps = 1000,
    [string]$Method = "dna",
    [string]$Out = "c:\lkwork\ht_work",
    [string]$PDBdir = "c:\lkwork"
)

$ErrorActionPreference = "Continue"
$engines = @{
  "exact"    = Join-Path $PDBdir "exact120\target\release\LKlight.exe"
  "cpu_grid" = Join-Path $PDBdir "grid120\target\release\LKlight.exe"
  "gpu"      = Join-Path $PDBdir "v120\target\release\LKlight.exe"
}

Write-Host "=== High-throughput three-engine GSO driver (5090) ==="
Write-Host "Glow=$Glow Steps=$Steps Method=$Method Out=$Out"

# -- 1) Stage lightdock_-prefixed PDBs (use grid engine setup to produce them, then copy) --
if (-not (Test-Path "$PDBdir\lightdock_1azp_receptor.pdb")) {
    $tmpSetup = "$PDBdir\setup_setuponly.json"
    @{seed=324324; anm_seed=324324; swarms=1; glowworms=2; starting_points_seed=324324;
      use_anm=$false; noh=$false; membrane=$false; noxt=$false; now=$false;
      receptor_pdb="1azp_receptor.pdb"; ligand_pdb="1azp_ligand.pdb"
      receptor_restraints=$null; ligand_restraints=$null; restraints=$null
      ftdock_file=$null; verbose_parser=$false; anm_rec=0; anm_lig=0} | ConvertTo-Json -Compress |
        Set-Content -Path $tmpSetup
    & $engines["cpu_grid"] setup $tmpSetup 2>&1 | Out-Null
    if (Test-Path "$PDBdir\lightdock_1azp_receptor.pdb") {
        Move-Item -Force "$PDBdir\lightdock_1azp_receptor.pdb" "$PDBdir\"
        Move-Item -Force "$PDBdir\lightdock_1azp_ligand.pdb"   "$PDBdir\"
    } else {
        Write-Error "Failed to produce lightdock_-prefixed PDBs"; exit 1
    }
    Remove-Item $tmpSetup -ErrorAction SilentlyContinue
}

# -- 2) Make work tree --
Remove-Item -Recurse -Force $Out -ErrorAction SilentlyContinue
New-Item -ItemType Directory -Path $Out | Out-Null

# -- 3) Write setup.json with the requested Glow count --
$setup = @{
  seed = 324324; anm_seed = 324324; ftdock_file = $null; noh = $false
  anm_rec = 0; anm_lig = 0; swarms = 1
  starting_points_seed = 324324; verbose_parser = $false; noxt = $false; now = $false
  restraints = $null; use_anm = $false; glowworms = $Glow; membrane = $false
  receptor_pdb = "1azp_receptor.pdb"; ligand_pdb = "1azp_ligand.pdb"
  receptor_restraints = $null; ligand_restraints = $null
} | ConvertTo-Json -Compress
Set-Content -Path (Join-Path $Out "setup.json") -Value $setup

# -- 4) initial_positions_0.dat: reuse pre-staged file if present, else generate via Python --
$initPath = Join-Path $Out "initial_positions_0.dat"
$stagedInit = Join-Path $PDBdir "initial_positions_0.dat"
if (Test-Path $stagedInit) {
    Copy-Item $stagedInit $initPath -Force
    Write-Host "reused pre-staged $stagedInit"
} else {
$py = @"
import random, math
rng = random.Random(20260904)
with open(r"$initPath","w") as f:
    for _ in range($Glow):
        mag = rng.uniform(6,20); th = rng.uniform(0,2*math.pi)
        ph = math.acos(rng.uniform(-1,1))
        tx = mag*math.sin(ph)*math.cos(th); ty = mag*math.sin(ph)*math.sin(th); tz = mag*math.cos(ph)
        a = [rng.uniform(-1,1) for _ in range(3)]
        L = math.sqrt(sum(x*x for x in a)) or 1
        ax,ay,az = [x/L for x in a]
        deg = rng.uniform(0,120); s = math.sin(math.radians(deg)/2); c = math.cos(math.radians(deg)/2)
        f.write(f"{tx:.6f} {ty:.6f} {tz:.6f} {c:.6f} {ax*s:.6f} {ay*s:.6f} {az*s:.6f}\n")
print("wrote",r"$initPath")
"@
$pyFile = Join-Path $Out "_geninit.py"
Set-Content -Path $pyFile -Value $py -Encoding ASCII
$pythonExe = (Get-Command python.exe -ErrorAction SilentlyContinue).Source
if (-not $pythonExe) { $pythonExe = (Get-Command py.exe -ErrorAction SilentlyContinue).Source }
if (-not $pythonExe) { Write-Error "Python not on PATH"; exit 2 }
& $pythonExe $pyFile
Remove-Item $pyFile -ErrorAction SilentlyContinue
}
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

$results = @()
foreach ($name in @("exact","cpu_grid","gpu")) {
    $bin = $engines[$name]
    if (-not (Test-Path $bin)) { Write-Warning "missing $bin"; continue }
    $dir = Join-Path $Out $name
    New-Item -ItemType Directory -Path $dir | Out-Null
    Copy-Item (Join-Path $Out "setup.json")              $dir
    Copy-Item $initPath                                    $dir
    Copy-Item (Join-Path $PDBdir "lightdock_1azp_receptor.pdb") $dir
    Copy-Item (Join-Path $PDBdir "lightdock_1azp_ligand.pdb")   $dir
    Write-Host ">>> $name  ($Method, ${Glow}glow x ${Steps}steps) ..."
    $sw = [System.Diagnostics.Stopwatch]::StartNew()
    $p  = Start-Process -FilePath $bin -ArgumentList @("run",(Join-Path $dir "setup.json"),(Join-Path $dir "initial_positions_0.dat"),"$Steps",$Method) `
           -WorkingDirectory $dir -RedirectStandardOutput (Join-Path $dir "stdout.log") -RedirectStandardError (Join-Path $dir "stderr.log") -NoNewWindow -Wait -PassThru
    $sw.Stop()
    $best = Best-Energy $dir
    $results += [PSCustomObject]@{ Engine=$name; RC=$p.ExitCode; WallSec=[math]::Round($sw.Elapsed.TotalSeconds,2); Best=$best }
    Write-Host ("    rc={0}  wall={1}s  best={2}" -f $p.ExitCode, [math]::Round($sw.Elapsed.TotalSeconds,2), $best)
}

# -- 6) Write summary CSV --
$csv = Join-Path $Out "ht_summary.csv"
$results | Export-Csv -Path $csv -NoTypeInformation -Encoding UTF8
Write-Host ""
Write-Host "=== summary ==="
$results | Format-Table | Out-String | Write-Host
Write-Host "ALL_DONE dir=$Out method=$Method glow=$Glow steps=$Steps"
