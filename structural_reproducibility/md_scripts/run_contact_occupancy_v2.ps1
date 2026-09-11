param([string]$ProjectRoot = '')

if ([string]::IsNullOrWhiteSpace($ProjectRoot)) {
    $ProjectRoot = Split-Path -Parent $PSScriptRoot
}

$ErrorActionPreference = 'Stop'
$gmx = 'D:\APP\gromacs\bin\gmx.exe'
$outRoot = Join-Path $ProjectRoot 'comparative_analysis_v1\contact_occupancy_v2'
$log = Join-Path $outRoot 'contact_occupancy_v2.log'
$stage = Join-Path $outRoot 'CURRENT_STAGE.txt'
New-Item -ItemType Directory -Force -Path $outRoot | Out-Null

function Write-Audit([string]$Message) {
    $line = "$(Get-Date -Format o)`t$Message"
    Add-Content -LiteralPath $log -Value $line
    Set-Content -LiteralPath $stage -Value $line
}

# Intersect with Protein a second time. Without this restriction, repeated
# residue numbers in solvent can be included by "same residue as".
$selection = 'res_com of (group "Protein" and same residue as (group "Protein" and within 0.4 of group "TPL"))'
$index = Join-Path $ProjectRoot 'systems\complex\build\index.ndx'

foreach ($rep in 1..3) {
    $repDir = Join-Path $ProjectRoot "systems\complex\rep$rep"
    $out = Join-Path $outRoot "complex_rep${rep}_protein_residue_contact_occupancy_0p4nm_20_100ns.xvg"
    $done = "$out.done"
    if (Test-Path -LiteralPath $done) {
        Write-Audit "SKIP complex/rep$rep corrected protein-residue occupancy (validated marker exists)"
        continue
    }

    Write-Audit "START complex/rep$rep corrected protein-residue occupancy, cutoff 0.4 nm, 20-100 ns"
    $runLog = Join-Path $outRoot "complex_rep${rep}_gmx_select.log"
    $ErrorActionPreference = 'SilentlyContinue'
    & $gmx select `
        -s (Join-Path $repDir 'md_100ns.tpr') `
        -f (Join-Path $repDir 'analysis\md_100ns_proteinfit.xtc') `
        -n $index `
        -select $selection `
        -of $out `
        -b 20000 -e 100000 -xvg none *> $runLog
    $gmxExit = $LASTEXITCODE
    $ErrorActionPreference = 'Stop'
    Add-Content -LiteralPath $log -Value (Get-Content -LiteralPath $runLog -Raw)
    if ($gmxExit -ne 0) {
        Write-Audit "FAIL complex/rep$rep corrected occupancy; gmx exit $gmxExit"
        exit $gmxExit
    }
    $rows = @(Get-Content -LiteralPath $out | Where-Object { $_ -match '^\s*\d' }).Count
    if ($rows -ne 186) {
        Write-Audit "FAIL complex/rep$rep corrected occupancy has $rows rows; expected 186 protein residues"
        exit 22
    }
    Set-Content -LiteralPath $done -Value "completed=$(Get-Date -Format o); rows=$rows; cutoff_nm=0.4; window_ps=20000-100000; protein_restricted=true"
    Write-Audit "PASS complex/rep$rep corrected protein-residue occupancy ($rows rows)"
}

Set-Content -LiteralPath (Join-Path $outRoot 'CONTACT_OCCUPANCY_V2_COMPLETE.txt') -Value "completed=$(Get-Date -Format o); replicas=3; residues=186; cutoff_nm=0.4; window_ps=20000-100000"
Write-Audit 'ALL_CORRECTED_CONTACT_OCCUPANCY_COMPLETE'
