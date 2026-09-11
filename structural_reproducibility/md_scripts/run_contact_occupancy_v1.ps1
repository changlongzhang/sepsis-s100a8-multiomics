param([string]$ProjectRoot = '')

if ([string]::IsNullOrWhiteSpace($ProjectRoot)) {
    $ProjectRoot = Split-Path -Parent $PSScriptRoot
}

$ErrorActionPreference = 'Stop'
$gmx = 'D:\APP\gromacs\bin\gmx.exe'
$outRoot = Join-Path $ProjectRoot 'comparative_analysis_v1\contact_occupancy'
$log = Join-Path $outRoot 'contact_occupancy.log'
$stage = Join-Path $outRoot 'CURRENT_STAGE.txt'
New-Item -ItemType Directory -Force -Path $outRoot | Out-Null

function Write-Audit([string]$Message) {
    $line = "$(Get-Date -Format o)`t$Message"
    Add-Content -LiteralPath $log -Value $line
    Set-Content -LiteralPath $stage -Value $line
}

$selection = 'res_com of (same residue as (group "Protein" and within 0.4 of group "TPL"))'
$index = Join-Path $ProjectRoot 'systems\complex\build\index.ndx'

foreach ($rep in 1..3) {
    $repDir = Join-Path $ProjectRoot "systems\complex\rep$rep"
    $out = Join-Path $outRoot "complex_rep${rep}_residue_contact_occupancy_0p4nm_20_100ns.xvg"
    $done = "$out.done"
    if (Test-Path -LiteralPath $done) {
        Write-Audit "SKIP complex/rep$rep residue contact occupancy (validated marker exists)"
        continue
    }

    Write-Audit "START complex/rep$rep residue contact occupancy, cutoff 0.4 nm, 20-100 ns"
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
        Write-Audit "FAIL complex/rep$rep residue contact occupancy; gmx exit $gmxExit"
        exit $gmxExit
    }
    $rows = @(Get-Content -LiteralPath $out | Where-Object { $_ -match '^\s*\d' }).Count
    if ($rows -lt 180) {
        Write-Audit "FAIL complex/rep$rep occupancy output has only $rows residue rows"
        exit 21
    }
    Set-Content -LiteralPath $done -Value "completed=$(Get-Date -Format o); rows=$rows; cutoff_nm=0.4; window_ps=20000-100000"
    Write-Audit "PASS complex/rep$rep residue contact occupancy ($rows residue rows)"
}

Set-Content -LiteralPath (Join-Path $outRoot 'CONTACT_OCCUPANCY_COMPLETE.txt') -Value "completed=$(Get-Date -Format o); replicas=3; cutoff_nm=0.4; window_ps=20000-100000"
Write-Audit 'ALL_CONTACT_OCCUPANCY_COMPLETE'
