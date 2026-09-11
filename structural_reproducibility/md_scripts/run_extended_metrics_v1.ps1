$ErrorActionPreference = 'Continue'
$root = Split-Path -Parent (Split-Path -Parent $MyInvocation.MyCommand.Path)
$gmx = 'D:\APP\gromacs\bin\gmx.exe'
$outRoot = Join-Path $root 'comparative_analysis_v1\per_trajectory'
$log = Join-Path $root 'comparative_analysis_v1\extended_metrics_supervisor.log'
$progress = Join-Path $root 'comparative_analysis_v1\CURRENT_STAGE.txt'
New-Item -ItemType Directory -Path $outRoot -Force | Out-Null
Add-Type @'
using System;
using System.Runtime.InteropServices;
public static class MetricsExecutionState {
    [DllImport("kernel32.dll", SetLastError=true)]
    public static extern uint SetThreadExecutionState(uint esFlags);
}
'@
function Log([string]$m) { Add-Content -LiteralPath $log -Value "$(Get-Date -Format o)`t$m" -Encoding UTF8 }
function Run-Metric([string]$label,[string]$marker,[scriptblock]$command,[string[]]$requiredOutputs) {
    if (Test-Path $marker) { Log "SKIP completed $label"; return }
    "$label`t$(Get-Date -Format o)" | Set-Content -LiteralPath $progress -Encoding UTF8
    Log "START $label"
    $o = (& $command 2>&1 | Out-String)
    $exit = $LASTEXITCODE
    Add-Content -LiteralPath $log -Value $o -Encoding UTF8
    if ($exit -ne 0) { throw "$label failed with exit code $exit" }
    foreach ($path in $requiredOutputs) { if (-not (Test-Path $path) -or (Get-Item $path).Length -eq 0) { throw "$label missing output $path" } }
    "PASS $(Get-Date -Format o)" | Set-Content -LiteralPath $marker -Encoding UTF8
    Log "PASS $label"
}
$jobs = @(
    [pscustomobject]@{system='complex';rep='rep1';tpr='md_100ns.tpr'},
    [pscustomobject]@{system='complex';rep='rep2';tpr='md_100ns.tpr'},
    [pscustomobject]@{system='complex';rep='rep3';tpr='md_100ns.tpr'},
    [pscustomobject]@{system='apo';rep='rep1';tpr='md_100ns_post_extension.tpr'},
    [pscustomobject]@{system='apo';rep='rep2';tpr='md_100ns_post_extension.tpr'},
    [pscustomobject]@{system='apo';rep='rep3';tpr='md_100ns_post_extension.tpr'}
)
try {
    [void][MetricsExecutionState]::SetThreadExecutionState(2147483649)
    foreach ($job in $jobs) {
        $repDir = Join-Path $root "systems\$($job.system)\$($job.rep)"
        $analysis = Join-Path $repDir 'analysis'
        $traj = Join-Path $analysis 'md_100ns_proteinfit.xtc'
        $tpr = Join-Path $repDir $job.tpr
        $index = Join-Path $root "systems\$($job.system)\build\index.ndx"
        $chainIndex = Join-Path $analysis 'chain_groups_v2.ndx'
        $out = Join-Path $outRoot "$($job.system)_$($job.rep)"
        New-Item -ItemType Directory -Path $out -Force | Out-Null
        if (-not (Test-Path $traj) -or -not (Test-Path $tpr) -or -not (Test-Path $index) -or -not (Test-Path $chainIndex)) { throw "Missing required input for $($job.system)/$($job.rep)" }
        $rmsf = Join-Path $out 'ca_rmsf_20_100ns.xvg';$avg = Join-Path $out 'ca_average_20_100ns.pdb'
        Run-Metric "$($job.system)/$($job.rep) C-alpha RMSF" (Join-Path $out 'RMSF.done') { @('3') | & $gmx rmsf -s $tpr -f $traj -n $index -o $rmsf -oq $avg -res -b 20000 -e 100000 } @($rmsf,$avg)
        $rg = Join-Path $out 'protein_rg_20_100ns.xvg'
        Run-Metric "$($job.system)/$($job.rep) protein Rg" (Join-Path $out 'RG.done') { @('1') | & $gmx gyrate -s $tpr -f $traj -n $index -o $rg -b 20000 -e 100000 } @($rg)
        $sasa = Join-Path $out 'protein_sasa_20_100ns.xvg';$sasaRes = Join-Path $out 'protein_sasa_residue_20_100ns.xvg'
        Run-Metric "$($job.system)/$($job.rep) protein SASA" (Join-Path $out 'SASA.done') { @('1','1') | & $gmx sasa -s $tpr -f $traj -n $index -o $sasa -or $sasaRes -b 20000 -e 100000 } @($sasa,$sasaRes)
        $interHb = Join-Path $out 'interchain_hbonds_20_100ns.xvg'
        Run-Metric "$($job.system)/$($job.rep) inter-chain hydrogen bonds" (Join-Path $out 'INTERCHAIN_HB.done') { @('0','1') | & $gmx hbond -s $tpr -f $traj -n $chainIndex -num $interHb -b 20000 -e 100000 } @($interHb)
        if ($job.system -eq 'complex') {
            $ligHb = Join-Path $out 'ligand_protein_hbonds_20_100ns.xvg'
            Run-Metric "$($job.system)/$($job.rep) ligand-protein hydrogen bonds" (Join-Path $out 'LIGAND_HB.done') { @('1','14') | & $gmx hbond -s $tpr -f $traj -n $index -num $ligHb -b 20000 -e 100000 } @($ligHb)
        }
    }
    "ALL_COMMON_METRICS_COMPLETE`t$(Get-Date -Format o)" | Set-Content -LiteralPath $progress -Encoding UTF8
    "COMPLETE $(Get-Date -Format o)" | Set-Content -LiteralPath (Join-Path $root 'comparative_analysis_v1\COMMON_METRICS_COMPLETE.txt') -Encoding UTF8
    Log 'All common per-trajectory metrics completed'
} catch {
    "FAILED`t$(Get-Date -Format o)`t$($_.Exception.Message)" | Set-Content -LiteralPath $progress -Encoding UTF8
    Log "FAILED $($_.Exception.Message)"
    throw
} finally {
    [void][MetricsExecutionState]::SetThreadExecutionState(2147483648)
}
