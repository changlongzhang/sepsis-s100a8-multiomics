$ErrorActionPreference = 'Continue'
$root = Split-Path -Parent (Split-Path -Parent $MyInvocation.MyCommand.Path)
$gmx = 'D:\APP\gromacs\bin\gmx.exe'
$repDir = Join-Path $root 'systems\apo\rep3'
$outRoot = Join-Path $root 'comparative_analysis_v1\chain_local_rmsf'
$tpr = Join-Path $repDir 'md_100ns_post_extension.tpr'
$traj = Join-Path $repDir 'analysis\md_100ns_proteinfit.xtc'
$ndx = Join-Path $repDir 'analysis\chain_groups_v2.ndx'
$log = Join-Path $outRoot 'apo3_parallel.log'

foreach ($chain in @([pscustomobject]@{name='B';group='2'}, [pscustomobject]@{name='D';group='3'})) {
    $out = Join-Path $outRoot "apo_rep3_chain$($chain.name)_backbone_localfit_rmsf_20_100ns.xvg"
    $marker = "$out.done"
    if (Test-Path -LiteralPath $marker) { continue }
    Add-Content -LiteralPath $log -Value "$(Get-Date -Format o)`tSTART apo_rep3_chain$($chain.name)" -Encoding UTF8
    $text = (@($chain.group) | & $gmx rmsf -s $tpr -f $traj -n $ndx -o $out -res -fit -b 20000 -e 100000 2>&1 | Out-String)
    $code = $LASTEXITCODE
    Add-Content -LiteralPath $log -Value $text -Encoding UTF8
    if ($code -ne 0) { throw "apo rep3 chain $($chain.name) failed with exit code $code" }
    $rows = @(Get-Content -LiteralPath $out | Where-Object { $_ -match '^\s*\d' }).Count
    if ($rows -ne 93) { throw "apo rep3 chain $($chain.name) produced $rows rows" }
    "completed=$(Get-Date -Format o); rows=93; fit=selected_chain_backbone; window=20-100ns" | Set-Content -LiteralPath $marker -Encoding UTF8
    Add-Content -LiteralPath $log -Value "$(Get-Date -Format o)`tPASS apo_rep3_chain$($chain.name)" -Encoding UTF8
}
