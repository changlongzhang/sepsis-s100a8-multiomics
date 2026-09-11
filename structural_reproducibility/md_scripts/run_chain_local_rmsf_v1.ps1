$ErrorActionPreference = 'Continue'
$root = Split-Path -Parent (Split-Path -Parent $MyInvocation.MyCommand.Path)
$gmx = 'D:\APP\gromacs\bin\gmx.exe'
$outRoot = Join-Path $root 'comparative_analysis_v1\chain_local_rmsf'
$log = Join-Path $outRoot 'chain_local_rmsf.log'
New-Item -ItemType Directory -Force -Path $outRoot | Out-Null

$jobs = @(
    [pscustomobject]@{system='complex';rep='rep1';tpr='md_100ns.tpr'},
    [pscustomobject]@{system='complex';rep='rep2';tpr='md_100ns.tpr'},
    [pscustomobject]@{system='complex';rep='rep3';tpr='md_100ns.tpr'},
    [pscustomobject]@{system='apo';rep='rep1';tpr='md_100ns_post_extension.tpr'},
    [pscustomobject]@{system='apo';rep='rep2';tpr='md_100ns_post_extension.tpr'},
    [pscustomobject]@{system='apo';rep='rep3';tpr='md_100ns_post_extension.tpr'}
)

foreach ($job in $jobs) {
    $repDir = Join-Path $root "systems\$($job.system)\$($job.rep)"
    $traj = Join-Path $repDir 'analysis\md_100ns_proteinfit.xtc'
    $tpr = Join-Path $repDir $job.tpr
    $ndx = Join-Path $repDir 'analysis\chain_groups_v2.ndx'
    foreach ($chain in @([pscustomobject]@{name='B';group='2'}, [pscustomobject]@{name='D';group='3'})) {
        $tag = "$($job.system)_$($job.rep)_chain$($chain.name)"
        $out = Join-Path $outRoot "${tag}_backbone_localfit_rmsf_20_100ns.xvg"
        $marker = "$out.done"
        if (Test-Path -LiteralPath $marker) { continue }
        Add-Content -LiteralPath $log -Value "$(Get-Date -Format o)`tSTART $tag" -Encoding UTF8
        $text = (@($chain.group) | & $gmx rmsf -s $tpr -f $traj -n $ndx -o $out -res -fit -b 20000 -e 100000 2>&1 | Out-String)
        $code = $LASTEXITCODE
        Add-Content -LiteralPath $log -Value $text -Encoding UTF8
        if ($code -ne 0) { throw "$tag failed with exit code $code" }
        $rows = @(Get-Content -LiteralPath $out | Where-Object { $_ -match '^\s*\d' }).Count
        if ($rows -ne 93) { throw "$tag produced $rows residue rows; expected 93" }
        "completed=$(Get-Date -Format o); rows=93; fit=selected_chain_backbone; window=20-100ns" | Set-Content -LiteralPath $marker -Encoding UTF8
        Add-Content -LiteralPath $log -Value "$(Get-Date -Format o)`tPASS $tag" -Encoding UTF8
    }
}

"completed=$(Get-Date -Format o); outputs=12; residues_per_chain=93; window=20-100ns" | Set-Content -LiteralPath (Join-Path $outRoot 'CHAIN_LOCAL_RMSF_COMPLETE.txt') -Encoding UTF8
