$ErrorActionPreference = 'Continue'
$root = Split-Path -Parent (Split-Path -Parent $MyInvocation.MyCommand.Path)
$gmx = 'D:\APP\gromacs\bin\gmx.exe'
$out = Join-Path $root 'comparative_analysis_v1\pooled_ca_pca_v1'
$log = Join-Path $out 'pca_supervisor.log'
$stage = Join-Path $out 'CURRENT_STAGE.txt'
New-Item -ItemType Directory -Force -Path $out | Out-Null

function Log([string]$m) { Add-Content -LiteralPath $log -Value "$(Get-Date -Format o)`t$m" -Encoding UTF8 }
function Run-Step([string]$label, [string]$marker, [scriptblock]$command, [string[]]$required) {
    if (Test-Path -LiteralPath $marker) { Log "SKIP $label"; return }
    "$label`t$(Get-Date -Format o)" | Set-Content -LiteralPath $stage -Encoding UTF8
    Log "START $label"
    $text = (& $command 2>&1 | Out-String)
    $code = $LASTEXITCODE
    Add-Content -LiteralPath $log -Value $text -Encoding UTF8
    if ($code -ne 0) { throw "$label failed with exit code $code" }
    foreach ($p in $required) {
        if (-not (Test-Path -LiteralPath $p) -or (Get-Item -LiteralPath $p).Length -eq 0) { throw "$label missing output $p" }
    }
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
    $refRep = Join-Path $root 'systems\complex\rep1'
    $ref = Join-Path $out 'common_ca_reference_20ns.gro'
    Run-Step 'extract common C-alpha reference at 20 ns' (Join-Path $out 'REFERENCE.done') {
        @('3') | & $gmx trjconv -s (Join-Path $refRep 'md_100ns.tpr') -f (Join-Path $refRep 'analysis\md_100ns_proteinfit.xtc') -n (Join-Path $root 'systems\complex\build\index.ndx') -o $ref -dump 20000
    } @($ref)

    $fitTrajs = @()
    foreach ($job in $jobs) {
        $tag = "$($job.system)_$($job.rep)"
        $repDir = Join-Path $root "systems\$($job.system)\$($job.rep)"
        $index = Join-Path $root "systems\$($job.system)\build\index.ndx"
        $extract = Join-Path $out "${tag}_ca_20_100ns_100ps.xtc"
        Run-Step "$tag extract C-alpha, 20-100 ns, 100 ps stride" (Join-Path $out "${tag}_EXTRACT.done") {
            @('3') | & $gmx trjconv -s (Join-Path $repDir $job.tpr) -f (Join-Path $repDir 'analysis\md_100ns_proteinfit.xtc') -n $index -o $extract -b 20000 -e 100000 -dt 100
        } @($extract)

        $fit = Join-Path $out "${tag}_ca_commonfit.xtc"
        Run-Step "$tag fit C-alpha to common reference" (Join-Path $out "${tag}_FIT.done") {
            @('0','0') | & $gmx trjconv -s $ref -f $extract -o $fit -fit rot+trans
        } @($fit)
        $fitTrajs += $fit
    }

    $pooled = Join-Path $out 'all_six_replicates_ca_commonfit_pooled.xtc'
    Run-Step 'concatenate all six common-fit C-alpha trajectories' (Join-Path $out 'POOLED.done') {
        & $gmx trjcat -f $fitTrajs -o $pooled -cat
    } @($pooled)

    $eval = Join-Path $out 'pooled_ca_eigenvalues.xvg'
    $evec = Join-Path $out 'pooled_ca_eigenvectors.trr'
    $avg = Join-Path $out 'pooled_ca_average.pdb'
    Run-Step 'pooled C-alpha covariance and eigendecomposition' (Join-Path $out 'COVAR.done') {
        @('0','0') | & $gmx covar -s $ref -f $pooled -o $eval -v $evec -av $avg -l (Join-Path $out 'pooled_ca_covar.log')
    } @($eval,$evec,$avg)

    foreach ($job in $jobs) {
        $tag = "$($job.system)_$($job.rep)"
        $fit = Join-Path $out "${tag}_ca_commonfit.xtc"
        $proj = Join-Path $out "${tag}_pc1_pc2_projection.xvg"
        Run-Step "$tag projection on pooled PC1-PC2 basis" (Join-Path $out "${tag}_PROJECTION.done") {
            @('0','0') | & $gmx anaeig -s $ref -f $fit -v $evec -first 1 -last 2 -proj $proj
        } @($proj)
    }

    "completed=$(Get-Date -Format o); atoms=186_C-alpha; window=20-100ns; stride=100ps; replicas=6; basis=pooled" | Set-Content -LiteralPath (Join-Path $out 'POOLED_CA_PCA_COMPLETE.txt') -Encoding UTF8
    "ALL_POOLED_CA_PCA_COMPLETE`t$(Get-Date -Format o)" | Set-Content -LiteralPath $stage -Encoding UTF8
    Log 'ALL_POOLED_CA_PCA_COMPLETE'
} catch {
    "FAILED`t$(Get-Date -Format o)`t$($_.Exception.Message)" | Set-Content -LiteralPath $stage -Encoding UTF8
    Log "FAILED $($_.Exception.Message)"
    throw
}
