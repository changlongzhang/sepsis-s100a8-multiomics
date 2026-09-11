$ErrorActionPreference = 'Continue'
$root = Split-Path -Parent (Split-Path -Parent $MyInvocation.MyCommand.Path)
$gmx = 'D:\APP\gromacs\bin\gmx.exe'
$source = Join-Path $root 'comparative_analysis_v1\pooled_ca_pca_v1'
$out = Join-Path $root 'comparative_analysis_v1\pooled_ca_pca_v2'
$log = Join-Path $out 'pca_v2_supervisor.log'
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

$tags = @('complex_rep1','complex_rep2','complex_rep3','apo_rep1','apo_rep2','apo_rep3')
$fits = @($tags | ForEach-Object { Join-Path $source "${_}_ca_commonfit.xtc" })
$ref = Join-Path $source 'common_ca_reference_20ns.gro'

try {
    foreach ($p in @($ref) + $fits) { if (-not (Test-Path -LiteralPath $p)) { throw "Missing validated v1 input $p" } }
    $pooled = Join-Path $out 'all_six_replicates_ca_commonfit_continuous_time.xtc'
    Run-Step 'concatenate common-fit trajectories with continuous 100 ps time' (Join-Path $out 'POOLED_CONTINUOUS.done') {
        @('0','80100','160200','240300','320400','400500') | & $gmx trjcat -f $fits -o $pooled -settime
    } @($pooled)

    $eval = Join-Path $out 'pooled_ca_eigenvalues.xvg'
    $evec = Join-Path $out 'pooled_ca_eigenvectors.trr'
    $avg = Join-Path $out 'pooled_ca_average.pdb'
    Run-Step 'continuous-time pooled C-alpha covariance' (Join-Path $out 'COVAR.done') {
        @('0','0') | & $gmx covar -s $ref -f $pooled -o $eval -v $evec -av $avg -l (Join-Path $out 'pooled_ca_covar.log')
    } @($eval,$evec,$avg)

    foreach ($tag in $tags) {
        $proj = Join-Path $out "${tag}_pc1_pc2_projection.xvg"
        Run-Step "$tag projection on continuous-time pooled PC1-PC2 basis" (Join-Path $out "${tag}_PROJECTION.done") {
            @('0','0') | & $gmx anaeig -s $ref -f (Join-Path $source "${tag}_ca_commonfit.xtc") -v $evec -first 1 -last 2 -proj $proj
        } @($proj)
    }

    "completed=$(Get-Date -Format o); frames=4806; timestep_ps=100; atoms=186_C-alpha; basis=pooled; time_discontinuities=false" | Set-Content -LiteralPath (Join-Path $out 'POOLED_CA_PCA_V2_COMPLETE.txt') -Encoding UTF8
    "ALL_POOLED_CA_PCA_V2_COMPLETE`t$(Get-Date -Format o)" | Set-Content -LiteralPath $stage -Encoding UTF8
    Log 'ALL_POOLED_CA_PCA_V2_COMPLETE'
} catch {
    "FAILED`t$(Get-Date -Format o)`t$($_.Exception.Message)" | Set-Content -LiteralPath $stage -Encoding UTF8
    Log "FAILED $($_.Exception.Message)"
    throw
}
