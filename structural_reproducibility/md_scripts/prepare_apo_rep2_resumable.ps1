$ErrorActionPreference = 'Continue'
try { (Get-Process -Id $PID).PriorityClass = 'AboveNormal' } catch {}
$root = Split-Path -Parent (Split-Path -Parent $MyInvocation.MyCommand.Path)
$gmx = 'D:\APP\gromacs\bin\gmx.exe'
$build = Join-Path $root 'systems\apo\build'
$rep = Join-Path $root 'systems\apo\rep2'
$mdp = Join-Path $root 'mdp'
New-Item -ItemType Directory -Path $rep -Force | Out-Null
$runLog = Join-Path $rep 'equilibration_supervisor.log'
function Log([string]$message) { Add-Content -LiteralPath $runLog -Value "$(Get-Date -Format o)`t$message" -Encoding UTF8 }
function Run-Gmx([string]$stage, [scriptblock]$command) {
    Log "START $stage"
    $output = (& $command 2>&1 | Out-String)
    Add-Content -LiteralPath $runLog -Value $output -Encoding UTF8
    if ($LASTEXITCODE -ne 0) { throw "$stage failed with exit code $LASTEXITCODE" }
    Log "PASS $stage"
}
function Assert-Checkpoint([string]$path, [double]$targetPs) {
    $output = (& $gmx check -f $path 2>&1 | Out-String)
    Add-Content -LiteralPath $runLog -Value $output -Encoding UTF8
    if ($LASTEXITCODE -ne 0 -or $output -notmatch "Last frame\s+-?\d+\s+time\s+$([regex]::Escape($targetPs.ToString('0.000')))" ) {
        throw "Checkpoint validation failed for $path at expected time $targetPs ps"
    }
}
try {
    if (-not (Test-Path (Join-Path $rep 'NVT_COMPLETE.txt'))) {
        Run-Gmx 'NVT grompp' { & $gmx grompp -f (Join-Path $mdp 'nvt_apo_rep2.mdp') -c (Join-Path $build 'em.gro') -r (Join-Path $build 'em.gro') -p (Join-Path $build 'topol.top') -n (Join-Path $build 'index.ndx') -o (Join-Path $rep 'nvt.tpr') -po (Join-Path $rep 'nvt_out.mdp') -maxwarn 0 }
        Run-Gmx 'NVT mdrun' { & $gmx mdrun -deffnm (Join-Path $rep 'nvt') -ntmpi 1 -ntomp 22 -nb gpu -pme gpu -bonded gpu -update cpu -pin on -cpt 2 }
        Assert-Checkpoint (Join-Path $rep 'nvt.cpt') 500.000
        "NVT_COMPLETE $(Get-Date -Format o) checkpoint_ps=500.000 gen_seed=507985 ld_seed=507985" | Set-Content -LiteralPath (Join-Path $rep 'NVT_COMPLETE.txt') -Encoding UTF8
    }
    if (-not (Test-Path (Join-Path $rep 'NPT1_COMPLETE.txt'))) {
        Run-Gmx 'NPT1 grompp' { & $gmx grompp -f (Join-Path $mdp 'npt_apo_rep2.mdp') -c (Join-Path $rep 'nvt.gro') -r (Join-Path $rep 'nvt.gro') -t (Join-Path $rep 'nvt.cpt') -p (Join-Path $build 'topol.top') -n (Join-Path $build 'index.ndx') -o (Join-Path $rep 'npt.tpr') -po (Join-Path $rep 'npt_out.mdp') -maxwarn 0 }
        Run-Gmx 'NPT1 mdrun' { & $gmx mdrun -deffnm (Join-Path $rep 'npt') -ntmpi 1 -ntomp 22 -nb gpu -pme gpu -bonded gpu -update cpu -pin on -cpt 2 }
        Assert-Checkpoint (Join-Path $rep 'npt.cpt') 1000.000
        "NPT1_COMPLETE $(Get-Date -Format o) checkpoint_ps=1000.000 ld_seed=507986" | Set-Content -LiteralPath (Join-Path $rep 'NPT1_COMPLETE.txt') -Encoding UTF8
    }
    if (-not (Test-Path (Join-Path $rep 'NPT_EXTENDED.txt'))) {
        Run-Gmx 'NPT extension grompp' { & $gmx grompp -f (Join-Path $mdp 'npt_apo_rep2_extend4ns.mdp') -c (Join-Path $rep 'npt.gro') -r (Join-Path $rep 'npt.gro') -t (Join-Path $rep 'npt.cpt') -p (Join-Path $build 'topol.top') -n (Join-Path $build 'index.ndx') -o (Join-Path $rep 'npt_ext4ns.tpr') -po (Join-Path $rep 'npt_ext4ns_out.mdp') -maxwarn 0 }
        Run-Gmx 'NPT extension mdrun' { & $gmx mdrun -deffnm (Join-Path $rep 'npt_ext4ns') -ntmpi 1 -ntomp 22 -nb gpu -pme gpu -bonded gpu -update cpu -pin on -cpt 2 }
        Assert-Checkpoint (Join-Path $rep 'npt_ext4ns.cpt') 4000.000
        "NPT_EXTENDED $(Get-Date -Format o) checkpoint_ps=4000.000 cumulative_npt_ps=5000 ld_seed=507986" | Set-Content -LiteralPath (Join-Path $rep 'NPT_EXTENDED.txt') -Encoding UTF8
    }
    "EQUILIBRATION_COMPLETE_PENDING_QC $(Get-Date -Format o); production not generated or released" | Set-Content -LiteralPath (Join-Path $rep 'EQUILIBRATION_COMPLETE_PENDING_QC.txt') -Encoding UTF8
    Log 'Equilibration completed; stopped at mandatory QC hold before production'
} catch {
    Log "FAILED $($_.Exception.Message)"
    throw
}
