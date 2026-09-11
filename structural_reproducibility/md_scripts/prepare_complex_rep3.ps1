$ErrorActionPreference = 'Continue'
try { (Get-Process -Id $PID).PriorityClass = 'AboveNormal' } catch {}

$root = Split-Path -Parent (Split-Path -Parent $MyInvocation.MyCommand.Path)
$gmx = 'D:\APP\gromacs\bin\gmx.exe'
$build = Join-Path $root 'systems\complex\build'
$rep = Join-Path $root 'systems\complex\rep3'
$mdp = Join-Path $root 'mdp'

if (Test-Path -LiteralPath $rep) { throw "Refusing to overwrite existing replicate directory: $rep" }
New-Item -ItemType Directory -Path $rep | Out-Null

function Run-Gmx([string]$stage, [scriptblock]$command) {
    $output = (& $command 2>&1 | Out-String)
    Add-Content -LiteralPath (Join-Path $rep 'preparation.log') -Value "`n===== $stage =====`n$output" -Encoding UTF8
    if ($LASTEXITCODE -ne 0) { throw "$stage failed" }
}

Run-Gmx 'NVT grompp' {
    & $gmx grompp -f (Join-Path $mdp 'nvt_complex_rep3.mdp') `
        -c (Join-Path $build 'em.gro') -r (Join-Path $build 'em.gro') `
        -p (Join-Path $build 'topol.top') -n (Join-Path $build 'index.ndx') `
        -o (Join-Path $rep 'nvt.tpr') -po (Join-Path $rep 'nvt_out.mdp') -maxwarn 0
}
Run-Gmx 'NVT mdrun' {
    & $gmx mdrun -deffnm (Join-Path $rep 'nvt') -ntmpi 1 -ntomp 22 `
        -nb gpu -pme gpu -bonded gpu -update cpu -pin on -cpt 2
}
Run-Gmx 'NPT grompp' {
    & $gmx grompp -f (Join-Path $mdp 'npt_complex_rep3.mdp') `
        -c (Join-Path $rep 'nvt.gro') -r (Join-Path $rep 'nvt.gro') `
        -t (Join-Path $rep 'nvt.cpt') -p (Join-Path $build 'topol.top') `
        -n (Join-Path $build 'index.ndx') -o (Join-Path $rep 'npt.tpr') `
        -po (Join-Path $rep 'npt_out.mdp') -maxwarn 0
}
Run-Gmx 'NPT mdrun' {
    & $gmx mdrun -deffnm (Join-Path $rep 'npt') -ntmpi 1 -ntomp 22 `
        -nb gpu -pme gpu -bonded gpu -update cpu -pin on -cpt 2
}
Run-Gmx 'Production grompp' {
    & $gmx grompp -f (Join-Path $mdp 'md_100ns_complex_rep3.mdp') `
        -c (Join-Path $rep 'npt.gro') -t (Join-Path $rep 'npt.cpt') `
        -p (Join-Path $build 'topol.top') -n (Join-Path $build 'index.ndx') `
        -o (Join-Path $rep 'md_100ns.tpr') -po (Join-Path $rep 'md_100ns_out.mdp') -maxwarn 0
}

$dump = (& $gmx dump -s (Join-Path $rep 'md_100ns.tpr') 2>&1 | Out-String)
Add-Content -LiteralPath (Join-Path $rep 'preparation.log') -Value "`n===== Production TPR audit =====`n$dump" -Encoding UTF8
if ($dump -notmatch 'nsteps\s+=\s+50000000') { throw 'Production TPR does not contain 50,000,000 steps' }
"PREPARED $(Get-Date -Format o) gen_seed=307985 nvt_ld_seed=307985 npt_ld_seed=307986 production_ld_seed=307987" |
    Set-Content -LiteralPath (Join-Path $rep 'PREPARED.txt') -Encoding UTF8
