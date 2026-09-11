$ErrorActionPreference = 'Continue'
$root = Split-Path -Parent (Split-Path -Parent $MyInvocation.MyCommand.Path)
$gmx = 'D:\APP\gromacs\bin\gmx.exe'
$build = Join-Path $root 'systems\apo\build'
$mdp = Join-Path $root 'mdp'

if (-not (Test-Path -LiteralPath (Join-Path $build 'apo_unsolvated.gro'))) {
    throw "Missing apo_unsolvated.gro"
}

function Run-Gmx([string]$stage, [scriptblock]$command) {
    $output = (& $command 2>&1 | Out-String)
    Add-Content -LiteralPath (Join-Path $build 'build.log') -Value "`n===== $stage =====`n$output" -Encoding UTF8
    if ($LASTEXITCODE -ne 0) { throw "$stage failed" }
}

if (-not (Test-Path -LiteralPath (Join-Path $build 'boxed.gro'))) {
    Run-Gmx 'Define dodecahedral box' {
        & $gmx editconf -f (Join-Path $build 'apo_unsolvated.gro') `
            -o (Join-Path $build 'boxed.gro') -bt dodecahedron -d 1.0
    }
}
if (-not (Test-Path -LiteralPath (Join-Path $build 'solvated.gro'))) {
    Run-Gmx 'Solvate' {
        & $gmx solvate -cp (Join-Path $build 'boxed.gro') -cs spc216.gro `
            -o (Join-Path $build 'solvated.gro') -p (Join-Path $build 'topol.top')
    }
}
if (-not (Test-Path -LiteralPath (Join-Path $build 'solv_ions.gro'))) {
    Run-Gmx 'Ion TPR' {
        & $gmx grompp -f (Join-Path $mdp 'ions.mdp') -c (Join-Path $build 'solvated.gro') `
            -p (Join-Path $build 'topol.top') -o (Join-Path $build 'ions.tpr') -maxwarn 1
    }
    Run-Gmx 'Add 0.15 M NaCl and neutralize' {
        'SOL' | & $gmx genion -s (Join-Path $build 'ions.tpr') `
            -o (Join-Path $build 'solv_ions.gro') -p (Join-Path $build 'topol.top') `
            -pname NA -nname CL -neutral -conc 0.15
    }
}
Run-Gmx 'Energy minimization TPR' {
    & $gmx grompp -f (Join-Path $mdp 'em.mdp') -c (Join-Path $build 'solv_ions.gro') `
        -p (Join-Path $build 'topol.top') -o (Join-Path $build 'em.tpr') -maxwarn 0
}
Run-Gmx 'Energy minimization' {
    & $gmx mdrun -deffnm (Join-Path $build 'em') -ntmpi 1 -ntomp 22 `
        -nb gpu -pme cpu -bonded cpu -update cpu -pin on
}
Run-Gmx 'Default index' {
    'q' | & $gmx make_ndx -f (Join-Path $build 'em.gro') -o (Join-Path $build 'index.ndx')
}

"BUILT $(Get-Date -Format o) source=complex/build/protein.gro calcium_count=4 ligand_count=0 solvent=TIP3P salt=0.15M" |
    Set-Content -LiteralPath (Join-Path $build 'BUILT.txt') -Encoding UTF8
