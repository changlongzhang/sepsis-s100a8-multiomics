$ErrorActionPreference = 'Stop'
$root = Split-Path -Parent (Split-Path -Parent $MyInvocation.MyCommand.Path)
$gmx = 'D:\APP\gromacs\bin\gmx.exe'
$build = Join-Path $root 'systems\complex\build'
$out = Join-Path $root 'systems\complex\rep1'
$mdp = Join-Path $root 'mdp\md_100ns.mdp'
New-Item -ItemType Directory -Path $out -Force | Out-Null

# Continue the validated 1 ns unrestrained benchmark to a cumulative 100 ns.
& $gmx convert-tpr -s (Join-Path $build 'md_1ns.tpr') -extend 99000 -o (Join-Path $out 'md_100ns.tpr')
if ($LASTEXITCODE -ne 0) { throw 'convert-tpr failed' }

& $gmx mdrun -s (Join-Path $out 'md_100ns.tpr') `
  -cpi (Join-Path $build 'md_1ns.cpt') `
  -deffnm (Join-Path $out 'md_100ns') `
  -ntmpi 1 -ntomp 8 -nb gpu -pme gpu -bonded gpu -update gpu -pin on -noappend
if ($LASTEXITCODE -ne 0) { throw 'production mdrun failed' }

"COMPLETED $(Get-Date -Format o)" | Set-Content -LiteralPath (Join-Path $out 'COMPLETED.txt')
