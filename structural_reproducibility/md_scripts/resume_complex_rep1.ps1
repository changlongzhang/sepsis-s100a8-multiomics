$ErrorActionPreference = 'Stop'
$root = Split-Path -Parent (Split-Path -Parent $MyInvocation.MyCommand.Path)
$gmx = 'D:\APP\gromacs\bin\gmx.exe'
$rep = Join-Path $root 'systems\complex\rep1'

& $gmx mdrun -s (Join-Path $rep 'md_100ns.tpr') `
  -cpi (Join-Path $rep 'md_100ns.cpt') `
  -deffnm (Join-Path $rep 'md_100ns') `
  -ntmpi 1 -ntomp 8 -nb gpu -pme gpu -bonded gpu -update cpu -pin on -noappend
if ($LASTEXITCODE -ne 0) { throw 'checkpoint resume failed' }

"COMPLETED $(Get-Date -Format o)" | Set-Content -LiteralPath (Join-Path $rep 'COMPLETED.txt')
