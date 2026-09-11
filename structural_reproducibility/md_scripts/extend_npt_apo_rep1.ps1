$ErrorActionPreference='Continue'
$root=Split-Path -Parent (Split-Path -Parent $MyInvocation.MyCommand.Path);$gmx='D:\APP\gromacs\bin\gmx.exe';$rep=Join-Path $root 'systems\apo\rep1';$build=Join-Path $root 'systems\apo\build';$mdp=Join-Path $root 'mdp'
if(Test-Path(Join-Path $rep 'NPT_EXTENDED.txt')){throw 'Apo rep1 NPT extension already completed'}
function Run([string]$stage,[scriptblock]$cmd){$o=(& $cmd 2>&1|Out-String);Add-Content -LiteralPath (Join-Path $rep 'npt_extension.log') -Value "`n===== $stage =====`n$o" -Encoding UTF8;if($LASTEXITCODE-ne 0){throw "$stage failed"}}
Run 'Extended NPT grompp' {& $gmx grompp -f (Join-Path $mdp 'npt_apo_rep1_extend4ns.mdp') -c (Join-Path $rep 'npt.gro') -r (Join-Path $rep 'npt.gro') -t (Join-Path $rep 'npt.cpt') -p (Join-Path $build 'topol.top') -n (Join-Path $build 'index.ndx') -o (Join-Path $rep 'npt_ext4ns.tpr') -po (Join-Path $rep 'npt_ext4ns_out.mdp') -maxwarn 0}
Run 'Extended NPT mdrun' {& $gmx mdrun -deffnm (Join-Path $rep 'npt_ext4ns') -ntmpi 1 -ntomp 22 -nb gpu -pme gpu -bonded gpu -update cpu -pin on -cpt 2}
"NPT_EXTENDED $(Get-Date -Format o) start_ps=1000 additional_ps=4000 target_total_ps=5000"|Set-Content -LiteralPath (Join-Path $rep 'NPT_EXTENDED.txt') -Encoding UTF8
