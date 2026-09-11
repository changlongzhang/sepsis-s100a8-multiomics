$ErrorActionPreference = 'Continue'
try { (Get-Process -Id $PID).PriorityClass = 'AboveNormal' } catch {}
$root=Split-Path -Parent (Split-Path -Parent $MyInvocation.MyCommand.Path)
$gmx='D:\APP\gromacs\bin\gmx.exe';$rep=Join-Path $root 'systems\apo\rep2'
$tpr=Join-Path $rep 'md_100ns_post_extension.tpr';$cpt=Join-Path $rep 'md_100ns.cpt';$deffnm=Join-Path $rep 'md_100ns'
$audit=Join-Path $rep 'guarded_segments.tsv';$log=Join-Path $rep 'guarded_supervisor.log';$targetPs=100000.0
function Log([string]$m){Add-Content -LiteralPath $log -Value "$(Get-Date -Format o)`t$m" -Encoding UTF8}
function Checkpoint-Time{if(-not(Test-Path $cpt)){return 0.0};$o=(& $gmx check -f $cpt 2>&1|Out-String);if($LASTEXITCODE-ne 0){throw "Checkpoint validation failed: $o"};$m=[regex]::Match($o,'Last frame\s+-?\d+\s+time\s+([0-9.]+)');if(-not $m.Success){throw "Cannot parse checkpoint time: $o"};return [double]$m.Groups[1].Value}
if(-not(Test-Path (Join-Path $rep 'PRODUCTION_RELEASED.txt'))){throw 'Production is on QC hold; missing PRODUCTION_RELEASED.txt'}
if(-not(Test-Path $tpr)){throw "Missing production TPR: $tpr"}
if(-not(Test-Path $audit)){"segment`tstart_iso`tend_iso`tbefore_ps`tafter_ps`texit_code`tcpt_sha256`tstatus"|Set-Content -LiteralPath $audit -Encoding UTF8}
$segment=0
while($true){
    $before=Checkpoint-Time
    if($before-ge($targetPs-0.1)){"COMPLETED $(Get-Date -Format o) checkpoint_ps=$before"|Set-Content -LiteralPath (Join-Path $rep 'COMPLETED.txt') -Encoding UTF8;Log "Target reached at $before ps";break}
    $segment++;$start=Get-Date -Format o;Log "Starting segment $segment from $before ps"
    $args=@('-s',$tpr,'-deffnm',$deffnm,'-ntmpi','1','-ntomp','22','-nb','gpu','-pme','gpu','-bonded','gpu','-update','cpu','-pin','on','-noappend','-cpt','2','-maxh','0.50')
    if(Test-Path $cpt){$args+=@('-cpi',$cpt)}
    $o=(& $gmx mdrun @args 2>&1|Out-String);$exit=$LASTEXITCODE;Add-Content -LiteralPath $log -Value $o -Encoding UTF8;$end=Get-Date -Format o
    if($exit-ne 0){"$segment`t$start`t$end`t$before`tNA`t$exit`t`tFAILED"|Add-Content -LiteralPath $audit -Encoding UTF8;throw "Production segment failed with exit code $exit"}
    $after=Checkpoint-Time;$hash=(Get-FileHash -LiteralPath $cpt -Algorithm SHA256).Hash
    if($after-le $before){throw "Checkpoint did not advance: $before -> $after"}
    "$segment`t$start`t$end`t$before`t$after`t$exit`t$hash`tOK"|Add-Content -LiteralPath $audit -Encoding UTF8
    Log "Segment $segment completed: $before -> $after ps; checkpoint_sha256=$hash"
}
