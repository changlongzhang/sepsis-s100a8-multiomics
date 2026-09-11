$ErrorActionPreference = 'SilentlyContinue'
$host.UI.RawUI.WindowTitle = 'S100A8 apo rep2 progress - display only'
$root = Split-Path -Parent (Split-Path -Parent $MyInvocation.MyCommand.Path)
$rep = Join-Path $root 'systems\apo\rep2'
while ($true) {
    $now = Get-Date
    $prod=Get-ChildItem -LiteralPath $rep -Filter 'md_100ns*.log' -File|Sort-Object LastWriteTime|Select-Object -Last 1
    if ($prod) { $stage='100 ns production'; $target=100.0; $log=$prod.FullName }
    elseif (Test-Path (Join-Path $rep 'EQUILIBRATION_COMPLETE_PENDING_QC.txt')) { $stage='QC hold'; $ns=5.0; $target=5.0 }
    elseif (Test-Path (Join-Path $rep 'npt_ext4ns.log')) { $stage='extended NPT (cumulative 1-5 ns)'; $target=4.0; $log=Join-Path $rep 'npt_ext4ns.log' }
    elseif (Test-Path (Join-Path $rep 'npt.log')) { $stage='initial NPT'; $target=1.0; $log=Join-Path $rep 'npt.log' }
    else { $stage='NVT'; $target=0.5; $log=Join-Path $rep 'nvt.log' }
    if ($null -eq $ns) { $ns=0.0; if (Test-Path $log) { $txt=(Get-Content -LiteralPath $log -Tail 350)-join "`n"; $m=[regex]::Matches($txt,'Step\s+Time\s+\d+\s+([0-9.]+)'); if($m.Count){$ns=[double]$m[$m.Count-1].Groups[1].Value/1000} } }
    $pct=[math]::Min(100,100*$ns/$target);$bar=('#'*[math]::Floor($pct/2.5)).PadRight(40,'-')
    $cpu=(Get-CimInstance Win32_Processor|Measure-Object LoadPercentage -Average).Average
    $gpu=(& nvidia-smi --query-gpu=utilization.gpu,memory.used,memory.total,temperature.gpu --format=csv,noheader,nounits|Select-Object -First 1)-split ','
    Clear-Host; Write-Host 'S100A8 apo replicate 2 equilibration' -ForegroundColor Cyan
    Write-Host ('Updated: {0}'-f$now.ToString('yyyy-MM-dd HH:mm:ss')); Write-Host ('Stage: {0}'-f$stage)
    Write-Host ('[{0}] {1:N2}%  {2:N3}/{3:N1} ns'-f$bar,$pct,$ns,$target)
    Write-Host ('CPU: {0}%  GPU: {1}%  VRAM: {2}/{3} MB  temp: {4} C'-f$cpu,$gpu[0].Trim(),$gpu[1].Trim(),$gpu[2].Trim(),$gpu[3].Trim())
    if(Get-Process -Name gmx -ErrorAction SilentlyContinue){Write-Host 'GROMACS RUNNING' -ForegroundColor Green}else{Write-Host 'GROMACS between stages / QC hold / stopped' -ForegroundColor Yellow}
    Write-Host 'Closing this display does not stop the simulation.' -ForegroundColor DarkGray
    $ns=$null; Start-Sleep 5
}
