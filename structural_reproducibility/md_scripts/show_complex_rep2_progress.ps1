$ErrorActionPreference = 'SilentlyContinue'
$host.UI.RawUI.WindowTitle = 'S100A8-triptolide complex rep2 progress - display only'
$root = Split-Path -Parent (Split-Path -Parent $MyInvocation.MyCommand.Path)
$rep = Join-Path $root 'systems\complex\rep2'
$lastNs = $null
$lastPoll = $null
$rate = $null

while ($true) {
    $now = Get-Date
    $prodLog = Get-ChildItem -LiteralPath $rep -Filter 'md_100ns*.log' -File | Sort-Object LastWriteTime | Select-Object -Last 1
    if ($prodLog) {
        $stage = 'complex rep2 production'
        $log = $prodLog
        $targetNs = 100.0
    } elseif (Test-Path -LiteralPath (Join-Path $rep 'npt.log')) {
        $stage = 'complex rep2 NPT equilibration'
        $log = Get-Item -LiteralPath (Join-Path $rep 'npt.log')
        $targetNs = $null
    } else {
        $stage = 'complex rep2 NVT equilibration'
        $log = Get-Item -LiteralPath (Join-Path $rep 'nvt.log')
        $targetNs = $null
    }

    $currentNs = 0.0
    if ($log) {
        $text = (Get-Content -LiteralPath $log.FullName -Tail 350) -join "`n"
        $m = [regex]::Matches($text, 'Step\s+Time\s+\d+\s+([0-9.]+)')
        if ($m.Count) { $currentNs = [double]$m[$m.Count-1].Groups[1].Value / 1000.0 }
    }

    if ($targetNs -and $lastNs -ne $null -and $currentNs -gt $lastNs) {
        $dh = ($now-$lastPoll).TotalHours
        if ($dh -gt 0) {
            $instant = ($currentNs-$lastNs)/$dh*24
            if ($rate -eq $null) {$rate=$instant} else {$rate=0.25*$instant+0.75*$rate}
        }
    }
    $lastNs=$currentNs; $lastPoll=$now
    $cpu=(Get-CimInstance Win32_Processor|Measure-Object LoadPercentage -Average).Average
    $os=Get-CimInstance Win32_OperatingSystem
    $free=[math]::Round($os.FreePhysicalMemory/1MB,2)
    $gpu=(& nvidia-smi --query-gpu=utilization.gpu,memory.used,memory.total,temperature.gpu,power.draw --format=csv,noheader,nounits|Select-Object -First 1)-split ','
    $gp=Get-Process -Name gmx -ErrorAction SilentlyContinue|Select-Object -First 1

    Clear-Host
    Write-Host 'S100A8-triptolide molecular dynamics' -ForegroundColor Cyan
    Write-Host ('Updated: {0}' -f $now.ToString('yyyy-MM-dd HH:mm:ss'))
    Write-Host ('Stage: {0}' -f $stage)
    if ($targetNs) {
        $pct=[math]::Min(100,100*$currentNs/$targetNs)
        $bar=('#'*[math]::Floor($pct/2.5)).PadRight(40,'-')
        Write-Host ('[{0}] {1:N2}%  {2:N3}/{3:N0} ns' -f $bar,$pct,$currentNs,$targetNs)
        if($rate -and $rate -gt 0){$eta=($targetNs-$currentNs)/$rate*24;Write-Host ('Speed: {0:N1} ns/day  ETA: {1:N1} h' -f $rate,$eta)}
    } else { Write-Host ('Current trajectory time: {0:N3} ns' -f $currentNs) }
    Write-Host ('CPU: {0}%  free RAM: {1} GB' -f $cpu,$free)
    Write-Host ('GPU: {0}%  VRAM: {1}/{2} MB  temp: {3} C  power: {4} W' -f $gpu[0].Trim(),$gpu[1].Trim(),$gpu[2].Trim(),$gpu[3].Trim(),$gpu[4].Trim())
    if($gp){Write-Host ('GROMACS RUNNING  PID {0}  priority {1}' -f $gp.Id,$gp.PriorityClass) -ForegroundColor Green}else{Write-Host 'GROMACS between stages / stopped' -ForegroundColor Yellow}
    Write-Host 'Closing this display does not stop the simulation.' -ForegroundColor DarkGray
    Start-Sleep -Seconds 5
}
