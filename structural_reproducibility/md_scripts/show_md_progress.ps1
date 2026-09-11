$ErrorActionPreference = 'SilentlyContinue'
$host.UI.RawUI.WindowTitle = 'S100A8-triptolide MD progress - display only'

$root = Split-Path -Parent (Split-Path -Parent $MyInvocation.MyCommand.Path)
$rep = Join-Path $root 'systems\complex\rep1'
$targetNs = 100.0
$lastNs = $null
$lastPoll = $null
$smoothedRate = $null

while ($true) {
    $now = Get-Date
    $latestLog = Get-ChildItem -LiteralPath $rep -Filter 'md_100ns.part*.log' -File |
        Sort-Object LastWriteTime | Select-Object -Last 1
    $currentNs = 0.0
    if ($latestLog) {
        $tail = Get-Content -LiteralPath $latestLog.FullName -Tail 350
        $text = $tail -join "`n"
        $matches = [regex]::Matches($text, 'Step\s+Time\s+\d+\s+([0-9.]+)')
        if ($matches.Count -gt 0) {
            $currentNs = [double]$matches[$matches.Count - 1].Groups[1].Value / 1000.0
        }
    }

    if ($lastNs -ne $null -and $lastPoll -ne $null -and $currentNs -gt $lastNs) {
        $deltaHours = ($now - $lastPoll).TotalHours
        if ($deltaHours -gt 0) {
            $instantRate = ($currentNs - $lastNs) / $deltaHours * 24.0
            if ($smoothedRate -eq $null) { $smoothedRate = $instantRate }
            else { $smoothedRate = 0.25 * $instantRate + 0.75 * $smoothedRate }
        }
    }
    $lastNs = $currentNs
    $lastPoll = $now

    $percent = [math]::Min(100.0, 100.0 * $currentNs / $targetNs)
    $filled = [math]::Floor($percent / 2.5)
    $bar = ('#' * $filled).PadRight(40, '-')
    $eta = 'calculating'
    if ($smoothedRate -and $smoothedRate -gt 0.1) {
        $hoursLeft = ($targetNs - $currentNs) / $smoothedRate * 24.0
        if ($hoursLeft -ge 0) { $eta = ('{0:N1} h' -f $hoursLeft) }
    }

    $cpuLoad = (Get-CimInstance Win32_Processor | Measure-Object LoadPercentage -Average).Average
    $os = Get-CimInstance Win32_OperatingSystem
    $freeGB = [math]::Round($os.FreePhysicalMemory / 1MB, 2)
    $gpu = (& nvidia-smi --query-gpu=utilization.gpu,memory.used,memory.total,temperature.gpu,power.draw --format=csv,noheader,nounits | Select-Object -First 1) -split ','
    $gmx = Get-Process -Name gmx -ErrorAction SilentlyContinue | Select-Object -First 1

    Clear-Host
    Write-Host 'S100A8-triptolide molecular dynamics' -ForegroundColor Cyan
    Write-Host ('Updated: {0}' -f $now.ToString('yyyy-MM-dd HH:mm:ss'))
    Write-Host ''
    Write-Host ('complex rep1  [{0}] {1:N2}%  {2:N3}/{3:N0} ns' -f $bar,$percent,$currentNs,$targetNs)
    if ($smoothedRate) { Write-Host ('Rolling speed: {0:N1} ns/day    ETA: {1}' -f $smoothedRate,$eta) }
    else { Write-Host ('Rolling speed: calculating    ETA: {0}' -f $eta) }
    Write-Host ''
    Write-Host ('CPU total: {0}%    free RAM: {1} GB' -f $cpuLoad,$freeGB)
    Write-Host ('GPU: {0}%    VRAM: {1}/{2} MB    temp: {3} C    power: {4} W' -f $gpu[0].Trim(),$gpu[1].Trim(),$gpu[2].Trim(),$gpu[3].Trim(),$gpu[4].Trim())
    if ($gmx) { Write-Host ('GROMACS: RUNNING    PID {0}    priority {1}' -f $gmx.Id,$gmx.PriorityClass) -ForegroundColor Green }
    else { Write-Host 'GROMACS: between segments / paused / stopped' -ForegroundColor Yellow }
    if ($latestLog) { Write-Host ('Latest part: {0}' -f $latestLog.BaseName) }
    Write-Host ''
    Write-Host 'Close this window or press Ctrl+C to close the display only; the simulation continues.' -ForegroundColor DarkGray
    Start-Sleep -Seconds 5
}
