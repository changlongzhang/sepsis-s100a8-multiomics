$ErrorActionPreference = 'SilentlyContinue'
$root = Split-Path -Parent (Split-Path -Parent $MyInvocation.MyCommand.Path)
$rep = Join-Path $root 'systems\apo\rep3'
$gmx = 'D:\APP\gromacs\bin\gmx.exe'
$runnerScript = Join-Path $root 'scripts\run_apo_rep3_guarded_highpower.ps1'
$cpt = Join-Path $rep 'md_100ns.cpt'
$monitorLog = Join-Path $rep 'persistent_monitor.log'
$targetPs = 100000.0
Add-Type @'
using System;
using System.Runtime.InteropServices;
public static class Apo3PersistentExecutionState {
    [DllImport("kernel32.dll", SetLastError=true)]
    public static extern uint SetThreadExecutionState(uint esFlags);
}
'@
function Log([string]$m) { Add-Content -LiteralPath $monitorLog -Value "$(Get-Date -Format o)`t$m" -Encoding UTF8 }
function Checkpoint-Time {
    if (-not (Test-Path $cpt)) { return 0.0 }
    $savedPreference = $ErrorActionPreference
    $ErrorActionPreference = 'Continue'
    $o = (& $gmx check -f $cpt 2>&1 | Out-String)
    $exitCode = $LASTEXITCODE
    $ErrorActionPreference = $savedPreference
    if ($exitCode -ne 0) { return -1.0 }
    $m = [regex]::Match($o, 'Last frame\s+-?\d+\s+time\s+([0-9.]+)')
    if (-not $m.Success) { return -1.0 }
    return [double]$m.Groups[1].Value
}
function Has-Blocking-Error {
    $logs = Get-ChildItem -LiteralPath $rep -Filter 'md_100ns.part*.log' -File
    foreach ($file in $logs) {
        $hits = Select-String -LiteralPath $file.FullName -Pattern '(?i)LINCS WARNING|Fatal error|CUDA error|segmentation fault|(^|[^A-Za-z])NaN([^A-Za-z]|$)'
        if ($hits) { Log "BLOCKING_ERROR file=$($file.Name) hits=$($hits.Count)"; return $true }
    }
    return $false
}
$continuous = 0x80000000
$systemRequired = 0x00000001
$lastTime = -1.0
$stalledRestarts = 0
Log 'Persistent monitor started'
powercfg /setactive 8c5e7fda-e8bf-4a96-9a85-a6e23a8c635c | Out-Null
while ($true) {
    [void][Apo3PersistentExecutionState]::SetThreadExecutionState($continuous -bor $systemRequired)
    $time = Checkpoint-Time
    if ($time -lt 0) { Log 'Checkpoint unreadable; monitor blocked'; 'CHECKPOINT_UNREADABLE' | Set-Content -LiteralPath (Join-Path $rep 'MONITOR_BLOCKED.txt') -Encoding UTF8; break }
    if ($time -ge ($targetPs - 0.1)) {
        "COMPLETED $(Get-Date -Format o) checkpoint_ps=$time" | Set-Content -LiteralPath (Join-Path $rep 'COMPLETED.txt') -Encoding UTF8
        Log "Target reached at $time ps; monitor stopping"
        powercfg /setactive 381b4222-f694-41f0-9685-ff5bb260df2e | Out-Null
        break
    }
    if (Has-Blocking-Error) { 'RUNTIME_ERROR_PATTERN' | Set-Content -LiteralPath (Join-Path $rep 'MONITOR_BLOCKED.txt') -Encoding UTF8; break }
    $gmxProcess = Get-Process -Name gmx -ErrorAction SilentlyContinue
    $runner = Get-CimInstance Win32_Process -Filter "Name='powershell.exe'" | Where-Object { $_.ProcessId -ne $PID -and $_.CommandLine -match '-File .*run_apo_rep3_guarded_highpower\.ps1' }
    if (-not $gmxProcess -and -not $runner) {
        if ($time -le ($lastTime + 0.01)) { $stalledRestarts++ } else { $stalledRestarts = 0 }
        if ($stalledRestarts -ge 3) { Log "Three restarts without checkpoint advance at $time ps; monitor blocked"; 'NO_CHECKPOINT_ADVANCE' | Set-Content -LiteralPath (Join-Path $rep 'MONITOR_BLOCKED.txt') -Encoding UTF8; break }
        Log "No worker found; launching supervisor from $time ps"
        Start-Process -FilePath 'powershell.exe' -ArgumentList @('-NoProfile','-ExecutionPolicy','Bypass','-File',$runnerScript) -WindowStyle Hidden | Out-Null
        $lastTime = $time
        Start-Sleep -Seconds 20
    }
    Start-Sleep -Seconds 30
}
[void][Apo3PersistentExecutionState]::SetThreadExecutionState($continuous)
