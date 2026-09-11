$ErrorActionPreference = 'SilentlyContinue'

Add-Type @'
using System;
using System.Runtime.InteropServices;
public static class MdExecutionState {
    [DllImport("kernel32.dll", SetLastError = true)]
    public static extern uint SetThreadExecutionState(uint esFlags);
}
'@

$ES_CONTINUOUS = 0x80000000
$ES_SYSTEM_REQUIRED = 0x00000001

while ($true) {
    $gmx = Get-Process -Name gmx -ErrorAction SilentlyContinue
    $supervisor = Get-CimInstance Win32_Process -Filter "Name='powershell.exe'" |
        Where-Object { $_.CommandLine -match 'run_complex_rep[123]_guarded_highpower\.ps1' }
    if (-not $gmx -and -not $supervisor) { break }
    [void][MdExecutionState]::SetThreadExecutionState($ES_CONTINUOUS -bor $ES_SYSTEM_REQUIRED)
    Start-Sleep -Seconds 30
}

[void][MdExecutionState]::SetThreadExecutionState($ES_CONTINUOUS)
