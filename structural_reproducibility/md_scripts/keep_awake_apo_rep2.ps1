$ErrorActionPreference = 'SilentlyContinue'
Add-Type @'
using System;
using System.Runtime.InteropServices;
public static class ApoMdExecutionState {
    [DllImport("kernel32.dll", SetLastError = true)]
    public static extern uint SetThreadExecutionState(uint esFlags);
}
'@
$continuous=0x80000000;$systemRequired=0x00000001
while ($true) {
    $gmx=Get-Process -Name gmx -ErrorAction SilentlyContinue
    $runner=Get-CimInstance Win32_Process -Filter "Name='powershell.exe'" | Where-Object {$_.CommandLine -match '(prepare_apo_rep2_resumable|run_apo_rep2_guarded_highpower)\.ps1'}
    if(-not $gmx -and -not $runner){break}
    [void][ApoMdExecutionState]::SetThreadExecutionState($continuous -bor $systemRequired)
    Start-Sleep 30
}
[void][ApoMdExecutionState]::SetThreadExecutionState($continuous)
