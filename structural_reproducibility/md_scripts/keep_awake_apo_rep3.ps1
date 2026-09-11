$ErrorActionPreference='SilentlyContinue'
Add-Type @'
using System;
using System.Runtime.InteropServices;
public static class Apo3ExecutionState { [DllImport("kernel32.dll", SetLastError=true)] public static extern uint SetThreadExecutionState(uint esFlags); }
'@
$c=0x80000000;$s=0x00000001
while($true){$g=Get-Process -Name gmx -ErrorAction SilentlyContinue;$r=Get-CimInstance Win32_Process -Filter "Name='powershell.exe'"|Where-Object{$_.CommandLine-match'(prepare_apo_rep3_resumable|run_apo_rep3_guarded_highpower)\.ps1'};if(-not$g-and-not$r){break};[void][Apo3ExecutionState]::SetThreadExecutionState($c-bor$s);Start-Sleep 30}
[void][Apo3ExecutionState]::SetThreadExecutionState($c)
