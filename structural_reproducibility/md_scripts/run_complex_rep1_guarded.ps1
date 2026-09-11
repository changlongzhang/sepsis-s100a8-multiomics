$ErrorActionPreference = 'Continue'

$root = Split-Path -Parent (Split-Path -Parent $MyInvocation.MyCommand.Path)
$gmx = 'D:\APP\gromacs\bin\gmx.exe'
$rep = Join-Path $root 'systems\complex\rep1'
$tpr = Join-Path $rep 'md_100ns.tpr'
$cpt = Join-Path $rep 'md_100ns.cpt'
$deffnm = Join-Path $rep 'md_100ns'
$audit = Join-Path $rep 'guarded_segments.tsv'
$consoleLog = Join-Path $rep 'guarded_supervisor.log'
$targetPs = 100000.0

function Write-SupervisorLog([string]$message) {
    $line = "$(Get-Date -Format o)`t$message"
    Add-Content -LiteralPath $consoleLog -Value $line -Encoding UTF8
}

function Get-CheckpointTimePs {
    $checkOutput = (& $gmx check -f $cpt 2>&1 | Out-String)
    if ($LASTEXITCODE -ne 0) {
        throw "gmx check failed for checkpoint: $checkOutput"
    }
    $match = [regex]::Match($checkOutput, 'Last frame\s+-?\d+\s+time\s+([0-9.]+)')
    if (-not $match.Success) {
        throw "Could not parse checkpoint time: $checkOutput"
    }
    return [double]$match.Groups[1].Value
}

if (-not (Test-Path -LiteralPath $gmx)) { throw "GROMACS executable missing: $gmx" }
if (-not (Test-Path -LiteralPath $tpr)) { throw "TPR missing: $tpr" }
if (-not (Test-Path -LiteralPath $cpt)) { throw "Checkpoint missing: $cpt" }

if (-not (Test-Path -LiteralPath $audit)) {
    "segment`tstart_iso`tend_iso`tbefore_ps`tafter_ps`texit_code`tcpt_sha256`tstatus" |
        Set-Content -LiteralPath $audit -Encoding UTF8
}

$segment = 0
while ($true) {
    $beforePs = Get-CheckpointTimePs
    if ($beforePs -ge ($targetPs - 0.1)) {
        "COMPLETED $(Get-Date -Format o) checkpoint_ps=$beforePs" |
            Set-Content -LiteralPath (Join-Path $rep 'COMPLETED.txt') -Encoding UTF8
        Write-SupervisorLog "Target reached at $beforePs ps"
        break
    }

    # Yield to interactive work when the machine is already heavily loaded.
    $cpuLoad = [double](Get-CimInstance Win32_Processor | Measure-Object -Property LoadPercentage -Average).Average
    $os = Get-CimInstance Win32_OperatingSystem
    $freeMemoryGB = [math]::Round($os.FreePhysicalMemory / 1MB, 2)
    $gpuFreeMB = 999999
    try {
        $gpuFreeMB = [double]((& nvidia-smi --query-gpu=memory.free --format=csv,noheader,nounits | Select-Object -First 1).Trim())
    } catch {
        Write-SupervisorLog "GPU availability query failed; stopping rather than changing compute backend"
        throw
    }

    if ($cpuLoad -ge 85 -or $freeMemoryGB -lt 2.5 -or $gpuFreeMB -lt 1500) {
        Write-SupervisorLog "Yielding: cpu_load=$cpuLoad free_memory_GB=$freeMemoryGB gpu_free_MB=$gpuFreeMB"
        Start-Sleep -Seconds 300
        continue
    }

    $segment++
    $startIso = Get-Date -Format o
    Write-SupervisorLog "Starting segment $segment from $beforePs ps; cpu_load=$cpuLoad free_memory_GB=$freeMemoryGB gpu_free_MB=$gpuFreeMB"

    $runOutput = (& $gmx mdrun -s $tpr -cpi $cpt -deffnm $deffnm `
        -ntmpi 1 -ntomp 4 -nb gpu -pme cpu -bonded cpu -update cpu `
        -pin off -noappend -cpt 2 -maxh 0.50 2>&1 | Out-String)
    $exitCode = $LASTEXITCODE
    Add-Content -LiteralPath $consoleLog -Value $runOutput -Encoding UTF8

    $endIso = Get-Date -Format o
    if ($exitCode -ne 0) {
        $hash = if (Test-Path -LiteralPath $cpt) { (Get-FileHash -LiteralPath $cpt -Algorithm SHA256).Hash } else { '' }
        "$segment`t$startIso`t$endIso`t$beforePs`tNA`t$exitCode`t$hash`tFAILED" |
            Add-Content -LiteralPath $audit -Encoding UTF8
        Write-SupervisorLog "Segment $segment failed with exit code $exitCode"
        throw "Guarded GROMACS segment failed with exit code $exitCode"
    }

    $afterPs = Get-CheckpointTimePs
    $hash = (Get-FileHash -LiteralPath $cpt -Algorithm SHA256).Hash
    if ($afterPs -le $beforePs) {
        "$segment`t$startIso`t$endIso`t$beforePs`t$afterPs`t$exitCode`t$hash`tNO_PROGRESS" |
            Add-Content -LiteralPath $audit -Encoding UTF8
        throw "Checkpoint did not advance: before=$beforePs after=$afterPs"
    }

    "$segment`t$startIso`t$endIso`t$beforePs`t$afterPs`t$exitCode`t$hash`tOK" |
        Add-Content -LiteralPath $audit -Encoding UTF8
    Write-SupervisorLog "Segment $segment completed: $beforePs -> $afterPs ps; checkpoint_sha256=$hash"
}
