$ErrorActionPreference = 'Continue'
$root = Split-Path -Parent (Split-Path -Parent $MyInvocation.MyCommand.Path)
$gmx = 'D:\APP\gromacs\bin\gmx.exe'
$outRoot = Join-Path $root 'comparative_analysis_v3'
$log = Join-Path $outRoot 'pbc_corrected_analysis_supervisor.log'
$stage = Join-Path $outRoot 'CURRENT_STAGE.txt'
New-Item -ItemType Directory -Force -Path $outRoot | Out-Null

Add-Type @'
using System;
using System.Runtime.InteropServices;
public static class PbcAnalysisExecutionState {
    [DllImport("kernel32.dll", SetLastError=true)]
    public static extern uint SetThreadExecutionState(uint esFlags);
}
'@

function Log([string]$m) { Add-Content -LiteralPath $log -Value "$(Get-Date -Format o)`t$m" -Encoding UTF8 }
function Run-Step([string]$label, [string]$marker, [scriptblock]$command, [string[]]$required) {
    if (Test-Path -LiteralPath $marker) { Log "SKIP $label"; return }
    "$label`t$(Get-Date -Format o)" | Set-Content -LiteralPath $stage -Encoding UTF8
    Log "START $label"
    $text = (& $command 2>&1 | Out-String)
    $code = $LASTEXITCODE
    Add-Content -LiteralPath $log -Value $text -Encoding UTF8
    if ($code -ne 0) { throw "$label failed with exit code $code" }
    foreach ($p in $required) {
        if (-not (Test-Path -LiteralPath $p) -or (Get-Item -LiteralPath $p).Length -eq 0) { throw "$label missing output $p" }
    }
    "PASS $(Get-Date -Format o)" | Set-Content -LiteralPath $marker -Encoding UTF8
    Log "PASS $label"
}
function Count-NumericRows([string]$path) {
    return @(Get-Content -LiteralPath $path | Where-Object { $_ -match '^\s*[-+0-9]' }).Count
}
function Require-Rows([string]$label, [string]$path, [int]$expected) {
    $rows = Count-NumericRows $path
    if ($rows -ne $expected) { throw "$label has $rows rows; expected $expected" }
}

$jobs = @(
    [pscustomobject]@{system='complex';rep='rep1';tpr='md_100ns.tpr'},
    [pscustomobject]@{system='complex';rep='rep2';tpr='md_100ns.tpr'},
    [pscustomobject]@{system='complex';rep='rep3';tpr='md_100ns.tpr'},
    [pscustomobject]@{system='apo';rep='rep1';tpr='md_100ns_post_extension.tpr'},
    [pscustomobject]@{system='apo';rep='rep2';tpr='md_100ns_post_extension.tpr'},
    [pscustomobject]@{system='apo';rep='rep3';tpr='md_100ns_post_extension.tpr'}
)

try {
    [void][PbcAnalysisExecutionState]::SetThreadExecutionState(2147483649)
    $originalComplexIndex = Join-Path $root 'systems\complex\build\index.ndx'
    $complexIndex = Join-Path $outRoot 'complex_analysis_v2.ndx'
    Run-Step 'build validated Protein_TPL index group' (Join-Path $outRoot 'COMPLEX_INDEX.done') {
        @('1 | 15','name 22 Protein_TPL','q') | & $gmx make_ndx -n $originalComplexIndex -o $complexIndex
    } @($complexIndex)
    if (-not (Select-String -LiteralPath $complexIndex -Pattern '^\[ Protein_TPL \]$' -Quiet)) { throw 'Protein_TPL group missing from corrected index' }

    foreach ($job in $jobs) {
        $tag = "$($job.system)_$($job.rep)"
        $repDir = Join-Path $root "systems\$($job.system)\$($job.rep)"
        $raw = Join-Path $repDir 'analysis\md_100ns_canonical_raw.xtc'
        $tpr = Join-Path $repDir $job.tpr
        $buildIndex = if ($job.system -eq 'complex') { $complexIndex } else { Join-Path $root 'systems\apo\build\index.ndx' }
        $chainIndex = Join-Path $repDir 'analysis\chain_groups_v2.ndx'
        $out = Join-Path $outRoot $tag
        New-Item -ItemType Directory -Force -Path $out | Out-Null
        foreach ($p in @($raw,$tpr,$buildIndex,$chainIndex)) { if (-not (Test-Path -LiteralPath $p)) { throw "$tag missing input $p" } }

        $clustered = Join-Path $out ($(if ($job.system -eq 'complex') {'clustered_nonwater.xtc'} else {'clustered_protein.xtc'}))
        if ($job.system -eq 'complex') {
            Run-Step "$tag protein+ligand PBC clustering" (Join-Path $out 'CLUSTER.done') {
                @('22','1','20') | & $gmx trjconv -s $tpr -f $raw -n $buildIndex -o $clustered -pbc cluster -center -ur compact
            } @($clustered)
        } else {
            Run-Step "$tag protein PBC clustering" (Join-Path $out 'CLUSTER.done') {
                @('1','1','1') | & $gmx trjconv -s $tpr -f $raw -n $buildIndex -o $clustered -pbc cluster -center -ur compact
            } @($clustered)
        }

        $checkMarker = Join-Path $out 'CLUSTER_CHECK.done'
        if (-not (Test-Path -LiteralPath $checkMarker)) {
            "$tag clustered trajectory integrity scan`t$(Get-Date -Format o)" | Set-Content -LiteralPath $stage -Encoding UTF8
            $checkText = (& $gmx check -f $clustered 2>&1 | Out-String)
            Add-Content -LiteralPath $log -Value $checkText -Encoding UTF8
            if ($LASTEXITCODE -ne 0 -or $checkText -notmatch 'Step\s+10001' -or $checkText -notmatch 'Coords\s+10001') { throw "$tag clustered trajectory failed 10001-frame integrity gate" }
            "PASS $(Get-Date -Format o); frames=10001; spacing_ps=10" | Set-Content -LiteralPath $checkMarker -Encoding UTF8
            Log "PASS $tag clustered trajectory integrity scan"
        }

        $subsetTpr = Join-Path $out 'analysis_subset.tpr'
        $subsetGroup = if ($job.system -eq 'complex') {'20'} else {'1'}
        Run-Step "$tag build topology matching clustered subset" (Join-Path $out 'SUBSET_TPR.done') {
            @($subsetGroup) | & $gmx convert-tpr -s $tpr -n $buildIndex -o $subsetTpr
        } @($subsetTpr)

        $subsetIndex = Join-Path $out 'analysis_subset.ndx'
        Run-Step "$tag build index matching clustered subset" (Join-Path $out 'SUBSET_INDEX.done') {
            @('q') | & $gmx make_ndx -f $subsetTpr -o $subsetIndex
        } @($subsetIndex)

        $nojump = Join-Path $out 'pbccluster_nojump.xtc'
        Run-Step "$tag temporal no-jump after molecular clustering" (Join-Path $out 'NOJUMP.done') {
            @('0') | & $gmx trjconv -s $subsetTpr -f $clustered -n $subsetIndex -o $nojump -pbc nojump
        } @($nojump)

        $fit = Join-Path $out 'pbccluster_proteinfit.xtc'
        Run-Step "$tag backbone fit after PBC clustering" (Join-Path $out 'FIT.done') {
            @('4','0') | & $gmx trjconv -s $subsetTpr -f $nojump -n $subsetIndex -o $fit -fit rot+trans
        } @($fit)

        $rmsd = Join-Path $out 'whole_backbone_rmsd_0_100ns.xvg'
        Run-Step "$tag whole-backbone RMSD" (Join-Path $out 'RMSD.done') {
            @('4','4') | & $gmx rms -s $subsetTpr -f $fit -n $subsetIndex -o $rmsd -tu ns -fit rot+trans
        } @($rmsd)
        Require-Rows "$tag RMSD" $rmsd 10001

        $chainB = Join-Path $out 'chainB_backbone_rmsd_0_100ns.xvg'
        Run-Step "$tag chain B backbone RMSD" (Join-Path $out 'CHAINB_RMSD.done') {
            @('2','2') | & $gmx rms -s $subsetTpr -f $fit -n $chainIndex -o $chainB -tu ns -fit rot+trans
        } @($chainB)
        Require-Rows "$tag chain B RMSD" $chainB 10001

        $chainD = Join-Path $out 'chainD_backbone_rmsd_0_100ns.xvg'
        Run-Step "$tag chain D backbone RMSD" (Join-Path $out 'CHAIND_RMSD.done') {
            @('3','3') | & $gmx rms -s $subsetTpr -f $fit -n $chainIndex -o $chainD -tu ns -fit rot+trans
        } @($chainD)
        Require-Rows "$tag chain D RMSD" $chainD 10001

        $com = Join-Path $out 'chainB_chainD_com_distance_pbc_nojump_0_100ns.xvg'
        Run-Step "$tag chain B-D COM distance on unrotated PBC trajectory" (Join-Path $out 'CHAIN_COM_PBC_NOJUMP.done') {
            & $gmx distance -s $subsetTpr -f $nojump -n $chainIndex -select 'com of group "chainB" plus com of group "chainD"' -oall $com -tu ns
        } @($com)
        Require-Rows "$tag chain COM distance" $com 10001

        $rg = Join-Path $out 'protein_rg_20_100ns.xvg'
        Run-Step "$tag protein Rg" (Join-Path $out 'RG.done') {
            @('1') | & $gmx gyrate -s $subsetTpr -f $fit -n $subsetIndex -o $rg -b 20000 -e 100000
        } @($rg)
        Require-Rows "$tag Rg" $rg 8001

        $sasa = Join-Path $out 'protein_sasa_20_100ns.xvg'
        Run-Step "$tag protein SASA" (Join-Path $out 'SASA.done') {
            @('1','1') | & $gmx sasa -s $subsetTpr -f $fit -n $subsetIndex -o $sasa -b 20000 -e 100000 -dt 100
        } @($sasa)
        Require-Rows "$tag SASA" $sasa 801

        $interHb = Join-Path $out 'interchain_hbonds_pbc_nojump_20_100ns.xvg'
        Run-Step "$tag inter-chain hydrogen bonds on unrotated PBC trajectory" (Join-Path $out 'INTERCHAIN_HB_PBC_NOJUMP.done') {
            @('0','1') | & $gmx hbond -s $subsetTpr -f $nojump -n $chainIndex -num $interHb -b 20000 -e 100000 -dt 100
        } @($interHb)
        Require-Rows "$tag inter-chain H-bonds" $interHb 801

        foreach ($chain in @([pscustomobject]@{name='B';group='2'}, [pscustomobject]@{name='D';group='3'})) {
            $local = Join-Path $out "chain$($chain.name)_backbone_localfit_rmsf_20_100ns.xvg"
            Run-Step "$tag chain $($chain.name) local-fit RMSF" (Join-Path $out "CHAIN$($chain.name)_LOCAL_RMSF.done") {
                @($chain.group) | & $gmx rmsf -s $subsetTpr -f $fit -n $chainIndex -o $local -res -fit -b 20000 -e 100000 -dt 100
            } @($local)
            Require-Rows "$tag chain $($chain.name) local RMSF" $local 93
        }

        if ($job.system -eq 'complex') {
            $pocketIndex = Join-Path $repDir 'analysis\pocket_groups.ndx'
            $ligRmsd = Join-Path $out 'ligand_rmsd_unwrapped_after_proteinfit_0_100ns.xvg'
            Run-Step "$tag unwrapped ligand RMSD after protein fit" (Join-Path $out 'LIGAND_RMSD.done') {
                @('15') | & $gmx rms -s $subsetTpr -f $fit -n $subsetIndex -o $ligRmsd -tu ns -fit none
            } @($ligRmsd)
            Require-Rows "$tag ligand RMSD" $ligRmsd 10001

            $pocket = Join-Path $out 'ligand_original_pocket_com_distance_unwrapped_0_100ns.xvg'
            Run-Step "$tag unwrapped ligand-original-pocket COM displacement" (Join-Path $out 'POCKET_DISTANCE.done') {
                & $gmx distance -s $subsetTpr -f $fit -n $pocketIndex -select 'com of group "ligand" plus com of group "pocket"' -oall $pocket -tu ns -nopbc
            } @($pocket)
            Require-Rows "$tag pocket distance" $pocket 10001

            $mindist = Join-Path $out 'ligand_protein_mindist_pbc_nojump_0_100ns.xvg'
            $contacts = Join-Path $out 'ligand_protein_contacts_pbc_nojump_0_100ns.xvg'
            Run-Step "$tag ligand-protein contacts on unrotated PBC trajectory" (Join-Path $out 'LIGAND_CONTACTS_PBC_NOJUMP.done') {
                @('1','15') | & $gmx mindist -s $subsetTpr -f $nojump -n $subsetIndex -od $mindist -on $contacts -d 0.45 -tu ns
            } @($mindist,$contacts)
            Require-Rows "$tag ligand minimum distance" $mindist 10001
            Require-Rows "$tag ligand contacts" $contacts 10001

            $ligHb = Join-Path $out 'ligand_protein_hbonds_pbc_nojump_20_100ns.xvg'
            Run-Step "$tag corrected ligand-protein hydrogen bonds on unrotated PBC trajectory" (Join-Path $out 'LIGAND_HB_PBC_NOJUMP.done') {
                @('1','15') | & $gmx hbond -s $subsetTpr -f $nojump -n $subsetIndex -num $ligHb -b 20000 -e 100000 -dt 100
            } @($ligHb)
            Require-Rows "$tag ligand-protein H-bonds" $ligHb 801

            $occ = Join-Path $out 'protein_residue_ligand_contact_occupancy_pbc_nojump_0p4nm_20_100ns.xvg'
            $selection = 'res_com of (group "Protein" and same residue as (group "Protein" and within 0.4 of group "TPL"))'
            Run-Step "$tag corrected residue contact occupancy on unrotated PBC trajectory" (Join-Path $out 'CONTACT_OCCUPANCY_PBC_NOJUMP.done') {
                & $gmx select -s $subsetTpr -f $nojump -n $subsetIndex -select $selection -of $occ -b 20000 -e 100000 -dt 100 -xvg none
            } @($occ)
            Require-Rows "$tag residue contact occupancy" $occ 186
        }

        "PASS $(Get-Date -Format o); pbc=cluster_then_nojump; center=Protein; frames=10001" | Set-Content -LiteralPath (Join-Path $out 'TRAJECTORY_METRICS_COMPLETE.txt') -Encoding UTF8
        Log "ALL METRICS COMPLETE $tag"
    }

    "completed=$(Get-Date -Format o); trajectories=6; pbc=cluster_then_nojump; canonical_inputs_preserved=true" | Set-Content -LiteralPath (Join-Path $outRoot 'PBC_CORRECTED_METRICS_V3_COMPLETE.txt') -Encoding UTF8
    "ALL_PBC_CORRECTED_METRICS_V3_COMPLETE`t$(Get-Date -Format o)" | Set-Content -LiteralPath $stage -Encoding UTF8
    Log 'ALL_PBC_CORRECTED_METRICS_V3_COMPLETE'
} catch {
    "FAILED`t$(Get-Date -Format o)`t$($_.Exception.Message)" | Set-Content -LiteralPath $stage -Encoding UTF8
    Log "FAILED $($_.Exception.Message)"
    throw
} finally {
    [void][PbcAnalysisExecutionState]::SetThreadExecutionState(2147483648)
}
