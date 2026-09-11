from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
import csv
import re
import subprocess


ROOT = Path(r"D:\桌面\sepsis\S100A8_triptolide_16OH_comparative_docking_topjournal")
TPL_ROOT = Path(r"D:\桌面\sepsis\S100A8_triptolide_docking_topjournal")
PREP = ROOT / "input" / "prepared"
TPL_PREP = TPL_ROOT / "input" / "prepared"
OUT = ROOT / "results" / "local_replicates"
OUT.mkdir(parents=True, exist_ok=True)
VINA = TPL_ROOT / "scripts" / "vina_1.2.7_win.exe"
LIGAND = PREP / "16-hydroxytriptolide_primary.pdbqt"
CENTERS = {"AC": [6.864, 11.296, -57.687], "BD": [15.878, 20.754, -89.453]}
SEEDS = [107985, 207985, 307985, 407985, 507985, 607985, 707985, 807985, 907985, 1007985]


def run(receptor: str, center: list[float], seed: int, scoring: str) -> dict:
    stem = f"{receptor}_{scoring}_seed{seed}"
    output = OUT / f"{stem}.pdbqt"
    command = [
        str(VINA), "--receptor", str(TPL_PREP / f"S100A8_{receptor}.pdbqt"),
        "--ligand", str(LIGAND),
        "--center_x", str(center[0]), "--center_y", str(center[1]), "--center_z", str(center[2]),
        "--size_x", "24", "--size_y", "24", "--size_z", "24",
        "--exhaustiveness", "64", "--num_modes", "20", "--energy_range", "5",
        "--seed", str(seed), "--cpu", "2", "--scoring", scoring, "--out", str(output),
    ]
    process = subprocess.run(command, capture_output=True, text=True, encoding="utf-8", errors="replace")
    (OUT / f"{stem}.log").write_text(
        "COMMAND\n" + subprocess.list2cmdline(command) + "\n\nSTDOUT\n" + process.stdout +
        "\nSTDERR\n" + process.stderr,
        encoding="utf-8",
    )
    scores = [float(value) for value in re.findall(r"^\s*\d+\s+(-?\d+(?:\.\d+)?)", process.stdout, re.MULTILINE)]
    return {
        "receptor": receptor, "scoring": scoring, "seed": seed,
        "best_score_kcal_mol": scores[0] if scores else None,
        "n_modes": len(scores), "returncode": process.returncode, "pose_file": output.name,
    }


jobs = [(receptor, center, seed, scoring)
        for receptor, center in CENTERS.items() for seed in SEEDS for scoring in ("vina", "vinardo")]
rows = []
with ThreadPoolExecutor(max_workers=4) as executor:
    futures = [executor.submit(run, *job) for job in jobs]
    for future in as_completed(futures):
        row = future.result()
        rows.append(row)
        print(row, flush=True)
rows.sort(key=lambda row: (row["receptor"], row["scoring"], row["seed"]))
if any(row["returncode"] != 0 or row["n_modes"] == 0 for row in rows):
    raise RuntimeError("At least one local docking job failed")
with (OUT / "local_replicates_summary.csv").open("w", newline="", encoding="utf-8-sig") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)
(OUT / "LOCAL_REPLICATES_COMPLETE.txt").write_text(
    f"completed; jobs={len(rows)}; receptors=2; scoring=2; seeds=10; exhaustiveness=64; failures=0\n",
    encoding="utf-8",
)
