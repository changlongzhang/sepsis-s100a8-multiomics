from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor, as_completed
from itertools import product
from pathlib import Path
import csv
import json
import re
import subprocess

import numpy as np


ROOT = Path(r"D:\桌面\sepsis\S100A8_triptolide_16OH_comparative_docking_topjournal")
TPL_ROOT = Path(r"D:\桌面\sepsis\S100A8_triptolide_docking_topjournal")
PREP = ROOT / "input" / "prepared"
TPL_PREP = TPL_ROOT / "input" / "prepared"
OUT = ROOT / "results" / "blind_tiling"
OUT.mkdir(parents=True, exist_ok=True)
VINA = TPL_ROOT / "scripts" / "vina_1.2.7_win.exe"
LIGAND = PREP / "16-hydroxytriptolide_primary.pdbqt"
extents = json.loads((TPL_PREP / "receptor_extents.json").read_text(encoding="utf-8"))


def centers(low: float, high: float, size: float = 28.0, pad: float = 4.0) -> list[float]:
    low, high = low - pad, high + pad
    count = max(1, int(np.ceil((high - low) / size)))
    if count == 1:
        return [(low + high) / 2]
    return np.linspace(low + size / 2, high - size / 2, count).tolist()


jobs = []
for receptor in ("AC", "BD"):
    extent = extents[receptor]
    axes = [centers(extent["min"][index], extent["max"][index]) for index in range(3)]
    for index, center in enumerate(product(*axes), start=1):
        jobs.append((receptor, f"{receptor}_tile{index:02d}", center))


def run(job: tuple) -> dict:
    receptor, stem, center = job
    output = OUT / f"{stem}.pdbqt"
    command = [
        str(VINA), "--receptor", str(TPL_PREP / f"S100A8_{receptor}.pdbqt"),
        "--ligand", str(LIGAND),
        "--center_x", str(center[0]), "--center_y", str(center[1]), "--center_z", str(center[2]),
        "--size_x", "28", "--size_y", "28", "--size_z", "28",
        "--exhaustiveness", "32", "--num_modes", "10", "--energy_range", "5",
        "--seed", str(107985 + len(stem) * 101 + int(abs(sum(center)) * 10)),
        "--cpu", "2", "--out", str(output),
    ]
    process = subprocess.run(command, capture_output=True, text=True, encoding="utf-8", errors="replace")
    (OUT / f"{stem}.log").write_text(
        "COMMAND\n" + subprocess.list2cmdline(command) + "\n\nSTDOUT\n" + process.stdout +
        "\nSTDERR\n" + process.stderr,
        encoding="utf-8",
    )
    match = re.search(r"^\s*1\s+(-?\d+(?:\.\d+)?)", process.stdout, re.MULTILINE)
    return {
        "receptor": receptor, "tile": stem,
        "center_x": center[0], "center_y": center[1], "center_z": center[2], "size": 28,
        "best_score_kcal_mol": float(match.group(1)) if match else None,
        "returncode": process.returncode,
    }


rows = []
with ThreadPoolExecutor(max_workers=4) as executor:
    futures = [executor.submit(run, job) for job in jobs]
    for future in as_completed(futures):
        row = future.result()
        rows.append(row)
        print(row, flush=True)
rows.sort(key=lambda row: (row["best_score_kcal_mol"] is None, row["best_score_kcal_mol"] or 999))
if any(row["returncode"] != 0 or row["best_score_kcal_mol"] is None for row in rows):
    raise RuntimeError("At least one blind docking job failed")
with (OUT / "blind_tiling_summary.csv").open("w", newline="", encoding="utf-8-sig") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)
(OUT / "BLIND_TILING_COMPLETE.txt").write_text(
    f"completed; jobs={len(rows)}; exhaustiveness=32; box=28A; failures=0\n", encoding="utf-8")
print("BEST", rows[:5])
