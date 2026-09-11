#!/usr/bin/env python3
"""Reassemble the two large RDS objects distributed across S3-S10 Files."""
from __future__ import annotations
import argparse
import hashlib
import json
import zipfile
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("--archives-dir", default=".", help="Directory containing S3-S10 zip files")
parser.add_argument("--output-dir", default="reassembled_large_rds")
args = parser.parse_args()

archives = Path(args.archives_dir)
output = Path(args.output_dir)
output.mkdir(parents=True, exist_ok=True)

with open(Path(__file__).with_name("large_rds_manifest.json"), "r", encoding="utf-8") as fh:
    manifest = json.load(fh)

for item in manifest["large_files"]:
    destination = output / item["original_path"]
    destination.parent.mkdir(parents=True, exist_ok=True)
    digest = hashlib.sha256()
    with destination.open("wb") as out:
        for part in item["parts"]:
            with zipfile.ZipFile(archives / part["archive"]) as zf:
                data = zf.read(part["member"])
            if hashlib.sha256(data).hexdigest() != part["sha256"]:
                raise RuntimeError(f"Checksum failure: {part['archive']}")
            out.write(data)
            digest.update(data)
    if digest.hexdigest() != item["sha256"]:
        raise RuntimeError(f"Reassembled checksum failure: {item['original_path']}")
    print(f"Verified: {destination}")
