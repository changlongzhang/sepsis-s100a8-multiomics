from pathlib import Path
import csv
import fitz
from PIL import Image
from datetime import datetime


OUT = Path(r"D:\桌面\sepsis\文章\返修投稿新主图")
STEMS = [f"Fig{i}" for i in range(1, 9)] + [f"S{i}_fig" for i in range(1, 9)]


def main():
    rows = []
    for stem in STEMS:
        ai, pdf, tif = (OUT / f"{stem}.{ext}" for ext in ("ai", "pdf", "tif"))
        with Image.open(tif) as im:
            width, height = im.size
            dpi = tuple(round(v) for v in im.info.get("dpi", (0, 0)))
            mode = im.mode
            compression = im.info.get("compression", "")
        doc = fitz.open(pdf)
        pages = len(doc)
        text_chars = sum(len(p.get_text()) for p in doc)
        doc.close()
        sizes = [p.stat().st_size for p in (ai, pdf, tif)]
        status = "PASS" if (
            all(p.exists() and p.stat().st_size > 0 for p in (ai, pdf, tif))
            and 789 <= width <= 2250 and height <= 2625
            and dpi[0] >= 299 and dpi[1] >= 299
            and mode in {"RGB", "RGBA"} and pages == 1
            and max(sizes) < 10 * 1024 * 1024
        ) else "CHECK"
        rows.append({
            "figure": stem, "width_px": width, "height_px": height,
            "dpi_x": dpi[0], "dpi_y": dpi[1], "mode": mode,
            "compression": compression, "pdf_pages": pages,
            "pdf_text_chars": text_chars,
            "ai_mb": round(sizes[0] / 1048576, 2),
            "pdf_mb": round(sizes[1] / 1048576, 2),
            "tif_mb": round(sizes[2] / 1048576, 2), "status": status,
        })
    with (OUT / "PLOS_figure_QA.csv").open("w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=rows[0].keys())
        w.writeheader(); w.writerows(rows)
    manifest = []
    for stem in STEMS:
        for ext in ("ai", "pdf", "tif"):
            path = OUT / f"{stem}.{ext}"
            manifest.append({"figure": stem, "format": ext.upper(), "filename": path.name,
                             "bytes": path.stat().st_size,
                             "modified": datetime.fromtimestamp(path.stat().st_mtime).isoformat(timespec="seconds")})
    with (OUT / "figure_export_manifest.csv").open("w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=manifest[0].keys())
        w.writeheader(); w.writerows(manifest)
    for row in rows:
        print(row)


if __name__ == "__main__":
    main()
