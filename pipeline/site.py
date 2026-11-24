"""Assemble the static site directory (web assets + data) for GitHub Pages.

    python -m pipeline.site --out site
"""
from __future__ import annotations

import argparse
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
ASSETS = ["index.html", "style.css", "colors.js", "wxlayer.js", "basemap.js", "ensstats.js", "viewer.js",
          "maplibre-gl.js", "maplibre-gl.css", "maplibre-gl.LICENSE.txt"]
ASSET_GLOBS = ["basemap_*.json"]


def build(out: Path) -> list[str]:
    out.mkdir(parents=True, exist_ok=True)
    copied = []
    names = list(ASSETS) + [p.name for g in ASSET_GLOBS for p in ROOT.glob(g)]
    for n in names:
        shutil.copy2(ROOT / n, out / n)
        copied.append(n)
    (out / ".nojekyll").write_text("")
    return copied


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="site")
    a = ap.parse_args()
    print(f"copied {len(build(Path(a.out)))} files to {a.out}")
