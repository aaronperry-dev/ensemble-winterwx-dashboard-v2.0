#!/usr/bin/env python3
"""Build the static basemap layers used by the viewer (run once; outputs are committed).

    python build_basemap.py [--ne-dir DIR] [--us-atlas counties-10m.json]

Inputs (downloaded automatically when --ne-dir is not given):
  * Natural Earth 1:10m vector data (public domain)  - coast, lakes, rivers, countries, states/provinces,
    populated places, marine names.     https://github.com/nvkelso/natural-earth-vector
  * us-atlas counties-10m.json (US Census cartographic boundaries, public domain), from the npm
    package `us-atlas`.                 https://github.com/topojson/us-atlas

Outputs (repo root, compact GeoJSON / JSON, coordinates rounded to 3 decimals):
  basemap_land.json  basemap_lakes.json  basemap_rivers.json  basemap_coast.json
  basemap_borders.json (international)   basemap_admin1.json (state/province polygons)
  basemap_counties.json (US counties)    basemap_labels.json (places, regions, water names)
"""
from __future__ import annotations

import argparse
import json
import sys
import urllib.request
from pathlib import Path

from shapely.geometry import box, mapping, shape
from shapely.ops import unary_union

BBOX = box(-175.0, 0.0, -30.0, 80.0)      # everything the viewer can pan to
NE = "https://raw.githubusercontent.com/nvkelso/natural-earth-vector/master/geojson/{}.geojson"
FILES = ["ne_10m_land", "ne_10m_lakes", "ne_10m_rivers_lake_centerlines", "ne_10m_coastline",
         "ne_10m_admin_0_boundary_lines_land", "ne_10m_admin_1_states_provinces",
         "ne_10m_admin_0_countries", "ne_10m_populated_places_simple", "ne_10m_geography_marine_polys"]
ROOT = Path(__file__).parent


def load(name: str, ne_dir: Path | None):
    if ne_dir and (ne_dir / f"{name}.geojson").exists():
        return json.loads((ne_dir / f"{name}.geojson").read_text())
    print("downloading", name, file=sys.stderr)
    with urllib.request.urlopen(NE.format(name), timeout=120) as r:
        return json.loads(r.read())


def rnd(obj, nd=3):
    if isinstance(obj, (list, tuple)):
        return [rnd(o, nd) for o in obj]
    if isinstance(obj, float):
        return round(obj, nd)
    return obj


def feat(geom, props=None):
    return {"type": "Feature", "properties": props or {}, "geometry": {"type": geom["type"], "coordinates": rnd(geom["coordinates"])}}


def clip_simplify(geom, tol):
    g = geom.intersection(BBOX)
    if g.is_empty:
        return None
    g = g.simplify(tol, preserve_topology=True)
    return None if g.is_empty else g


def dump(name: str, fc: dict):
    p = ROOT / name
    p.write_text(json.dumps(fc, separators=(",", ":"), ensure_ascii=False))
    print(f"{name:26s} {p.stat().st_size/1e6:6.2f} MB  {len(fc.get('features', [])) if isinstance(fc, dict) and 'features' in fc else ''}")


def fc(features):
    return {"type": "FeatureCollection", "features": features}


# ----------------------------------------------------------------------------- us-atlas decode
def topo_arcs(topo):
    sx, sy = topo["transform"]["scale"]
    tx, ty = topo["transform"]["translate"]
    arcs = []
    for a in topo["arcs"]:
        x = y = 0
        pts = []
        for dx, dy in a:
            x += dx
            y += dy
            pts.append((x * sx + tx, y * sy + ty))
        arcs.append(pts)
    return arcs


def topo_ring(arc_ids, arcs):
    ring = []
    for i in arc_ids:
        pts = arcs[i] if i >= 0 else arcs[~i][::-1]
        ring.extend(pts if not ring else pts[1:])
    return ring


def topo_polygons(geom, arcs):
    if geom["type"] == "Polygon":
        return [[topo_ring(r, arcs) for r in geom["arcs"]]]
    return [[topo_ring(r, arcs) for r in poly] for poly in geom["arcs"]]


def county_label(fips: str, name: str) -> str:
    st, co = fips[:2], int(fips[2:])
    if st == "22":
        return f"{name} Parish"
    if st == "02":
        return name
    if (st in ("51", "24", "29", "32") and co >= 500) or fips == "11001":
        return name if name.lower().endswith("city") else f"{name} (city)"
    return f"{name} County"


def build_counties(us_atlas: Path, postal_by_fips: dict[str, str]):
    topo = json.loads(us_atlas.read_text())
    arcs = topo_arcs(topo)
    states = {g["id"]: g["properties"]["name"] for g in topo["objects"]["states"]["geometries"]}
    feats = []
    for g in topo["objects"]["counties"]["geometries"]:
        fips = str(g["id"]).zfill(5)
        polys = topo_polygons(g, arcs)
        geom = {"type": "MultiPolygon" if len(polys) > 1 else "Polygon",
                "coordinates": polys if len(polys) > 1 else polys[0]}
        name = g.get("properties", {}).get("name", fips)
        feats.append(feat(geom, {"n": county_label(fips, name), "s": postal_by_fips.get(fips[:2], ""),
                                 "f": fips}))
    return fc(feats)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ne-dir", type=Path)
    ap.add_argument("--us-atlas", type=Path, default=Path("counties-10m.json"))
    a = ap.parse_args()
    d = {n: load(n, a.ne_dir) for n in FILES}

    # land ------------------------------------------------------------------
    land = unary_union([shape(f["geometry"]) for f in d["ne_10m_land"]["features"]])
    g = clip_simplify(land, 0.008)
    dump("basemap_land.json", fc([feat(mapping(g))]))

    # lakes (large ones only) ----------------------------------------------------
    lakes, lake_labels = [], []
    for f in d["ne_10m_lakes"]["features"]:
        s = clip_simplify(shape(f["geometry"]), 0.008)
        if s is None or s.area < 0.02:
            continue
        lakes.append(feat(mapping(s)))
        p = f["properties"]
        if p.get("name") and p.get("scalerank", 99) <= 3:
            rp = s.representative_point()
            lake_labels.append([p["name"], round(rp.x, 2), round(rp.y, 2), p.get("min_label", 4)])
    dump("basemap_lakes.json", fc(lakes))

    # rivers --------------------------------------------------------------------
    riv = []
    for f in d["ne_10m_rivers_lake_centerlines"]["features"]:
        p = f["properties"]
        if p.get("featurecla") != "River" or p.get("scalerank", 99) > 5:
            continue
        s = clip_simplify(shape(f["geometry"]), 0.02)
        if s is not None:
            riv.append(feat(mapping(s), {"n": p.get("name", "")}))
    dump("basemap_rivers.json", fc(riv))

    # coastline ----------------------------------------------------------------
    co = unary_union([shape(f["geometry"]) for f in d["ne_10m_coastline"]["features"]])
    g = clip_simplify(co, 0.008)
    dump("basemap_coast.json", fc([feat(mapping(g))]))

    # international borders --------------------------------------------------------
    b = []
    for f in d["ne_10m_admin_0_boundary_lines_land"]["features"]:
        s = clip_simplify(shape(f["geometry"]), 0.01)
        if s is not None:
            b.append(feat(mapping(s)))
    dump("basemap_borders.json", fc(b))

    # admin-1 (US / Canada / Mexico) -------------------------------------------------
    adm1, adm1_labels, postal_by_fips = [], [], {}
    for f in d["ne_10m_admin_1_states_provinces"]["features"]:
        p = f["properties"]
        if p["adm0_a3"] not in ("USA", "CAN", "MEX"):
            continue
        s = clip_simplify(shape(f["geometry"]), 0.015)
        if s is None:
            continue
        if p["adm0_a3"] == "USA" and p.get("fips", "").startswith("US"):
            postal_by_fips[p["fips"][2:]] = p.get("postal", "")
        adm1.append(feat(mapping(s), {"n": p["name"], "c": p["adm0_a3"], "p": p.get("postal") or ""}))
        adm1_labels.append([p["name"], round(p["longitude"], 2), round(p["latitude"], 2), p.get("min_label", 5),
                            p["adm0_a3"]])
    dump("basemap_admin1.json", fc(adm1))

    # counties ------------------------------------------------------------------------
    if a.us_atlas.exists():
        dump("basemap_counties.json", build_counties(a.us_atlas, postal_by_fips))
    else:
        print("WARNING: us-atlas counties-10m.json not found; skipping counties", file=sys.stderr)

    # labels -----------------------------------------------------------------------------
    inb = lambda x, y: -170 <= x <= -45 and 8 <= y <= 74
    countries = []
    for f in d["ne_10m_admin_0_countries"]["features"]:
        p = f["properties"]
        x, y = p.get("LABEL_X"), p.get("LABEL_Y")
        if x is not None and inb(x, y):
            countries.append([p["NAME"], round(x, 2), round(y, 2), p.get("LABELRANK", 5)])
    places = []
    for f in d["ne_10m_populated_places_simple"]["features"]:
        p = f["properties"]
        x, y = p["longitude"], p["latitude"]
        if not inb(x, y) or p["adm0_a3"] not in ("USA", "CAN", "MEX", "CUB", "BHS", "GTM", "BLZ", "HND", "SLV", "JAM", "HTI", "DOM", "GRL", "SPM"):
            continue
        if p["min_zoom"] > 7.6 and p["pop_max"] < 15000:
            continue
        cap = 2 if p["adm0cap"] else (1 if p["featurecla"] == "Admin-1 capital" else 0)
        places.append([p["name"], round(x, 3), round(y, 3), p["pop_max"], round(p["min_zoom"], 1), cap])
    marine = []
    for f in d["ne_10m_geography_marine_polys"]["features"]:
        p = f["properties"]
        if p.get("featurecla") not in ("ocean", "sea", "gulf", "bay", "strait", "sound") or p.get("min_label", 9) > 5.5:
            continue
        s = shape(f["geometry"]).intersection(BBOX)
        if s.is_empty:
            continue
        rp = s.representative_point()
        if inb(rp.x, rp.y):
            marine.append([p["name"].title() if p["name"].isupper() else p["name"], round(rp.x, 2), round(rp.y, 2),
                           p.get("min_label", 3), p["featurecla"]])
    lab = {"countries": countries, "admin1": adm1_labels, "places": places, "marine": marine, "lakes": lake_labels}
    (ROOT / "basemap_labels.json").write_text(json.dumps(lab, separators=(",", ":"), ensure_ascii=False))
    print("basemap_labels.json         ", {k: len(v) for k, v in lab.items()})


if __name__ == "__main__":
    main()
