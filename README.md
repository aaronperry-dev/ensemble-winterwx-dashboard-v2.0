# GEFS Winter Weather

A 31-member GEFS ensemble viewer for **QPF, snowfall and precipitation type** (Revised Bourgouin method,
Birk et al. 2021, *Weather and Forecasting*). A Python pipeline runs on GitHub Actions four times a day, writes
small gzipped binary files, and deploys them with a static MapLibre/WebGL viewer to GitHub Pages. No data is
committed to the repo.

## Set up (once)

1. Create the repo and push this folder to `main`.
2. **Settings → Pages → Build and deployment → Source: GitHub Actions.**
3. **Actions → "GEFS winter weather" → Run workflow** (use `max_fhr = 48` for a quick first test).
4. Your site is at `https://<user>.github.io/<repo>/`.

The schedule is `40 4,10,16,22 * * *` UTC, about 4.5 hours after each cycle. The job waits up to 100 minutes
for files that are not posted yet and publishes whatever is complete (the manifest says `complete: false` while a
run is still filling in).

## What the pipeline does

| Step | Detail |
|---|---|
| Members | `c00` + `p01…p30` (31), 0.5° `pgrb2a` + `pgrb2b`, 3-hourly to F240 |
| Source order | AWS S3 `noaa-gefs-pds` → Google Cloud `gfs-ensemble-forecast-system` → NOMADS GRIB filter → NOMADS HTTP |
| Byte ranges | S3, GCS and NOMADS HTTP: read the `.idx`, then request only the messages needed (byte-adjacent ones merged). About 20 requests and 5 MB per member-step |
| NOMADS GRIB filter | A server-side CGI, so there is no file to byte-range. It subsets by variable, level and lat/lon box instead |
| Rate limits | Token bucket per host, job-wide totals divided across worker processes (aws 80/s, gcs 60/s, NOMADS 1.5/s with 2 connections). `Retry-After` is honored, NOMADS has a 404 budget (8/min) so missing files cannot get the runner IP blocked, and each source has a circuit breaker |
| Precip type | Wet-bulb profile → melting/refreezing energy → ProbIce, ProbSN, ProbPL, ProbFZRA (Eqs. 2, 7–9). Each member's 3-hr QPF is split among types by probability share |
| Snowfall | QPF × snow share × SLR / 10 (cm). SLR is **Kuchera** by default (clamped 5–30), or fixed with `--slr-mode fixed:12.5` / the workflow input |
| Heights | Hypsometric integration of virtual temperature, so no height download is needed |

Run it locally:

```bash
pip install -r requirements-dev.txt
python -m pytest -q tests
python -m pipeline.run --out site/data --max-fhr 48 --workers 6
python -m pipeline.site --out site
cd site && python -m http.server 8000
```

### Note on Eq. 7

The freezing-rain equation reads `ProbFZRA_i = −2.1·RE + 0.2·ME + 458`. Text extraction of the paper PDF turns the
minus sign into a "2" and the decimal point into a ":", which reads as `−22.1`. The code uses −2.1, and a regression test
(`test_shallow_cold_layer_is_freezing_rain`) guards it.

## Data format (`data/<yyyymmdd_cc>/…`, `data/latest.json`)

`latest.json` lists the run, domain, scales and every step. Each step has two files:

* `cum/fHHH.bin.gz`: **cumulative-since-init** QPF (`u16`, 0.05 mm per count) then snowfall (`u16`, 0.1 cm per count), each `[member][y][x]`.
* `ptype/fHHH.bin.gz`: `u8 [4][member][y][x]` = % RA, SN, FZRA, PL (0 where the member's 3-hr QPF is below 0.05 mm).

Cumulative storage means any rolling 3/6/12/24/48-hour total is just `cum(t) − cum(t−w)`, two small files fetched
by the browser. A full run is about 22 MB. Domain: 60–20°N, 135–50°W, 0.5° (81 × 171). The 64-byte `GWX1` header is
documented in `pipeline/encode.py`.

## Viewer

* MapLibre GL JS 4.7.1 (vendored), a self-hosted Natural Earth + US Census basemap (countries, states/provinces, counties, rivers, lakes), a graticule, and a canvas label layer with collision handling. There is no tile or font server.
* A custom WebGL2 layer paints the data with crisp class edges at any zoom, using the NWS-style ramps from the three reference images.
* Hover (desktop) or tap (phone) for a tooltip with mean, median, min, max, 10th/90th percentile, chance of exceeding thresholds, a dot strip of all 31 members and, for type, member counts and mean Bourgouin probabilities. Click to pin.
* Statistics: mean, median, percentiles, min/max, probability of at least *x*, and any single member. Units can be switched between inches and mm/cm.
* State is kept in the URL hash, so links are shareable.

## Layout

```
pipeline/        data pipeline (python -m pipeline.run, pipeline.site)
tests/           pytest suite, no network needed
index.html style.css viewer.js basemap.js wxlayer.js colors.js   viewer
basemap_*.json   boundaries and labels (build_basemap.py regenerates them)
maplibre-gl.*    vendored MapLibre
.github/workflows/gefs.yml
```

## Caveats

* The Google Cloud source could not be reached from the development sandbox, so its fallback is covered only by the mock-server tests. The AWS and NOMADS filter paths were run against real data.
* The 2-D GEFS fields are 0.5°, so precipitation type is smooth and not suitable for sub-county detail.
* Basemap: Natural Earth (public domain), US Census cartographic boundary data via us-atlas. Ensemble data: NOAA/NCEP GEFS (public domain).
