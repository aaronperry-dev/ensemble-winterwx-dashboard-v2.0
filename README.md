# Ensemble Winter Weather: GEFS · ECMWF IFS ENS · GEPS · blend

An ensemble viewer for **QPF, snowfall, freezing rain, sleet and precipitation type** (Revised Bourgouin method,
Birk et al. 2021, *Weather and Forecasting*) from three ensemble families, shown separately or **blended into one
product**: NOAA **GEFS** (31 members), ECMWF **IFS ENS** (50 members) and Environment Canada's **GEPS** (21 members).
A Python pipeline runs on GitHub Actions, writes small gzipped binary files per model, and deploys them with a
static MapLibre/WebGL viewer to GitHub Pages. No data is committed to the repo.

## The three models and the blend

| | GEFS | IFS ENS | GEPS |
|---|---|---|---|
| Members | 31 (control + 30) | 50 (perturbed only; the control is not in the open data) | 21 (control + 20) |
| Cycles | 00/06/12/18Z | 00/12Z | 00/12Z |
| Hours | 3-hourly to F240 | 3-hourly to F144, then 6-hourly | 3-hourly to F192, then 6-hourly |
| Source | AWS S3 → GCS → NOMADS | `data.ecmwf.int` (AWS copy as fallback) | MSC datamart `dd.weather.gc.ca` |
| Grid used | 0.5° | 0.25°, sub-sampled to the 0.5° points (exact, no interpolation) | 0.5° |
| Pressure levels in the profile | 15 (1000–400 hPa) | 5 (1000, 925, 850, 700, 500) | 5 (same) |
| Extra field | `CPOFP` (percent frozen), used for the freezing-rain cap | none | none |
| Fetch | byte ranges per member (`.idx`) | byte ranges per member (JSON `.index`), 14 fields each | whole files, one per field, all members inside |
| Credit | NOAA (public domain) | ECMWF open data, CC BY 4.0 | ECCC, Open Government Licence – Canada |

All three go through the same decision tree (below) and are written in the same file format on the same 0.5° grid,
so the viewer treats them alike.

**The blend is made in the browser**, not in the pipeline: it pools the members of every model that shares the run
time (up to 102 members) and applies the same statistics. Nothing is duplicated on the server.

* **Weighting.** Default *equal by model*: each model carries the same total weight (a third each with all three), so
  IFS's 50 members do not drown out GEPS's 21. The alternative, *equal by member*, counts every member the same. Mean,
  probabilities and the type shares are weighted averages; percentiles are the weighted analogue of linear interpolation
  between ranked members (identical to the usual one for equal weights); PMM and LPMM match amplitudes by cumulative
  weight instead of by count. Single-model views use the original unweighted code, unchanged.
* **Which run time.** IFS and GEPS exist only at 00/12Z. For the blend, GEFS also publishes the 00/12Z cycle that IFS and
  GEPS have (`synoptic.json`, picked by `pipeline.run --synoptic-too`); on a 06/18Z run that is a second, mostly reused,
  GEFS build. The blend uses the run time shared by the most models. A model still on an older run is left out and the
  panel says so, so the blend never mixes initialisation times.
* **Hours.** The blend has the hours every included model has: 3-hourly to F144, then 6-hourly. A "step" is the gap to
  the previous shared hour (6 h beyond F144), and the accumulation windows stay exact.
* The tooltip colours each member's dot by model and lists each model's mean next to the blend's statistics.

**Fewer pressure levels.** IFS and GEPS open data carry 5 levels in the 1000–500 hPa range that matter (versus 15 for
GEFS). Re-running the Bourgouin calculation on one GEFS cycle (8 members, 16 member-steps, 92,000 wet cells) with only
those 5 levels gave the same dominant type as the full profile in 97.4 % of wet cells (7 levels: 98.3 %), but about 45 %
more freezing-rain cells (3,177 vs 2,211 with 5 levels), and interpolating the coarse profile back to 15 levels did not
remove that. So treat IFS and GEPS freezing rain as somewhat generous compared with GEFS, and remember it is the one
type the ASOS check found hard to forecast anyway. The blend dilutes it by the model weights.

**IFS and GEPS differ from GEFS in two more ways.** The percent-frozen cap needs `CPOFP`, which only GEFS publishes, so it is
not applied to them. IFS gives 2 m dew point instead of humidity (converted with the Bolton formula), and IFS/GEPS
precipitation accumulates from the start of the run, so each hour's total is the difference of two cumulative values.

## Set up (once)

1. Create the repo and push this folder to `main`.
2. **Settings → Pages → Build and deployment → Source: GitHub Actions.**
3. **Actions → "Ensemble winter weather (GEFS, IFS ENS, GEPS)" → Run workflow** (use `max_fhr = 48` for a quick first test).
4. Your site is at `https://<user>.github.io/<repo>/`.

The workflow has four jobs: `test`; `build`, a matrix of **three parallel jobs** (`gefs`, `ifs`, `geps`), each with its
own runner, timeout and reuse; `assemble`, which collects whatever the three produced into one site; and `deploy`.
One model being late or down does not stop the others: the site is deployed with the models that finished, and the
viewer shows only those (the Blend button appears when at least two share a run time).

The schedule is `40 4,8,10,16,20,22 * * *` UTC: the 04:40, 10:40, 16:40 and 22:40 runs pick up the GEFS cycles (about 4.5 hours
after each), and the 08:40 and 20:40 runs pick up the new IFS cycle, which ECMWF posts about 8 hours after 00/12Z (GEPS is out
after about 5.5). A run that finds a model's newest cycle not posted yet uses that model's previous cycle, which is
already published and is simply reused. Each job waits up to 100 minutes for files that are not posted yet and publishes
whatever is complete (the manifest says `complete: false` while a run is still filling in).

Rough cost of a full 10-day run: GEFS about 31 members × 80 steps; IFS 50 × 64 steps × 14 messages (about 16 GB of ranged reads from
ECMWF, ~45,000 requests held under 28/s); GEPS 72 steps × 14 whole files (about 4.5 GB). The published data is about 100 MB in
total, since the files are mostly zeros and compress well.

## What the pipeline does

| Step | Detail |
|---|---|
| Members / range (GEFS) | `c00` + `p01…p30` (31), 0.5° `pgrb2a` + `pgrb2b`. **10 days: 80 three-hourly steps, F003 to F240** (`--max-fhr` accepts up to 384; steps beyond F240 are 6-hourly) |
| Source order | AWS S3 `noaa-gefs-pds` → Google Cloud `gfs-ensemble-forecast-system` → NOMADS GRIB filter → NOMADS HTTP |
| Byte ranges | S3, GCS and NOMADS HTTP: read the `.idx`, then request only the messages needed (byte-adjacent ones merged). About 20 requests and 5 MB per member-step |
| NOMADS GRIB filter | A server-side CGI, so there is no file to byte-range. It subsets by variable, level and lat/lon box instead |
| Rate limits | Token bucket per host, job-wide totals divided across worker processes (aws 80/s, gcs 60/s, NOMADS 1.5/s with 2 connections). `Retry-After` is honored, NOMADS has a 404 budget (8/min) so missing files cannot get the runner IP blocked, and each source has a circuit breaker |
| Precip type | Wet-bulb profile → melting/refreezing energy → ProbIce, ProbSN, ProbPL, ProbFZRA (Eqs. 2, 7–9). Each member's 3-hr QPF is split among types by probability share |
| Freezing rain, sleet | Liquid-equivalent accumulations: each step's QPF × the member's freezing-rain (or sleet) share, summed. Not radial ice accretion and not pellet depth |
| Snowfall | QPF × snow share × SLR / 10 (cm). SLR is **Kuchera** by default (clamped 5–30, 12:1 baseline), or fixed with `--slr-mode fixed:12.5` / the workflow input. The files also store the snow liquid equivalent, so the viewer's **Snow-to-liquid ratio** selector (Kuchera, 13:1, 12:1, 10:1) needs no new data |
| 2 m temperature gates | Freezing rain becomes rain where 2 m temperature is above +1 °C, and ice pellets become rain above +3 °C (`FZRA_MAX_T2_C`, `PL_MAX_T2_C` in `pipeline/config.py`). 2 m temperature is already downloaded, so there is no extra cost. On one cycle (F024 and F060, 31 members) it removed 13 % of Bourgouin's freezing-rain cells; in the ASOS check it left freezing-rain detection unchanged and improved rain accuracy |
| Freezing-rain cap | GEFS publishes `CPOFP` (percent frozen precipitation). **Only where the model says at least 90 % is frozen** (`CPOFP_CAP_MIN_PCT`), a member's freezing-rain share is capped at (100 − CPOFP)% of its QPF; the excess goes to snow and sleet in proportion to their shares (all to snow if both are zero). Rain is never turned into snow. Below 90 the field did not help when checked against ASOS (see Verification), and a full cap halved freezing-rain detection. If `CPOFP` is missing for a member-hour, no cap is applied |
| Step-averaged shares | Accumulations (snow, freezing rain, sleet, snow liquid) use the mean of the type shares at the start and end of each 3-hr step, because the QPF fell over the whole step. The precip-type map keeps the end-of-step shares. Where only one end has a type, that end is used |
| Heights | Hypsometric integration of virtual temperature, so no height download is needed |

### Precipitation-type decision tree

Applied to every member at every output step, in this order, for all three models. For GEFS only fields in the 0.5° files are used (the 0.25° GEFS files have no pressure-level profiles); IFS uses its 0.25° fields sub-sampled to 0.5°.

1. **Is it precipitating?** Member 3-hr QPF ≥ 0.05 mm. Otherwise no type.
2. **Profile-based type (Revised Bourgouin).** Wet-bulb profile from the 2 m level and pressure levels → melting and refreezing energy, ice-nucleation probability → % rain, snow, freezing rain, sleet. Surface wet-bulb above 0 °C turns any freezing rain into rain (paper, Appendix c).
3. **2 m temperature gates** (the last step of the RAP/HRRR tree). Freezing rain needs a surface near freezing: above **+1 °C** its share becomes rain. Ice pellets need cold air: above **+3 °C** their share becomes rain. Snow is not gated, because the melting-energy test already handles it and only about 0.5 % of snow cells are above 3 °C. Where 2 m temperature is missing, no gate applies.
4. **Model liquid-fraction cap** (GEFS only, since only GEFS publishes `CPOFP`), only where the model says at least 90 % of the precipitation is frozen (`CPOFP` ≥ 90): freezing rain cannot exceed the liquid fraction it reports (100 − CPOFP). The excess goes to snow and sleet in proportion to their shares. Rain is never converted to snow. Below 90, or where CPOFP is missing, no cap applies.
5. **Accumulations** (snow, freezing rain, sleet, snow liquid) use the mean of the type shares at the start and end of each step. The type map shows the end-of-step shares.
6. **NCEP's own dominant-type flags are not used.** They were evaluated against ASOS (below) and Bourgouin scored slightly better for freezing rain (CSI 0.22 vs 0.20; AUC 0.77 vs 0.72) and overall (0.873 vs 0.864). Bourgouin also gives each member fractional type shares, which the QPF split and the type totals need, and it comes from the profile that is already downloaded.

Known limits: Bourgouin's supercooled-liquid branch (cloud tops warmer than about −10 °C, no warm layer) still produces freezing rain at very cold surface temperatures where snow is usually observed (model ZR below −5 °C was observed as snow about 80 % of the time). A simple cold-surface rule fixed that but cost more freezing-rain detection than it saved, so it is left out.

### Verification against ASOS present weather

Checked with archived GEFS fields (8 members, 3–24 h leads, 00Z cycles) on the 14 days with the most freezing-rain and sleet reports in the last three winters (2023-24, 2024-25, 2025-26; 36 states, 1,622 ASOS stations, METAR codes ±1 h of the valid time). Sample: station-steps where precipitation was observed, the model was wet and NCEP had a type (323,527). Freezing rain (ZR) POD = share of observed ZR the method called; FAR = share of its ZR calls with no ZR observed; CSI combines both.

| Method | Type accuracy | ZR POD | ZR FAR | ZR CSI |
|---|---|---|---|---|
| Raw Bourgouin | 0.859 | 0.50 | 0.72 | 0.22 |
| + 2 m gates (+1 °C) + full CPOFP cap (an earlier version) | 0.876 | 0.24 | 0.70 | 0.16 |
| **+ 2 m gates + cap only where CPOFP ≥ 90 (current)** | **0.873** | **0.43** | **0.69** | **0.22** |
| NCEP dominant-type flags (not used) | 0.864 | 0.33 | 0.67 | 0.20 |

* The same ordering holds in the held-out 2025-26 winter (current ZR CSI 0.22 vs NCEP 0.19; full cap 0.15). Bootstrapping over days puts the raw ZR POD at 0.43–0.57 and the full-cap POD at 0.19–0.32 (90 % intervals), so the cap's loss is not noise.
* Where Bourgouin says freezing rain, ZR was observed 30–33 % of the time for CPOFP below 90 but 17 % at 90 and above, which is why the cap is limited to that range.
* The +1 °C freezing-rain gate leaves ZR POD unchanged (0.50) and slightly improves rain accuracy: model ZR with 2 m temperature of 1–2 °C was observed as rain 85 % of the time, and above 2 °C 94 %.
* A rule turning freezing rain into snow below about −5 °C raised overall accuracy but lowered ZR detection without raising ZR CSI, so it is not used.
* Sleet is poorly forecast by every method (POD 0.1–0.5, FAR above 0.9) and sparsely reported, so the sleet numbers are too noisy to tune against.

Caveats: the days were picked for freezing-rain activity, so scores describe event days, not climatology. Many ASOS sites have no freezing-rain sensor and report rain or snow instead, which inflates every method's FAR (POD is less affected). The model grid is 0.5°, compared against point observations. The type is a snapshot at the valid time against observations within an hour.

Run it locally:

```bash
pip install -r requirements-dev.txt
python -m pytest -q tests && node --test tests/ensstats.test.js
python -m pipeline.run --out site/data --max-fhr 48 --workers 6      # add --min-members 20 --straggler-minutes 15 to tune
python -m pipeline.site --out site
cd site && python -m http.server 8000
```

### Note on Eq. 7

The freezing-rain equation reads `ProbFZRA_i = −2.1·RE + 0.2·ME + 458`. Text extraction of the paper PDF turns the
minus sign into a "2" and the decimal point into a ":", which reads as `−22.1`. The code uses −2.1, and a regression test
(`test_shallow_cold_layer_is_freezing_rain`) guards it.

## When members or data are missing

The pipeline keeps producing every product for every forecast hour that it can. Nothing is all-or-nothing.

| Situation | What happens |
|---|---|
| A file or request fails | Retried with backoff (`Retry-After` honored), then the next source in the chain (AWS → GCS → NOMADS filter → NOMADS HTTP) |
| A file is not posted yet | Polled every 45 s. This is not an error and never counts toward the 3 strikes. AWS and GCS are asked first; NOMADS (NCEP's own servers, with a strict request budget) is asked only for the hour being waited on, and only after the mirrors have lacked the file for `--origin-after-seconds` (default 120). A mirror that errors, rather than saying "not posted", skips that wait. When NOMADS' 404 budget is used up the pipeline simply does not ask, which is also not a failure and does not trip the circuit breaker |
| A member is still missing after the other members of that hour arrived (`--straggler-minutes`, default 15) | **Dropped from that hour onward**; the run continues with the rest |
| A member fails 3 times with a real error (corrupt GRIB, a missing field) | Dropped the same way, instead of ending the run |
| Nothing arrives for `--wait-minutes` at one hour | If at least `--min-members` (default 20) have posted, the stragglers are dropped. Otherwise the run stops there and publishes the hours it has (`complete: false`) |
| Fewer than `--min-members` would remain | The run stops there and publishes the hours it has |

A dropped member is never skipped for one hour only, because GEFS accumulations are chained per member and a gap would break its totals. Each file's header holds a member bitmask and the manifest records `n`, `mask` and `dropped` for every step, so the viewer
works on exactly the members that exist for the selected window (a window needs the member in both its end and start files), and shows a warning with the member names. Statistics, probabilities, PMM/LPMM and the single-member list all use the available members. The decision logic is the pure function `step_decision()` in `pipeline/run.py`, covered by `tests/test_policy.py` (scenario tests with a fake model: member missing from one hour on, repeated errors, a late result that must not revive a dropped member, too few members, nothing at all).

## Not downloading the same cycle twice

Every deploy builds a fresh site, so the only copy of a finished cycle is the one on the live Pages site. At the start of a run the pipeline reads `data/<model>/latest.json` (GEFS also `synoptic.json`) from that site (`--reuse-from`, set by the workflow) and, if it is the **same cycle** made by the **same calculation code** (`code_rev` in the manifest) with the same members, grid and snow-ratio mode, it:

* copies the published hours back from your own site (never from NOAA) and checks every file;
* if the cycle is complete, stops there with no NOAA requests at all, so a push that only changes the viewer costs seconds;
* if the cycle is partial, restores the running accumulators from the published files and downloads **only the missing hours**. Members dropped earlier stay dropped. Restored totals carry at most one file-quantization step of rounding (0.05 mm QPF, 0.1 cm snow);
* falls back to a normal download when the cycle, code, members or snow mode differ, when no site exists yet, or from the first hour whose file is missing or damaged.

Run the workflow by hand with `force_refresh` = true (or pass `--force`) to download everything again, for example to pick up members that were dropped on the first attempt. `code_rev` covers `bourgouin.py`, `config.py`, `encode.py`, `families.py`, `geps.py`, `grib.py`, `ifs.py`, `member.py`, `products.py` and `run.py`; changes to the source/transport code or to the viewer do not invalidate published data.

## Data format (`data/<model>/<yyyymmdd_cc>/…`, `data/<model>/latest.json`)

`<model>` is `gefs`, `ifs` or `geps`. `latest.json` lists the run, `family`, `members`, pressure `levels`, domain, scales and every step (GEFS also has `synoptic.json`, the 00/12Z cycle that matches IFS/GEPS). Each step has two files:

* `cum/fHHH.bin.gz`: **cumulative-since-init** `u16 [member][y][x]` arrays, in order: QPF (0.05 mm per count), snowfall (0.1 cm), freezing rain (0.01 mm liquid equivalent), sleet (0.01 mm liquid equivalent), snow liquid equivalent (0.02 mm, the fifth array). The header's `narr` byte gives the array count, so older 2- and 4-array files still read; the fifth scale sits at header offset 52.
* `ptype/fHHH.bin.gz`: `u8 [4][member][y][x]` = % RA, SN, FZRA, PL (0 where the member's 3-hr QPF is below 0.05 mm).

Cumulative storage means any rolling 3/6/12/24/48-hour total is just `cum(t) − cum(t−w)`, two small files fetched
by the browser. A full GEFS run is about 25 MB (IFS about 40 MB, GEPS about 20 MB). Domain: 60–20°N, 135–50°W, 0.5° (81 × 171). The 64-byte `GWX1` header is
documented in `pipeline/encode.py`; its member bitmask spans 64 members (bits 0–31 at offset 48, 32–63 at offset 56) because IFS has 50.

## Viewer

* MapLibre GL JS 4.7.1 (vendored), a self-hosted Natural Earth + US Census basemap (countries, states/provinces, counties, rivers, lakes), a graticule, and a canvas label layer with collision handling. There is no tile or font server.
* A custom WebGL2 layer paints the data with crisp class edges at any zoom, using the NWS-style ramps from the three reference images.
* Hover (desktop) or tap (phone) for a tooltip with mean, median, min, max, 10th/90th percentile, chance of exceeding thresholds, a dot strip of all members (coloured by model in the blend, with each model's mean) and, for type, member counts and mean Bourgouin probabilities. Click to pin.
* A **Model** selector (Blend, GEFS, IFS ENS, GEPS) and, for the blend, a weighting selector; both are kept in the URL (`m=`, `w=`).
* Statistics (computed in the browser from the member totals of the selected window; weighted in the blend as described above):
  * **Min, P5, P10, P25, P50, Mean, P75, P90, P95, Max.** Percentiles interpolate linearly between ranked members (31 members, so P5 sits between the 2nd and 3rd lowest).
  * **PMM**, the probability-matched mean (Ebert 2001). The ensemble mean supplies the spatial pattern. All member values pooled over the whole map are ranked, and every 31st value is dealt back onto the mean field by rank. Amplitudes are no longer smoothed away, and the map-wide total stays close to the mean's.
  * **LPMM**, the local probability-matched mean (Clark 2017; Snook et al. 2019/2020). The same matching is done independently in a square neighborhood around each grid cell, which avoids mixing very different regimes (a coastal low and a dry interior). The neighborhood is selectable (±1° to ±4°; default ±2°, about 140 miles). Published LPM settings are tuned for 3-km convection-allowing ensembles, so treat the default as a starting point for 0.5° GEFS. Cells whose ensemble mean is zero stay zero.
  * Precip type has two views: the most likely type (QPF-weighted) and the per-type probability (% of members). The tooltip lists member counts for each type and the mean Bourgouin probabilities.
  * Snowfall has a ratio selector (Kuchera, 13:1, 12:1, 10:1) that rescales from the stored snow liquid equivalent and is kept in the URL (`slr=`). Files from before this change lack that array, and the selector stays hidden.
  * Probability of at least *x*, and any single member. Everything in the viewer is shown in imperial units: inches for QPF, snowfall, freezing rain and sleet, and miles/feet on the scale bar. (The data files keep mm/cm internally; conversion happens in the browser.)
  * The PMM/LPMM math is in `ensstats.js` and tested with `node --test tests/ensstats.test.js` (distribution matching, rank order, equality with PMM for a map-sized window, brute-force comparison).
* Speed: QPF, snowfall, freezing rain and sleet share the same cumulative files, so switching among them reuses data already in memory (about 15 ms). After each view loads, the browser quietly prefetches, one file at a time at idle priority, the precip-type file, the start files for the other accumulation periods, and the next three hours in whichever direction you are stepping. It skips this when Data Saver or a 2G connection is on. PMM and LPMM run in a Web Worker (`statworker.js`) so the map stays responsive, and finished fields are memoized. Decoded files are held under a memory budget (90 MB on low-memory phones, 220 MB otherwise).
* State is kept in the URL hash, so links are shareable.

## Layout

```
pipeline/        data pipeline (python -m pipeline.run --family gefs|ifs|geps, pipeline.site); pipeline/reuse.py skips hours already published
                 families.py (the three models), ifs.py and geps.py (their fetchers), run.py (shared orchestration)
tests/           pytest suite, no network needed
index.html style.css viewer.js basemap.js wxlayer.js colors.js ensstats.js statworker.js   viewer
basemap_*.json   boundaries and labels (build_basemap.py regenerates them)
maplibre-gl.*    vendored MapLibre
.github/workflows/gefs.yml
```

## Caveats

* The Google Cloud source could not be reached from the development sandbox, so its fallback is covered only by the mock-server tests. The AWS and NOMADS filter paths were run against real data.
* The fields are 0.5°, so precipitation type is smooth and not suitable for sub-county detail.
* IFS ENS open data has no control member, so the IFS ensemble is 50 perturbed members. GEPS was only checked against the one live cycle available during development; its datamart keeps files for a limited time, so a run started long after the cycle may find them gone.
* The ASOS verification above covered GEFS only; IFS and GEPS are processed the same way but have not been verified separately.
* A model that fails to build for a cycle is missing from that deploy (the site is rebuilt from the job outputs each time), until its next successful run.
* Basemap: Natural Earth (public domain), US Census cartographic boundary data via us-atlas. Ensemble data: NOAA/NCEP GEFS (public domain), ECMWF IFS ENS open data (CC BY 4.0, attribution required), ECCC GEPS (Open Government Licence – Canada).
