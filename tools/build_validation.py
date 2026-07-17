#!/usr/bin/env python3
"""Build the bundled lab-validation dataset for the autocurator, from GrowthDB
(experimental growth/uptake/secretion rates) + the Media DB (BiGG-exchange media).

Inputs (build-time, kept outside the repo in /data/modelseed_cache alongside seed tsvs):
  curl -o /data/modelseed_cache/gdb_growth_records.json \\
    https://raw.githubusercontent.com/omidard/GrowthDB/main/data/growth_records.json
  Media exchange compositions come from the local Media repo working copy at
  /data/media_curate/data/media/*.json (per-medium files with BiGG exchange bounds).

Outputs:
  docs/data/growthdb.json : {species:[...], records:[{sp,org,mu,dt,up,sec,med,cond,cit,doi}]}
  docs/data/media_ex.json : {media_id: {name, aerobic, ex:[[exchange, lb, ub], ...]}}
      only the media that GrowthDB records reference (674 of 13k).
"""
import json, os, glob, re

GDB = "/data/modelseed_cache/gdb_growth_records.json"
MEDIA_DIR = "/data/media_curate/data/media"
OUT = "/data/gem_autocurator/docs/data"

# ---- standard strain id extraction (mirrored in autocurator.js: strainStd) ----
_CC = r"ATCC|DSMZ|DSM|NCTC|CCUG|JCM|NBRC|IFO|CECT|LMG|CIP|NCIMB|BCRC|KCTC|NRRL|CGMCC|MCCC|PCC|UTEX|CBS|KACC|VPI|NCDO|NCFB|BCCM|KCCM"

def strain_std(text):
    """Normalise a free-text strain string to a comparable standard token, or None.
    Priority: culture-collection accession > str./substr. designation > lab designation > K-12."""
    if not text:
        return None
    t = str(text)
    m = re.search(r"\b(" + _CC + r")\s*[-: ]?\s*(\d+[A-Za-z]?)\b", t, re.I)
    if m:
        return (m.group(1) + m.group(2)).upper()
    m = re.search(r"substr\.?\s+([A-Za-z0-9][A-Za-z0-9\-]{1,})", t, re.I) or re.search(r"(?:str\.?|strain)\s+([A-Za-z0-9][A-Za-z0-9\-]{1,})", t, re.I)
    if m and re.search(r"\d", m.group(1)):
        return re.sub(r"[^A-Z0-9]", "", m.group(1).upper())
    m = re.search(r"\b([A-Z]{1,4}\d{2,6}[A-Za-z]?)\b", t)
    if m:
        return m.group(1).upper()
    if re.search(r"\bK-?12\b", t, re.I):
        return "K12"
    return None

def strain_display(text):
    if not text:
        return None
    # first clause before an explanatory parenthesis/comma, trimmed
    s = re.split(r"\s*[(,]", str(text).strip())[0].strip()
    return (s or str(text).strip())[:48]

gr = json.load(open(GDB))

def trim_rates(lst):
    out = []
    for x in (lst or []):
        ex = x.get("exchange"); r = x.get("rate")
        if not ex or r is None:
            continue
        out.append({"ex": ex, "met": x.get("bigg_metabolite"), "r": round(float(r), 4), "u": x.get("units") or ""})
    return out

records, species, referenced, tax2sp = [], {}, set(), {}
for r in gr:
    mu = r.get("growth_rate_per_h")
    up = trim_rates(r.get("uptake_rates"))
    sec = trim_rates(r.get("secretion_rates"))
    if mu is None and not up and not sec:
        continue                                   # nothing to validate against
    med = r.get("medium") or {}
    mid = med.get("media_id")
    if mid:
        referenced.add(mid)
    tax = r.get("ncbi_tax_id")
    sp_key = r.get("gtdb_species") or r.get("species") or r.get("organism")
    if tax and sp_key and str(tax) not in tax2sp:
        tax2sp[str(tax)] = sp_key
    cond = r.get("conditions") or {}
    prov = r.get("provenance") or {}
    sp = r.get("gtdb_species") or r.get("species") or r.get("organism") or "unknown"
    strain_raw = r.get("strain")
    # fall back to a strain designation embedded in the organism name (e.g. "... str. K-12 substr. MG1655")
    sstd = strain_std(strain_raw) or strain_std(r.get("organism"))
    rec = {
        "sp": sp, "org": r.get("organism"),
        "strain": strain_display(strain_raw), "sstd": sstd,
        "mu": round(float(mu), 4) if mu is not None else None,
        "dt": r.get("doubling_time_h"),
        "up": up, "sec": sec,
        "med": {"id": mid, "key": med.get("canonical_key"), "name": med.get("canonical_name") or med.get("description")},
        "cond": {"o2": cond.get("oxygen"), "t": cond.get("temperature_C"), "pH": cond.get("pH"),
                 "mode": cond.get("culture_mode"), "D": cond.get("dilution_rate_per_h")},
        "cit": (prov.get("citation") or "")[:180], "doi": prov.get("doi"),
    }
    records.append(rec)
    species[sp] = species.get(sp, 0) + 1

species_list = sorted(species.keys())
json.dump({"species": species_list, "records": records, "tax2sp": tax2sp}, open(OUT + "/growthdb.json", "w"), separators=(",", ":"))

# media exchange compositions for referenced media
media_ex = {}
files = {os.path.basename(f)[:-5]: f for f in glob.glob(MEDIA_DIR + "/*.json")}
for mid in referenced:
    f = files.get(mid)
    if not f:
        continue
    try:
        d = json.load(open(f))
    except Exception:
        continue
    ex = []
    for comp in d.get("components", []):
        e = comp.get("exchange")
        if not e:
            continue
        lb = comp.get("lower_bound"); ub = comp.get("upper_bound")
        ex.append([e, -1000.0 if lb is None else float(lb), 1000.0 if ub is None else float(ub)])
    if ex:
        media_ex[mid] = {"name": d.get("name"), "aerobic": d.get("aerobic"), "ex": ex}
json.dump(media_ex, open(OUT + "/media_ex.json", "w"), separators=(",", ":"))

n_mu = sum(1 for r in records if r["mu"] is not None)
n_up = sum(1 for r in records if r["up"])
n_sec = sum(1 for r in records if r["sec"])
print("growthdb.json: %d records (%d with μ, %d uptake, %d secretion), %d species" % (len(records), n_mu, n_up, n_sec, len(species_list)))
print("media_ex.json: %d media compositions" % len(media_ex))
for f in ("growthdb.json", "media_ex.json"):
    print(" ", f, os.path.getsize(OUT + "/" + f), "bytes")
