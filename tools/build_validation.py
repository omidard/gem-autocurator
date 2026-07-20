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
    m = re.search(r"substr\.?\s+([A-Za-z0-9][A-Za-z0-9\-]{1,})", t, re.I) or re.search(r"(?:str\.?|strain)\s+([A-Za-z0-9][A-Za-z0-9\-]{1,})", t)
    if m:
        tok = m.group(1)
        if not re.match(r"^(wild|type|unknown|unspecified|isolate|sp|strain|and|the|not|clinical|reference|derivative|derivatives|mutant|parent|parental)$", tok, re.I) and (re.search(r"\d", tok) or re.match(r"^[A-Z]", tok)):
            return re.sub(r"[^A-Z0-9]", "", tok.upper())
    m = re.search(r"\b([A-Z]{1,4}\d{1,6}[A-Za-z]?)\b", t)
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
        # fu = flux_usable: the rate is biomass-specific (mmol/gDW/h) and can be a real FBA bound
        out.append({"ex": ex, "met": x.get("bigg_metabolite"), "r": round(float(r), 4),
                    "u": x.get("units") or "", "fu": bool(x.get("flux_usable"))})
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
        "mu_ok": r.get("mu_usable"), "mu_qc": r.get("mu_qc"),
        "dt": r.get("doubling_time_h"),
        "up": up, "sec": sec,
        "med": {"id": mid, "key": med.get("canonical_key"), "name": med.get("canonical_name") or med.get("description"),
                # exchanges formulated from the paper's own recipe (when the medium isn't a linked Media DB id)
                "ex": [[e["exchange"], e["lb"]] for e in (med.get("exchanges") or [])] or None,
                "fmt": med.get("formulation"),
                # medium_type: defined media give a real FBA validation; in_vivo/environmental/complex_undefined
                # cannot yield an exact exchange set, so the autocurator must say so rather than fake it
                "mt": med.get("medium_type"), "formulable": med.get("formulable")},
        "cond": {"o2": cond.get("oxygen"), "t": cond.get("temperature_C"), "pH": cond.get("pH"),
                 "mode": cond.get("culture_mode"), "D": cond.get("dilution_rate_per_h")},
        "cit": (prov.get("citation") or "")[:180], "doi": prov.get("doi"),
    }
    records.append(rec)
    species[sp] = species.get(sp, 0) + 1

species_list = sorted(species.keys())

# GrowthDB-fitted maintenance (NGAM) + biomass yield (Yxs), from qS-vs-µ Pirt fits (r²>=0.6 only)
maint = {}
SPIDX = "/data/modelseed_cache/gdb_species_index.json"
if os.path.exists(SPIDX):
    sidx = json.load(open(SPIDX))
    rows = sidx["species"] if isinstance(sidx, dict) and "species" in sidx else (sidx if isinstance(sidx, list) else list(sidx.values()))
    for r in rows:
        if isinstance(r, dict) and (r.get("ngam") or r.get("yxs")):
            maint[r["s"]] = {"ngam": r.get("ngam"), "yxs": r.get("yxs")}

# provenance/version stamp so the modeller knows how current the validation bundle is
import datetime, subprocess
def _git_rev(repo):
    try:
        return subprocess.check_output(["git", "-C", repo, "rev-parse", "--short", "HEAD"], stderr=subprocess.DEVNULL).decode().strip()
    except Exception:
        return None
meta = {"built": datetime.date.today().isoformat(),
        "n_records": len(records), "n_species": len(species_list),
        "n_mu": sum(1 for r in records if r["mu"] is not None),
        "growthdb_rev": _git_rev("/data/GrowthDB_work") or _git_rev("/data/GrowthDB"),
        "media_rev": _git_rev("/data/media_curate")}

json.dump({"meta": meta, "species": species_list, "records": records, "tax2sp": tax2sp, "maint": maint}, open(OUT + "/growthdb.json", "w"), separators=(",", ":"))
print("bundle stamp:", meta)
print("maintenance fits (NGAM/Yxs): %d species" % len(maint))

# ---- substrate-utilisation SPECTRUM per species (for grows-on-X confusion-matrix validation) ----
# positive-precedence: an exchange with ANY growth evidence (utilizable / carbon_growth / growth-positive
# phenotype) is 'grows'; negatives are growth-negative phenotypes not seen positive anywhere.
SPECIES_DIR = "/data/GrowthDB_work/data/species"
spectrum = {}
if os.path.isdir(SPECIES_DIR):
    for f in glob.glob(SPECIES_DIR + "/*.json"):
        try:
            d = json.load(open(f))
        except Exception:
            continue
        sp = d.get("species")
        if not sp:
            continue
        dv = d.get("derived") or {}
        pos, neg = set(), set()
        by_strain = {}
        for ex in (dv.get("utilizable_exchanges") or []):
            pos.add(ex)
        for r in (d.get("carbon_growth") or []):
            ex = r.get("exchange")
            if ex:
                pos.add(ex)
                st = r.get("strain")
                if st:
                    by_strain.setdefault(st, set()).add(ex)
        for p in (d.get("phenotypes") or []):
            if p.get("ptype") in ("growth", "growth+acid") and p.get("exchange"):
                (pos if p.get("phenotype") == "positive" else neg).add(p["exchange"])
        neg -= pos                                    # positive evidence wins a conflict
        if pos or neg:
            spectrum[sp] = {"p": sorted(pos), "n": sorted(neg),
                            "s": {st: sorted(ex) for st, ex in by_strain.items() if ex}}
with open(OUT + "/spectrum.json", "w") as fh:
    json.dump(spectrum, fh, separators=(",", ":"))
n_pos = sum(len(v["p"]) for v in spectrum.values())
n_neg = sum(len(v["n"]) for v in spectrum.values())
print("spectrum: %d species | %d grows-on + %d no-grow exchanges | %d bytes"
      % (len(spectrum), n_pos, n_neg, os.path.getsize(OUT + "/spectrum.json")))

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
