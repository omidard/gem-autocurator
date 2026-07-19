#!/usr/bin/env python3
"""Build the COMPREHENSIVE, structure-clustered identifier backbone for the autocurator.

MetaNetX (MNXref) already clusters every compound/reaction from BiGG, KEGG, ModelSEED,
ChEBI, HMDB, MetaCyc, Rhea, SABIO-RK … by structure (one MNXM/MNXR per InChIKey cluster).
We ingest that, keep the GEM-relevant clusters (those carrying a BiGG/KEGG/ModelSEED/
MetaCyc/Rhea id), and give EACH cluster ONE canonical id:
  * its BiGG id if the cluster has one           (biggr = True)
  * otherwise ONE stable BiGG-like id, shared by ALL of the cluster's database ids and
    FROZEN in tools/backbone_freeze.json so it never changes between builds.
Then every database id (and InChIKey / synonym) in the cluster maps to that one id — so a
model using glc__D, C00031, cpd00027 or CHEBI:4167 all resolve to the same standard id and
nobody can introduce a discrepancy with their own id logic.

Inputs (build-time, /data/mnx_tmp, from https://www.metanetx.org/ftp/4.4/):
  chem_prop.tsv  chem_xref.tsv  reac_prop.tsv  reac_xref.tsv
Outputs (docs/data/, gzipped — the app fetches + gunzips them):
  metabolite_map.json.gz  reaction_map.json.gz  bigg_met_props.json.gz  backbone_coverage.json
"""
import json, os, re, gzip, hashlib

MNX = "/data/mnx_tmp"
OUT = "/data/gem_autocurator/docs/data"
FREEZE = "/data/gem_autocurator/tools/backbone_freeze.json"

# ---- namespace normalisation for chem_xref / reac_xref source prefixes ----
MET_NS = {"bigg.metabolite": "bigg", "biggm": "old", "kegg.compound": "kegg", "keggc": "kegg",
          "seed.compound": "seed", "seedm": "seed", "chebi": "chebi", "metacyc.compound": "biocyc",
          "metacycm": "biocyc", "hmdb": "hmdb", "mnx": "mnx"}
RXN_NS = {"bigg.reaction": "bigg", "biggr": "old", "kegg.reaction": "kegg", "keggr": "kegg",
          "seed.reaction": "seed", "seedr": "seed", "rhea": "rhea", "rhear": "rhea",
          "metacyc.reaction": "biocyc", "metacycr": "biocyc", "mnx": "mnx"}
MET_METABOLIC = {"bigg", "old", "kegg", "seed", "biocyc"}
RXN_METABOLIC = {"bigg", "old", "kegg", "seed", "rhea", "biocyc"}

GREEK = {"alpha": "a", "beta": "b", "gamma": "g", "delta": "d", "epsilon": "e", "omega": "o"}
STOP = {"acid", "ion", "the", "of", "and", "a", "an"}

def slugify(name):
    s = (name or "").lower()
    s = re.sub(r"\([rs+\-]\)-?|(^|[^a-z])[dl]-", " ", s)
    for g, v in GREEK.items():
        s = s.replace(g, v)
    words = [w for w in re.split(r"[^a-z0-9]+", s) if w and w not in STOP]
    joined = "".join(words)
    if len(joined) <= 10:
        base = joined
    elif len(words) >= 2:
        base = "".join(w if len(w) <= 4 else w[0] + re.sub(r"[aeiou]", "", w[1:])[:3] for w in words)[:12]
    else:
        base = re.sub(r"[aeiou]", "", joined)[:8] or joined[:8]
    return re.sub(r"[^a-z0-9]", "", base)

def h4(s):
    return hashlib.md5(s.encode()).hexdigest()[:4]

def load_freeze():
    try:
        return json.load(open(FREEZE))
    except Exception:
        return {"met": {}, "rxn": {}}

def parse_prop_names(path, wanted):
    """MNXM/MNXR -> (name, formula, charge, inchikey14) for wanted ids (chem_prop)."""
    props = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 2 or p[0] not in wanted:
                continue
            name = p[1] if len(p) > 1 else ""
            formula = p[3] if len(p) > 3 else ""
            charge = p[4] if len(p) > 4 else ""
            ik = p[7].replace("InChIKey=", "") if len(p) > 7 and p[7] else ""
            props[p[0]] = (name, formula, charge, ik)
    return props

def collect(xref_path, nsmap, metabolic_ns):
    """Two passes over an xref file -> {mnx: {ns: [ids], 'names': set}} for metabolic clusters."""
    metabolic = set()
    with open(xref_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 2 or ":" not in p[0]:
                continue
            ns = nsmap.get(p[0].split(":", 1)[0].lower())
            if ns in metabolic_ns:
                metabolic.add(p[1])
    clusters = {}
    with open(xref_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 2 or p[1] not in metabolic:
                continue
            mnx = p[1]
            c = clusters.setdefault(mnx, {})
            if ":" not in p[0]:
                continue
            db, ext = p[0].split(":", 1)
            ns = nsmap.get(db.lower())
            if not ns:
                continue
            desc = p[2] if len(p) > 2 else ""
            if ns == "bigg":
                prim = "secondary/obsolete/fantasy" not in desc
                c.setdefault("bigg", [])
                if not any(e == ext for e, _ in c["bigg"]):
                    c["bigg"].append((ext, prim))
                if prim and desc:
                    c.setdefault("_name", desc.split("||")[0])
            else:
                c.setdefault(ns, [])
                if ext not in c[ns]:
                    c[ns].append(ext)
    return metabolic, clusters

def clean_bigg(b):
    return re.sub(r"^[RM]_", "", b)

def build(kind, xref_path, prop_path, nsmap, metabolic_ns, freeze, out_name):
    metabolic, clusters = collect(xref_path, nsmap, metabolic_ns)
    props = parse_prop_names(prop_path, metabolic)
    frz = freeze[kind]
    # seed taken with the frozen bigg-like ids AND every real BiGG id, so a bigg-like never collides with a BiGG id
    taken = set(frz.values())
    for c in clusters.values():
        for b, _p in c.get("bigg", []):
            taken.add(clean_bigg(b))
    idmap = {}
    n_bigg = n_like = 0
    outprops = {}
    def put(ns, ext, canon, biggr):
        if not ext:
            return
        k = ns + ":" + ext
        if k not in idmap:
            idmap[k] = {"bigg": canon, "biggr": 1 if biggr else 0}
    for mnx in sorted(clusters):
        c = clusters[mnx]
        name, formula, charge, ik = props.get(mnx, ("", "", "", ""))
        name = c.get("_name") or name
        bigg_entries = c.get("bigg", [])
        primset = {clean_bigg(b) for b, p in bigg_entries if p}
        allbigg = {clean_bigg(b) for b, _ in bigg_entries}
        covered = bool(bigg_entries)
        if covered:
            canon = min(primset or allbigg, key=lambda x: (len(x), x))   # a real, primary BiGG id
            n_bigg += 1
        else:
            n_like += 1
            if mnx in frz:
                canon = frz[mnx]
            else:
                if kind == "met":
                    base = slugify(name) or ("cpd" + h4(mnx))
                else:
                    base = re.sub(r"[^a-z0-9]", "", ("kr" + c["kegg"][0].lower()) if c.get("kegg") else ("xr" + h4(mnx)))
                canon = base if base not in taken else base + h4(mnx)
                taken.add(canon)
                frz[mnx] = canon
        # bigg ids: a PRIMARY id maps to itself (it is a valid BiGG id); secondary/obsolete -> the canonical primary
        for b, p in bigg_entries:
            cb = clean_bigg(b)
            target = cb if cb in primset else canon
            put("bigg", b, target, True)
        # every other database id in the cluster -> the one canonical id
        for ns, vals in c.items():
            if ns in ("bigg", "_name"):
                continue
            for ext in vals:
                put(ns, ext, canon, covered)
        put("mnx", mnx, canon, covered)
        if ik:
            put("inchikey", ik, canon, covered)
            put("inchikey14", ik.split("-")[0], canon, covered)
        if kind == "met" and name:
            put("name", re.sub(r"[^a-z0-9]", "", name.lower()), canon, covered)
        outprops[canon] = {"f": formula, "c": (None if charge == "" else _num(charge)), "n": name[:50]}
    _dump_gz(idmap, os.path.join(OUT, out_name))
    print(f"{kind}: {len(idmap)} keys | {n_bigg} BiGG-covered clusters + {n_like} BiGG-like ({len(metabolic)} clusters)")
    return idmap, outprops, {"clusters": len(metabolic), "bigg": n_bigg, "bigglike": n_like, "keys": len(idmap)}

def _num(s):
    try:
        return int(s)
    except ValueError:
        try:
            return float(s)
        except ValueError:
            return None

def _dump_gz(obj, path):
    with gzip.open(path + ".gz", "wt", encoding="utf-8") as fh:
        json.dump(obj, fh, separators=(",", ":"))

freeze = load_freeze()
print("METABOLITES …")
mmap, mprops, mcov = build("met", MNX + "/chem_xref.tsv", MNX + "/chem_prop.tsv", MET_NS, MET_METABOLIC, freeze, "metabolite_map.json")
print("REACTIONS …")
rmap, rprops, rcov = build("rxn", MNX + "/reac_xref.tsv", MNX + "/reac_prop.tsv", RXN_NS, RXN_METABOLIC, freeze, "reaction_map.json")
_dump_gz(mprops, os.path.join(OUT, "bigg_met_props.json"))
json.dump(freeze, open(FREEZE, "w"), separators=(",", ":"))
json.dump({"metabolites": mcov, "reactions": rcov}, open(os.path.join(OUT, "backbone_coverage.json"), "w"), indent=1)
print("coverage:", json.dumps({"met": mcov, "rxn": rcov}))
for f in ("metabolite_map.json.gz", "reaction_map.json.gz", "bigg_met_props.json.gz"):
    print(" ", f, os.path.getsize(os.path.join(OUT, f)), "bytes")
