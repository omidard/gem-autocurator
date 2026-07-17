#!/usr/bin/env python3
"""Build the ID-mapping backbone for the GEM autocurator.

Produces two lookup tables the client uses to canonicalise any model's metabolite
and reaction identifiers:
  data/metabolite_map.json : "<ns>:<id>" -> {bigg, biggr(bool), name?}   (ns: bigg,mnx,seed,kegg,chebi,inchikey,inchikey14,name)
  data/reaction_map.json   : "<ns>:<id>" -> {bigg, biggr(bool)}          (ns: bigg,mnx,kegg,seed,rhea,old)

Strategy: BiGG is the curated target namespace. Every metabolite/reaction that BiGG
has (via its own xrefs) maps to its BiGG id (biggr=True). Compounds/reactions that
other major DBs (KEGG, ModelSEED/SEED, MetaNetX, ChEBI, RHEA) have but BiGG does NOT
get a deterministic **BiGG-like** id (biggr=False) so downstream models can still be
canonicalised consistently.
"""
import json, re, os, csv, sys

OUT = "/data/gem_autocurator/docs/data"
os.makedirs(OUT, exist_ok=True)
BIGG_MET = json.load(open("/data/media_curate/tools/bigg_metabolite_dict.json"))
XF = json.load(open("/data/media_curate/tools/xref_fallback.json"))
BIGG_RXN_FILE = "/data/reactome/ESKAPE_AMR/kegg/bigg_models_reactions.txt"
REAC_XREF = "/data/bioconversion/thermo/reac_xref.tsv"

def norm_name(s):
    return re.sub(r"[^a-z0-9]", "", (s or "").lower())

# ---------------- METABOLITES ----------------
met_map = {}
def mput(ns, idv, bigg, biggr, name=None):
    if not idv:
        return
    k = ns + ":" + str(idv)
    if k not in met_map:
        met_map[k] = {"bigg": bigg, "biggr": biggr}
        if name:
            met_map[k]["name"] = name[:60]

# 1) everything BiGG has -> its BiGG id (biggr=True)
mnx2bigg = {}
for b, rec in BIGG_MET.items():
    x = rec.get("xrefs", {})
    nm = (rec.get("name") or b)
    mput("bigg", b, b, True, nm)
    if x.get("mnx"):
        mnx2bigg[x["mnx"]] = b
        mput("mnx", x["mnx"], b, True, nm)
    mput("seed", x.get("seed"), b, True, nm)
    mput("kegg", x.get("kegg"), b, True, nm)
    if x.get("chebi"):
        mput("chebi", str(x["chebi"]), b, True, nm)
    ik = x.get("inchikey")
    if ik:
        mput("inchikey", ik, b, True, nm)
        mput("inchikey14", ik.split("-")[0], b, True, nm)
    for syn in rec.get("synonyms", [])[:2]:
        mput("name", norm_name(syn), b, True, nm)

# 2) BiGG-like ids for compounds other DBs have but BiGG does not (via MetaNetX MNX hubs)
#    BiGG-like id = "x_<seed|kegg|mnx id>" (deterministic, unique, namespaced).
used_like = set()
def biggify(prefix, srcid):
    base = "x_" + re.sub(r"[^a-z0-9]", "", (prefix + srcid).lower())
    return base

# KEGG compounds not in BiGG
n_like = 0
for kegg, mnx in XF["kegg2mnx"].items():
    if mnx in mnx2bigg:
        continue                      # BiGG already covers it
    like = biggify("kc", kegg)
    seed = XF["mnx2seed"].get(mnx)
    for ns, idv in (("kegg", kegg), ("mnx", mnx), ("seed", seed)):
        if idv and (ns + ":" + idv) not in met_map:
            met_map[ns + ":" + idv] = {"bigg": like, "biggr": False}
    n_like += 1

n_bigg_targets = len({v["bigg"] for v in met_map.values() if v["biggr"]})
print("METABOLITES: map keys", len(met_map), "| BiGG targets", n_bigg_targets, "| BiGG-like from KEGG", n_like)
json.dump(met_map, open(OUT + "/metabolite_map.json", "w"))

# ---------------- REACTIONS ----------------
rxn_map = {}
def rput(ns, idv, bigg, biggr):
    if not idv:
        return
    k = ns + ":" + str(idv)
    if k not in rxn_map:
        rxn_map[k] = {"bigg": bigg, "biggr": biggr}

bigg_rxn_ids = set()
with open(BIGG_RXN_FILE) as fh:
    r = csv.DictReader(fh, delimiter="\t")
    for row in r:
        b = row["bigg_id"]
        bigg_rxn_ids.add(b)
        rput("bigg", b, b, True)
        for old in (row.get("old_bigg_ids") or "").split(";"):
            rput("old", old.strip(), b, True)
        # database_links: extract RHEA / MetaNetX ids
        dbl = row.get("database_links") or ""
        for rid in re.findall(r"rhea/(\d+)", dbl):
            rput("rhea", rid, b, True)
        for mnx in re.findall(r"metanetx\.reaction/(MNXR\d+)", dbl):
            rput("mnx", mnx, b, True)
print("BiGG reactions:", len(bigg_rxn_ids))

# MetaNetX reac_xref: source_id \t mnx \t desc  -> link kegg/seed reactions to bigg via mnx
mnxr2bigg = {k.split(":", 1)[1]: v["bigg"] for k, v in rxn_map.items() if k.startswith("mnx:")}
# also build bigg.reaction:X -> its mnx (from reac_xref where source is bigg.reaction)
n_r_bigg = n_r_like = 0
try:
    with open(REAC_XREF) as fh:
        # first pass: bigg.reaction:X -> mnx (fills mnxr2bigg for bigg rxns not caught above)
        bigg_src2mnx = {}
        rows = []
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            src, mnx = parts[0], parts[1]
            rows.append((src, mnx))
            if src.startswith("bigg.reaction:") and mnx.startswith("MNXR"):
                bb = src.split(":", 1)[1]
                if bb in bigg_rxn_ids:
                    mnxr2bigg.setdefault(mnx, bb)
        for src, mnx in rows:
            if ":" not in src or not mnx.startswith("MNXR"):
                continue
            ns, idv = src.split(":", 1)
            nsk = {"kegg.reaction": "kegg", "seed.reaction": "seed", "metacyc.reaction": "biocyc",
                   "rhea": "rhea", "bigg.reaction": "bigg"}.get(ns)
            if not nsk:
                continue
            if mnx in mnxr2bigg:
                rput(nsk, idv, mnxr2bigg[mnx], True); n_r_bigg += 1
            else:
                like = "xr_" + re.sub(r"[^a-z0-9]", "", idv.lower())
                rput(nsk, idv, like, False)
                rput("mnx", mnx, like, False)
                n_r_like += 1
except FileNotFoundError:
    print("reac_xref not found; reactions limited to BiGG universe")
print("REACTIONS: map keys", len(rxn_map), "| linked-to-BiGG", n_r_bigg, "| BiGG-like", n_r_like)
json.dump(rxn_map, open(OUT + "/reaction_map.json", "w"))

# a compact metabolite property table (formula + charge) for mass/charge QC
prop = {}
for b, rec in BIGG_MET.items():
    x = rec.get("xrefs", {})
    prop[b] = {"f": x.get("formula") or rec.get("formula") or "", "c": rec.get("charge"),
               "n": (rec.get("name") or b)[:50]}
json.dump(prop, open(OUT + "/bigg_met_props.json", "w"))
print("met props:", len(prop))
print("sizes:", {f: os.path.getsize(OUT + "/" + f) for f in ("metabolite_map.json", "reaction_map.json", "bigg_met_props.json")})
