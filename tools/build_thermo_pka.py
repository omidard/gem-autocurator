#!/usr/bin/env python3
"""Build bundled pKa + thermodynamics tables for the autocurator, from the
ModelSEED biochemistry database (Biochemistry/compounds.tsv, reactions.tsv).

Inputs (build-time only, not tracked in the repo — kept in /data/modelseed_cache):
  curl -o /data/modelseed_cache/compounds.tsv \\
    https://raw.githubusercontent.com/ModelSEED/ModelSEEDDatabase/master/Biochemistry/compounds.tsv
  curl -o /data/modelseed_cache/reactions.tsv \\
    https://raw.githubusercontent.com/ModelSEED/ModelSEEDDatabase/master/Biochemistry/reactions.tsv

Outputs:
  docs/data/pka_table.json  : "<ns>:<id>" -> {a:[acid pKa...], b:[base conj-acid pKa...]}
      ns in {bigg, kegg}. Used to recompute metabolite charge at an arbitrary pH by
      Henderson-Hasselbalch, anchored to the model's own reference charge (pH 7).
  docs/data/thermo_rxn.json : "<ns>:<id>" -> {rev, dg, dge}
      ns in {bigg}. rev in {>,<,=}  (ModelSEED thermodynamic reversibility, a consensus
      of group-contribution + eQuilibrator). dg = drG'm (kJ/mol), dge = uncertainty.

pKa convention (verified against glycine/lysine/glutamate/histidine):
  pka col  = acid sites HA<->A- + H+          : f_deprot(pH)=1/(1+10^(pKa-pH)), charge -f
  pkb col  = pKa of the conjugate acid BH+<->B : f_prot(pH) =1/(1+10^(pH-pKa)), charge +f
"""
import json, csv, os, re, sys

SEED = "/data/modelseed_cache"   # raw ModelSEED TSVs live outside the repo (build-time inputs)
OUT = "/data/gem_autocurator/docs/data"
NULL = 10000000  # ModelSEED sentinel for unknown deltag

def parse_alias(field, want):
    """Return list of ids for db `want` from a ModelSEED aliases string."""
    out = []
    for seg in (field or "").split("|"):
        seg = seg.strip()
        if seg.lower().startswith(want.lower() + ":"):
            ids = seg.split(":", 1)[1]
            for tok in re.split(r"[;\s]+", ids.strip()):
                if tok:
                    out.append(tok)
    return out

def parse_pk(s):
    vals = []
    for entry in (s or "").split(";"):
        entry = entry.strip()
        if not entry:
            continue
        parts = entry.split(":")
        try:
            v = float(parts[-1])
        except ValueError:
            continue
        if -10 <= v <= 25:            # keep plausible pKa; drop junk
            vals.append(round(v, 2))
    return vals

# ---------------- COMPOUNDS -> pKa table ----------------
pka = {}
n_c = 0
with open(os.path.join(SEED, "compounds.tsv")) as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        if row.get("is_obsolete") == "1":
            continue
        a = parse_pk(row.get("pka"))
        b = parse_pk(row.get("pkb"))
        if not a and not b:
            continue
        rec = {}
        if a: rec["a"] = a
        if b: rec["b"] = b
        biggs = parse_alias(row.get("aliases"), "BiGG")
        keggs = parse_alias(row.get("aliases"), "KEGG")
        wrote = False
        for bid in biggs:
            pka.setdefault("bigg:" + bid, rec); wrote = True
        for kid in keggs:
            if kid.startswith("C"):
                pka.setdefault("kegg:" + kid, rec); wrote = True
        if wrote:
            n_c += 1
json.dump(pka, open(os.path.join(OUT, "pka_table.json"), "w"), separators=(",", ":"))

# ---------------- REACTIONS -> thermo table ----------------
thermo = {}
n_r = 0
with open(os.path.join(SEED, "reactions.tsv")) as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        if row.get("is_obsolete") == "1" or row.get("status", "").startswith("CPD"):
            continue
        biggs = parse_alias(row.get("aliases"), "BiGG")
        if not biggs:
            continue
        rev = row.get("reversibility") or "?"
        rec = {"rev": rev}
        try:
            dg = float(row.get("deltag"))
            if abs(dg) < NULL:
                rec["dg"] = round(dg, 2)
                dge = float(row.get("deltagerr"))
                if abs(dge) < NULL:
                    rec["dge"] = round(dge, 2)
        except (TypeError, ValueError):
            pass
        for bid in biggs:
            thermo.setdefault("bigg:" + bid, rec)
        n_r += 1
json.dump(thermo, open(os.path.join(OUT, "thermo_rxn.json"), "w"), separators=(",", ":"))

print("pka_table:", len(pka), "keys from", n_c, "compounds")
print("thermo_rxn:", len(thermo), "keys from", n_r, "reactions")
for f in ("pka_table.json", "thermo_rxn.json"):
    print(" ", f, os.path.getsize(os.path.join(OUT, f)), "bytes")
