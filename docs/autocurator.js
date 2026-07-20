/* GEM Autocurator — client-side curation engine. Everything runs in the browser. */
import GLPK from './vendor/glpk.esm.js';
const $ = (id) => document.getElementById(id);
const el = (t, c, h) => { const e = document.createElement(t); if (c) e.className = c; if (h != null) e.innerHTML = h; return e; };
const esc = (s) => String(s == null ? '' : s).replace(/[&<>"]/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));
const fmt = (n) => typeof n === 'number' ? n.toLocaleString() : n;

const STAGES = ['upload', 'ids', 'discrep', 'structure', 'charge', 'thermo', 'validate', 'report'];
const CRUMB = { upload: 'Upload model', ids: 'Identifier mapping', discrep: 'Discrepancies', structure: 'Structural QC', charge: 'Mass & charge @ pH', thermo: 'Thermodynamics & cycles', validate: 'Lab validation', report: 'Report & export' };
let REF = null;        // {met, rxn, props}
let MODEL = null;      // parsed model
let RESULT = null;     // curation result {ids, discrep, structure, charge}
const APPROVED = {};   // issueId -> true/false (undefined = pending)
try { Object.defineProperty(window, '__ac', { get: () => ({ MODEL, RESULT, REF, APPROVED, exportMAT, exportJson, applyApproved, GDB, MEDIA, SPECTRUM, conditionFromRecord, validationCoverage, spectrumFor }) }); } catch (e) {}  // debug/headless hook

/* ---------------- reference maps ---------------- */
async function fetchJsonGz(url) {   // GitHub Pages serves .gz raw (no Content-Encoding) -> decompress client-side
  const r = await fetch(url); if (!r.ok) throw new Error(url);
  const buf = new Uint8Array(await r.arrayBuffer());
  if (buf[0] === 0x1f && buf[1] === 0x8b && typeof DecompressionStream === 'function') {   // gzip magic -> raw bytes, decompress
    const s = new Response(buf).body.pipeThrough(new DecompressionStream('gzip'));
    return JSON.parse(await new Response(s).text());
  }
  return JSON.parse(new TextDecoder().decode(buf));   // server already decompressed (Content-Encoding) or plain json
}
async function loadRef() {
  if (REF) return REF;
  $('ac-load').style.display = 'flex';
  const setMsg = m => { const e = $('ac-load-msg'); if (e) e.textContent = m; };
  setMsg('Loading the identifier backbone (BiGG · KEGG · ModelSEED · MetaNetX · ChEBI · Rhea)…');
  const [met, rxn, props, pka, thermo, cov] = await Promise.all([
    fetchJsonGz('data/metabolite_map.json.gz'),
    fetchJsonGz('data/reaction_map.json.gz'),
    fetchJsonGz('data/bigg_met_props.json.gz'),
    fetch('data/pka_table.json').then(r => r.json()).catch(() => ({})),
    fetch('data/thermo_rxn.json').then(r => r.json()).catch(() => ({})),
    fetch('data/backbone_coverage.json').then(r => r.json()).catch(() => null),
  ]);
  REF = { met, rxn, props, pka, thermo, cov };
  // canonical-id sets: every id the backbone HANDS OUT (real BiGG + clustered BiGG-like). A model that
  // already uses a canonical id must map to itself, so re-curating a curated model changes nothing
  // (idempotency) — without this, a synthesised BiGG-like id (e.g. cycldi35gnyl) isn't a lookup key and
  // gets re-synthesised from a leftover DB-id name, flipping cycldi35gnyl -> c16463 on reload.
  REF.metCanon = new Set(); for (const v of Object.values(met)) if (v && v.bigg) REF.metCanon.add(v.bigg);
  REF.rxnCanon = new Set(); for (const v of Object.values(rxn)) if (v && v.bigg) REF.rxnCanon.add(v.bigg);
  $('ac-load').style.display = 'none';
  return REF;
}

/* ---------------- parsing ---------------- */
function parseModel(text, filename) {
  const t = text.trimStart();
  if (t[0] === '{') return parseCobraJson(JSON.parse(text), filename);
  if (t[0] === '<') return parseSBML(text, filename);
  throw new Error('Unrecognised format. Provide a COBRA .json or an SBML .xml/.sbml file.');
}
function parseCobraJson(j, filename) {
  const mets = (j.metabolites || []).map(m => ({
    id: m.id, name: m.name || m.id, formula: m.formula || '', charge: (m.charge == null ? null : +m.charge),
    compartment: m.compartment || compOf(m.id), anno: xrefFromAnno(m.annotation),
  }));
  const rxns = (j.reactions || []).map(r => ({
    id: r.id, name: r.name || r.id, s: r.metabolites || {}, lb: r.lower_bound, ub: r.upper_bound,
    gpr: r.gene_reaction_rule || '', anno: xrefFromAnno(r.annotation),
  }));
  return { id: j.id || filename, name: j.name || j.id || filename, filename, format: 'COBRA JSON',
    mets, rxns, genes: (j.genes || []).length, raw: j };
}
function tagText(node, ns, local) {
  const els = node.getElementsByTagName('*');
  return els;
}
function parseSBML(text, filename) {
  const doc = new DOMParser().parseFromString(text, 'application/xml');
  if (doc.getElementsByTagName('parsererror').length) throw new Error('SBML failed to parse (invalid XML).');
  const all = t => Array.from(doc.getElementsByTagNameNS('*', t));
  const attr = (n, name) => n.getAttribute(name) ?? n.getAttributeNS('http://www.sbml.org/sbml/level3/version1/fbc/version2', name);
  // fbc bounds are parameters
  const params = {};
  all('parameter').forEach(p => { params[p.getAttribute('id')] = parseFloat(p.getAttribute('value')); });
  const mets = all('species').map(s => ({
    id: s.getAttribute('id'), name: s.getAttribute('name') || s.getAttribute('id'),
    formula: attr(s, 'chemicalFormula') || '', charge: s.hasAttribute('fbc:charge') || attr(s, 'charge') != null ? parseInt(attr(s, 'charge')) : null,
    compartment: s.getAttribute('compartment') || compOf(s.getAttribute('id')), anno: xrefFromCV(s),
  }));
  const rxns = all('reaction').map(r => {
    const s = {};
    const side = (tag, sign) => Array.from(r.getElementsByTagNameNS('*', tag)).forEach(ref => {
      const sp = ref.getAttribute('species'); const st = parseFloat(ref.getAttribute('stoichiometry') || '1');
      if (sp) s[sp] = (s[sp] || 0) + sign * st;
    });
    // listOfReactants/Products contain speciesReference; select by parent
    const lr = r.getElementsByTagNameNS('*', 'listOfReactants')[0];
    const lp = r.getElementsByTagNameNS('*', 'listOfProducts')[0];
    const collect = (list, sign) => { if (!list) return; Array.from(list.getElementsByTagNameNS('*', 'speciesReference')).forEach(ref => { const sp = ref.getAttribute('species'); const st = parseFloat(ref.getAttribute('stoichiometry') || '1'); if (sp) s[sp] = (s[sp] || 0) + sign * st; }); };
    collect(lr, -1); collect(lp, +1);
    const rev = r.getAttribute('reversible') === 'true';
    const lbP = attr(r, 'lowerFluxBound'), ubP = attr(r, 'upperFluxBound');
    let lb = lbP != null ? params[lbP] : (rev ? -1000 : 0);
    let ub = ubP != null ? params[ubP] : 1000;
    return { id: r.getAttribute('id'), name: r.getAttribute('name') || r.getAttribute('id'), s, lb, ub,
      gpr: gprOf(r), anno: xrefFromCV(r) };
  });
  const modelEl = all('model')[0] || {};
  return { id: (modelEl.getAttribute && modelEl.getAttribute('id')) || filename, name: (modelEl.getAttribute && modelEl.getAttribute('name')) || filename,
    filename, format: 'SBML', mets, rxns, genes: all('geneProduct').length, raw: text };
}
function gprOf(r) {
  const g = r.getElementsByTagNameNS('*', 'geneProductRef');
  return Array.from(g).map(x => x.getAttribute('fbc:geneProduct') || x.getAttribute('geneProduct')).filter(Boolean).join(' or ');
}
function compOf(id) { const m = /_(\w)$|\[(\w)\]$/.exec(id || ''); return m ? (m[1] || m[2]) : 'c'; }
function baseId(id) { return (id || '').replace(/^M_|^R_/, '').replace(/_[a-z]\d?$/, '').replace(/\[[a-z]\d?\]$/, ''); }
function xrefFromAnno(a) {
  const out = {};
  if (!a) return out;
  const list = Array.isArray(a) ? a : Object.entries(a).flatMap(([k, v]) => (Array.isArray(v) ? v : [v]).map(x => [k, x]));
  for (const item of list) {
    let k, v;
    if (Array.isArray(item)) { [k, v] = item; } else { continue; }
    const key = String(k).toLowerCase(); const val = String(v);
    if (/kegg/.test(key)) out.kegg = val.replace(/.*[:/]/, '');
    else if (/seed|modelseed/.test(key)) out.seed = val.replace(/.*[:/]/, '');
    else if (/chebi/.test(key)) out.chebi = val.replace(/.*CHEBI:?/i, '');
    else if (/inchikey/.test(key)) out.inchikey = val.replace(/.*[:/]/, '');
    else if (/metanetx|mnx/.test(key)) out.mnx = val.replace(/.*[:/]/, '');
    else if (/rhea/.test(key)) out.rhea = val.replace(/.*[:/]/, '');
  }
  return out;
}
function xrefFromCV(node) {
  const out = {};
  Array.from(node.getElementsByTagNameNS('*', 'li')).forEach(li => {
    const r = li.getAttribute('rdf:resource') || li.getAttributeNS('http://www.w3.org/1999/02/22-rdf-syntax-ns#', 'resource') || '';
    if (/kegg\.compound|kegg\.reaction|\/C\d|\/R\d/.test(r)) out[/reaction/.test(r) ? 'kegg' : 'kegg'] = r.replace(/.*\//, '');
    else if (/seed|modelseed/i.test(r)) out.seed = r.replace(/.*\//, '');
    else if (/chebi/i.test(r)) out.chebi = r.replace(/.*CHEBI:?/i, '');
    else if (/inchikey/i.test(r)) out.inchikey = r.replace(/.*\//, '');
    else if (/metanetx/i.test(r)) out.mnx = r.replace(/.*\//, '');
    else if (/rhea/i.test(r)) out.rhea = r.replace(/.*\//, '');
  });
  return out;
}

/* ---------------- curation engine ---------------- */
// recognise when an identifier IS a database id, so we probe the right cross-ref namespace
const DBID_MET = [[/^C\d{5}$/i, 'kegg'], [/^cpd\d{5,}$/i, 'seed'], [/^MNXM\d+$/i, 'mnx'], [/^CHEBI[:_]?(\d+)$/i, 'chebi'], [/^HMDB\d+$/i, 'hmdb'], [/^([A-Z]{14}-[A-Z]{10}-[A-Z])$/i, 'inchikey']];
const DBID_RXN = [[/^R\d{5}$/i, 'kegg'], [/^rxn\d{5,}$/i, 'seed'], [/^MNXR\d+$/i, 'mnx'], [/^(\d{5,})$/, 'rhea']];
const DBNS_NAME = { kegg: 'KEGG', seed: 'ModelSEED', mnx: 'MetaNetX', chebi: 'ChEBI', rhea: 'Rhea', hmdb: 'HMDB', inchikey: 'InChIKey' };
function dbIdOf(base, table) {
  for (const [re, ns] of table) {
    const mm = base.match(re);
    if (mm) {
      let id = mm[1] || base;
      if (ns === 'kegg' || ns === 'mnx' || ns === 'hmdb' || ns === 'inchikey') id = id.toUpperCase();   // KEGG/MNX/HMDB ids are uppercase in the backbone
      else if (ns === 'seed') id = id.toLowerCase();                                                     // ModelSEED ids are lowercase
      return { ns, id };
    }
  }
  return null;
}
// normalise a namespaced id to the case the backbone uses (KEGG/MNX/HMDB upper, SEED lower, ChEBI numeric)
function nsKey(ns, v) { v = String(v); if (ns === 'kegg' || ns === 'mnx' || ns === 'hmdb' || ns === 'inchikey') v = v.toUpperCase(); else if (ns === 'seed') v = v.toLowerCase(); else if (ns === 'chebi') v = v.replace(/^CHEBI:?/i, ''); return ns + ':' + v; }
function canonMet(m) {
  const M = REF.met; const base = baseId(m.id);
  let hit = M['bigg:' + base];                                       // direct BiGG id
  if (!hit && REF.metCanon.has(base)) return { bigg: base, biggr: false, generated: false };   // already a canonical BiGG-like id -> map to itself (idempotent)
  if (!hit) { const d = dbIdOf(base, DBID_MET); if (d) hit = M[d.ns + ':' + d.id]; }   // id is itself a DB id (KEGG/SEED/MNX/ChEBI/InChIKey)
  if (!hit) for (const [ns, v] of Object.entries(m.anno || {})) { hit = M[nsKey(ns, v)]; if (hit) break; }
  if (!hit && m.anno && m.anno.inchikey) hit = M['inchikey:' + m.anno.inchikey.toUpperCase()];   // full InChIKey = exact chemical identity
  if (!hit) { const cp = base.replace(/_copy\d*$/i, ''); if (cp !== base) hit = M['bigg:' + cp]; }   // "_copy2" duplicate suffix
  if (!hit && m.name) hit = M['name:' + m.name.toLowerCase().replace(/[^a-z0-9]/g, '')];
  return hit || null;
}
function canonRxn(r) {
  const R = REF.rxn; const base = baseId(r.id);
  let hit = R['bigg:' + base];
  if (!hit && REF.rxnCanon.has(base)) return { bigg: base, biggr: false, generated: false };   // already a canonical BiGG-like id -> map to itself (idempotent)
  if (!hit) hit = R['old:' + base];
  if (!hit) { const d = dbIdOf(base, DBID_RXN); if (d) hit = R[d.ns + ':' + d.id]; }
  if (!hit) for (const [ns, v] of Object.entries(r.anno || {})) { hit = R[nsKey(ns, v)]; if (hit) break; }
  if (!hit) { const cp = base.replace(/_copy\d*$/i, ''); if (cp !== base) hit = R['bigg:' + cp] || R['old:' + cp]; }   // "_copy2" duplicate suffix
  return hit || null;
}
// strongest available chemical identity for merging duplicates (compartment handled separately)
function identityKey(m) {
  if (m.anno && m.anno.inchikey) return 'ik:' + m.anno.inchikey;         // exact structure incl. stereo
  if (m.canon && !m.canon.generated) return 'bg:' + m.canon.bigg;        // resolved BiGG / BiGG-like
  if (m.anno && m.anno.chebi) return 'chebi:' + String(m.anno.chebi).replace(/^CHEBI:?/i, '');
  return null;                                                            // identity unknown -> never merge
}

/* ---------------- robust BiGG-like id synthesis ----------------
   For metabolites/reactions absent from BiGG and our cross-ref tables, mint a
   deterministic, readable, BiGG-styled id from the name, following BiGG conventions
   (lowercase compound slugs, greek→latin, drop stereodescriptors, compartment suffix).
   Uniqueness is enforced against real BiGG ids and across generated ids. */
const GREEK = { alpha: 'a', beta: 'b', gamma: 'g', delta: 'd', epsilon: 'e', omega: 'o', 'α': 'a', 'β': 'b', 'γ': 'g', 'δ': 'd', 'ε': 'e', 'ω': 'o' };
const STOP = new Set(['acid', 'ion', 'the', 'of', 'and', 'a', 'an']);
function hash36(s) { let h = 0; for (let i = 0; i < s.length; i++) h = (h * 31 + s.charCodeAt(i)) >>> 0; return h.toString(36); }
function normName(name) {
  let s = (name || '').toLowerCase();
  s = s.replace(/\([rs+\-]\)-?|(^|[^a-z])[dl]-/g, ' ');        // (R)- (S)- (+)- D- L- stereodescriptors
  Object.keys(GREEK).forEach(g => { s = s.replace(new RegExp(g, 'g'), GREEK[g]); });
  return s;
}
function slugify(name, kind) {
  let s = normName(name);
  const words = s.split(/[^a-z0-9]+/).filter(w => w && !STOP.has(w));
  let base;
  const joined = words.join('');
  if (joined.length <= 10) base = joined;
  else if (words.length >= 2) base = words.map(w => w.length <= 4 ? w : w[0] + w.slice(1).replace(/[aeiou]/g, '').slice(0, 3)).join('').slice(0, 12);
  else base = joined.replace(/[aeiou]/g, m => '').slice(0, 8) || joined.slice(0, 8);
  base = base.replace(/[^a-z0-9]/g, '');
  if (!base) base = 'cpd' + hash36(name || 'x').slice(0, 4);
  return kind === 'rxn' ? base.toUpperCase() : base;
}
function mintBiggLike(name, fallback, kind, taken, realHas) {
  let base = slugify(name || fallback || '', kind);
  // never silently claim a real BiGG id, and never collide with another minted id (unless same source name)
  const key0 = base;
  if (taken.has(base) && taken.get(base) !== (name || fallback)) base = base + hash36(name || fallback).slice(0, 3);
  else if (realHas(base)) base = base + hash36(name || fallback).slice(0, 3);
  taken.set(base, name || fallback);
  return base;
}
const EXRE = /^(EX_|DM_|SK_|sink_|R_EX_|demand)/i;
function isExchange(r) { return EXRE.test(r.id) || Object.keys(r.s).length === 1; }
function isBiomass(r) { return /biomass|_bio\b|^BIO_|objective/i.test(r.id) || /biomass/i.test(r.name || ''); }

/* ---------------- LP layer (GLPK-WASM) ---------------- */
let _glpk = null;
async function glpk() { if (!_glpk) _glpk = await GLPK(); return _glpk; }
function buildLP(g, model, ov, extra) {
  const rows = {}; model.mets.forEach(m => { rows[m.id] = { name: m.id, vars: [], bnds: { type: g.GLP_FX, lb: 0, ub: 0 } }; });
  const bounds = [];
  model.rxns.forEach(r => {
    let lb = r.lb == null ? -1000 : r.lb, ub = r.ub == null ? 1000 : r.ub;
    if (ov && ov[r.id]) { if (ov[r.id].lb != null) lb = ov[r.id].lb; if (ov[r.id].ub != null) ub = ov[r.id].ub; }
    let type = g.GLP_DB; if (lb === ub) type = g.GLP_FX; else if (lb <= -1e30 && ub >= 1e30) type = g.GLP_FR; else if (lb <= -1e30) type = g.GLP_UP; else if (ub >= 1e30) type = g.GLP_LO;
    bounds.push({ name: r.id, type, lb, ub });
    Object.entries(r.s).forEach(([mid, c]) => { if (rows[mid]) rows[mid].vars.push({ name: r.id, coef: c }); });
  });
  if (extra) extra.forEach(e => { bounds.push({ name: e.id, type: g.GLP_DB, lb: e.lb, ub: e.ub }); Object.entries(e.s).forEach(([mid, c]) => { if (rows[mid]) rows[mid].vars.push({ name: e.id, coef: c }); }); });
  const subjectTo = Object.values(rows).filter(rw => rw.vars.length).map(rw => ({ name: rw.name, vars: rw.vars, bnds: rw.bnds }));
  return { name: 'lp', objective: { direction: g.GLP_MAX, name: 'z', vars: [] }, subjectTo, bounds };
}
async function optOf(g, lp, rid, dir) { lp.objective = { direction: dir, name: 'z', vars: [{ name: rid, coef: 1 }] }; const r = await g.solve(lp, { msglev: g.GLP_MSG_OFF, presol: true }); return (r.result && r.result.z) || 0; }
async function blockedReactions(model, onProg) {
  const g = await glpk(); const lp = buildLP(g, model);
  const internal = model.rxns.filter(r => !isExchange(r));
  const blocked = []; let i = 0;
  for (const r of internal) {
    const mx = await optOf(g, lp, r.id, g.GLP_MAX); const mn = await optOf(g, lp, r.id, g.GLP_MIN);
    if (Math.abs(mx) < 1e-7 && Math.abs(mn) < 1e-7) blocked.push(r.id);
    if (onProg && (++i % 20 === 0 || i === internal.length)) { onProg(i, internal.length); await new Promise(z => setTimeout(z, 0)); }
  }
  return blocked;
}
async function detectEGC(model) {
  const g = await glpk();
  const idOf = (bigg) => { const m = model.mets.find(x => x.id === bigg + '_c' || x.id === bigg || (x.canon && x.canon.bigg === bigg && x.compartment === 'c')); return m ? m.id : null; };
  const atp = idOf('atp'), adp = idOf('adp'), pi = idOf('pi'), h = idOf('h'), h2o = idOf('h2o');
  if (!atp || !adp || !pi) return { tested: false, reason: 'no cytosolic ATP/ADP/Pi found' };
  const ov = {}; model.rxns.forEach(r => { if (isExchange(r)) ov[r.id] = { lb: 0 }; });   // close all uptake
  const s = {}; s[atp] = -1; if (h2o) s[h2o] = -1; s[adp] = 1; s[pi] = 1; if (h) s[h] = 1;
  const lp = buildLP(g, model, ov, [{ id: '__EGC_ATP__', s, lb: 0, ub: 1000 }]);
  lp.objective = { direction: g.GLP_MAX, name: 'z', vars: [{ name: '__EGC_ATP__', coef: 1 }] };
  const r = await g.solve(lp, { msglev: g.GLP_MSG_OFF, presol: true });
  const val = (r.result && r.result.z) || 0;
  const carriers = val > 1e-6 && r.result.vars ? Object.entries(r.result.vars).filter(([k, v]) => k !== '__EGC_ATP__' && Math.abs(v) > 1e-6).map(([k, v]) => ({ id: k, flux: +v.toFixed(2) })).sort((a, b) => Math.abs(b.flux) - Math.abs(a.flux)).slice(0, 40) : [];
  return { tested: true, egc: val > 1e-6, atpFlux: +val.toFixed(2), carriers };
}
/* Flux-balance analysis: maximise objId, return objective + full flux vector. */
async function fba(model, objId, ov) {
  const g = await glpk(); const lp = buildLP(g, model, ov);
  lp.objective = { direction: g.GLP_MAX, name: 'z', vars: [{ name: objId, coef: 1 }] };
  const r = await g.solve(lp, { msglev: g.GLP_MSG_OFF, presol: true });
  return { obj: (r.result && r.result.z) || 0, vars: (r.result && r.result.vars) || {}, status: r.result ? r.result.status : 0 };
}
function findBiomass(model) {
  const bios = model.rxns.filter(isBiomass);
  const pool = (bios.length ? bios : model.rxns.filter(r => !isExchange(r)));
  return pool.slice().sort((a, b) => Object.keys(b.s).length - Object.keys(a.s).length)[0] || null;
}
// the ATP-maintenance (NGAM) reaction: by id/name, else the atp+h2o->adp+pi(+h) hydrolysis reaction
function findMaintenance(model) {
  let r = model.rxns.find(x => /^(R_)?ATPM$|_ATPM$|\bNGAM\b|maintenance/i.test(x.id) || /atp\s*maintenance|non-?growth/i.test(x.name || ''));
  if (r) return r;
  const bg = (m) => { const x = model.mets.find(z => z.id === m); return x ? (x.canon ? x.canon.bigg : baseId(x.id)) : null; };
  return model.rxns.find(x => { if (isExchange(x) || isBiomass(x)) return false; const b = Object.keys(x.s).map(bg); const s = new Set(b);
    return s.has('atp') && s.has('h2o') && s.has('adp') && s.has('pi') && b.length <= 5; }) || null;
}
/* proton/water balancing: can the reaction be balanced by adding x H+ and y H2O? returns fix or null */
function protonWaterFix(bal, cbal) {
  const H = bal.H || 0, O = bal.O || 0;                 // products - reactants imbalance
  const x = -(cbal || 0) || 0;                           // protons to add to products (H+ charge +1)
  const y = -O || 0;                                     // waters to add to products (normalise -0)
  const others = Object.entries(bal).filter(([e, v]) => e !== 'H' && e !== 'O' && Math.abs(v) > 1e-6);
  if (others.length) return null;                        // non-H/O imbalance -> manual
  if (Math.abs(H + x + 2 * y) > 1e-6) return null;       // protons+water can't reconcile H -> manual
  if (!x && !y) return null;
  return { h: x, h2o: y };                               // add x H+ (+ y H2O) to product side (neg = reactant side)
}
function fixStr(fix, comp) {
  const parts = [];
  if (fix.h) parts.push(`${Math.abs(fix.h)} H⁺ (<code>h_${comp}</code>) to the ${fix.h > 0 ? 'product' : 'reactant'} side`);
  if (fix.h2o) parts.push(`${Math.abs(fix.h2o)} H₂O (<code>h2o_${comp}</code>) to the ${fix.h2o > 0 ? 'product' : 'reactant'} side`);
  return parts.join(' and ');
}

/* ---------------- pH-dependent charge (Henderson-Hasselbalch, ModelSEED pKa) ----------------
   The metabolite's model charge is the reference at pH 7. We apply the titration slope from
   the compound's acid/base pKa sites, anchored so the recomputed charge equals the model
   charge at pH 7 (keeps balanced reactions balanced at reference, flags real pH-induced shifts). */
const DEFAULT_PH = 7.2, REF_PH = 7.2;   // BiGG charges are defined at pH 7.2; anchor there
function pkaRec(met) {
  const P = REF.pka; if (!P) return null;
  const b = met.canon ? met.canon.bigg : baseId(met.id);
  let rec = P['bigg:' + b];
  if (!rec && met.anno) { const k = met.anno.kegg || met.anno['kegg.compound']; if (k) rec = P['kegg:' + k]; }
  return rec || null;
}
function zSites(rec, pH) {                                   // absolute charge from pKa sites at pH
  let z = 0;
  if (rec.a) rec.a.forEach(a => z -= 1 / (1 + Math.pow(10, a - pH)));   // acid: deprotonates -> -1
  if (rec.b) rec.b.forEach(bp => z += 1 / (1 + Math.pow(10, pH - bp))); // base: protonates -> +1
  return z;
}
function chargeAtPH(met, refCharge, pH) {
  if (refCharge == null) return null;
  const rec = pkaRec(met);
  if (!rec || Math.abs(pH - REF_PH) < 1e-9) return refCharge;
  return refCharge + (zSites(rec, pH) - zSites(rec, REF_PH));
}

function computeMassCharge(model, pH) {
  const massIss = [], chargeIss = [];
  const byId = {}; model.mets.forEach(m => { byId[m.id] = m; });
  const F = (mid) => { const m = byId[mid]; let f = (m && m.formula) || ''; if (!f && m && m.canon && REF.props[m.canon.bigg]) f = REF.props[m.canon.bigg].f; return parseFormula(f); };
  const refCH = (m) => { if (m && m.charge != null) return m.charge; if (m && m.canon && REF.props[m.canon.bigg] && REF.props[m.canon.bigg].c != null) return REF.props[m.canon.bigg].c; return null; };
  const CH = (mid) => { const m = byId[mid]; const q = refCH(m); return q == null ? null : chargeAtPH(m, q, pH); };
  const compById = {}; model.mets.forEach(m => { compById[m.id] = m.compartment; });
  const hIn = {}, h2oIn = {}; model.mets.forEach(m => { const b = m.canon ? m.canon.bigg : baseId(m.id); if (b === 'h') hIn[m.compartment] = m.id; if (b === 'h2o') h2oIn[m.compartment] = m.id; });
  const titrated = Math.abs(pH - REF_PH) > 0.05;
  const phTag = titrated ? ` at pH ${(+pH).toFixed(1)}` : '';
  const chgTol = titrated ? 0.5 : 1e-3;   // at reference pH catch integer bugs; when titrated flag only ≥½-proton shifts
  model.rxns.forEach(r => { if (isExchange(r) || isBiomass(r)) return;
    const bal = {}; let known = true;
    Object.entries(r.s).forEach(([mid, c]) => { const f = F(mid); if (!f || !Object.keys(f).length) known = false; Object.entries(f).forEach(([e, n]) => bal[e] = (bal[e] || 0) + c * n); });
    let cbal = 0, ck = true; Object.entries(r.s).forEach(([mid, c]) => { const q = CH(mid); if (q == null) ck = false; else cbal += c * q; });
    const massOff = known ? Object.entries(bal).filter(([e, v]) => Math.abs(v) > 1e-6) : [];
    const chgOff = ck && Math.abs(cbal) > chgTol;
    if (!massOff.length && !chgOff) return;
    const comp = compById[Object.keys(r.s)[0]] || 'c';
    const fix = (known && ck) ? protonWaterFix(bal, cbal) : null;
    const canApply = fix && (!fix.h || hIn[comp]) && (!fix.h2o || h2oIn[comp]);
    const fixTxt = fix ? (canApply ? ` — <b>fix:</b> add ${fixStr(fix, comp)}` : ' — a proton/water fix balances it, but the model lacks h/h2o in this compartment') : ' — not a proton/water imbalance; needs manual curation';
    const applyObj = canApply ? { rxn: r.id, proton: fix, h: hIn[comp], h2o: h2oIn[comp] } : null;
    if (massOff.length) massIss.push({ id: 'mass_' + r.id, cat: 'charge', sub: 'mass', sev: 'bad', kind: 'mass',
      title: `Mass imbalance in <code>${esc(r.id)}</code>`, note: 'unbalanced: ' + massOff.map(([e, v]) => `${e}${v > 0 ? '+' : ''}${(+v.toFixed(2))}`).join(', ') + fixTxt, apply: applyObj });
    else chargeIss.push({ id: 'chg_' + r.id, cat: 'charge', sub: 'charge', sev: 'warn', kind: 'charge',
      title: `Charge imbalance in <code>${esc(r.id)}</code>${phTag}`, note: `net charge ${cbal > 0 ? '+' : ''}${(+cbal.toFixed(2))}` + fixTxt, apply: applyObj });
  });
  return { massIss, chargeIss };
}

function curate(model) {
  // --- identifiers ---
  const realHasMet = (b) => !!REF.met['bigg:' + b], realHasRxn = (b) => !!REF.rxn['bigg:' + b];
  const takenMet = new Map(), takenRxn = new Map();
  const unresolved = [];   // recognised DB ids not in the backbone -> given a placeholder BiGG-like id (never left as-is)
  let mBigg = 0, mLike = 0, mGen = 0, mDb = 0;
  model.mets.forEach(m => { const c = canonMet(m); const base = baseId(m.id);
    if (c) { m.canon = c; if (c.biggr) mBigg++; else mLike++; return; }
    const db = dbIdOf(base, DBID_MET);
    const like = mintBiggLike(m.name, base, 'met', takenMet, realHasMet);   // always mint a BiGG-like id
    m.canon = { bigg: like, biggr: false, generated: db ? 'dbmiss' : 'name', dbns: db ? db.ns : null };
    if (db) { mDb++; unresolved.push({ id: m.id, base, ns: db.ns, kind: 'met', like }); } else mGen++; });
  let rBigg = 0, rLike = 0, rGen = 0, rDb = 0;
  model.rxns.forEach(r => { if (isExchange(r)) { r.canon = { bigg: baseId(r.id), biggr: true, exch: true }; return; } const c = canonRxn(r); const base = baseId(r.id);
    if (c) { r.canon = c; if (c.biggr) rBigg++; else rLike++; return; }
    const db = dbIdOf(base, DBID_RXN);
    const like = mintBiggLike(r.name, base, 'rxn', takenRxn, realHasRxn);   // always mint a BiGG-like id
    r.canon = { bigg: like, biggr: false, generated: db ? 'dbmiss' : 'name', dbns: db ? db.ns : null };
    if (db) { rDb++; unresolved.push({ id: r.id, base, ns: db.ns, kind: 'rxn', like }); } else rGen++; });
  const idIssues = [];
  const dbn = (c) => DBNS_NAME[c.dbns] || c.dbns;
  const noteMet = (c) => c.biggr ? 'canonical BiGG id' : c.generated === 'dbmiss' ? `BiGG-like id — the id is a valid ${dbn(c)} id but that compound is not in MetaNetX, so it was given a stable placeholder BiGG-like id` : c.generated === 'name' ? 'BiGG-like id synthesised from the metabolite name (compound not in the backbone)' : 'canonical BiGG-like id — one stable id shared across this compound\'s database ids via MetaNetX (BiGG has none)';
  const noteRxn = (c) => c.biggr ? 'canonical BiGG reaction id' : c.generated === 'dbmiss' ? `BiGG-like id — the id is a valid ${dbn(c)} id but that reaction is not in MetaNetX, so it was given a stable placeholder BiGG-like id` : c.generated === 'name' ? 'BiGG-like id synthesised from the reaction name (reaction not in the backbone)' : 'canonical BiGG-like id — one stable id shared across this reaction\'s database ids via MetaNetX (BiGG has none)';
  // every metabolite/reaction has a canonical id; show a rename wherever it differs from the model's id
  model.mets.forEach(m => { if (m.canon && !m.canon.exch && m.canon.bigg !== baseId(m.id)) idIssues.push({ id: 'mid_' + m.id, cat: 'ids', sev: m.canon.biggr ? 'info' : 'warn', kind: 'met', generated: !!m.canon.generated,
    title: `Metabolite <code>${esc(m.id)}</code>`, from: baseId(m.id), to: m.canon.bigg + '_' + m.compartment, biggr: m.canon.biggr,
    note: noteMet(m.canon), apply: { met: m.id, newId: m.canon.bigg + '_' + m.compartment } }); });
  model.rxns.forEach(r => { if (r.canon && !r.canon.exch && r.canon.bigg !== baseId(r.id)) idIssues.push({ id: 'rid_' + r.id, cat: 'ids', sev: r.canon.biggr ? 'info' : 'warn', kind: 'rxn', generated: !!r.canon.generated,
    title: `Reaction <code>${esc(r.id)}</code>`, from: baseId(r.id), to: r.canon.bigg, biggr: r.canon.biggr,
    note: noteRxn(r.canon), apply: { rxn: r.id, newId: r.canon.bigg } }); });

  // --- discrepancies: distinct model ids that are the SAME chemical entity IN THE SAME COMPARTMENT ---
  // Identity is decided by InChIKey > resolved BiGG id > ChEBI (compartment is NOT an identity — _c/_e/_h stay separate).
  const dup = {};
  model.mets.forEach(m => { const ik = identityKey(m); if (!ik) return; const key = ik + '@@' + m.compartment; (dup[key] = dup[key] || []).push(m); });
  const discrep = [];
  const byId = {}; model.mets.forEach(m => { byId[m.id] = m; });
  Object.values(dup).forEach(arr => { if (arr.length < 2) return;
    const comp = arr[0].compartment;
    const target = (arr.map(m => m.canon).find(c => c && c.biggr) || arr[0].canon).bigg;   // prefer a real BiGG id as the merge target
    const into = target + '_' + comp;
    const forms = new Set(arr.map(m => (m.formula || '').replace(/\s/g, '')).filter(Boolean));
    const charges = new Set(arr.map(m => m.charge).filter(c => c != null));
    const conflict = forms.size > 1 || charges.size > 1;
    const caveat = conflict ? ` ⚠ Their formula/charge annotations disagree (${[...forms].join(' vs ') || 'formula'}${charges.size > 1 ? '; charge ' + [...charges].join(' vs ') : ''}) — the merge keeps one; verify it is the correct compound.` : '';
    discrep.push({ id: 'dupm_' + into, cat: 'discrep', sev: 'bad', kind: 'dupmet',
      title: `${arr.length} metabolites are the same compound`, ids: arr.map(m => m.id), to: into,
      note: `<code>${arr.map(m => esc(m.id)).join('</code>, <code>')}</code> are the same species in compartment <b>${esc(comp)}</b> (matched by ${arr.some(m => m.anno && m.anno.inchikey) ? 'InChIKey' : 'shared identifier'}) → merge to <span class="to">${esc(into)}</span>.${caveat}`,
      apply: { merge: arr.map(m => m.id), into } });
  });
  // reaction signature over CANONICAL metabolite ids (+ compartment); direction-normalised
  const rxnSig = (r, sign) => Object.entries(r.s).map(([mid, c]) => { const m = byId[mid]; const cm = (m && m.canon) ? m.canon.bigg + '_' + m.compartment : baseId(mid); return cm + ':' + (sign * c); }).sort().join('|');
  const rdup = {};
  model.rxns.forEach(r => { if (!r.canon || r.canon.exch || r.canon.generated) return; (rdup[r.canon.bigg] = rdup[r.canon.bigg] || []).push(r); });
  Object.entries(rdup).forEach(([b, arr]) => { if (arr.length < 2) return;
    const idsCode = '<code>' + arr.map(r => esc(r.id)).join('</code>, <code>') + '</code>';
    const s0 = rxnSig(arr[0], 1);
    const identical = arr.every(r => rxnSig(r, 1) === s0 || rxnSig(r, -1) === s0);   // same stoichiometry (allowing reversed direction)
    if (identical) {
      discrep.push({ id: 'dupr_' + b, cat: 'discrep', sev: 'bad', kind: 'duprxn',
        title: `${arr.length} identical reactions → one`, ids: arr.map(r => r.id), to: b,
        note: `${idsCode} resolve to <span class="to">${esc(b)}</span> and have <b>identical stoichiometry</b> (same substrates → products) — true duplicates, safe to merge.`,
        apply: { mergeRxn: arr.map(r => r.id), into: b } });
    } else {
      discrep.push({ id: 'rxc_' + b, cat: 'discrep', sev: 'warn', kind: 'rxnconflict', advisory: true, ids: arr.map(r => r.id), to: b,
        title: `${arr.length} reactions share id <code>${esc(b)}</code> but their formulas differ`,
        note: `${idsCode} map to the same BiGG id <b>${esc(b)}</b> but their reaction equations are <b>different</b> (different metabolites or stoichiometry) — these are <b>not</b> duplicates, just a name/id collision. Review each; they are <b>not</b> auto-merged.` });
    }
  });

  // --- structural QC ---
  const prod = {}, cons = {}, deg = {};
  model.mets.forEach(m => { prod[m.id] = 0; cons[m.id] = 0; deg[m.id] = 0; });
  model.rxns.forEach(r => { const rev = r.lb < 0 && r.ub > 0; Object.entries(r.s).forEach(([mid, c]) => { if (deg[mid] == null) { deg[mid] = 0; prod[mid] = 0; cons[mid] = 0; } deg[mid]++; if (c > 0 || rev) prod[mid]++; if (c < 0 || rev) cons[mid]++; }); });
  const exMet = new Set(); model.rxns.forEach(r => { if (isExchange(r)) Object.keys(r.s).forEach(x => exMet.add(x)); });
  const structure = [];
  model.mets.forEach(m => { if (exMet.has(m.id)) return; if (deg[m.id] === 0) return;
    if (prod[m.id] === 0 || cons[m.id] === 0) structure.push({ id: 'dead_' + m.id, cat: 'structure', sev: 'warn', kind: 'deadend',
      title: `Dead-end metabolite <code>${esc(m.id)}</code>`, note: `only ${prod[m.id] === 0 ? 'consumed, never produced' : 'produced, never consumed'} — cannot carry steady-state flux; blocks its reactions.`, apply: { flag: m.id } }); });
  const orphan = model.rxns.filter(r => !isExchange(r) && !(r.gpr && r.gpr.trim())).length;

  // --- mass & charge imbalance (charge is pH-dependent) ---
  const { massIss, chargeIss } = computeMassCharge(model, DEFAULT_PH);

  return { idIssues, discrep, structure, orphan, massIss, chargeIss, unresolved, ph: DEFAULT_PH,
    counts: { mBigg, mLike, mGen, mDb, rBigg, rLike, rGen, rDb, mets: model.mets.length, rxns: model.rxns.length, genes: model.genes,
      exch: model.rxns.filter(isExchange).length, dead: structure.length, mass: massIss.length, charge: chargeIss.length,
      dupMet: discrep.filter(d => d.kind === 'dupmet').length, dupRxn: discrep.filter(d => d.kind === 'duprxn').length, rxnConflict: discrep.filter(d => d.kind === 'rxnconflict').length } };
}
function parseFormula(f) {
  const out = {}; if (!f || /[^A-Za-z0-9().]/.test(f.replace(/\s/g, ''))) { }
  const re = /([A-Z][a-z]?)(\d*)/g; let m;
  while ((m = re.exec(f || ''))) { if (!m[1]) continue; out[m[1]] = (out[m[1]] || 0) + (m[2] ? +m[2] : 1); }
  return out;
}

/* ---------------- rendering ---------------- */
function goStage(stage) {
  STAGES.forEach(s => { const b = document.querySelector(`.ac-nav-item[data-stage="${s}"]`); if (b) b.classList.toggle('active', s === stage); $('stage-' + s).classList.toggle('active', s === stage); });
  $('ac-crumbs').textContent = CRUMB[stage];
  const r = { ids: renderIds, discrep: renderDiscrep, structure: renderStructure, charge: renderCharge, thermo: renderThermo, validate: renderValidate, report: renderReport }[stage];
  if (r) r();
}
function enableNav() { STAGES.forEach(s => { const b = document.querySelector(`.ac-nav-item[data-stage="${s}"]`); if (b && s !== 'upload') b.disabled = false; }); }
function navCount(stage, n, cls) { const b = document.querySelector(`.ac-nav-item[data-stage="${stage}"]`); if (!b) return; let c = b.querySelector('.cnt'); if (!c) { c = el('span', 'cnt'); b.appendChild(c); } c.className = 'cnt ' + cls; c.textContent = n; }

function issueRow(iss) {
  const row = el('div', 'ac-issue sev-' + iss.sev);
  row.dataset.iid = iss.id;
  const txt = el('div', 'txt');
  txt.appendChild(el('div', 'ttl', iss.title));
  let fix = iss.note || '';
  if (iss.from && iss.to) fix = `<span class="from">${esc(iss.from)}</span> → <span class="to">${esc(iss.to)}</span> — ${esc(iss.note || '')}`;
  txt.appendChild(el('div', 'fix', fix));
  if (iss.advisory) { row.append(txt); return row; }   // review-only finding: no approve/reject / no auto-apply
  const acts = el('div', 'ac-acts');
  const yes = el('button', 'yes', '✓'); const no = el('button', 'no', '✕');
  yes.title = 'Approve fix'; no.title = 'Reject';
  const sync = () => { yes.classList.toggle('on', APPROVED[iss.id] === true); no.classList.toggle('on', APPROVED[iss.id] === false); row.classList.toggle('approved', APPROVED[iss.id] === true); row.classList.toggle('rejected', APPROVED[iss.id] === false); };
  yes.onclick = () => { APPROVED[iss.id] = APPROVED[iss.id] === true ? undefined : true; sync(); };
  no.onclick = () => { APPROVED[iss.id] = APPROVED[iss.id] === false ? undefined : false; sync(); };
  acts.append(yes, no); row.append(txt, acts); sync();
  return row;
}
function batchBar(issues) {
  const bar = el('div', 'ac-batchbar');
  const info = el('span', '', `<b>${issues.length}</b> item${issues.length === 1 ? '' : 's'} · <span style="color:var(--ac-ok)">${issues.filter(i => APPROVED[i.id] === true).length} approved</span> · <span style="color:var(--ac-bad)">${issues.filter(i => APPROVED[i.id] === false).length} rejected</span>`);
  info.style.cssText = 'font-size:12.5px;color:var(--ink-2);margin-right:auto';
  const all = el('button', 'ac-btn ok', '✓ Approve all'); all.onclick = () => { issues.forEach(i => APPROVED[i.id] = true); goStageRefresh(); };
  const none = el('button', 'ac-btn bad', '✕ Reject all'); none.onclick = () => { issues.forEach(i => APPROVED[i.id] = false); goStageRefresh(); };
  bar.append(info, all, none); return bar;
}
let _curStage = 'ids';
function goStageRefresh() { const cur = document.querySelector('.ac-nav-item.active'); if (cur) goStage(cur.dataset.stage); }
function issueList(container, issues, emptyMsg) {
  if (!issues.length) { container.appendChild(el('div', 'ac-empty', `<span class="big">✓</span>${emptyMsg}`)); return; }
  container.appendChild(batchBar(issues));
  const list = el('div', 'ac-issues'); issues.forEach(i => list.appendChild(issueRow(i))); container.appendChild(list);
}
function kpi(v, l, cls) { return `<div class="ac-kpi ${cls || ''}"><div class="bar"></div><div class="v tabular">${fmt(v)}</div><div class="l">${esc(l)}</div></div>`; }
function drawVizIds(c) {
  if (!window.Plotly) return;
  const dark = document.documentElement.getAttribute('data-theme') === 'dark';
  const ink = dark ? '#C2CFE0' : '#334155';
  const base = { paper_bgcolor: 'rgba(0,0,0,0)', plot_bgcolor: 'rgba(0,0,0,0)', font: { family: 'Fira Sans, sans-serif', size: 11, color: ink }, margin: { l: 8, r: 8, t: 26, b: 8 } };
  const cfg = { responsive: true, displayModeBar: false };
  const donut = (id, title, vals) => window.Plotly.newPlot(id, [{ type: 'pie', hole: .62, labels: ['BiGG', 'BiGG-like (clustered)', 'BiGG-like (from name)', 'BiGG-like (out-of-backbone id)'], values: vals, marker: { colors: ['#15803D', '#2563EB', '#7C3AED', '#B45309'] }, textinfo: 'percent', textfont: { size: 10 }, hovertemplate: '%{label}: %{value}<extra></extra>', sort: false }], Object.assign({}, base, { title: { text: title, font: { size: 12 } }, showlegend: false, annotations: [{ text: fmt(vals[0] + vals[1] + vals[2] + vals[3]), x: .5, y: .5, showarrow: false, font: { size: 15, color: ink } }] }), cfg);
  donut('viz-met', 'Metabolite ids', [c.mBigg, c.mLike, c.mGen, c.mDb]);
  donut('viz-rxn', 'Reaction ids', [c.rBigg, c.rLike, c.rGen, c.rDb]);
  const labels = ['Renames', 'Dup metabolites', 'Dup reactions', 'Dead-ends', 'Mass imbalance', 'Charge imbalance'];
  const vals = [RESULT.idIssues.length, c.dupMet, c.dupRxn, c.dead, c.mass, c.charge];
  const cols = ['#2563EB', '#DC2626', '#DC2626', '#B45309', '#DC2626', '#B45309'];
  window.Plotly.newPlot('viz-iss', [{ type: 'bar', orientation: 'h', y: labels.slice().reverse(), x: vals.slice().reverse(), marker: { color: cols.slice().reverse() }, text: vals.slice().reverse().map(v => v || ''), textposition: 'outside', textfont: { size: 10 }, hovertemplate: '%{y}: %{x}<extra></extra>' }], Object.assign({}, base, { title: { text: 'Findings', font: { size: 12 } }, margin: { l: 110, r: 24, t: 26, b: 24 }, xaxis: { gridcolor: dark ? 'rgba(255,255,255,.08)' : '#EEF2F8', zeroline: false }, yaxis: { tickfont: { size: 10 } } }), cfg);
}

function renderIds() {
  const s = $('stage-ids'); s.innerHTML = ''; const c = RESULT.counts;
  const cov = REF.cov;
  s.appendChild(el('div', 'ac-sh', `<h2>Identifier mapping</h2><p>Every identifier is resolved against a <b>structure-clustered backbone</b>: MetaNetX groups every compound and reaction from BiGG · KEGG · ModelSEED · ChEBI · MetaCyc · Rhea · HMDB by InChIKey, so each entity has <b>one</b> canonical id — its BiGG id, or, when BiGG has none, <b>one stable BiGG-like id shared across all its database ids</b>. Two models using <code>glc__D</code>, <code>C00031</code> or <code>cpd00027</code> all resolve to the same id.${cov ? ` The backbone spans <b>${fmt(cov.metabolites.clusters)}</b> metabolite clusters (<b>${fmt(cov.metabolites.bigg)}</b> in BiGG, <b>${fmt(cov.metabolites.bigglike)}</b> BiGG-like) and <b>${fmt(cov.reactions.clusters)}</b> reaction clusters.` : ''}</p>`));
  s.appendChild(el('div', 'ac-kpis', kpi(c.mBigg, 'metabolites → BiGG', 'ok') + kpi(c.mLike, 'metabolites → BiGG-like', 'info') + kpi(c.mGen + c.mDb, 'metabolites → generated id', (c.mGen + c.mDb) ? 'warn' : 'ok') + kpi(c.rBigg, 'reactions → BiGG', 'ok') + kpi(c.rLike, 'reactions → BiGG-like', 'info') + kpi(c.rGen + c.rDb, 'reactions → generated id', (c.rGen + c.rDb) ? 'warn' : 'ok')));
  const pct = Math.round(100 * (c.mBigg + c.mLike) / Math.max(1, c.mets));
  // --- dashboard visualisation: coverage donuts + issue overview ---
  const viz = el('div', 'ac-card'); viz.appendChild(el('h3', '', 'Curation overview'));
  viz.appendChild(el('div', 'sub', 'Reference coverage and every finding across the pipeline, at a glance.'));
  const grid = el('div', ''); grid.style.cssText = 'display:grid;grid-template-columns:1fr 1fr 1.4fr;gap:14px';
  const d1 = el('div', 'ac-plot'); d1.id = 'viz-met'; d1.style.height = '210px';
  const d2 = el('div', 'ac-plot'); d2.id = 'viz-rxn'; d2.style.height = '210px';
  const d3 = el('div', 'ac-plot'); d3.id = 'viz-iss'; d3.style.height = '210px';
  grid.append(d1, d2, d3); viz.appendChild(grid); s.appendChild(viz);
  setTimeout(() => drawVizIds(c), 30);
  const card = el('div', 'ac-card'); card.appendChild(el('h3', '', 'Proposed identifier renames'));
  card.appendChild(el('div', 'sub', 'Approve to rename in the exported model; reject to keep the original id. Exchange/demand reactions are left as-is.'));
  issueList(card, RESULT.idIssues, 'Every identifier already matches its canonical form. Nothing to rename.');
  s.appendChild(card);
  // recognised DB ids not in the backbone -> given a placeholder BiGG-like id (still fully resolved, just flagged)
  if (RESULT.unresolved && RESULT.unresolved.length) {
    const uc = el('div', 'ac-card'); uc.appendChild(el('h3', '', 'Database ids not in the backbone — placeholder ids assigned'));
    uc.appendChild(el('div', 'sub', 'These are valid identifiers in another database (KEGG / ModelSEED / MetaNetX / ChEBI / Rhea) but that specific entity is not in MetaNetX 4.4, so it could not be clustered. Each was still given a <b>stable BiGG-like id</b> (so nothing is left un-standardised) — verify against the source if you need its exact chemistry.'));
    const ul = el('div', 'ac-issues');
    RESULT.unresolved.slice(0, 200).forEach(u => ul.appendChild(el('div', 'ac-issue sev-info', `<div class="txt" style="flex:1"><div class="ttl"><code>${esc(u.id)}</code> → <span class="to">${esc(u.like)}</span></div><div class="note" style="font-size:12px">valid <b>${esc(DBNS_NAME[u.ns] || u.ns)}</b> ${u.kind === 'met' ? 'compound' : 'reaction'} id, not in MetaNetX — assigned a placeholder BiGG-like id</div></div>`)));
    uc.appendChild(ul);
    if (RESULT.unresolved.length > 200) uc.appendChild(el('div', 'sub', `… and ${fmt(RESULT.unresolved.length - 200)} more.`));
    s.appendChild(uc);
  }
  const gtot = c.mGen + c.mDb + c.rGen + c.rDb;
  s.appendChild(el('div', 'ac-interp', `<b>Interpretation.</b> Every metabolite and reaction now carries a canonical identifier. ${pct}% of metabolites resolved to a real BiGG id (${fmt(c.mBigg)} covered), ${fmt(c.mLike)} to a clustered BiGG-like id. ${gtot ? `${fmt(c.mGen + c.mDb)} metabolites and ${fmt(c.rGen + c.rDb)} reactions were not in the backbone and were given a generated BiGG-like id${(c.mDb + c.rDb) ? ` (of which ${fmt(c.mDb + c.rDb)} are valid database ids simply missing from MetaNetX)` : ''} — nothing is left as a raw, un-standardised id.` : 'Nothing needed a generated id — the whole model maps to the backbone.'}`));
  navCount('ids', RESULT.idIssues.length, RESULT.idIssues.length ? 'warn' : 'ok');
}
function renderDiscrep() {
  const s = $('stage-discrep'); s.innerHTML = '';
  s.appendChild(el('div', 'ac-sh', `<h2>Discrepancies</h2><p>The most consequential curation: when two different identifiers in your model refer to the <b>same</b> compound or reaction. Left uncurated, these split flux, break mass balance and inflate the network. Each is resolved to one canonical entity.</p>`));
  const c = RESULT.counts;
  s.appendChild(el('div', 'ac-kpis', kpi(c.dupMet, 'duplicate metabolites', c.dupMet ? 'bad' : 'ok') + kpi(c.dupRxn, 'identical reactions', c.dupRxn ? 'bad' : 'ok') + kpi(c.rxnConflict || 0, 'id conflicts (review)', c.rxnConflict ? 'warn' : 'ok') + kpi(c.mDb + c.rDb, 'out-of-backbone ids', (c.mDb + c.rDb) ? 'warn' : 'ok')));
  const card = el('div', 'ac-card'); card.appendChild(el('h3', '', 'Duplicate entities to merge'));
  card.appendChild(el('div', 'sub', 'Approve to merge the listed ids into the single canonical id in the exported model.'));
  issueList(card, RESULT.discrep, 'No two identifiers collapse to the same entity — your namespace is clean.');
  s.appendChild(card);
  navCount('discrep', RESULT.discrep.length, RESULT.discrep.length ? 'bad' : 'ok');
}
function renderStructure() {
  const s = $('stage-structure'); s.innerHTML = ''; const c = RESULT.counts;
  s.appendChild(el('div', 'ac-sh', `<h2>Structural QC</h2><p>Network-topology problems that silently block flux: dead-end metabolites (produced or consumed but never both), and reactions with no gene association.</p>`));
  s.appendChild(el('div', 'ac-kpis', kpi(c.dead, 'dead-end metabolites', c.dead ? 'warn' : 'ok') + kpi(RESULT.orphan, 'orphan reactions (no GPR)', RESULT.orphan ? 'info' : 'ok') + kpi(c.exch, 'exchange/demand rxns', 'info')));
  const card = el('div', 'ac-card'); card.appendChild(el('h3', '', 'Dead-end metabolites'));
  card.appendChild(el('div', 'sub', 'Approve to flag these for gap-filling / removal in the report. Blocked-reaction detection (flux-based) is in the thermodynamics stage.'));
  issueList(card, RESULT.structure, 'No dead-end metabolites — every internal metabolite can be both produced and consumed.');
  s.appendChild(card);
  navCount('structure', c.dead, c.dead ? 'warn' : 'ok');
}
function renderCharge() {
  const s = $('stage-charge'); s.innerHTML = ''; const c = RESULT.counts;
  const covered = MODEL.mets.filter(m => pkaRec(m)).length, pct = MODEL.mets.length ? Math.round(100 * covered / MODEL.mets.length) : 0;
  s.appendChild(el('div', 'ac-sh', `<h2>Mass &amp; charge balance</h2><p>Every non-exchange reaction is checked for element and charge conservation. Metabolite charges are recomputed at your simulation pH by Henderson-Hasselbalch over ModelSEED pKa sites, anchored to the model's charge at pH 7.2.</p>`));
  const phbar = el('div', 'ac-card');
  phbar.innerHTML = `<h3>Simulation pH</h3><div class="sub">Charges titrate live as you move the slider: <b>${fmt(covered)}</b> of ${fmt(MODEL.mets.length)} metabolites (${pct}%) carry pKa data and re-charge; the rest keep their reference charge. Charge balance and the proton fixes below are recomputed at the chosen pH.</div>
    <div style="display:flex;align-items:center;gap:14px;margin-top:10px"><span style="font-size:12px;color:var(--ink-2)">4</span><input type="range" id="ph-slider" min="4" max="10" step="0.1" value="${RESULT.ph}" style="flex:1"><span style="font-size:12px;color:var(--ink-2)">10</span><span id="ph-val" class="tabular" style="font-size:20px;font-weight:700;color:var(--primary-2);min-width:56px;text-align:right">${(+RESULT.ph).toFixed(1)}</span></div>`;
  s.appendChild(phbar);
  const body = el('div', ''); s.appendChild(body);
  const draw = () => {
    body.innerHTML = ''; const cc0 = RESULT.counts;
    body.appendChild(el('div', 'ac-kpis', kpi(cc0.mass, 'mass-imbalanced reactions', cc0.mass ? 'bad' : 'ok') + kpi(cc0.charge, 'charge-imbalanced reactions', cc0.charge ? 'warn' : 'ok') + kpi((+RESULT.ph).toFixed(1), 'simulation pH', 'info')));
    const mc = el('div', 'ac-card'); mc.appendChild(el('h3', '', 'Mass imbalance')); mc.appendChild(el('div', 'sub', 'Elements do not conserve across the reaction (using model or BiGG formulas). pH-independent.'));
    issueList(mc, RESULT.massIss, 'Every reaction with known formulas conserves mass.'); body.appendChild(mc);
    const cc = el('div', 'ac-card'); cc.appendChild(el('h3', '', 'Charge imbalance')); cc.appendChild(el('div', 'sub', `Net charge is non-zero at pH ${(+RESULT.ph).toFixed(1)} — usually a missing proton.`));
    issueList(cc, RESULT.chargeIss, 'Charge conserves across every reaction with known charges at this pH.'); body.appendChild(cc);
    navCount('charge', cc0.mass + cc0.charge, cc0.mass ? 'bad' : (cc0.charge ? 'warn' : 'ok'));
  };
  const recompute = (pH) => {
    const { massIss, chargeIss } = computeMassCharge(MODEL, pH);
    RESULT.massIss = massIss; RESULT.chargeIss = chargeIss; RESULT.ph = pH;
    RESULT.counts.mass = massIss.length; RESULT.counts.charge = chargeIss.length;
    draw();
  };
  const sl = phbar.querySelector('#ph-slider'); let t = null;
  sl.oninput = () => { phbar.querySelector('#ph-val').textContent = (+sl.value).toFixed(1); clearTimeout(t); t = setTimeout(() => recompute(+sl.value), 120); };
  draw();
  return;
}
function renderThermo() {
  const s = $('stage-thermo'); s.innerHTML = '';
  s.appendChild(el('div', 'ac-sh', `<h2>Thermodynamics &amp; cycles</h2><p>Flux-based checks solved live with the bundled GLPK-WASM LP solver: reactions that can never carry flux, and thermodynamically infeasible energy-generating cycles that let the model make ATP from nothing.</p>`));
  if (!RESULT.thermo) {
    const internal = MODEL.rxns.filter(r => !isExchange(r)).length;
    const run = el('div', 'ac-card');
    run.innerHTML = `<h3>Run the flux scan</h3><div class="sub">Energy-generating-cycle detection is one LP; blocked-reaction detection runs flux-variability over all ${fmt(internal)} internal reactions (${fmt(internal * 2)} LPs)${internal > 700 ? ' — this may take up to a minute on a large model' : ''}.</div>`;
    const btn = el('button', 'ac-btn primary', '∮ Run thermodynamic scan'); const prog = el('div', '', ''); prog.style.cssText = 'font-size:12.5px;color:var(--ink-2);margin-top:10px';
    btn.onclick = async () => { btn.disabled = true; btn.textContent = 'Solving…';
      const egc = await detectEGC(MODEL); prog.textContent = 'Energy-generating cycles checked. Scanning for blocked reactions…';
      const blocked = await blockedReactions(MODEL, (i, n) => { prog.innerHTML = `Blocked-reaction FVA: <b>${i}/${n}</b> reactions…`; });
      RESULT.thermo = { egc, blocked }; renderThermo(); };
    run.append(btn, prog); s.appendChild(run); navCount('thermo', '?', 'info'); return;
  }
  const { egc, blocked } = RESULT.thermo;
  s.appendChild(el('div', 'ac-kpis', kpi(blocked.length, 'blocked reactions', blocked.length ? 'warn' : 'ok') + kpi(egc.tested ? (egc.egc ? egc.atpFlux : 0) : '—', 'ATP from nothing (EGC flux)', egc.egc ? 'bad' : 'ok') + kpi(egc.tested && egc.egc ? egc.carriers.length : 0, 'reactions in the EGC loop', egc.egc ? 'bad' : 'ok')));
  const ec = el('div', 'ac-card'); ec.appendChild(el('h3', '', 'Energy-generating cycle (EGC)'));
  ec.appendChild(el('div', 'sub', 'All uptake is closed and ATP hydrolysis is maximised. Any positive flux means the network can make ATP from nothing — a thermodynamically infeasible loop.'));
  if (!egc.tested) ec.appendChild(el('div', 'ac-empty', 'ATP/ADP/Pi not found in the cytosol — EGC test skipped.'));
  else if (!egc.egc) ec.appendChild(el('div', 'ac-empty', '<span class="big">✓</span>No energy-generating cycle — the model cannot make ATP without an energy source.'));
  else { ec.appendChild(el('div', 'ac-interp', `<b>Infeasible cycle found.</b> The model produces <b>${egc.atpFlux}</b> units of ATP with no substrate uptake, carried by ${egc.carriers.length} reactions. Constrain the directionality of one of these to break it:`));
    const list = el('div', 'ac-issues'); egc.carriers.forEach(cw => { const r = MODEL.rxns.find(x => x.id === cw.id); list.appendChild(issueRow({ id: 'egc_' + cw.id, sev: 'bad', title: `<code>${esc(cw.id)}</code> carries flux ${cw.flux}`, note: (r && r.name && r.name !== cw.id ? esc(r.name) + ' — ' : '') + 'candidate to make irreversible (constrain lb or ub to 0).', apply: { flag: cw.id } })); }); ec.appendChild(list); }
  s.appendChild(ec);
  const bc = el('div', 'ac-card'); bc.appendChild(el('h3', '', 'Blocked reactions'));
  bc.appendChild(el('div', 'sub', 'Reactions that cannot carry flux in any feasible state (FVA min = max = 0) — usually caused by a dead-end metabolite or a gap. Approve to flag for gap-filling / removal.'));
  const bIss = blocked.map(id => { const r = MODEL.rxns.find(x => x.id === id); return { id: 'blk_' + id, sev: 'warn', title: `Blocked reaction <code>${esc(id)}</code>`, note: (r && r.name && r.name !== id ? esc(r.name) + ' — ' : '') + 'cannot carry flux; likely a gap or a dead-end substrate.', apply: { flag: id } }; });
  issueList(bc, bIss, 'No blocked reactions — every reaction can carry flux somewhere in the solution space.'); s.appendChild(bc);
  const dirIss = directionalityIssues(MODEL);
  const dc = el('div', 'ac-card'); dc.appendChild(el('h3', '', 'Reaction directionality (ΔG consensus)'));
  dc.appendChild(el('div', 'sub', `Bounds are reconciled with ModelSEED's thermodynamic reversibility (a group-contribution + eQuilibrator consensus) and its Δ<sub>r</sub>G′<sup>m</sup>. ${fmt(RESULT.thermo.dirScanned)} of ${fmt(MODEL.rxns.filter(r=>!isExchange(r)).length)} internal reactions carry ΔG data. Only the actionable case is flagged below — the model allows both directions but ΔG makes the reaction one-way. A further <b>${fmt(RESULT.thermo.dirStricter)}</b> reaction${RESULT.thermo.dirStricter===1?' is':'s are'} more restrictive than ModelSEED's default reversibility, which is usually intentional curation and left as-is.`));
  issueList(dc, dirIss, 'No reaction is more permissive than its ΔG allows — every reversible reaction is thermodynamically bidirectional.'); s.appendChild(dc);
  navCount('thermo', blocked.length + (egc.egc ? egc.carriers.length : 0) + dirIss.length, (egc.egc || blocked.length || dirIss.length) ? 'warn' : 'ok');
}
/* directionality: compare model bounds to ModelSEED thermodynamic reversibility + ΔrG'm.
   Orientation-independent (compares the reversible/irreversible CLASS, not the sign).
   Only the high-confidence, actionable direction is flagged: the model allows both ways but
   ΔG makes the reaction effectively one-way (an infeasible-cycle risk). The opposite case —
   the model is stricter than ModelSEED's permissive default — is usually intentional curation,
   so it is reported as a count, not flagged. */
function directionalityIssues(model) {
  const T = REF.thermo; if (!T) { RESULT.thermo.dirScanned = 0; RESULT.thermo.dirStricter = 0; return []; }
  const out = []; let scanned = 0, stricter = 0;
  model.rxns.forEach(r => {
    if (isExchange(r) || isBiomass(r)) return;
    const canon = canonRxn(r); if (!canon) return;
    const t = T['bigg:' + canon.bigg]; if (!t) return;
    scanned++;
    const lb = r.lb == null ? -1000 : r.lb, ub = r.ub == null ? 1000 : r.ub;
    const modelReversible = lb < -1e-9 && ub > 1e-9;
    const dgTxt = t.dg != null ? ` Δ<sub>r</sub>G′<sup>m</sup> = ${t.dg > 0 ? '+' : ''}${t.dg}${t.dge != null ? ' ± ' + t.dge : ''} kJ/mol.` : '';
    const nm = r.name && r.name !== r.id ? esc(r.name) + ' — ' : '';
    const confIrrev = t.dg != null && Math.abs(t.dg) > 2 * (t.dge || 2) + 10;   // |ΔG| well outside uncertainty & sizeable
    if ((t.rev === '>' || t.rev === '<') && confIrrev && modelReversible) {
      out.push({ id: 'dir_' + r.id, sev: 'warn', title: `<code>${esc(r.id)}</code> is reversible but ΔG makes it one-way`,
        note: nm + `The model allows both directions (bounds ${(+lb)}…${(+ub)}), yet ModelSEED and its Δ<sub>r</sub>G make this reaction effectively irreversible — a route to infeasible cycles.${dgTxt} Constrain the uphill direction to 0 (verify the reaction's orientation first).`, apply: { flag: r.id } });
    } else if (t.rev === '=' && !modelReversible) {
      stricter++;
    }
  });
  RESULT.thermo.dirScanned = scanned; RESULT.thermo.dirStricter = stricter;
  return out;
}
/* ---------------- lab validation (GrowthDB × MediaDB × FBA) ---------------- */
let GDB = null, MEDIA = null, SPECTRUM = null;   // lazy-loaded bundles
const _valState = { sp: null, query: '' };
async function loadValidation() {
  if (GDB) return;
  const [g, m, s] = await Promise.all([
    fetch('data/growthdb.json').then(r => r.json()).catch(() => ({ species: [], records: [] })),
    fetch('data/media_ex.json').then(r => r.json()).catch(() => ({})),
    fetch('data/spectrum.json').then(r => r.json()).catch(() => ({})),
  ]);
  GDB = g; MEDIA = m; SPECTRUM = s;
  // group records by species for fast lookup
  GDB.bySpecies = {}; GDB.records.forEach((r, i) => { (GDB.bySpecies[r.sp] = GDB.bySpecies[r.sp] || []).push(i); });
}
// lightweight fuzzy score: normalized token/substring match (typo-tolerant via edit distance on words)
function editDist(a, b) { const m = a.length, n = b.length; if (!m || !n) return Math.max(m, n); const d = Array.from({ length: m + 1 }, (_, i) => [i, ...Array(n).fill(0)]); for (let j = 0; j <= n; j++) d[0][j] = j; for (let i = 1; i <= m; i++) for (let j = 1; j <= n; j++) d[i][j] = Math.min(d[i - 1][j] + 1, d[i][j - 1] + 1, d[i - 1][j - 1] + (a[i - 1] === b[j - 1] ? 0 : 1)); return d[m][n]; }
function matchSpecies(query, limit) {
  const q = query.toLowerCase().trim(); if (!q) return [];
  const qt = q.split(/\s+/);
  return GDB.species.map(sp => {
    const s = sp.toLowerCase();
    let score = 0;
    if (s === q) score = 100; else if (s.startsWith(q)) score = 80; else if (s.includes(q)) score = 60;
    else { const st = s.split(/\s+/); let hit = 0; qt.forEach(w => { if (st.some(x => x === w || (w.length > 3 && x.startsWith(w)) || editDist(w, x) <= 1)) hit++; }); score = hit ? 30 + 12 * hit - Math.min(20, editDist(q, s)) : 0; }
    return { sp, score, n: (GDB.bySpecies[sp] || []).length };
  }).filter(x => x.score > 0).sort((a, b) => b.score - a.score || b.n - a.n).slice(0, limit || 8);
}
/* ---- organism + strain auto-detection (strain name / GCF-GCA accession / BV-BRC id / taxon) ---- */
const _SCC = 'ATCC|DSMZ|DSM|NCTC|CCUG|JCM|NBRC|IFO|CECT|LMG|CIP|NCIMB|BCRC|KCTC|NRRL|CGMCC|MCCC|PCC|UTEX|CBS|KACC|VPI|NCDO|NCFB|BCCM|KCCM';
const _STRAIN_STOP = /^(wild|type|unknown|unspecified|isolate|sp|strain|and|the|not|clinical|reference|derivative|derivatives|mutant|parent|parental)$/i;
function strainStd(text) {   // mirror of build_validation.py strain_std: free-text strain -> comparable token
  if (!text) return null; const t = String(text);
  let m = t.match(new RegExp('\\b(' + _SCC + ')\\s*[-: ]?\\s*(\\d+[A-Za-z]?)\\b', 'i'));
  if (m) return (m[1] + m[2]).toUpperCase();
  // designation after substr./str./strain — a digit OR a capitalised proper name (e.g. "str. Sakai")
  m = t.match(/substr\.?\s+([A-Za-z0-9][A-Za-z0-9\-]+)/i) || t.match(/(?:str\.?|strain)\s+([A-Za-z0-9][A-Za-z0-9\-]+)/);
  if (m && !_STRAIN_STOP.test(m[1]) && (/\d/.test(m[1]) || /^[A-Z]/.test(m[1]))) return m[1].toUpperCase().replace(/[^A-Z0-9]/g, '');
  m = t.match(/\b([A-Z]{1,4}\d{1,6}[A-Za-z]?)\b/);
  if (m) return m[1].toUpperCase();
  if (/\bK-?12\b/i.test(t)) return 'K12';
  return null;
}
function strainFromName(name) {   // infraspecific remainder of an NCBI organism name = the strain designation
  if (!name) return null;
  const sp = speciesNorm(name); const rest = String(name).slice(String(name).toLowerCase().indexOf(sp.toLowerCase()) + sp.length).trim();
  return rest ? rest.replace(/^(str\.?|substr\.?)\s+/i, '') : null;
}
function strainMatch(recSstd, recStrain, detTok) {
  if (!detTok) return false;
  if (recSstd && recSstd === detTok) return true;
  const rn = (recStrain || '').toUpperCase().replace(/[^A-Z0-9]/g, '');
  if (rn && rn.includes(detTok)) return true;
  if (recSstd && (recSstd.includes(detTok) || detTok.includes(recSstd)) && detTok.length >= 4 && recSstd.length >= 4) return true;
  return false;
}
async function ncbiJson(url) { try { const r = await fetch(url); if (!r.ok) return null; return await r.json(); } catch (e) { return null; } }
async function ncbiOrgFromAccession(acc) { const d = await ncbiJson('https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/' + encodeURIComponent(acc) + '/dataset_report'); const rep = d && d.reports && d.reports[0]; const o = rep && rep.organism; if (!o) return null; const inf = o.infraspecific_names || {}; return { name: o.organism_name, taxid: o.tax_id, strain: inf.strain || inf.isolate || null }; }
async function ncbiOrgFromTaxon(tid) { const d = await ncbiJson('https://api.ncbi.nlm.nih.gov/datasets/v2/taxonomy/taxon/' + encodeURIComponent(tid)); const n = d && d.taxonomy_nodes && d.taxonomy_nodes[0] && d.taxonomy_nodes[0].taxonomy; return n ? { name: n.organism_name, taxid: n.tax_id } : null; }
function withStrain(res, srcText, ncbiStrain, ncbiName) {
  res.strain = ncbiStrain || strainFromName(ncbiName) || null;   // human-readable strain designation
  res.sstd = strainStd(ncbiStrain || '') || strainStd(ncbiName || '') || strainStd(srcText || '') || null;
  return res;
}
async function resolveOrganism(text) {
  if (!text) return null; const t = String(text).trim();
  let m = t.match(/\bGC[AF]_?\d{6,}(?:\.\d+)?/i);
  if (m) { let acc = m[0].toUpperCase().replace(/^(GC[AF])(\d)/, '$1_$2'); if (!/\.\d+$/.test(acc)) acc += '.1'; const o = await ncbiOrgFromAccession(acc); if (o) return withStrain({ name: o.name, taxid: o.taxid, via: 'assembly ' + acc + ' · NCBI', conf: true }, t, o.strain, o.name); }
  m = t.match(/\b\d{9}(?:\.\d+)?\b/);
  if (m) { const ver = /\.\d+$/.test(m[0]) ? m[0] : m[0] + '.1'; for (const pre of ['GCF_', 'GCA_']) { const o = await ncbiOrgFromAccession(pre + ver); if (o) return withStrain({ name: o.name, taxid: o.taxid, via: 'assembly ' + pre + ver + ' · NCBI', conf: true }, t, o.strain, o.name); } }
  m = t.match(/\b(\d{3,7})\.\d{1,3}\b/);
  if (m) { const tid = m[1]; const o = await ncbiOrgFromTaxon(tid); const sp = (GDB.tax2sp && GDB.tax2sp[tid]) || (o && o.name); if (sp) return withStrain({ name: sp, taxid: tid, via: 'BV-BRC id ' + m[0] + ' · taxon ' + tid + (o ? ' · NCBI' : ''), conf: true }, t, null, o && o.name); }
  m = t.match(/^\d{3,7}$/);
  if (m) { const tid = m[0]; const o = await ncbiOrgFromTaxon(tid); const sp = (GDB.tax2sp && GDB.tax2sp[tid]) || (o && o.name); if (sp) return withStrain({ name: sp, taxid: tid, via: 'taxon ' + tid + (o ? ' · NCBI' : ''), conf: true }, t, null, o && o.name); }
  const norm = t.replace(/[_\-]+/g, ' ').replace(/\.(xml|json|mat|sbml)$/i, '');
  if (/[a-z]{3,}/i.test(norm)) { const best = matchSpecies(norm, 1)[0]; if (best) return withStrain({ name: best.sp, via: 'model name', conf: best.score >= 60 }, t, null, t); }
  return null;
}
async function detectOrganism() {
  for (const c of [MODEL.name, MODEL.id]) { const r = await resolveOrganism(c); if (r) return r; }
  return null;
}
function speciesNorm(name) {   // strain/organism name -> "Genus species"
  let s = String(name || '').replace(/[\[\]]/g, '').replace(/^Candidatus\s+/i, '').trim();
  const parts = s.split(/\s+/);
  if (parts.length >= 2 && /^[A-Za-z]/.test(parts[0])) return parts[0] + ' ' + parts[1].replace(/[^A-Za-z0-9.\-]/g, '');
  return s;
}
function speciesToGdb(name) {   // resolved organism name -> best GrowthDB species (with data) or null
  if (!name) return null;
  if (GDB.bySpecies[name]) return name;
  const norm = speciesNorm(name);
  if (GDB.bySpecies[norm]) return norm;
  const best = matchSpecies(norm, 1)[0]; return best && best.score >= 55 ? best.sp : null;
}

function renderValidate() {
  const s = $('stage-validate'); s.innerHTML = '';
  s.appendChild(el('div', 'ac-sh', `<h2>Lab validation</h2><p>The organism is auto-detected from the model — a strain name, a GCF/GCA assembly accession (with or without prefix), or a BV-BRC id. Its measured growth from <b>GrowthDB</b> and medium from the <b>Media DB</b> are loaded into an editable condition; you tune the exchange fluxes and simulate, and the model's predicted growth and secretion are compared to the experiment.</p>`));
  const banner = el('div', 'ac-card'); banner.id = 'val-banner'; banner.innerHTML = `<div class="ac-load" style="display:flex;gap:10px;align-items:center"><div class="ac-spin"></div><span>Loading GrowthDB & Media DB, detecting organism…</span></div>`;
  s.appendChild(banner);
  const box = el('div', 'ac-card'); box.innerHTML = `<h3>Organism</h3><div class="sub">Auto-detected from the model id/name. Override here — matching is typo-tolerant.</div>
    <input id="val-q" type="text" placeholder="e.g. Escherichia coli, B. subtilis, GCF_000005845.2, 511145.12…" style="width:100%;margin-top:10px;padding:9px 12px;border:1px solid var(--line);border-radius:8px;background:var(--surface-2);color:var(--ink);font-size:14px">
    <div id="val-hits" style="margin-top:10px"></div>
    <div style="font-size:11px;color:var(--ink-3);margin-top:8px">Accession / id resolution sends only the identifier string to the NCBI Datasets API — your model never leaves the page.</div>`;
  s.appendChild(box);
  const results = el('div', ''); results.id = 'val-results'; s.appendChild(results);
  loadValidation().then(async () => {
    const mm = GDB.meta || {};
    box.querySelector('.sub').innerHTML += ` <b>${fmt(GDB.records.length)}</b> records / <b>${fmt(GDB.species.length)}</b> species loaded.` +
      (mm.built ? ` <span style="color:var(--ink-3);font-size:11.5px">· validation bundle built ${esc(mm.built)}${mm.growthdb_rev ? ' · GrowthDB @' + esc(mm.growthdb_rev) : ''}${mm.media_rev ? ' · Media @' + esc(mm.media_rev) : ''}</span>` : '');
    const q = $('val-q'), hits = $('val-hits');
    const showHits = () => {
      const m = matchSpecies(q.value, 8); hits.innerHTML = '';
      if (!q.value.trim()) return;
      if (!m.length) { hits.innerHTML = `<div class="ac-empty" style="padding:10px">No GrowthDB species matches “${esc(q.value)}”. <button class="ac-btn" id="val-manual2" style="margin-left:8px">Enter my own data →</button></div>`; const mb = $('val-manual2'); if (mb) mb.onclick = () => showNoData(q.value.trim()); return; }
      m.forEach(h => { const b = el('button', 'ac-chip'); b.style.cssText = 'margin:3px 6px 3px 0;padding:6px 11px;border:1px solid var(--line);border-radius:16px;background:var(--surface-2);color:var(--ink);cursor:pointer;font-size:13px';
        b.innerHTML = `<i>${esc(h.sp)}</i> <span style="color:var(--primary-2);font-weight:600">${h.n}</span>`;
        b.onclick = () => { _valState.sp = h.sp; showRecords(h.sp, strainStd(q.value), null); }; hits.appendChild(b); });
    };
    q.oninput = () => { _valState.query = q.value; showHits(); };
    // auto-detect
    const det = await detectOrganism();
    const bn = $('val-banner');
    const strainLine = (d) => d.sstd ? ` Target strain <b>${esc(d.strain || d.sstd)}</b> (<code>${esc(d.sstd)}</code>) — records are matched to it below.` : ` No strain designation in the model, so matching is at species level.`;
    if (!det) { bn.innerHTML = `<h3>Organism not auto-detected</h3><div class="sub">Couldn't read a strain name, assembly accession or BV-BRC id from <code>${esc(MODEL.id || '')}</code>${MODEL.name && MODEL.name !== MODEL.id ? ' / <code>' + esc(MODEL.name) + '</code>' : ''}. Type your strain in the box below.</div>`; }
    else {
      const gsp = det.conf ? speciesToGdb(det.name) : null;
      if (det.conf && gsp) { bn.innerHTML = `<h3>✓ <i>${esc(gsp)}</i>${det.sstd ? ' · strain <span style="color:var(--primary-2)">' + esc(det.strain || det.sstd) + '</span>' : ''}</h3><div class="sub">Detected via ${esc(det.via)}. GrowthDB has <b>${fmt((GDB.bySpecies[gsp] || []).length)}</b> records for this species.${strainLine(det)}</div>`; _valState.sp = gsp; q.value = det.strain ? gsp + ' ' + det.strain : gsp; showRecords(gsp, det.sstd, det.strain); }
      else if (det.conf) { bn.innerHTML = `<h3>Detected <i>${esc(det.name)}</i>${det.sstd ? ' · strain ' + esc(det.strain || det.sstd) : ''}</h3><div class="sub">Via ${esc(det.via)}. <b>GrowthDB has no measured growth data for this species</b>, so a database-backed validation isn't possible. You can still validate by entering your own measurements and medium.</div>`; const b = el('button', 'ac-btn primary', 'Enter my own measurements →'); b.style.marginTop = '10px'; b.onclick = () => showNoData(det.name); bn.appendChild(b); q.value = det.name; }
      else { bn.innerHTML = `<h3>Best guess: <i>${esc(det.name)}</i></h3><div class="sub">Read from the model name (low confidence). Confirm below or search for a different organism.</div>`; q.value = det.name; showHits(); }
    }
    navCount('validate', '✓', 'info');
  });
  navCount('validate', '…', 'info');
}
function o2Aerobic(rec) { const o = ((rec.cond && rec.cond.o2) || '').toLowerCase(); return /aerob|oxic/.test(o) && !/anaerob|anoxic|micro/.test(o); }
function recordRow(i, tier) {
  const r = GDB.records[i];
  const row = el('div', 'ac-issue sev-info'); row.style.cursor = 'default';
  const subs = r.up.map(u => esc(u.met || u.ex)).join(', '), prods = r.sec.map(u => esc(u.met || u.ex)).join(', ');
  const cond = [r.cond.o2, r.cond.mode, r.cond.t != null ? r.cond.t + '°C' : null, r.cond.pH != null ? 'pH ' + r.cond.pH : null].filter(Boolean).join(' · ');
  const strainBadge = r.strain ? `<span style="display:inline-block;padding:1px 8px;border-radius:10px;font-size:10.5px;font-weight:600;color:#fff;background:${tier === 'same' ? '#15803D' : tier === 'other' ? '#64748B' : '#94A3B8'};margin-left:6px">${esc(r.strain)}${r.sstd && r.sstd !== r.strain ? '' : ''}</span>` : `<span style="display:inline-block;padding:1px 8px;border-radius:10px;font-size:10.5px;color:var(--ink-3);background:var(--surface-2);margin-left:6px">strain n/a</span>`;
  const txt = el('div', 'txt'); txt.style.flex = '1';
  txt.innerHTML = `<div class="ttl">${r.mu != null ? 'μ = <b>' + r.mu + '</b> h⁻¹' : (r.dt ? 'doubling ' + r.dt + ' h' : 'rates only')}${strainBadge}${r.med.name ? ' · <span style="color:var(--ink-2)">' + esc(r.med.name) + '</span>' : ''}</div>
    <div class="note" style="font-size:12px">${cond ? esc(cond) + ' — ' : ''}${subs ? 'uptake: ' + subs + '. ' : ''}${prods ? 'secretes: ' + prods + '. ' : ''}${r.cit ? '<span style="color:var(--ink-3)">' + esc(r.cit) + '</span>' : ''}</div>`;
  const act = el('div', ''); act.style.cssText = 'flex:none;align-self:center'; const sim = el('button', 'ac-btn primary', 'Load & edit →'); sim.style.whiteSpace = 'nowrap'; sim.onclick = () => renderCondition(conditionFromRecord(i)); act.appendChild(sim);
  row.append(txt, act); return row;
}
function recText(i) { const r = GDB.records[i]; return [r.strain || '', (r.med && r.med.name) || '', r.up.map(u => u.met).join(' '), r.sec.map(u => u.met).join(' '), r.cit || '', r.cond.o2 || '', r.cond.mode || ''].join(' ').toLowerCase(); }
const PAGE = 25;
function showRecords(sp, detTok, detStrainName) {
  const results = $('val-results'); results.innerHTML = '';
  const idxs = GDB.bySpecies[sp] || [];
  if (!idxs.length) { showNoData(sp); return; }
  const same = [], other = [], nostrain = [];
  idxs.forEach(i => { const r = GDB.records[i]; if (detTok && strainMatch(r.sstd, r.strain, detTok)) same.push(i); else if (r.strain || r.sstd) other.push(i); else nostrain.push(i); });
  const head = el('div', 'ac-card');
  head.innerHTML = `<h3><i>${esc(sp)}</i> — ${idxs.length} record${idxs.length === 1 ? '' : 's'}</h3>
    <div class="sub">Pick a condition to load into the editable simulator.${detTok ? ` Records for strain <code>${esc(detTok)}</code> are grouped first; other strains are collapsed below.` : ''}</div>
    <input id="val-recq" type="text" placeholder="🔎 Filter by strain, medium, substrate or citation…" style="width:100%;margin-top:10px;padding:9px 12px;border:1px solid var(--line);border-radius:8px;background:var(--surface-2);color:var(--ink);font-size:13px">`;
  if (detTok && !same.length) head.appendChild(el('div', 'ac-interp', `<b>No strain-matched experimental data.</b> GrowthDB has ${idxs.length} record${idxs.length === 1 ? '' : 's'} for <i>${esc(sp)}</i> but none for strain <b>${esc(detStrainName || detTok)}</b> — a strain-exact validation isn't possible. Use the other-strain records below as a species-level guide, or enter your own.`));
  const mb = el('button', 'ac-btn', '✎ Enter my own condition instead'); mb.style.marginTop = '10px'; mb.onclick = () => showNoData(detStrainName ? sp + ' ' + detStrainName : sp); head.appendChild(mb);
  results.appendChild(head);
  renderCoveragePanel(results, sp, detTok, detStrainName);        // per-strain validation-resource map
  renderSpectrumPanel(results, sp, detTok);                       // grows-on-X confusion-matrix validation
  const wrap = el('div', ''); results.appendChild(wrap);
  const state = { q: '' }; const renderers = [];
  const addGroup = (label, arr, tier, color, open) => {
    if (!arr.length) return;
    const card = el('div', 'ac-card'); card.style.marginTop = '12px';
    const hdr = el('div', ''); hdr.style.cssText = 'display:flex;align-items:center;gap:8px;cursor:pointer;user-select:none';
    const chev = el('span', '', '▸'); chev.style.cssText = 'color:var(--ink-3);font-size:12px;width:12px';
    const title = el('div', ''); title.style.cssText = `font-size:13px;font-weight:700;color:${color};flex:1;letter-spacing:.02em`;
    hdr.append(chev, title); card.appendChild(hdr);
    const body = el('div', ''); body.style.marginTop = '10px'; card.appendChild(body);
    let expanded = open, page = PAGE;
    const render = () => {
      const q = state.q; const filt = q ? arr.filter(i => recText(i).includes(q)) : arr;
      const eff = q ? filt.length > 0 : expanded;   // filtering auto-reveals matching groups
      title.innerHTML = `${esc(label)} <span style="color:var(--ink-3);font-weight:500">(${filt.length}${q && filt.length !== arr.length ? ' of ' + arr.length : ''})</span>`;
      chev.textContent = eff ? '▾' : '▸';
      if (!eff) { body.style.display = 'none'; return; }
      body.style.display = ''; body.innerHTML = '';
      const shown = filt.slice(0, page);
      const l = el('div', 'ac-issues'); shown.forEach(i => l.appendChild(recordRow(i, tier))); body.appendChild(l);
      if (filt.length > page) { const more = el('button', 'ac-btn', `Show ${Math.min(PAGE, filt.length - page)} more · ${filt.length - page} remaining`); more.style.marginTop = '10px'; more.onclick = () => { page += PAGE; render(); }; body.appendChild(more); }
      if (!filt.length) body.appendChild(el('div', 'ac-empty', 'Nothing matches the filter here.'));
    };
    hdr.onclick = () => { if (state.q) return; expanded = !expanded; page = PAGE; render(); };
    renderers.push(render); wrap.appendChild(card); render();
  };
  if (detTok) {
    addGroup('THIS STRAIN — ' + (detStrainName || detTok), same, 'same', '#15803D', true);
    addGroup('OTHER STRAINS OF ' + sp.toUpperCase(), other, 'other', 'var(--ink-2)', false);
    addGroup('STRAIN UNSPECIFIED', nostrain, 'na', 'var(--ink-3)', false);
  } else {
    addGroup('ALL RECORDS', other.concat(same, nostrain), 'other', 'var(--ink-2)', true);
  }
  const q = $('val-recq'); let t = null;
  q.oninput = () => { clearTimeout(t); t = setTimeout(() => { state.q = q.value.trim().toLowerCase(); renderers.forEach(g => g()); }, 150); };
}
function showNoData(name) {
  const results = $('val-results'); results.innerHTML = '';
  renderCondition(blankCondition(name || 'your strain'), true);
}
function exIndexByMet(model) {
  const idx = {};
  model.rxns.forEach(r => { if (!isExchange(r)) return; const mid = Object.keys(r.s)[0]; if (!mid) return; const m = model.mets.find(x => x.id === mid); const b = m ? (m.canon ? m.canon.bigg : baseId(m.id)) : baseId(mid); if (b && !(b in idx)) idx[b] = r.id; });
  return idx;
}
const INORGANIC = ['h2o', 'h', 'pi', 'nh4', 'so4', 'k', 'na1', 'mg2', 'ca2', 'fe2', 'fe3', 'cl', 'co2', 'mobd', 'cu2', 'mn2', 'zn2', 'ni2', 'cobalt2', 'cbl1', 'sel', 'slnt', 'tungs'];
const INORG_SET = new Set(INORGANIC.concat(['o2']));
const CARBON_UPTAKE = 10;
function metOfExId(exId) { return exId.replace(/^R_/, '').replace(/^EX_/i, '').replace(/_[a-z0-9]+$/i, ''); }
function metName(model, bigg, exId) { const mid = Object.keys((model.rxns.find(r => r.id === exId) || { s: {} }).s)[0]; const m = mid && model.mets.find(x => x.id === mid); return (m && m.name && m.name !== m.id) ? m.name : bigg; }
// build an editable condition from a GrowthDB record: inorganics + medium + measured (experimental) uptakes
function conditionFromRecord(idx) {
  const rec = GDB.records[idx]; const exByMet = exIndexByMet(MODEL);
  const rows = [], seen = new Set(); const miss = [], medMiss = [];
  const add = (met, lb, src, meas) => { const ex = exByMet[met]; if (!ex) { if (src === 'exp') miss.push(met); return null; } if (seen.has(ex)) return rows.find(r => r.ex === ex); seen.add(ex); const row = { ex, met, name: metName(MODEL, met, ex), lb, ub: 1000, src, meas: meas || null }; rows.push(row); return row; };
  INORGANIC.forEach(b => add(b, -1000, 'inorg'));
  if (o2Aerobic(rec)) add('o2', -1000, 'inorg');
  // medium components: a linked Media DB medium, else the exchanges GrowthDB formulated from the paper's recipe.
  // a medium component with no exchange in the model can't be supplied -> track it (organic OR mineral) so a
  // no-growth verdict isn't blamed on metabolism when it's really a missing exchange.
  const medEx = (rec.med.id && MEDIA[rec.med.id]) ? MEDIA[rec.med.id].ex : (rec.med.ex || []);
  let medTot = 0;
  medEx.forEach(([e]) => { const b = metOfExId(e); const isInorg = INORG_SET.has(b); medTot++;
    const row = add(b, isInorg ? -1000 : -CARBON_UPTAKE, isInorg ? 'inorg' : 'medium');
    if (!row && !medMiss.includes(b)) medMiss.push(b); });
  rec.up.forEach(u => { if (!u.met) return;
    // a flux-usable (biomass-specific mmol/gDW/h) uptake rate is a REAL bound; otherwise normalise to -10
    const usable = u.fu && u.r < 0;
    const lb = INORG_SET.has(u.met) ? -1000 : (usable ? u.r : -CARBON_UPTAKE);
    const row = add(u.met, lb, 'exp', { r: u.r, u: u.u, fu: !!u.fu });
    if (row) { row.src = 'exp'; row.meas = { r: u.r, u: u.u, fu: !!u.fu }; if (!INORG_SET.has(u.met)) row.lb = usable ? u.r : -CARBON_UPTAKE; } });
  const mf = maintFor(rec.sp);
  return { species: rec.strain ? rec.sp + ' · ' + rec.strain : rec.sp, mu: rec.mu, mu_ok: rec.mu_ok, mu_qc: rec.mu_qc, rows, sec: rec.sec.map(x => ({ met: x.met, r: x.r, u: x.u })), miss, medMiss, medHave: medTot - medMiss.length, cit: rec.cit, doi: rec.doi, medName: rec.med.name, manual: false, ngam: mf ? mf.ngam : null, yxs: mf ? mf.yxs : null, ngamSrc: mf ? 'growthdb' : null };
}
function maintFor(sp) { return (GDB.maint && (GDB.maint[sp] || GDB.maint[speciesNorm(sp)])) || null; }
function blankCondition(species) {
  const exByMet = exIndexByMet(MODEL); const rows = [], seen = new Set();
  const add = (met, lb, src) => { const ex = exByMet[met]; if (!ex || seen.has(ex)) return; seen.add(ex); rows.push({ ex, met, name: metName(MODEL, met, ex), lb, ub: 1000, src, meas: null }); };
  INORGANIC.forEach(b => add(b, -1000, 'inorg')); add('o2', -1000, 'inorg');
  const mf = maintFor(species);
  return { species, mu: null, rows, sec: [], miss: [], cit: null, doi: null, medName: null, manual: true, ngam: mf ? mf.ngam : null, yxs: mf ? mf.yxs : null, ngamSrc: mf ? 'growthdb' : null };
}
const SRC_BADGE = { exp: ['experimental', '#15803D'], medium: ['medium', '#2563EB'], inorg: ['mineral', '#64748B'], manual: ['manual', '#7C3AED'] };
function renderCondition(cond, isManual) {
  const results = $('val-results');
  let card = document.getElementById('val-cond'); if (card) card.remove();
  card = el('div', 'ac-card'); card.id = 'val-cond';
  const exOpts = Object.entries(exIndexByMet(MODEL));
  card.innerHTML = `<h3>${cond.manual ? 'Your condition' : 'Condition'} — <i>${esc(cond.species)}</i>${cond.medName ? ' · <span style="color:var(--ink-2);font-weight:400">' + esc(cond.medName) + '</span>' : ''}</h3>
    <div class="sub">${cond.manual ? 'No GrowthDB data — enter your measured growth rate and formulate the medium by choosing exchanges and fluxes.' : 'Edit any flux before simulating. Negative = uptake (mmol gDW⁻¹ h⁻¹). Rows from measured GrowthDB rates are tagged <b>experimental</b>; medium/mineral rows are <b>pre-set (no experimental backup)</b>.'}</div>
    <div style="display:flex;align-items:center;gap:10px;margin:12px 0 6px"><label style="font-size:13px;color:var(--ink-2)">Measured μ (h⁻¹):</label><input id="val-mu" type="number" step="0.01" min="0" placeholder="unknown" value="${cond.mu == null ? '' : cond.mu}" style="width:110px;padding:6px 9px;border:1px solid var(--line);border-radius:7px;background:var(--surface-2);color:var(--ink)">${cond.mu == null ? '<span style="font-size:11.5px;color:var(--ink-3)">leave blank for a qualitative growth check</span>' : (cond.mu_ok === false ? `<span style="font-size:11.5px;color:#B45309">⚠ GrowthDB flags this μ as <b>${esc(cond.mu_qc || 'suspect')}</b> — verify before trusting</span>` : '<span style="font-size:11.5px;color:#15803D">experimental</span>')}</div>`;
  // ATP maintenance (NGAM) — from GrowthDB fit or user-entered, applied to the model's ATPM reaction
  const maintRxn = findMaintenance(MODEL);
  const mbar = el('div', ''); mbar.style.cssText = 'display:flex;align-items:center;gap:10px;margin:2px 0 8px;flex-wrap:wrap';
  mbar.innerHTML = `<label style="font-size:13px;color:var(--ink-2)">ATP maintenance / NGAM (mmol gDW⁻¹ h⁻¹):</label>
    <input id="val-ngam" type="number" step="0.01" min="0" value="${cond.ngam == null ? '' : cond.ngam}" placeholder="${maintRxn ? '0' : 'n/a'}" ${maintRxn ? '' : 'disabled'} style="width:100px;padding:6px 9px;border:1px solid var(--line);border-radius:7px;background:var(--surface-2);color:var(--ink)">
    <span style="font-size:11.5px;color:var(--ink-3)">${!maintRxn ? '— model has no ATP-maintenance reaction to constrain' : (cond.ngamSrc === 'growthdb' ? '<b style="color:#15803D">GrowthDB-fitted</b>' + (cond.yxs ? ` · biomass yield Yxs=${cond.yxs} gDW mmol⁻¹` : '') : 'no GrowthDB maintenance fit for this species — enter your own or leave blank') + ` · applies to <code>${esc(maintRxn.id)}</code>`}</span>`;
  card.appendChild(mbar);
  const tbl = el('div', ''); tbl.id = 'val-tbl'; card.appendChild(tbl);
  // add-exchange control
  const addbar = el('div', ''); addbar.style.cssText = 'display:flex;gap:8px;margin-top:10px;flex-wrap:wrap;align-items:center';
  const sel = el('select', ''); sel.style.cssText = 'flex:1;min-width:200px;max-width:340px;padding:7px;border:1px solid var(--line);border-radius:7px;background:var(--surface-2);color:var(--ink);font-size:13px';
  sel.innerHTML = '<option value="">+ add exchange…</option>' + exOpts.map(([met, ex]) => `<option value="${esc(ex)}|${esc(met)}">${esc(met)} — ${esc(ex)}</option>`).join('');
  const fluxIn = el('input', ''); fluxIn.type = 'number'; fluxIn.step = 'any'; fluxIn.value = '-10'; fluxIn.style.cssText = 'width:90px;padding:7px;border:1px solid var(--line);border-radius:7px;background:var(--surface-2);color:var(--ink)';
  const addBtn = el('button', 'ac-btn', 'Add'); addBtn.onclick = () => { if (!sel.value) return; const [ex, met] = sel.value.split('|'); if (cond.rows.some(r => r.ex === ex)) { cond.rows.find(r => r.ex === ex).lb = +fluxIn.value; } else cond.rows.push({ ex, met, name: metName(MODEL, met, ex), lb: +fluxIn.value, ub: 1000, src: 'manual', meas: null }); drawTbl(); sel.value = ''; };
  addbar.append(sel, fluxIn, addBtn); card.appendChild(addbar);
  // run
  const run = el('button', 'ac-btn primary', '▶ Run simulation'); run.style.marginTop = '14px'; run.onclick = () => { const mv = $('val-mu').value; cond.mu = mv === '' ? null : +mv; const nv = $('val-ngam'); cond.ngam = (nv && nv.value !== '') ? +nv.value : null; runCondition(cond); };
  card.appendChild(run);
  const out = el('div', ''); out.id = 'val-out-wrap'; card.appendChild(out);
  results.appendChild(card);
  const drawTbl = () => {
    tbl.innerHTML = `<div style="display:grid;grid-template-columns:1.5fr 96px 1fr 26px;gap:8px;font-size:11px;color:var(--ink-3);font-weight:600;padding:2px 0;border-bottom:1px solid var(--line)"><div>EXCHANGE</div><div>FLUX (lb)</div><div>SOURCE</div><div></div></div>`;
    cond.rows.forEach((row, i) => {
      const r = el('div', ''); r.style.cssText = 'display:grid;grid-template-columns:1.5fr 96px 1fr 26px;gap:8px;align-items:center;padding:5px 0;border-bottom:1px solid var(--line)';
      const badge = SRC_BADGE[row.src] || SRC_BADGE.manual;
      const measTxt = row.src === 'exp' && row.meas ? (row.meas.fu ? `measured ${row.meas.r} ${esc(row.meas.u || '')} — specific rate, used as the bound` : `measured ${row.meas.r} ${esc(row.meas.u || '')} — not gDW-specific, normalised to −10`) : (row.src === 'inorg' ? 'unlimited mineral' : 'pre-set, no experimental backup');
      r.innerHTML = `<div style="font-size:12.5px"><code>${esc(row.met)}</code> <span style="color:var(--ink-3);font-size:11px">${esc(row.name && row.name !== row.met ? row.name : row.ex)}</span></div>
        <input type="number" step="any" value="${row.lb}" data-i="${i}" style="width:88px;padding:5px 7px;border:1px solid var(--line);border-radius:6px;background:var(--surface-2);color:var(--ink)">
        <div style="font-size:11.5px"><span style="display:inline-block;padding:2px 8px;border-radius:10px;color:#fff;background:${badge[1]};font-size:10.5px;font-weight:600">${badge[0]}</span> <span style="color:var(--ink-3)">${measTxt}</span></div>
        <button data-del="${i}" style="background:none;border:none;color:var(--ink-3);cursor:pointer;font-size:14px">✕</button>`;
      tbl.appendChild(r);
    });
    tbl.querySelectorAll('input[data-i]').forEach(inp => inp.oninput = () => { cond.rows[+inp.dataset.i].lb = +inp.value; });
    tbl.querySelectorAll('button[data-del]').forEach(b => b.onclick = () => { cond.rows.splice(+b.dataset.del, 1); drawTbl(); });
  };
  drawTbl();
  if (cond.miss && cond.miss.length) card.insertBefore(el('div', '', `<div style="font-size:11.5px;color:#B45309;margin:8px 0">Note: measured substrate${cond.miss.length > 1 ? 's' : ''} <b>${cond.miss.map(esc).join(', ')}</b> ${cond.miss.length > 1 ? 'have' : 'has'} no exchange in this model and could not be bound.</div>`), tbl);
  // medium exchange-completeness: components the medium provides but the model can't take up (a no-growth verdict here may be an exchange gap, not a real metabolic gap)
  if (cond.medMiss && cond.medMiss.length) card.insertBefore(el('div', '', `<div style="font-size:11.5px;color:#B45309;margin:8px 0"><b>Medium incompletely represented:</b> ${cond.medHave} of ${cond.medHave + cond.medMiss.length} medium components have an exchange; the model has <b>no exchange</b> for <b>${cond.medMiss.map(esc).join(', ')}</b>. If this run predicts no growth, add the missing exchange(s)/transport before concluding it is a metabolic gap.</div>`), tbl);
}
async function runCondition(cond) {
  const wrap = $('val-out-wrap'); wrap.innerHTML = `<div class="ac-load" style="display:flex;gap:10px;align-items:center;margin-top:14px"><div class="ac-spin"></div><span>Solving FBA…</span></div>`;
  await new Promise(r => setTimeout(r, 30));
  const bio = findBiomass(MODEL);
  if (!bio) { wrap.innerHTML = '<div class="ac-empty" style="margin-top:14px">No biomass/objective reaction found in this model — cannot predict growth.</div>'; return; }
  const ov = {}; MODEL.rxns.forEach(r => { if (isExchange(r)) ov[r.id] = { lb: 0, ub: 1000 }; });
  let carbon = false; cond.rows.forEach(row => { if (!row.ex) return; ov[row.ex] = { lb: row.lb, ub: row.ub == null ? 1000 : row.ub }; if (row.lb < -1e-9 && !INORG_SET.has(row.met)) carbon = true; });
  const boundN = cond.rows.filter(r => r.lb < -1e-9).length;
  const maintRxn = findMaintenance(MODEL);
  let ngamInfo = null;
  if (cond.ngam != null && cond.ngam > 0 && maintRxn) { ov[maintRxn.id] = { lb: cond.ngam }; ngamInfo = { rxn: maintRxn.id, value: cond.ngam, src: cond.ngamSrc, yxs: cond.yxs }; }
  const card = el('div', 'ac-card'); card.id = 'val-out'; card.style.marginTop = '14px';
  const exByMet = exIndexByMet(MODEL);
  if (!carbon) { wrap.innerHTML = ''; wrap.appendChild(card); renderValidationResult(card, cond, { predMu: 0, mediaBound: boundN, bio, secCheck: [], feasible: false, noCarbon: true }); return; }
  const res = await fba(MODEL, bio.id, ov);
  const predMu = Math.max(0, res.obj);
  const secCheck = cond.sec.map(p => { const rid = exByMet[p.met]; const flux = rid ? (res.vars[rid] || 0) : null; return { met: p.met, ex: rid, measured: p.r, u: p.u, predFlux: flux == null ? null : +flux.toFixed(3), secreted: flux != null && flux > 1e-6 }; });
  wrap.innerHTML = ''; wrap.appendChild(card);
  renderValidationResult(card, cond, { predMu, mediaBound: boundN, bio, secCheck, feasible: res.obj > 1e-6, ngam: ngamInfo });
  navCount('validate', '✓', predMu > 1e-6 ? 'ok' : 'warn');
}
function renderValidationResult(card, rec, R) {
  const mu = rec.mu, pred = +R.predMu.toFixed(4);
  const spName = rec.species || rec.sp;
  let verdict, vcls, vtext;
  if (R.noCarbon) { card.innerHTML = ''; card.appendChild(el('h3', '', `Prediction vs experiment — <i>${esc(spName)}</i>`));
    card.appendChild(el('div', 'ac-empty', `<span class="big">—</span>No carbon/energy source with a negative flux is bound to an exchange that exists in <b>${esc(MODEL.id || 'this model')}</b>, so growth cannot be predicted. Add an organic exchange (uptake, negative flux) to the table above.`));
    return; }
  if (mu == null) { verdict = R.feasible ? 'Grows' : 'No growth'; vcls = R.feasible ? 'ok' : 'bad';
    vtext = R.feasible ? `The model grows on this medium (μ<sub>max</sub> = ${pred} h⁻¹). No measured μ in this record to compare against — the qualitative growth call is the check.` : `The model cannot grow on this medium — check that the reported carbon source has an exchange and pathway in the model.`; }
  else if (!R.feasible) { verdict = 'False negative'; vcls = 'bad';
    vtext = `The organism grew at μ = ${mu} h⁻¹ in the lab, but the model predicts <b>no growth</b> on this medium. Likely a gap: a missing transporter for the carbon source or a blocked biosynthetic route. See the Structural QC and Thermodynamics stages.`; }
  else { const ratio = pred / mu; const pctErr = Math.abs(pred - mu) / mu * 100;
    if (ratio >= 0.7 && ratio <= 1.45) { verdict = 'Consistent'; vcls = 'ok'; vtext = `Predicted μ<sub>max</sub> = ${pred} h⁻¹ is within ${Math.round(pctErr)}% of the measured ${mu} h⁻¹ — the model reproduces the observed growth on this medium.`; }
    else if (ratio > 1.45) { verdict = 'Over-predicts'; vcls = 'warn'; vtext = `The model's μ<sub>max</sub> (${pred} h⁻¹) exceeds the measured ${mu} h⁻¹ by ${Math.round(ratio * 100 - 100)}%. FBA gives the theoretical maximum; the gap is expected but a large excess can flag a missing maintenance (NGAM/ATPM) constraint or an over-open medium.`; }
    else { verdict = 'Under-predicts'; vcls = 'bad'; vtext = `The model's μ<sub>max</sub> (${pred} h⁻¹) is well below the measured ${mu} h⁻¹ — a likely gap in a biosynthetic pathway, or too tight a medium/uptake bound (uptake was normalised to 10 mmol gDW⁻¹ h⁻¹).`; }
  }
  const secGood = R.secCheck.filter(x => x.secreted).length, secTot = R.secCheck.length;
  card.innerHTML = '';
  card.appendChild(el('h3', '', `Prediction vs experiment — <i>${esc(spName)}</i>`));
  card.appendChild(el('div', 'ac-kpis', kpi(mu == null ? '—' : mu, 'measured μ (h⁻¹)', 'info') + kpi(pred, 'predicted μ_max (h⁻¹)', vcls) + kpi(`${secGood}/${secTot}`, 'secretion products reproduced', secTot ? (secGood === secTot ? 'ok' : 'warn') : 'info') + kpi(R.mediaBound, 'exchanges taking up flux', 'info')));
  const bcol = { ok: '#15803D', warn: '#B45309', bad: '#DC2626', info: '#2563EB' }[vcls] || '#2563EB';
  const chip = el('div', ''); chip.style.cssText = 'margin:4px 0 12px'; chip.innerHTML = `<span style="display:inline-block;padding:5px 14px;border-radius:14px;font-weight:600;font-size:13px;color:#fff;background:${bcol}">${verdict}</span>`;
  card.appendChild(chip);
  const plot = el('div', 'ac-plot'); plot.id = 'val-plot'; plot.style.height = '230px'; card.appendChild(plot);
  card.appendChild(el('div', 'ac-interp', vtext));
  if (R.ngam) card.appendChild(el('div', '', `<div style="font-size:12px;color:var(--ink-2);margin-top:8px;padding:8px 10px;background:var(--surface-2);border-radius:8px">⚙ <b>ATP maintenance applied.</b> NGAM constrained to <b>${R.ngam.value}</b> mmol gDW⁻¹ h⁻¹ on <code>${esc(R.ngam.rxn)}</code> (${R.ngam.src === 'growthdb' ? 'GrowthDB-fitted' : 'your value'}${R.ngam.yxs ? `, biomass yield Yxs=${R.ngam.yxs}` : ''}). This burns ATP that would otherwise drive growth, so the predicted μ is the maintenance-corrected maximum — closer to the real rate than the unconstrained theoretical maximum.</div>`));
  if (secTot) {
    const st = el('div', ''); st.style.marginTop = '10px';
    st.innerHTML = '<div style="font-size:12.5px;color:var(--ink-2);margin-bottom:6px"><b>Secretion pattern</b> (does the model route flux to each measured by-product?):</div>' +
      R.secCheck.map(x => `<span style="display:inline-block;margin:0 6px 6px 0;padding:4px 10px;border-radius:12px;font-size:12.5px;background:var(--surface-2);border:1px solid var(--line)">${x.secreted ? '✓' : (x.ex ? '✕' : '—')} <b>${esc(x.met)}</b> <span style="color:var(--ink-3)">meas ${x.measured} ${esc(x.u || '')}${x.ex ? '; model ' + (x.predFlux != null ? x.predFlux : 'n/a') : '; no exchange in model'}</span></span>`).join('');
    card.appendChild(st);
  }
  if (rec.cit) card.appendChild(el('div', '', `<div style="font-size:11.5px;color:var(--ink-3);margin-top:12px;border-top:1px solid var(--line);padding-top:8px">Source: ${esc(rec.cit)}${rec.doi ? ' · <a href="https://doi.org/' + esc(rec.doi) + '" target="_blank" rel="noopener" style="color:var(--primary-2)">doi:' + esc(rec.doi) + '</a>' : ''}</div>`));
  setTimeout(() => drawValPlot(rec, R, verdict, vcls), 30);
}
function drawValPlot(rec, R, verdict, vcls) {
  if (!window.Plotly) return;
  const dark = document.documentElement.getAttribute('data-theme') === 'dark';
  const ink = dark ? '#C2CFE0' : '#334155';
  const col = { ok: '#15803D', warn: '#B45309', bad: '#DC2626', info: '#2563EB' }[vcls] || '#2563EB';
  const traces = [], hasMu = rec.mu != null;
  const cats = [], meas = [], predv = [];
  if (hasMu) { cats.push('growth μ (h⁻¹)'); meas.push(rec.mu); predv.push(+R.predMu.toFixed(4)); }
  window.Plotly.newPlot('val-plot', [
    { type: 'bar', name: 'measured', x: cats, y: meas, marker: { color: '#94A3B8' }, text: meas.map(v => v), textposition: 'outside', textfont: { size: 11 } },
    { type: 'bar', name: 'predicted', x: cats, y: predv, marker: { color: col }, text: predv.map(v => v), textposition: 'outside', textfont: { size: 11 } },
  ], { barmode: 'group', paper_bgcolor: 'rgba(0,0,0,0)', plot_bgcolor: 'rgba(0,0,0,0)', font: { family: 'Fira Sans, sans-serif', size: 12, color: ink }, margin: { l: 48, r: 16, t: 30, b: 30 }, title: { text: 'Measured vs predicted growth rate', font: { size: 13 } }, legend: { orientation: 'h', y: 1.15, x: .5, xanchor: 'center' }, yaxis: { gridcolor: dark ? 'rgba(255,255,255,.08)' : '#EEF2F8', zeroline: false, rangemode: 'tozero' } }, { responsive: true, displayModeBar: false });
}

/* ---------------- substrate-utilisation spectrum validation (grows-on-X confusion matrix) ---------------- */
function spectrumFor(sp, strainTok) {
  const s = SPECTRUM && SPECTRUM[sp]; if (!s) return null;
  const pos = new Set(s.p || []), neg = new Set(s.n || []);
  let strainN = 0;
  if (strainTok && s.s) for (const [st, exs] of Object.entries(s.s)) {
    if (strainMatch(strainStd(st), st, strainTok)) { strainN += exs.length; exs.forEach(e => { pos.add(e); neg.delete(e); }); }
  }
  return { pos, neg, strainN };
}
// pure substrate sweep -> confusion counts (no DOM). onProg(i,n) for progress; cap limits solves.
async function computeSpectrum(sp, strainTok, aerobic, cap, onProg) {
  const spec = spectrumFor(sp, strainTok); const bio = findBiomass(MODEL);
  if (!spec || !bio) return null;
  const exByMet = exIndexByMet(MODEL);
  const obsMap = new Map();
  spec.pos.forEach(ex => obsMap.set(metOfExId(ex), 'pos'));
  spec.neg.forEach(ex => { const b = metOfExId(ex); if (!obsMap.has(b)) obsMap.set(b, 'neg'); });
  let items = []; obsMap.forEach((obs, b) => { const rid = exByMet[b]; if (rid) items.push({ b, rid, obs }); });
  const notInModel = obsMap.size - items.length;
  // cap must KEEP CLASS BALANCE: positives are inserted before negatives, so a naive first-N slice
  // takes only grows-on calls -> zero negatives -> TN=FP=0 -> MCC collapses to 0. Sample both classes.
  const capped = cap && items.length > cap;
  if (capped) {
    const pos = items.filter(x => x.obs === 'pos'), neg = items.filter(x => x.obs === 'neg');
    const nNeg = Math.min(neg.length, Math.max(neg.length ? 1 : 0, Math.round(cap * neg.length / items.length)));
    const nPos = Math.min(pos.length, cap - nNeg);
    items = pos.slice(0, nPos).concat(neg.slice(0, nNeg));
  }
  const NPS = ['nh4', 'pi', 'so4'];
  let TP = 0, TN = 0, FP = 0, FN = 0; const fp = [], fn = [];
  for (let i = 0; i < items.length; i++) {
    const it = items[i];
    const ov = {}; MODEL.rxns.forEach(r => { if (isExchange(r)) ov[r.id] = { lb: 0, ub: 1000 }; });
    INORGANIC.concat(NPS).forEach(b => { const rid = exByMet[b]; if (rid) ov[rid] = { lb: -1000, ub: 1000 }; });
    if (aerobic) { const o = exByMet['o2']; if (o) ov[o] = { lb: -1000, ub: 1000 }; }
    ov[it.rid] = { lb: -10, ub: 1000 };                         // the one carbon source under test
    const res = await fba(MODEL, bio.id, ov);
    const pred = res.obj > 1e-4, obs = it.obs === 'pos';
    if (pred && obs) TP++; else if (!pred && !obs) TN++;
    else if (pred && !obs) { FP++; fp.push(it); } else { FN++; fn.push(it); }
    if (onProg && (i % 15 === 0 || i === items.length - 1)) { onProg(i + 1, items.length); await new Promise(z => setTimeout(z, 0)); }
  }
  return { TP, TN, FP, FN, fp, fn, tested: items.length, notInModel, obsTotal: obsMap.size, capped, strainN: spec.strainN };
}
function spectrumMCC(R) { const den = Math.sqrt((R.TP + R.FP) * (R.TP + R.FN) * (R.TN + R.FP) * (R.TN + R.FN)); return den ? (R.TP * R.TN - R.FP * R.FN) / den : 0; }
async function runSpectrumValidation(sp, strainTok, host, aerobic) {
  const spec = spectrumFor(sp, strainTok); const bio = findBiomass(MODEL);
  host.innerHTML = `<div class="ac-load" style="display:flex;gap:10px;align-items:center"><div class="ac-spin"></div><span>Sweeping substrates…</span></div>`;
  if (!bio) { host.innerHTML = '<div class="ac-empty">No biomass reaction — cannot run the spectrum.</div>'; return; }
  if (!spec) { host.innerHTML = `<div class="ac-empty">No substrate spectrum in GrowthDB for <i>${esc(sp)}</i>.</div>`; return; }
  const prog = el('div', ''); prog.style.cssText = 'font-size:12.5px;color:var(--ink-2);margin-top:8px'; host.innerHTML = ''; host.appendChild(prog);
  const R = await computeSpectrum(sp, strainTok, aerobic, 0, (i, n) => { prog.innerHTML = `Testing substrate <b>${i}/${n}</b>…`; });
  if (!R.tested) { host.innerHTML = `<div class="ac-empty">None of the ${fmt(R.obsTotal)} substrates in GrowthDB's spectrum for <i>${esc(sp)}</i> have an exchange in this model — nothing to test.</div>`; return; }
  renderSpectrumResult(host, sp, { ...R, strainTok, aerobic });
}
function renderSpectrumResult(host, sp, R) {
  const n = R.TP + R.TN + R.FP + R.FN;
  const acc = n ? (R.TP + R.TN) / n : 0;
  const sens = (R.TP + R.FN) ? R.TP / (R.TP + R.FN) : 0, spc = (R.TN + R.FP) ? R.TN / (R.TN + R.FP) : 0;
  const den = Math.sqrt((R.TP + R.FP) * (R.TP + R.FN) * (R.TN + R.FP) * (R.TN + R.FN));
  const mcc = den ? (R.TP * R.TN - R.FP * R.FN) / den : 0;
  const mcls = mcc >= 0.6 ? 'ok' : mcc >= 0.3 ? 'warn' : 'bad';
  host.innerHTML = '';
  host.appendChild(el('div', 'ac-kpis', kpi((100 * acc).toFixed(0) + '%', 'accuracy', acc >= 0.8 ? 'ok' : 'warn') + kpi(mcc.toFixed(2), 'MCC', mcls) + kpi((100 * sens).toFixed(0) + '%', 'sensitivity (recall)', 'info') + kpi((100 * spc).toFixed(0) + '%', 'specificity', 'info') + kpi(n, 'substrates tested', 'info')));
  // 2x2 confusion matrix
  const cm = el('div', ''); cm.style.cssText = 'display:grid;grid-template-columns:auto 1fr 1fr;gap:2px;max-width:420px;margin:6px 0 12px;font-size:12.5px';
  const cell = (t, bg, cl) => `<div style="padding:8px 10px;background:${bg};color:${cl || 'var(--ink)'};border-radius:6px;text-align:center">${t}</div>`;
  cm.innerHTML = cell('', 'transparent') + cell('<b>obs: grows</b>', 'var(--surface-2)') + cell('<b>obs: no-grow</b>', 'var(--surface-2)') +
    cell('<b>model: grows</b>', 'var(--surface-2)') + cell(`<b>${R.TP}</b><br>true +`, '#dcfce7', '#15803D') + cell(`<b>${R.FP}</b><br>false +`, '#fef3c7', '#B45309') +
    cell('<b>model: no-grow</b>', 'var(--surface-2)') + cell(`<b>${R.FN}</b><br>false −`, '#fee2e2', '#DC2626') + cell(`<b>${R.TN}</b><br>true −`, '#dcfce7', '#15803D');
  host.appendChild(cm);
  host.appendChild(el('div', 'ac-interp', `Tested <b>${R.tested}</b> substrates that have an exchange in the model${R.notInModel ? ` (${R.notInModel} more in GrowthDB's spectrum have no exchange here — a coverage gap)` : ''}. ${R.strainN ? `<b>${R.strainN}</b> calls are specific to strain <code>${esc(R.strainTok)}</code>; the rest are the species consensus.` : 'Species-level spectrum (no strain-specific calls matched).'} O₂ ${R.aerobic ? 'open (aerobic test)' : 'closed (anaerobic test)'}.`));
  const list = (title, arr, note, color) => { if (!arr.length) return; const c = el('div', 'ac-card'); c.style.marginTop = '8px'; c.appendChild(el('h3', '', title)); c.appendChild(el('div', 'sub', note));
    const l = el('div', ''); l.style.cssText = 'display:flex;flex-wrap:wrap;gap:6px'; arr.slice(0, 60).forEach(x => { const s = el('span', ''); s.style.cssText = `padding:3px 9px;border-radius:12px;font-size:12px;background:var(--surface-2);border:1px solid ${color}`; s.innerHTML = `<code>${esc(x.b)}</code>`; l.appendChild(s); }); c.appendChild(l);
    if (arr.length > 60) c.appendChild(el('div', 'sub', `…and ${arr.length - 60} more.`)); host.appendChild(c); };
  list(`False negatives — model can't grow but the organism does (${R.fn.length})`, R.fn, 'The strongest curation targets: a missing transporter or biosynthetic/catabolic pathway blocks growth on these. Gap-fill candidates.', '#DC2626');
  list(`False positives — model grows but the organism doesn't (${R.fp.length})`, R.fp, 'The model is over-permissive on these: a missing regulatory constraint, an erroneous reaction, or a gap-fill that shouldn\'t be there.', '#B45309');
}
function validationCoverage(sp, detTok) {
  const idxs = GDB.bySpecies[sp] || [];
  let mu = 0, muStrain = 0, flux = 0, media = 0, strainRecs = 0;
  idxs.forEach(i => { const r = GDB.records[i]; const hit = detTok && strainMatch(r.sstd, r.strain, detTok);
    if (hit) strainRecs++;
    if (r.mu != null && r.mu_ok !== false) { mu++; if (hit) muStrain++; }
    if ((r.up || []).concat(r.sec || []).some(x => x.fu)) flux++;
    if (r.med && (r.med.id || (r.med.ex && r.med.ex.length))) media++;
  });
  const spec = spectrumFor(sp, detTok); const maint = maintFor(sp);
  return { cond: idxs.length, strainRecs, mu, muStrain, flux, media,
    specPos: spec ? spec.pos.size : 0, specNeg: spec ? spec.neg.size : 0, specStrain: spec ? spec.strainN : 0, maint: !!maint };
}
function renderCoveragePanel(host, sp, detTok, detStrainName) {
  const c = validationCoverage(sp, detTok);
  const card = el('div', 'ac-card'); card.style.cssText = 'margin-top:12px';
  card.appendChild(el('h3', '', `Validation resources — <i>${esc(sp)}</i>${detTok ? ' · strain <span style="color:var(--primary-2)">' + esc(detStrainName || detTok) + '</span>' : ''}`));
  card.appendChild(el('div', 'sub', 'What GrowthDB can validate for this organism. Coverage differs by strain — strain-specific evidence is the most reliable; species-level fills the gaps at lower confidence.'));
  card.appendChild(el('div', 'ac-kpis',
    kpi(c.mu, 'growth-rate (µ) conditions', c.mu ? 'ok' : 'warn') +
    kpi(detTok ? c.muStrain : '—', 'of them strain-specific', c.muStrain ? 'ok' : 'info') +
    kpi(c.flux, 'flux-usable rate conditions', c.flux ? 'ok' : 'info') +
    kpi(c.media, 'GEM-ready media', c.media ? 'ok' : 'warn') +
    kpi(c.specPos + c.specNeg, 'substrate spectrum calls', (c.specPos + c.specNeg) ? 'ok' : 'warn') +
    kpi(detTok ? c.specStrain : '—', 'strain-specific substrates', c.specStrain ? 'ok' : 'info') +
    kpi(c.maint ? 'yes' : 'no', 'GrowthDB NGAM/yield fit', c.maint ? 'ok' : 'info')));
  const avail = [];
  if (c.mu) avail.push('<b>growth-rate</b> (µ vs FBA)'); if (c.flux) avail.push('<b>uptake/secretion flux</b>');
  if (c.specPos + c.specNeg) avail.push('<b>substrate-utilisation spectrum</b> (confusion matrix)'); if (c.maint) avail.push('<b>maintenance/yield</b>');
  card.appendChild(el('div', 'ac-interp', `${avail.length ? 'Runnable validations: ' + avail.join(' · ') + '.' : 'No quantitative validation data for this organism.'} ${detTok && !c.strainRecs ? `<b>No growth records for strain ${esc(detStrainName || detTok)}</b> — the µ/rate/media checks fall back to the species (lower confidence); only the ${c.specStrain} strain-specific substrate calls are strain-exact.` : ''}`));
  // one-click aggregate scorecard across every validation type
  if (avail.length) {
    const bar = el('div', ''); bar.style.cssText = 'display:flex;gap:10px;align-items:center;margin-top:12px;flex-wrap:wrap';
    const btn = el('button', 'ac-btn primary', '◆ Score this GEM against GrowthDB');
    const aer = el('label', ''); aer.style.cssText = 'font-size:12.5px;color:var(--ink-2);display:flex;align-items:center;gap:5px';
    aer.innerHTML = '<input type="checkbox" id="score-aer" checked> assume aerobic';
    const sc = el('div', ''); sc.id = 'score-out'; sc.style.marginTop = '12px';
    btn.onclick = () => { btn.disabled = true; runScorecard(sp, detTok, sc, document.getElementById('score-aer').checked).finally(() => btn.disabled = false); };
    bar.append(btn, aer); card.appendChild(bar); card.appendChild(sc);
  }
  host.appendChild(card);
}
// pure FBA of an editable condition (no DOM) — used by the aggregate scorecard
async function solveCondition(cond) {
  const bio = findBiomass(MODEL); if (!bio) return null;
  const ov = {}; MODEL.rxns.forEach(r => { if (isExchange(r)) ov[r.id] = { lb: 0, ub: 1000 }; });
  let carbon = false; cond.rows.forEach(row => { if (!row.ex) return; ov[row.ex] = { lb: row.lb, ub: row.ub == null ? 1000 : row.ub }; if (row.lb < -1e-9 && !INORG_SET.has(row.met)) carbon = true; });
  const maintRxn = findMaintenance(MODEL);
  if (cond.ngam != null && cond.ngam > 0 && maintRxn) ov[maintRxn.id] = { lb: cond.ngam };
  if (!carbon) return { predMu: 0, feasible: false, noCarbon: true, res: null };
  const res = await fba(MODEL, bio.id, ov);
  return { predMu: Math.max(0, res.obj), feasible: res.obj > 1e-6, res };
}
async function runScorecard(sp, detTok, host, aerobic) {
  host.innerHTML = `<div class="ac-load" style="display:flex;gap:10px;align-items:center"><div class="ac-spin"></div><span>Scoring the GEM across all validation types…</span></div>`;
  await new Promise(r => setTimeout(r, 20));
  const bio = findBiomass(MODEL);
  if (!bio) { host.innerHTML = '<div class="ac-empty">No biomass reaction — cannot score this model.</div>'; return; }
  const exByMet = exIndexByMet(MODEL);
  const idxs = GDB.bySpecies[sp] || [];
  // conditions with a usable medium (linked or formulated); strain-specific first, then cap for responsiveness
  const withMed = idxs.filter(i => { const r = GDB.records[i]; return r.med && ((r.med.id && MEDIA[r.med.id]) || (r.med.ex && r.med.ex.length)); });
  withMed.sort((a, b) => (detTok && strainMatch(GDB.records[b].sstd, GDB.records[b].strain, detTok) ? 1 : 0) - (detTok && strainMatch(GDB.records[a].sstd, GDB.records[a].strain, detTok) ? 1 : 0));
  const CAP = 40; const sample = withMed.slice(0, CAP);
  const prog = el('div', ''); prog.style.cssText = 'font-size:12.5px;color:var(--ink-2)'; host.innerHTML = ''; host.appendChild(prog);
  let feasN = 0, grow = 0; const muPairs = []; let secHit = 0, secTot = 0; let strainCond = 0;
  for (let k = 0; k < sample.length; k++) {
    const i = sample[k]; const cond = conditionFromRecord(i);
    if (detTok && strainMatch(GDB.records[i].sstd, GDB.records[i].strain, detTok)) strainCond++;
    const r = await solveCondition(cond);
    if (r && !r.noCarbon) { feasN++; if (r.feasible) grow++;
      if (cond.mu != null && cond.mu_ok !== false && r.feasible) muPairs.push([cond.mu, +r.predMu.toFixed(4)]);
      cond.sec.forEach(p => { const rid = exByMet[p.met]; if (!rid) return; secTot++; const f = r.res ? (r.res.vars[rid] || 0) : 0; if (f > 1e-6) secHit++; }); }
    if (k % 4 === 0 || k === sample.length - 1) { prog.innerHTML = `Simulating growth condition <b>${k + 1}/${sample.length}</b>…`; await new Promise(z => setTimeout(z, 0)); }
  }
  prog.innerHTML = 'Sweeping the substrate spectrum…'; await new Promise(z => setTimeout(z, 0));
  const spec = await computeSpectrum(sp, detTok, aerobic, 90, (i, n) => { prog.innerHTML = `Sweeping substrate <b>${i}/${n}</b>…`; });
  renderScorecard(host, sp, detTok, {
    feasN, grow, feasAcc: feasN ? grow / feasN : null, sampled: sample.length, withMed: withMed.length, strainCond,
    muPairs, secHit, secTot, secRecall: secTot ? secHit / secTot : null, spec, aerobic
  });
}
function renderScorecard(host, sp, detTok, S) {
  host.innerHTML = '';
  const dims = [];   // {label, val, cls, detail}
  if (S.feasN) { const a = S.feasAcc; dims.push({ label: 'Growth feasibility', val: (100 * a).toFixed(0) + '%', cls: a >= 0.9 ? 'ok' : a >= 0.7 ? 'warn' : 'bad', detail: `${S.grow}/${S.feasN} should-grow conditions produce growth on their reported medium`, w: a }); }
  if (S.muPairs.length) {
    const ratios = S.muPairs.map(([m, p]) => m > 0 ? p / m : null).filter(x => x != null && isFinite(x));
    ratios.sort((x, y) => x - y); const med = ratios.length ? ratios[Math.floor(ratios.length / 2)] : null;
    const within = S.muPairs.filter(([m, p]) => m > 0 && p / m >= 0.5 && p / m <= 2).length; const frac = within / S.muPairs.length;
    dims.push({ label: 'Growth-rate agreement', val: (100 * frac).toFixed(0) + '%', cls: frac >= 0.7 ? 'ok' : frac >= 0.4 ? 'warn' : 'bad', detail: `${within}/${S.muPairs.length} predicted µ within 2× of measured (median pred/meas ${med ? med.toFixed(2) : 'n/a'}×). FBA maximises growth, so it sets an upper bound and tends to over-predict.`, w: frac });
  }
  if (S.spec && S.spec.tested) { const mcc = spectrumMCC(S.spec); const norm = (mcc + 1) / 2;
    dims.push({ label: 'Substrate spectrum (MCC)', val: mcc.toFixed(2), cls: mcc >= 0.6 ? 'ok' : mcc >= 0.3 ? 'warn' : 'bad', detail: `${S.spec.TP + S.spec.TN}/${S.spec.tested} correct on grows-on/no-grow${S.spec.capped ? ` (${S.spec.tested}-substrate balanced sample)` : ''}; ${S.spec.FN} false-neg (gap-fill), ${S.spec.FP} false-pos (over-permissive)`, w: norm }); }
  if (S.secTot) { const a = S.secRecall; dims.push({ label: 'Secretion recall', val: (100 * a).toFixed(0) + '%', cls: a >= 0.7 ? 'ok' : a >= 0.4 ? 'warn' : 'bad', detail: `${S.secHit}/${S.secTot} measured secretion products the model can produce under the reported condition`, w: a }); }
  if (!dims.length) { host.appendChild(el('div', 'ac-empty', 'No validation dimension could be scored (no simulable medium and no spectrum).')); return; }
  const overall = dims.reduce((s, d) => s + d.w, 0) / dims.length;
  const grade = overall >= 0.85 ? ['A', 'ok'] : overall >= 0.7 ? ['B', 'ok'] : overall >= 0.55 ? ['C', 'warn'] : overall >= 0.4 ? ['D', 'warn'] : ['F', 'bad'];
  const gcol = { ok: '#15803D', warn: '#B45309', bad: '#DC2626' }[grade[1]];
  const card = el('div', 'ac-card'); card.style.cssText = 'border:2px solid ' + gcol + '33';
  card.appendChild(el('h3', '', `GEM validation scorecard — <i>${esc(sp)}</i>${detTok ? ' · ' + esc(detTok) : ''}`));
  const head = el('div', ''); head.style.cssText = 'display:flex;align-items:center;gap:16px;margin:8px 0 4px;flex-wrap:wrap';
  head.innerHTML = `<div style="font-size:44px;font-weight:800;line-height:1;color:${gcol}">${grade[0]}</div><div style="font-size:13px;color:var(--ink-2)">Composite <b>${(100 * overall).toFixed(0)}%</b> across ${dims.length} validation dimension${dims.length > 1 ? 's' : ''}, from <b>${S.sampled}</b>${S.withMed > S.sampled ? ' of ' + S.withMed : ''} growth condition${S.sampled > 1 ? 's' : ''}${S.strainCond ? ` (${S.strainCond} strain-specific)` : ''} + the substrate spectrum. Each dimension is unweighted; the letter reflects the mean.</div>`;
  card.appendChild(head);
  dims.forEach(d => { const cls = { ok: '#15803D', warn: '#B45309', bad: '#DC2626' }[d.cls];
    const row = el('div', ''); row.style.cssText = 'display:grid;grid-template-columns:180px 70px 1fr;gap:12px;align-items:baseline;padding:8px 0;border-top:1px solid var(--line)';
    row.innerHTML = `<div style="font-size:13px;font-weight:600">${esc(d.label)}</div><div style="font-size:16px;font-weight:700;color:${cls}">${d.val}</div><div style="font-size:12px;color:var(--ink-3)">${d.detail}</div>`;
    card.appendChild(row); });
  card.appendChild(el('div', 'ac-interp', `A single number hides detail — open each validation below for the per-condition table, confusion matrix, and the exact false-positive / false-negative substrate lists (the actionable curation targets). ${S.withMed > S.sampled ? `Only the first ${S.sampled} conditions with a simulable medium were scored for responsiveness; ` : ''}growth-feasibility counts a condition as passing if the model grows at all on the reported medium.`));
  host.appendChild(card);
}
function renderSpectrumPanel(host, sp, detTok) {
  const spec = spectrumFor(sp, detTok);
  const card = el('div', 'ac-card'); card.style.cssText = 'margin-top:12px;border:2px solid #dfe6f5';
  if (!spec || (!spec.pos.size && !spec.neg.size)) { card.innerHTML = `<h3>Substrate-utilisation validation</h3><div class="sub">No grows-on/no-grow spectrum in GrowthDB for <i>${esc(sp)}</i> — can't run the substrate confusion matrix for this organism.</div>`; host.appendChild(card); return; }
  card.innerHTML = `<h3>Substrate-utilisation validation <span class="note">— the standard GEM Biolog check</span></h3>
    <div class="sub">Sweeps every substrate GrowthDB knows <i>${esc(sp)}</i> ${spec.strainN ? '(and strain <code>' + esc(detTok) + '</code>) ' : ''}grows / doesn't grow on (<b>${fmt(spec.pos.size)}</b> grows-on, <b>${fmt(spec.neg.size)}</b> no-grow): sets it as the sole carbon source, runs FBA, and builds a confusion matrix. False negatives = gap-fill targets; false positives = over-permissive reactions.</div>
    <label style="display:inline-flex;align-items:center;gap:6px;font-size:12.5px;color:var(--ink-2);margin:10px 0"><input type="checkbox" id="spec-aer" checked> aerobic (open O₂)</label>`;
  const btn = el('button', 'ac-btn primary', '▶ Run substrate-spectrum validation'); btn.style.marginLeft = '12px';
  const out = el('div', ''); out.style.marginTop = '12px';
  btn.onclick = () => { btn.disabled = true; btn.textContent = 'Running…'; runSpectrumValidation(sp, detTok, out, card.querySelector('#spec-aer').checked).then(() => { btn.disabled = false; btn.textContent = '↻ Re-run'; }); };
  card.appendChild(btn); card.appendChild(out); host.appendChild(card);
}

function renderReport() {
  const s = $('stage-report'); s.innerHTML = ''; const c = RESULT.counts;
  const allIss = [...RESULT.idIssues, ...RESULT.discrep, ...RESULT.structure, ...RESULT.massIss, ...RESULT.chargeIss];
  const appr = allIss.filter(i => APPROVED[i.id] === true).length, rej = allIss.filter(i => APPROVED[i.id] === false).length, pend = allIss.length - appr - rej;
  s.appendChild(el('div', 'ac-sh', `<h2>Report &amp; export</h2><p>Your supervised curation summary. Approved fixes are written into the exported model; rejected and pending items are left unchanged but recorded in the report.</p>`));
  s.appendChild(el('div', 'ac-kpis', kpi(allIss.length, 'total findings', 'info') + kpi(appr, 'approved', 'ok') + kpi(rej, 'rejected', 'bad') + kpi(pend, 'pending', pend ? 'warn' : 'ok')));
  const card = el('div', 'ac-card'); card.appendChild(el('h3', '', 'Download the curated model & report'));
  card.appendChild(el('div', 'sub', `${appr} approved fix${appr === 1 ? '' : 'es'} applied. Choose a format:`));
  const exp = el('div', 'ac-export');
  const dl = (label, sub, tag, fn) => { const d = el('div', 'ac-dl'); d.innerHTML = `<div class="ic">${tag}</div><div><b>${label}</b><span>${sub}</span></div>`; d.onclick = fn; return d; };
  exp.append(
    dl('Curated COBRA JSON', 'cobrapy-ready model', 'JSON', () => download(exportJson(), MODEL.id + '_curated.json', 'application/json')),
    dl('Curated SBML', 'SBML L3 FBC', 'XML', () => download(exportSBML(), MODEL.id + '_curated.xml', 'application/xml')),
    dl('Curated MATLAB (.mat)', 'COBRA Toolbox model struct', 'MAT', () => download(exportMAT(), MODEL.id + '_curated.mat', 'application/octet-stream')),
    dl('Curation report', 'full findings + decisions', 'MD', () => download(exportReport(allIss), MODEL.id + '_curation_report.md', 'text/markdown')),
  );
  card.appendChild(exp); s.appendChild(card);
  navCount('report', appr, appr ? 'ok' : 'warn');
}

/* ---------------- export ---------------- */
function applyApproved(model) {
  const m = JSON.parse(JSON.stringify({ mets: model.mets.map(x => ({ id: x.id, name: x.name, formula: x.formula, charge: x.charge, compartment: x.compartment })), rxns: model.rxns.map(r => ({ id: r.id, name: r.name, s: r.s, lb: r.lb, ub: r.ub, gpr: r.gpr })), id: model.id, name: model.name, genes: model.genes }));
  const log = { renames: 0, merges: 0, protons: 0 };
  const rxnById = {}; m.rxns.forEach(r => { rxnById[r.id] = r; });
  // 1) proton/water balancing fixes (original ids still intact)
  [...RESULT.massIss, ...RESULT.chargeIss].forEach(i => { if (APPROVED[i.id] === true && i.apply && i.apply.proton) {
    const r = rxnById[i.apply.rxn]; if (!r) return; const p = i.apply.proton;
    if (p.h && i.apply.h) { r.s[i.apply.h] = (r.s[i.apply.h] || 0) + p.h; if (!r.s[i.apply.h]) delete r.s[i.apply.h]; }
    if (p.h2o && i.apply.h2o) { r.s[i.apply.h2o] = (r.s[i.apply.h2o] || 0) + p.h2o; if (!r.s[i.apply.h2o]) delete r.s[i.apply.h2o]; }
    log.protons++;
  } });
  // 2) metabolite merges (duplicate species -> one canonical id)
  const metRe = {};
  RESULT.discrep.filter(d => d.kind === 'dupmet').forEach(d => { if (APPROVED[d.id] !== true) return; d.apply.merge.forEach(mid => metRe[mid] = d.apply.into); log.merges++; });
  if (Object.keys(metRe).length) {
    const seen = new Set(); m.mets = m.mets.filter(x => { const nid = metRe[x.id] || x.id; if (seen.has(nid)) return false; seen.add(nid); x.id = nid; return true; });
    m.rxns.forEach(r => { const ns = {}; Object.entries(r.s).forEach(([mid, c]) => { const nid = metRe[mid] || mid; ns[nid] = (ns[nid] || 0) + c; }); r.s = ns; });
  }
  // 3) reaction merges (keep first, drop the rest)
  const drop = new Set();
  RESULT.discrep.filter(d => d.kind === 'duprxn').forEach(d => { if (APPROVED[d.id] !== true) return; d.apply.mergeRxn.slice(1).forEach(rid => drop.add(rid)); log.merges++; });
  if (drop.size) m.rxns = m.rxns.filter(r => !drop.has(r.id));
  // 4) identifier renames
  const ren = {};
  RESULT.idIssues.forEach(i => { if (APPROVED[i.id] === true && i.apply) { if (i.apply.met && !metRe[i.apply.met]) ren['m:' + i.apply.met] = i.apply.newId; if (i.apply.rxn && !drop.has(i.apply.rxn)) ren['r:' + i.apply.rxn] = i.apply.newId; log.renames++; } });
  const seenM = new Set(); m.mets = m.mets.filter(x => { const n = ren['m:' + x.id] || x.id; if (seenM.has(n)) return false; seenM.add(n); x.id = n; return true; });
  m.rxns.forEach(r => { const n = ren['r:' + r.id]; if (n) r.id = n; const ns = {}; Object.entries(r.s).forEach(([mid, c]) => { ns[ren['m:' + mid] || mid] = c; }); r.s = ns; });
  m._log = log;
  return m;
}
function exportJson() {
  const m = applyApproved(MODEL);
  const out = { id: m.id, name: m.name, version: '1', metabolites: m.mets.map(x => ({ id: x.id, name: x.name, formula: x.formula || '', charge: x.charge == null ? 0 : x.charge, compartment: x.compartment })),
    reactions: m.rxns.map(r => ({ id: r.id, name: r.name, metabolites: r.s, lower_bound: r.lb == null ? -1000 : r.lb, upper_bound: r.ub == null ? 1000 : r.ub, gene_reaction_rule: r.gpr || '' })),
    genes: [], compartments: {}, notes: { curated_by: 'GEM Autocurator', date: new Date().toISOString().slice(0, 10) } };
  return JSON.stringify(out, null, 1);
}
function exportSBML() {
  const m = applyApproved(MODEL);
  const bounds = new Map(); const bid = v => { const k = 'B' + String(v).replace(/[.\-]/g, '_'); bounds.set(k, v); return k; };
  const rx = m.rxns.map(r => { const lb = bid(r.lb == null ? -1000 : r.lb), ub = bid(r.ub == null ? 1000 : r.ub);
    const reac = Object.entries(r.s).filter(([, c]) => c < 0).map(([mid, c]) => `<speciesReference species="${esc(mid)}" stoichiometry="${Math.abs(c)}" constant="true"/>`).join('');
    const prod = Object.entries(r.s).filter(([, c]) => c > 0).map(([mid, c]) => `<speciesReference species="${esc(mid)}" stoichiometry="${c}" constant="true"/>`).join('');
    return `<reaction id="${esc(r.id)}" name="${esc(r.name)}" reversible="${r.lb < 0}" fast="false" fbc:lowerFluxBound="${lb}" fbc:upperFluxBound="${ub}"><listOfReactants>${reac}</listOfReactants><listOfProducts>${prod}</listOfProducts></reaction>`; }).join('\n');
  const sp = m.mets.map(x => `<species id="${esc(x.id)}" name="${esc(x.name)}" compartment="${esc(x.compartment)}" hasOnlySubstanceUnits="false" boundaryCondition="false" constant="false" fbc:charge="${x.charge == null ? 0 : x.charge}" fbc:chemicalFormula="${esc(x.formula || '')}"/>`).join('\n');
  const comps = [...new Set(m.mets.map(x => x.compartment))].map(c => `<compartment id="${esc(c)}" constant="true"/>`).join('');
  const params = [...bounds].map(([k, v]) => `<parameter id="${k}" value="${v}" constant="true"/>`).join('');
  return `<?xml version="1.0" encoding="UTF-8"?>\n<sbml xmlns="http://www.sbml.org/sbml/level3/version1/core" xmlns:fbc="http://www.sbml.org/sbml/level3/version1/fbc/version2" level="3" version="1" fbc:required="false">\n<model id="${esc(m.id)}" name="${esc(m.name)} (autocurated)" fbc:strict="true">\n<listOfCompartments>${comps}</listOfCompartments>\n<listOfParameters>${params}</listOfParameters>\n<listOfSpecies>\n${sp}\n</listOfSpecies>\n<listOfReactions>\n${rx}\n</listOfReactions>\n</model>\n</sbml>`;
}
function exportReport(allIss) {
  const c = RESULT.counts; const d = new Date().toISOString().slice(0, 10);
  const sect = (name, arr) => { if (!arr.length) return `\n## ${name}\n\n_none_\n`; return `\n## ${name} (${arr.length})\n\n` + arr.map(i => `- [${APPROVED[i.id] === true ? 'x' : ' '}] **${i.title.replace(/<[^>]+>/g, '')}** — ${(i.note || '').replace(/<[^>]+>/g, '')}${i.from ? ` (\`${i.from}\` → \`${i.to}\`)` : ''}`).join('\n') + '\n'; };
  return `# GEM Autocuration report — ${MODEL.name}\n\nGenerated ${d} · ${MODEL.format} · ${c.mets} metabolites · ${c.rxns} reactions · ${c.genes} genes\n\n## Summary\n\n| Check | Finding |\n|---|---|\n| Metabolites → BiGG / BiGG-like / unmapped | ${c.mBigg} / ${c.mLike} / ${c.mUn} |\n| Reactions → BiGG / BiGG-like / unmapped | ${c.rBigg} / ${c.rLike} / ${c.rUn} |\n| Duplicate metabolites / reactions | ${c.dupMet} / ${c.dupRxn} |\n| Dead-end metabolites | ${c.dead} |\n| Mass-imbalanced reactions | ${c.mass} |\n| Charge-imbalanced reactions | ${c.charge} |\n| Findings approved / rejected / pending | ${allIss.filter(i => APPROVED[i.id] === true).length} / ${allIss.filter(i => APPROVED[i.id] === false).length} / ${allIss.filter(i => APPROVED[i.id] == null).length} |\n${sect('Identifier renames', RESULT.idIssues)}${sect('Discrepancies (duplicates)', RESULT.discrep)}${sect('Dead-end metabolites', RESULT.structure)}${sect('Mass imbalance', RESULT.massIss)}${sect('Charge imbalance', RESULT.chargeIss)}\n---\nCurated with the GEM Autocurator against BiGG + KEGG + ModelSEED + MetaNetX + ChEBI + RHEA.\n`;
}
function download(content, name, type) { const b = content instanceof Blob ? content : new Blob([content], { type }); const u = URL.createObjectURL(b); const a = el('a'); a.href = u; a.download = name; a.click(); setTimeout(() => URL.revokeObjectURL(u), 2000); }

/* ---- MATLAB v5 (.mat) COBRA-model writer ---- */
function _cat(...a) { let n = 0; a.forEach(x => n += x.length); const o = new Uint8Array(n); let p = 0; a.forEach(x => { o.set(x, p); p += x.length; }); return o; }
function _pad8(a) { const r = a.length % 8; return r ? _cat(a, new Uint8Array(8 - r)) : a; }
function _tag(type, n) { const b = new Uint8Array(8); const d = new DataView(b.buffer); d.setUint32(0, type, true); d.setUint32(4, n, true); return b; }
function _elem(type, data) { return _cat(_tag(type, data.length), _pad8(data)); }
function _i32(nums) { const b = new Uint8Array(nums.length * 4); const d = new DataView(b.buffer); nums.forEach((n, i) => d.setInt32(i * 4, n, true)); return b; }
function _dbl(nums) { const b = new Uint8Array(nums.length * 8); const d = new DataView(b.buffer); nums.forEach((n, i) => d.setFloat64(i * 8, n, true)); return b; }
function _str(s) { const b = new Uint8Array(s.length); for (let i = 0; i < s.length; i++) b[i] = s.charCodeAt(i) & 0xff; return b; }
function _u16(codes) { const b = new Uint8Array(codes.length * 2); const d = new DataView(b.buffer); codes.forEach((c, i) => d.setUint16(i * 2, c, true)); return b; }
function _mx(cls, dims, name, content) { const flags = _elem(6, _i32([cls, 0])); const dim = _elem(5, _i32(dims)); const nm = _elem(1, _str(name)); const body = _cat(flags, dim, nm, content); return _cat(_tag(14, body.length), body); }
function mDbl(name, arr, rows, cols) { return _mx(6, [rows, cols], name, _elem(9, _dbl(Array.from(arr)))); }
function mChar(name, str) { return _mx(4, [str ? 1 : 0, str.length], name, _elem(4, _u16(Array.from(str).map(c => c.charCodeAt(0))))); }
function mCell(name, strs) { const cells = strs.map(s => mChar('', s)); return _mx(1, [strs.length, 1], name, cells.length ? _cat(...cells) : new Uint8Array(0)); }
function mStruct(name, fields) { const L = 32; const fnLen = _elem(5, _i32([L])); const names = new Uint8Array(fields.length * L); fields.forEach((f, i) => { const nb = _str(f[0]); names.set(nb.subarray(0, L - 1), i * L); }); const namesEl = _elem(1, names); const vals = fields.length ? _cat(...fields.map(f => f[1])) : new Uint8Array(0); return _mx(2, [1, 1], name, _cat(fnLen, namesEl, vals)); }
function exportMAT() {
  const m = applyApproved(MODEL);
  const rxns = m.rxns.map(r => r.id), mets = m.mets.map(x => x.id), metIdx = {}; mets.forEach((id, i) => metIdx[id] = i);
  const nM = mets.length, nR = rxns.length;
  const S = new Float64Array(nM * nR); m.rxns.forEach((r, j) => Object.entries(r.s).forEach(([mid, c]) => { const i = metIdx[mid]; if (i != null) S[i + j * nM] = c; }));
  const fields = [
    ['rxns', mCell('', rxns)], ['mets', mCell('', mets)], ['S', mDbl('', S, nM, nR)],
    ['lb', mDbl('', m.rxns.map(r => r.lb == null ? -1000 : r.lb), nR, 1)], ['ub', mDbl('', m.rxns.map(r => r.ub == null ? 1000 : r.ub), nR, 1)],
    ['c', mDbl('', m.rxns.map(r => (/biomass/i.test(r.id) || /biomass/i.test(r.name || '')) ? 1 : 0), nR, 1)], ['b', mDbl('', mets.map(() => 0), nM, 1)],
    ['rev', mDbl('', m.rxns.map(r => r.lb < 0 ? 1 : 0), nR, 1)],
    ['rxnNames', mCell('', m.rxns.map(r => r.name || r.id))], ['metNames', mCell('', m.mets.map(x => x.name || x.id))],
    ['metFormulas', mCell('', m.mets.map(x => x.formula || ''))], ['grRules', mCell('', m.rxns.map(r => r.gpr || ''))],
    ['description', mChar('', (m.name || m.id || 'model') + ' — autocurated')],
  ];
  const struct = mStruct('model', fields);
  const head = new Uint8Array(128); const desc = 'MATLAB 5.0 MAT-file, Platform: GEM Autocurator, Created: ' + new Date().toISOString().slice(0, 10);
  for (let i = 0; i < 116; i++) head[i] = i < desc.length ? desc.charCodeAt(i) : 0x20;
  new DataView(head.buffer).setUint16(124, 0x0100, true); head[126] = 0x49; head[127] = 0x4D;
  return new Blob([head, struct], { type: 'application/octet-stream' });
}

/* ---------------- flow ---------------- */
async function ingest(text, filename) {
  await loadRef();
  $('ac-load').style.display = 'flex'; $('ac-load-msg').textContent = 'Parsing & curating…';
  await new Promise(r => setTimeout(r, 20));
  try {
    MODEL = parseModel(text, filename);
    RESULT = curate(MODEL);
  } catch (e) { $('ac-load').style.display = 'none'; alert('Could not curate this model:\n' + e.message); console.error(e); return; }
  $('ac-load').style.display = 'none';
  const chip = $('ac-modelchip'); chip.style.display = ''; chip.innerHTML = `<b>${esc(MODEL.name)}</b> · ${MODEL.format} · ${fmt(MODEL.mets.length)} mets · ${fmt(MODEL.rxns.length)} rxns`;
  enableNav();
  navCount('discrep', RESULT.discrep.length, RESULT.discrep.length ? 'bad' : 'ok');
  goStage('ids');
}
function initUpload() {
  const drop = $('ac-drop'), file = $('ac-file');
  drop.onclick = () => file.click();
  file.onchange = () => { const f = file.files[0]; if (f) f.text().then(t => ingest(t, f.name.replace(/\.[^.]+$/, ''))); };
  ['dragover', 'dragenter'].forEach(ev => drop.addEventListener(ev, e => { e.preventDefault(); drop.classList.add('drag'); }));
  ['dragleave', 'drop'].forEach(ev => drop.addEventListener(ev, e => { e.preventDefault(); drop.classList.remove('drag'); }));
  drop.addEventListener('drop', e => { const f = e.dataTransfer.files[0]; if (f) f.text().then(t => ingest(t, f.name.replace(/\.[^.]+$/, ''))); });
  document.querySelectorAll('.ac-ex').forEach(b => b.onclick = () => fetch('data/e_coli_core.json').then(r => r.text()).then(t => ingest(t, 'e_coli_core')).catch(() => alert('Demo model not available yet.')));
  document.querySelectorAll('.ac-nav-item').forEach(b => b.onclick = () => { if (!b.disabled) goStage(b.dataset.stage); });
  $('ac-theme').onclick = () => { const t = document.documentElement.getAttribute('data-theme') === 'dark' ? 'light' : 'dark'; document.documentElement.setAttribute('data-theme', t); try { localStorage.setItem('gem-browser-theme', t); } catch (e) {} };
}
initUpload();
// deep-link / headless demo auto-load
if (/[?#]demo/.test(location.search + location.hash)) {
  const st = (/stage=(\w+)/.exec(location.search + location.hash) || [])[1];
  window.addEventListener('load', () => fetch('data/e_coli_core.json').then(r => r.text()).then(t => ingest(t, 'e_coli_core')).then(() => { if (st) goStage(st); }));
}
