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
try { Object.defineProperty(window, '__ac', { get: () => ({ MODEL, RESULT, REF, APPROVED, exportMAT, applyApproved }) }); } catch (e) {}  // debug/headless hook

/* ---------------- reference maps ---------------- */
async function loadRef() {
  if (REF) return REF;
  $('ac-load').style.display = 'flex';
  const setMsg = m => { const e = $('ac-load-msg'); if (e) e.textContent = m; };
  setMsg('Loading identifier maps (BiGG + KEGG + ModelSEED + MetaNetX)…');
  const [met, rxn, props] = await Promise.all([
    fetch('data/metabolite_map.json').then(r => r.json()),
    fetch('data/reaction_map.json').then(r => r.json()),
    fetch('data/bigg_met_props.json').then(r => r.json()),
  ]);
  REF = { met, rxn, props };
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
function canonMet(m) {
  const M = REF.met;
  let hit = M['bigg:' + baseId(m.id)];                       // direct BiGG id
  if (!hit) for (const [ns, v] of Object.entries(m.anno || {})) { hit = M[ns + ':' + v]; if (hit) break; }
  if (!hit && m.name) hit = M['name:' + m.name.toLowerCase().replace(/[^a-z0-9]/g, '')];
  return hit || null;
}
function canonRxn(r) {
  const R = REF.rxn;
  let hit = R['bigg:' + baseId(r.id)];
  if (!hit) hit = R['old:' + baseId(r.id)];
  if (!hit) for (const [ns, v] of Object.entries(r.anno || {})) { hit = R[ns + ':' + v]; if (hit) break; }
  return hit || null;
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

function curate(model) {
  // --- identifiers ---
  let mBigg = 0, mLike = 0, mUn = 0;
  model.mets.forEach(m => { const c = canonMet(m); m.canon = c; if (!c) mUn++; else if (c.biggr) mBigg++; else mLike++; });
  let rBigg = 0, rLike = 0, rUn = 0;
  model.rxns.forEach(r => { if (isExchange(r)) { r.canon = { bigg: baseId(r.id), biggr: true, exch: true }; return; } const c = canonRxn(r); r.canon = c; if (!c) rUn++; else if (c.biggr) rBigg++; else rLike++; });
  const idIssues = [];
  model.mets.forEach(m => { if (m.canon && m.canon.bigg !== baseId(m.id)) idIssues.push({ id: 'mid_' + m.id, cat: 'ids', sev: m.canon.biggr ? 'info' : 'warn', kind: 'met',
    title: `Metabolite <code>${esc(m.id)}</code>`, from: baseId(m.id), to: m.canon.bigg + '_' + m.compartment, biggr: m.canon.biggr,
    note: m.canon.biggr ? 'canonical BiGG id' : 'BiGG-like id (compound not in BiGG; canonical across KEGG/ModelSEED/MetaNetX)', apply: { met: m.id, newId: m.canon.bigg + '_' + m.compartment } }); });
  model.rxns.forEach(r => { if (r.canon && !r.canon.exch && r.canon.bigg !== baseId(r.id)) idIssues.push({ id: 'rid_' + r.id, cat: 'ids', sev: r.canon.biggr ? 'info' : 'warn', kind: 'rxn',
    title: `Reaction <code>${esc(r.id)}</code>`, from: baseId(r.id), to: r.canon.bigg, biggr: r.canon.biggr,
    note: r.canon.biggr ? 'canonical BiGG reaction id' : 'BiGG-like id (reaction not in BiGG; canonical across KEGG/ModelSEED/RHEA)', apply: { rxn: r.id, newId: r.canon.bigg } }); });

  // --- discrepancies: distinct model ids that canonicalise to the SAME BiGG id ---
  const dup = {};
  model.mets.forEach(m => { if (!m.canon) return; const key = m.canon.bigg + '@' + m.compartment; (dup[key] = dup[key] || []).push(m); });
  const discrep = [];
  Object.entries(dup).forEach(([key, arr]) => { if (arr.length > 1) discrep.push({ id: 'dupm_' + key, cat: 'discrep', sev: 'bad', kind: 'dupmet',
    title: `${arr.length} metabolites collapse to one`, ids: arr.map(m => m.id), to: key.replace('@', '_'),
    note: `<code>${arr.map(m => esc(m.id)).join('</code>, <code>')}</code> all resolve to <span class="to">${esc(key.replace('@', '_'))}</span> — likely duplicate species with different names/ids.`,
    apply: { merge: arr.map(m => m.id), into: key.split('@')[0] + '_' + key.split('@')[1] } }); });
  const rdup = {};
  model.rxns.forEach(r => { if (!r.canon || r.canon.exch) return; (rdup[r.canon.bigg] = rdup[r.canon.bigg] || []).push(r); });
  Object.entries(rdup).forEach(([b, arr]) => { if (arr.length > 1) discrep.push({ id: 'dupr_' + b, cat: 'discrep', sev: 'bad', kind: 'duprxn',
    title: `${arr.length} reactions collapse to one`, ids: arr.map(r => r.id), to: b,
    note: `<code>${arr.map(r => esc(r.id)).join('</code>, <code>')}</code> all resolve to <span class="to">${esc(b)}</span> — duplicate reactions.`, apply: { mergeRxn: arr.map(r => r.id), into: b } }); });
  // name conflicts: same canon, different display names
  const nmeBy = {};
  model.mets.forEach(m => { if (!m.canon) return; (nmeBy[m.canon.bigg] = nmeBy[m.canon.bigg] || new Set()).add((m.name || '').trim()); });

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

  // --- mass & charge imbalance ---
  const massIss = [], chargeIss = [];
  const byId = {}; model.mets.forEach(m => { byId[m.id] = m; });
  const F = (mid) => { const m = byId[mid]; let f = (m && m.formula) || ''; if (!f && m && m.canon && REF.props[m.canon.bigg]) f = REF.props[m.canon.bigg].f; return parseFormula(f); };
  const CH = (mid) => { const m = byId[mid]; if (m && m.charge != null) return m.charge; if (m && m.canon && REF.props[m.canon.bigg] && REF.props[m.canon.bigg].c != null) return REF.props[m.canon.bigg].c; return null; };
  const compById = {}; model.mets.forEach(m => { compById[m.id] = m.compartment; });
  const hIn = {}, h2oIn = {}; model.mets.forEach(m => { const b = m.canon ? m.canon.bigg : baseId(m.id); if (b === 'h') hIn[m.compartment] = m.id; if (b === 'h2o') h2oIn[m.compartment] = m.id; });
  model.rxns.forEach(r => { if (isExchange(r) || isBiomass(r)) return;
    const bal = {}; let known = true;
    Object.entries(r.s).forEach(([mid, c]) => { const f = F(mid); if (!f || !Object.keys(f).length) known = false; Object.entries(f).forEach(([e, n]) => bal[e] = (bal[e] || 0) + c * n); });
    let cbal = 0, ck = true; Object.entries(r.s).forEach(([mid, c]) => { const q = CH(mid); if (q == null) ck = false; else cbal += c * q; });
    const massOff = known ? Object.entries(bal).filter(([e, v]) => Math.abs(v) > 1e-6) : [];
    const chgOff = ck && Math.abs(cbal) > 1e-6;
    if (!massOff.length && !chgOff) return;
    const comp = compById[Object.keys(r.s)[0]] || 'c';
    const fix = (known && ck) ? protonWaterFix(bal, cbal) : null;
    const canApply = fix && (!fix.h || hIn[comp]) && (!fix.h2o || h2oIn[comp]);
    const fixTxt = fix ? (canApply ? ` — <b>fix:</b> add ${fixStr(fix, comp)}` : ' — a proton/water fix balances it, but the model lacks h/h2o in this compartment') : ' — not a proton/water imbalance; needs manual curation';
    const applyObj = canApply ? { rxn: r.id, proton: fix, h: hIn[comp], h2o: h2oIn[comp] } : null;
    if (massOff.length) massIss.push({ id: 'mass_' + r.id, cat: 'charge', sub: 'mass', sev: 'bad', kind: 'mass',
      title: `Mass imbalance in <code>${esc(r.id)}</code>`, note: 'unbalanced: ' + massOff.map(([e, v]) => `${e}${v > 0 ? '+' : ''}${(+v.toFixed(2))}`).join(', ') + fixTxt, apply: applyObj });
    else chargeIss.push({ id: 'chg_' + r.id, cat: 'charge', sub: 'charge', sev: 'warn', kind: 'charge',
      title: `Charge imbalance in <code>${esc(r.id)}</code>`, note: `net charge ${cbal > 0 ? '+' : ''}${(+cbal.toFixed(2))}` + fixTxt, apply: applyObj });
  });

  return { idIssues, discrep, structure, orphan, massIss, chargeIss,
    counts: { mBigg, mLike, mUn, rBigg, rLike, rUn, mets: model.mets.length, rxns: model.rxns.length, genes: model.genes,
      exch: model.rxns.filter(isExchange).length, dead: structure.length, mass: massIss.length, charge: chargeIss.length,
      dupMet: discrep.filter(d => d.kind === 'dupmet').length, dupRxn: discrep.filter(d => d.kind === 'duprxn').length } };
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
  const donut = (id, title, vals) => window.Plotly.newPlot(id, [{ type: 'pie', hole: .62, labels: ['BiGG', 'BiGG-like', 'unmapped'], values: vals, marker: { colors: ['#15803D', '#2563EB', '#B45309'] }, textinfo: 'percent', textfont: { size: 10 }, hovertemplate: '%{label}: %{value}<extra></extra>', sort: false }], Object.assign({}, base, { title: { text: title, font: { size: 12 } }, showlegend: false, annotations: [{ text: fmt(vals[0] + vals[1] + vals[2]), x: .5, y: .5, showarrow: false, font: { size: 15, color: ink } }] }), cfg);
  donut('viz-met', 'Metabolite ids', [c.mBigg, c.mLike, c.mUn]);
  donut('viz-rxn', 'Reaction ids', [c.rBigg, c.rLike, c.rUn]);
  const labels = ['Renames', 'Dup metabolites', 'Dup reactions', 'Dead-ends', 'Mass imbalance', 'Charge imbalance'];
  const vals = [RESULT.idIssues.length, c.dupMet, c.dupRxn, c.dead, c.mass, c.charge];
  const cols = ['#2563EB', '#DC2626', '#DC2626', '#B45309', '#DC2626', '#B45309'];
  window.Plotly.newPlot('viz-iss', [{ type: 'bar', orientation: 'h', y: labels.slice().reverse(), x: vals.slice().reverse(), marker: { color: cols.slice().reverse() }, text: vals.slice().reverse().map(v => v || ''), textposition: 'outside', textfont: { size: 10 }, hovertemplate: '%{y}: %{x}<extra></extra>' }], Object.assign({}, base, { title: { text: 'Findings', font: { size: 12 } }, margin: { l: 110, r: 24, t: 26, b: 24 }, xaxis: { gridcolor: dark ? 'rgba(255,255,255,.08)' : '#EEF2F8', zeroline: false }, yaxis: { tickfont: { size: 10 } } }), cfg);
}

function renderIds() {
  const s = $('stage-ids'); s.innerHTML = ''; const c = RESULT.counts;
  s.appendChild(el('div', 'ac-sh', `<h2>Identifier mapping</h2><p>Every metabolite and reaction identifier is resolved against <b>BiGG</b> — and, when a compound or reaction exists only in KEGG / ModelSEED / MetaNetX / ChEBI / RHEA, a deterministic <b>BiGG-like</b> id so your model stays internally consistent and portable.</p>`));
  s.appendChild(el('div', 'ac-kpis', kpi(c.mBigg, 'metabolites → BiGG', 'ok') + kpi(c.mLike, 'metabolites → BiGG-like', 'info') + kpi(c.mUn, 'metabolites unmapped', c.mUn ? 'warn' : 'ok') + kpi(c.rBigg, 'reactions → BiGG', 'ok') + kpi(c.rLike, 'reactions → BiGG-like', 'info') + kpi(c.rUn, 'reactions unmapped', c.rUn ? 'warn' : 'ok')));
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
  s.appendChild(el('div', 'ac-interp', `<b>Interpretation.</b> ${pct}% of metabolites resolved to a canonical id (${c.mBigg} exact BiGG, ${c.mLike} BiGG-like). ${c.mUn ? `<b>${c.mUn}</b> metabolites and <b>${c.rUn}</b> reactions could not be resolved from their id or annotations — they keep their original ids and are listed under Discrepancies for review.` : 'Full coverage — every id maps to a curated reference.'}`));
  navCount('ids', RESULT.idIssues.length, RESULT.idIssues.length ? 'warn' : 'ok');
}
function renderDiscrep() {
  const s = $('stage-discrep'); s.innerHTML = '';
  s.appendChild(el('div', 'ac-sh', `<h2>Discrepancies</h2><p>The most consequential curation: when two different identifiers in your model refer to the <b>same</b> compound or reaction. Left uncurated, these split flux, break mass balance and inflate the network. Each is resolved to one canonical entity.</p>`));
  const c = RESULT.counts;
  s.appendChild(el('div', 'ac-kpis', kpi(c.dupMet, 'duplicate metabolites', c.dupMet ? 'bad' : 'ok') + kpi(c.dupRxn, 'duplicate reactions', c.dupRxn ? 'bad' : 'ok') + kpi(c.mUn + c.rUn, 'unresolved ids', (c.mUn + c.rUn) ? 'warn' : 'ok')));
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
  s.appendChild(el('div', 'ac-sh', `<h2>Mass &amp; charge balance</h2><p>Every non-exchange reaction is checked for element and charge conservation. Charges follow the protonation state at your simulation pH.</p>`));
  const phbar = el('div', 'ac-card');
  phbar.innerHTML = `<h3>Simulation pH</h3><div class="sub">BiGG charges are defined at pH 7.2. Set the pH you will simulate at; a full pKa-microspecies recharge is on the roadmap — for now charge balance uses the reference protonation state.</div>
    <div style="display:flex;align-items:center;gap:14px"><input type="range" id="ph-slider" min="4" max="9" step="0.1" value="7.2" style="flex:1"><span id="ph-val" class="tabular" style="font-size:20px;font-weight:700;color:var(--primary-2);min-width:52px">7.2</span></div>`;
  s.appendChild(phbar);
  const sl = phbar.querySelector('#ph-slider'); sl.oninput = () => { phbar.querySelector('#ph-val').textContent = (+sl.value).toFixed(1); };
  s.appendChild(el('div', 'ac-kpis', kpi(c.mass, 'mass-imbalanced reactions', c.mass ? 'bad' : 'ok') + kpi(c.charge, 'charge-imbalanced reactions', c.charge ? 'warn' : 'ok')));
  const mc = el('div', 'ac-card'); mc.appendChild(el('h3', '', 'Mass imbalance')); mc.appendChild(el('div', 'sub', 'Elements do not conserve across the reaction (using model or BiGG formulas).'));
  issueList(mc, RESULT.massIss, 'Every reaction with known formulas conserves mass.'); s.appendChild(mc);
  const cc = el('div', 'ac-card'); cc.appendChild(el('h3', '', 'Charge imbalance')); cc.appendChild(el('div', 'sub', 'Net charge is non-zero — usually a missing proton at this pH.'));
  issueList(cc, RESULT.chargeIss, 'Charge conserves across every reaction with known charges.'); s.appendChild(cc);
  navCount('charge', c.mass + c.charge, (c.mass) ? 'bad' : (c.charge ? 'warn' : 'ok'));
}
function roadmap(stage, title, items) {
  const s = $('stage-' + stage); s.innerHTML = '';
  s.appendChild(el('div', 'ac-sh', `<h2>${esc(title)}</h2>`));
  const r = el('div', 'ac-roadmap'); r.innerHTML = `<h3>${esc(title)} <span class="tag">v2 — engine wiring in progress</span></h3><p style="font-size:13px;color:var(--ink-2);margin:6px 0 0">This stage uses the in-browser LP solver and our GrowthDB / MediaDB resources. The reference data and algorithm are specified; the interactive panel lands next.</p><ul>${items.map(i => `<li>${i}</li>`).join('')}</ul>`;
  s.appendChild(r);
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
  const dir = el('div', 'ac-roadmap'); dir.innerHTML = `<h3>Reaction directionality <span class="tag">ΔG dataset — v2</span></h3><p style="font-size:12.8px;color:var(--ink-2);margin:6px 0 0">Next: reconcile each reaction's reversibility with ΔG°′ (component-contribution / eQuilibrator) and a consensus of KEGG · MetaCyc · ModelSEED · Rhea, then re-set bounds. Bundling the ΔG reference is the remaining step.</p>`;
  s.appendChild(dir);
  navCount('thermo', blocked.length + (egc.egc ? egc.carriers.length : 0), (egc.egc || blocked.length) ? 'warn' : 'ok');
}
function renderValidate() { roadmap('validate', 'Lab validation', [
  '<b>Strain match</b> — type your strain; we fuzzy-match it against <b>GrowthDB</b> (typo-tolerant, with close hits shown) and import its measured growth / uptake / secretion rates.',
  '<b>Media</b> — the exact growth medium of each measurement is resolved from <b>MediaDB</b> (13k media) and bound onto the model.',
  '<b>GAM / NGAM</b> — growth- and non-growth-associated maintenance are fitted from the GrowthDB rate series.',
  '<b>Prediction vs experiment</b> — every condition with lab data is simulated and compared, with a designed figure + interpretation.']); navCount('validate', '—', 'info'); }

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
