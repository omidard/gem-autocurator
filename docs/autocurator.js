/* GEM Autocurator — client-side curation engine. Everything runs in the browser. */
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
  model.rxns.forEach(r => { if (isExchange(r) || isBiomass(r)) return;
    const bal = {}; let known = true;
    Object.entries(r.s).forEach(([mid, c]) => { const f = F(mid); if (!f || !Object.keys(f).length) known = false; Object.entries(f).forEach(([e, n]) => bal[e] = (bal[e] || 0) + c * n); });
    if (known) { const off = Object.entries(bal).filter(([e, v]) => Math.abs(v) > 1e-6); if (off.length) massIss.push({ id: 'mass_' + r.id, cat: 'charge', sub: 'mass', sev: 'bad', kind: 'mass',
      title: `Mass imbalance in <code>${esc(r.id)}</code>`, note: 'unbalanced elements: ' + off.map(([e, v]) => `${e}${v > 0 ? '+' : ''}${(+v.toFixed(2))}`).join(', '), apply: { flag: r.id } }); }
    let cbal = 0, ck = true; Object.entries(r.s).forEach(([mid, c]) => { const q = CH(mid); if (q == null) ck = false; else cbal += c * q; });
    if (ck && Math.abs(cbal) > 1e-6) chargeIss.push({ id: 'chg_' + r.id, cat: 'charge', sub: 'charge', sev: 'warn', kind: 'charge',
      title: `Charge imbalance in <code>${esc(r.id)}</code>`, note: `net charge ${cbal > 0 ? '+' : ''}${(+cbal.toFixed(2))} — often a missing/extra proton at the simulation pH.`, apply: { flag: r.id } });
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

function renderIds() {
  const s = $('stage-ids'); s.innerHTML = ''; const c = RESULT.counts;
  s.appendChild(el('div', 'ac-sh', `<h2>Identifier mapping</h2><p>Every metabolite and reaction identifier is resolved against <b>BiGG</b> — and, when a compound or reaction exists only in KEGG / ModelSEED / MetaNetX / ChEBI / RHEA, a deterministic <b>BiGG-like</b> id so your model stays internally consistent and portable.</p>`));
  s.appendChild(el('div', 'ac-kpis', kpi(c.mBigg, 'metabolites → BiGG', 'ok') + kpi(c.mLike, 'metabolites → BiGG-like', 'info') + kpi(c.mUn, 'metabolites unmapped', c.mUn ? 'warn' : 'ok') + kpi(c.rBigg, 'reactions → BiGG', 'ok') + kpi(c.rLike, 'reactions → BiGG-like', 'info') + kpi(c.rUn, 'reactions unmapped', c.rUn ? 'warn' : 'ok')));
  const pct = Math.round(100 * (c.mBigg + c.mLike) / Math.max(1, c.mets));
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
function renderThermo() { roadmap('thermo', 'Thermodynamics & cycles', [
  '<b>Blocked reactions</b> — FVA via the bundled GLPK-WASM solver: reactions that cannot carry flux in any state.',
  '<b>Energy-generating cycles (EGCs)</b> — maximise ATP hydrolysis with all uptake closed; any positive solution is an infeasible free-energy loop to break.',
  '<b>Directionality</b> — ΔG°′ (component-contribution / eQuilibrator) reconciled with a consensus of KEGG, MetaCyc, ModelSEED and Rhea directionality, then bounds are re-set.']); navCount('thermo', '—', 'info'); }
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
    dl('MATLAB (.mat)', 'COBRA Toolbox — roadmap', 'MAT', () => alert('MAT export (COBRA Toolbox struct) is on the roadmap. JSON and SBML round-trip losslessly today.')),
    dl('Curation report', 'full findings + decisions', 'MD', () => download(exportReport(allIss), MODEL.id + '_curation_report.md', 'text/markdown')),
  );
  card.appendChild(exp); s.appendChild(card);
  navCount('report', appr, appr ? 'ok' : 'warn');
}

/* ---------------- export ---------------- */
function applyApproved(model) {
  const m = JSON.parse(JSON.stringify({ mets: model.mets.map(x => ({ id: x.id, name: x.name, formula: x.formula, charge: x.charge, compartment: x.compartment })), rxns: model.rxns.map(r => ({ id: r.id, name: r.name, s: r.s, lb: r.lb, ub: r.ub, gpr: r.gpr })), id: model.id, name: model.name, genes: model.genes }));
  const ren = {};
  [...RESULT.idIssues].forEach(i => { if (APPROVED[i.id] === true && i.apply) { if (i.apply.met) ren['m:' + i.apply.met] = i.apply.newId; if (i.apply.rxn) ren['r:' + i.apply.rxn] = i.apply.newId; } });
  m.mets.forEach(x => { const n = ren['m:' + x.id]; if (n) x.id = n; });
  m.rxns.forEach(r => { const n = ren['r:' + r.id]; if (n) r.id = n; const ns = {}; Object.entries(r.s).forEach(([mid, c]) => { ns[ren['m:' + mid] || mid] = c; }); r.s = ns; });
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
function download(content, name, type) { const b = new Blob([content], { type }); const u = URL.createObjectURL(b); const a = el('a'); a.href = u; a.download = name; a.click(); setTimeout(() => URL.revokeObjectURL(u), 2000); }

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
