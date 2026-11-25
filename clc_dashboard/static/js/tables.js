/* =========================================================
   Tables: visor tipo Excel de TODOS los componentes del modelo (/modelo del H5).

   - Hojas: cada componente abierto es una hoja (como en Excel) con su propio estado
     (vista, filtros, orden, pagina, fila seleccionada). Clic en la lista = abrir en la hoja
     activa; clic en un enlace de la formula = abrir en una hoja nueva (o ir a la que ya lo
     tiene); Ctrl+clic = pestaña nueva del navegador.
   - Formula: codigo de la regla (con enlaces que, si hay una fila seleccionada, apuntan al
     elemento exacto: los argumentos de la regla se resuelven con los indices de esa fila),
     expresion expandida del elemento (exacta, enlaces a las Var) y LaTeX cuando existe.
   - "Appears in": que reglas mencionan este componente (sale del codigo de las reglas).

   Los numeros llegan sin redondear; el redondeo es solo de presentacion (Digits). Copiar y
   CSV van siempre con precision completa.
   ========================================================= */

const $ = id => document.getElementById(id);
const TYPES = ['Var', 'Expression', 'Param', 'Constraint', 'Objective'];
const VIOL_TOL = 1e-6;
const PAGE = 500;
const MAX_TABS = 12;
const NOT_MODEL = new Set(['pyo', 'np', 'math', 'd', 'Constraint', 'Expression', 'Var', 'Param']);
const KW = /\b(def|return|if|elif|else|for|in|and|or|not|lambda|sum|None|True|False)\b/g;

const G = {
  run: null, comps: [], byKey: new Map(), sets: {}, refsIn: new Map(),
  type: 'all', stage: 'all', query: '',
  cstat: 'all',   // restricciones en la lista: 'all' | 'viol' (algun elemento violado) | 'ok'
  values: 'both', idxMode: 'phys', digits: '6', formulaCollapsed: false,
  tabs: [], active: -1, seq: 0,
  cmpRun: null,   // corrida B con la que se compara la tabla (el mismo componente)
};
const DIFF_RTOL = 1e-9;  // |B − A| por encima de esto (relativo a max(1,|A|)) = "difiere"
const FIELD_LABEL = { valor: 'value (Pyomo)', valor_fisico: 'value (physical)', lb: 'lb', ub: 'ub', cuerpo: 'body',
                      inferior: 'lower', superior: 'upper', residuo: 'residual', dual: 'dual' };
const cur = () => G.tabs[G.active] || null;

const list = new RunList({
  multi: false,
  onPick: name => { if (name !== G.run) loadRun(name); },
  isSelected: n => n === G.run,
  filter: r => r.has_model,
  emptyMessage: 'No run has the full model export yet. Re-export a run with the new notebook cell.',
});

/** Desglosable = escalar con regla (la funcion objetivo y sus Expressions de costo). */
const canBreakdown = c => !!c && c.n === 1 && !!c.regla && (c.tipo === 'Objective' || c.tipo === 'Expression');

function newTab(key) {
  const comp = G.byKey.get(key);
  return { key, comp, data: null, view: comp && comp.tipo === 'Objective' && canBreakdown(comp) ? 'breakdown' : 'long',
           sort: null, dir: 1, page: 0, bd: null, bdB: null, bdClosed: new Set(),
           filters: {}, rowStat: 'all', pivot: { row: 0, col: -1, field: 'valor', fixed: {} },
           sel: null, fview: null, pending: null,
           cmp: null, cmpField: null, onlyDiff: false };
}

// ---------------------------------------------------------------------------
// Formato
// ---------------------------------------------------------------------------

function num(v) {
  if (v === null || v === undefined) return '—';
  if (typeof v !== 'number') return String(v);
  if (G.digits === 'full') return String(v);
  const d = +G.digits;
  if (v === 0) return '0';
  const a = Math.abs(v);
  if (a < 1e-4 || a >= 1e15) return v.toExponential(d - 1).replace(/\.?0+e/, 'e');
  // mas digitos enteros que los pedidos: el entero completo (con miles), nunca ceros inventados
  // (toPrecision(6) de 1898144.88 daria 1898140)
  if (a >= 10 ** d) return Math.round(v).toLocaleString('en-US');
  return String(Number(v.toPrecision(d)));
}

function cellHTML(v, cls = '', title = '') {
  const t = title ? ` title="${esc(title)}"` : '';
  if (v === null || v === undefined) return `<td class="null ${cls}"${t}>—</td>`;
  if (typeof v === 'boolean') return v ? `<td class="bool-yes ${cls}"${t}>✓</td>` : `<td class="bool-no ${cls}"${t}>·</td>`;
  if (typeof v === 'string') return `<td class="${cls}"${t}>${esc(v)}</td>`;
  return `<td class="${cls}"${t}>${num(v)}</td>`;
}

// ---------------------------------------------------------------------------
// Indices y columnas
// ---------------------------------------------------------------------------

function dimLabel(tab, i) {
  const dm = tab.comp.dims[i];
  const phys = G.idxMode === 'phys' && tab.data.idx_fis[i];
  if (!phys) return { main: dm.set, sub: dm.fisico ? 'dimensionless' : '' };
  const sym = (dm.regla.split('=')[0] || dm.set).trim();
  return { main: `${sym} (${dm.unidades})`, sub: `${dm.set} · ${dm.regla}` };
}

function idxValue(tab, i, r) {
  const fis = G.idxMode === 'phys' && tab.data.idx_fis[i];
  return fis ? fis[r] : tab.data.idx[i][r];
}

/** El indice i de la fila r coincide con `val` (texto)? Numerico con tolerancia: Python
 *  escribe 1.0 donde JS escribe 1. */
function idxMatch(tab, i, r, val) {
  const v = tab.data.idx[i][r];
  if (typeof v === 'number') { const x = +val; return isFinite(x) && Math.abs(v - x) <= 1e-9 * Math.max(1, Math.abs(v)); }
  return String(v) === String(val);
}

function uniqueVals(tab, i) {
  const seen = new Map();
  tab.data.idx[i].forEach((v, r) => { const k = String(v); if (!seen.has(k)) seen.set(k, r); });
  const keys = [...seen.keys()];
  if (typeof tab.data.idx[i][0] === 'number') keys.sort((a, b) => +a - +b);
  return keys.map(k => ({ key: k, label: labelOfIdx(tab, i, seen.get(k)) }));
}

function labelOfIdx(tab, i, r) {
  const v = idxValue(tab, i, r);
  return typeof v === 'number' ? num(v) : String(v);
}

/** Columnas de valores de la vista larga. Cada una trae get(r) -> valor de la fila r. */
function valueColumns(tab) {
  if (cmpActive(tab)) return compareColumns(tab);
  const c = tab.comp, F = tab.data.fields, cols = [];
  const unit = c.unidades ? ` [${c.unidades}]` : '';
  const col = (key, main, sub) => ({ key, main, sub, get: r => F[key][r] });
  if (c.tipo === 'Constraint') {
    [['cuerpo', 'body'], ['inferior', 'lower'], ['superior', 'upper'], ['residuo', 'residual'], ['dual', 'dual'], ['activa', 'active']]
      .forEach(([k, l]) => { if (F[k]) cols.push(col(k, l, '')); });
    return cols;
  }
  const hasPhys = !!F.valor_fisico;
  if (G.values !== 'phys' || !hasPhys) cols.push(col('valor', 'value', hasPhys ? 'Pyomo (scaled)' : (c.escala_conocida ? 'no scaling' : 'Pyomo · scale unknown')));
  if (G.values !== 'pyomo' && hasPhys) cols.push(col('valor_fisico', 'value', `physical${unit}`));
  if (c.tipo === 'Var') {
    [['lb', 'lb', 'Pyomo'], ['ub', 'ub', 'Pyomo'], ['fijo', 'fixed', '']].forEach(([k, l, s]) => { if (F[k]) cols.push(col(k, l, s)); });
  }
  return cols;
}

// ---------------------------------------------------------------------------
// Comparacion con otra corrida: el MISMO componente, elemento por elemento
// ---------------------------------------------------------------------------

const cmpActive = tab => !!(G.cmpRun && tab.cmp && tab.cmp.run === G.cmpRun && tab.cmp.data);

/** Campo que se compara por defecto: el valor que se esta mirando (fisico o Pyomo); en una
 *  restriccion, el residuo. */
function baseField(tab) {
  if (tab.comp.tipo === 'Constraint') return 'residuo';
  return tab.data.fields.valor_fisico && G.values !== 'pyomo' ? 'valor_fisico' : 'valor';
}

function cmpFields(tab) {
  const FA = tab.data.fields, FB = tab.cmp.data.fields;
  return Object.keys(FIELD_LABEL).filter(k => FA[k] && FB[k]);
}

function cmpField(tab) {
  const ok = cmpFields(tab);
  return ok.includes(tab.cmpField) ? tab.cmpField : (ok.includes(baseField(tab)) ? baseField(tab) : ok[0]);
}

/** Clave de un elemento por el valor de sus indices: asi se emparejan A y B aunque el orden
 *  de los elementos difiera (o la malla sea distinta: lo que no empareja queda vacio). */
function idxKey(data, r) {
  return data.idx.map(col => typeof col[r] === 'number' ? Number(col[r]).toPrecision(12) : String(col[r])).join('\u0001');
}

function buildRowMap(tab) {
  const A = tab.data, B = tab.cmp.data;
  if (A.idx.length !== B.idx.length) { tab.cmp.error = 'different number of indices in B'; tab.cmp.data = null; return; }
  const scalar = A.idx.length === 0;
  const nB = scalar ? (Object.values(B.fields)[0] || []).length : B.idx[0].length;
  const posB = new Map();
  if (!scalar) for (let r = 0; r < nB; r++) posB.set(idxKey(B, r), r);
  const nA = tab.comp.n;
  const rowOf = new Array(nA);
  let matched = 0;
  for (let r = 0; r < nA; r++) {
    const b = scalar ? (nB ? 0 : undefined) : posB.get(idxKey(A, r));
    rowOf[r] = b === undefined ? -1 : b;
    if (b !== undefined) matched++;
  }
  Object.assign(tab.cmp, { rowOf, matched, nB });
}

function compareColumns(tab) {
  const f = cmpField(tab), FA = tab.data.fields, FB = tab.cmp.data.fields, rowOf = tab.cmp.rowOf;
  const A = r => FA[f][r];
  const B = r => rowOf[r] < 0 ? null : FB[f][rowOf[r]];
  const D = r => { const a = A(r), b = B(r); return typeof a === 'number' && typeof b === 'number' ? b - a : null; };
  const P = r => { const a = A(r), d = D(r); return d === null ? null : a === 0 ? (d === 0 ? 0 : null) : 100 * d / Math.abs(a); };
  const unit = f === 'valor_fisico' && tab.comp.unidades ? ` [${tab.comp.unidades}]` : '';
  const lab = FIELD_LABEL[f] + unit;
  return [
    { key: 'cmpA', main: 'A', sub: lab, get: A },
    { key: 'cmpB', main: 'B', sub: lab, get: B, cmp: true },
    { key: 'cmpD', main: 'Δ = B − A', sub: lab, get: D, diff: true },
    { key: 'cmpP', main: 'Δ %', sub: 'relative to |A|', get: P, diff: true },
  ];
}

/** El elemento r existe en A y en B y su valor cambia. Lo que no esta en B (otra malla) no
 *  cuenta como diferencia: se reporta aparte como "matched". */
function differs(tab, r) {
  if (tab.cmp.rowOf[r] < 0) return false;
  const cols = compareColumns(tab), a = cols[0].get(r), d = cols[2].get(r);
  if (d === null) return cols[0].get(r) !== cols[1].get(r);  // vacio en uno solo de los dos
  return Math.abs(d) > DIFF_RTOL * Math.max(1, Math.abs(a));
}

async function ensureCmp(tab) {
  const run = G.cmpRun;
  tab.cmp = { run, data: null, loading: true };
  try {
    const d = await api(`/api/model/${encodeURIComponent(run)}/component?key=${encodeURIComponent(tab.key)}`);
    if (!tab.cmp || tab.cmp.run !== run) return;
    tab.cmp.data = d;
    buildRowMap(tab);
  } catch (e) {
    if (!tab.cmp || tab.cmp.run !== run) return;
    tab.cmp.error = /not found/i.test(e.message) ? 'this component does not exist in B' : e.message;
  }
  tab.cmp.loading = false;
  if (cur() === tab) render();
}

function fillCmpSelect() {
  const sel = $('selCmp');
  const cands = list.runs.filter(r => r.has_model && r.name !== G.run);
  if (G.cmpRun && !cands.some(r => r.name === G.cmpRun)) G.cmpRun = null;
  sel.innerHTML = `<option value="">— none —</option>` + cands.map(r =>
    `<option value="${esc(r.name)}" ${r.name === G.cmpRun ? 'selected' : ''}>B: ${fmtDate(r.date)} · ${esc(r.label)}</option>`).join('');
  sel.classList.toggle('on', !!G.cmpRun);
}

function atBound(tab, r) {
  const F = tab.data.fields;
  if (tab.comp.tipo !== 'Var' || !F.lb || (F.fijo && F.fijo[r])) return null;
  const v = F.valor[r];
  if (v === null) return null;
  const near = b => b !== null && Math.abs(v - b) <= 1e-7 * Math.max(1, Math.abs(b));
  return near(F.lb[r]) ? 'lb' : near(F.ub[r]) ? 'ub' : null;
}

function cmp(a, b) {
  if (a === b) return 0;
  if (a === null || a === undefined) return 1;
  if (b === null || b === undefined) return -1;
  if (typeof a === 'string' || typeof b === 'string') return String(a).localeCompare(String(b));
  return a < b ? -1 : 1;
}

function longRows(tab) {
  const F = tab.data.fields;
  let rows = [...Array(tab.comp.n).keys()];
  for (const [i, val] of Object.entries(tab.filters)) {
    if (val === '' || val === undefined) continue;
    rows = rows.filter(r => idxMatch(tab, +i, r, val));
  }
  if (tab.rowStat !== 'all' && F.residuo && !cmpActive(tab)) {
    const viol = r => F.residuo[r] !== null && F.residuo[r] > VIOL_TOL;
    rows = rows.filter(tab.rowStat === 'viol' ? viol : r => !viol(r));
  }
  if (tab.onlyDiff && cmpActive(tab)) rows = rows.filter(r => differs(tab, r));
  if (tab.sort !== null) {
    // '@i' = indice i; cualquier otra clave = una columna de valores
    const col = tab.sort.startsWith('@') ? null : valueColumns(tab).find(c => c.key === tab.sort);
    if (col || tab.sort.startsWith('@')) {
      const get = col ? col.get : r => idxValue(tab, +tab.sort.slice(1), r);
      rows.sort((a, b) => tab.dir * cmp(get(a), get(b)));
    }
  }
  return rows;
}

// ---------------------------------------------------------------------------
// Render principal
// ---------------------------------------------------------------------------

function render() {
  renderTabs();
  const tab = cur();
  if (!tab) { showTableEmpty(G.run ? 'Pick a component from the list.' : 'Pick a run, then a component.'); return; }
  if (!tab.data) { showTableEmpty('Loading…'); return; }
  if (G.cmpRun && (!tab.cmp || tab.cmp.run !== G.cmpRun)) ensureCmp(tab);  // pide B en segundo plano
  if (!G.cmpRun) tab.cmp = null;
  renderInfo(tab);
  renderFormula(tab);
  const pivotOK = tab.comp.dims.length >= 1, bdOK = canBreakdown(tab.comp);
  if (tab.view === 'breakdown' && !bdOK) tab.view = 'long';
  document.querySelectorAll('#segView .cfg-tab').forEach(b => {
    b.classList.toggle('active', b.dataset.v === tab.view);
    if (b.dataset.v === 'pivot') b.disabled = !pivotOK;
    if (b.dataset.v === 'breakdown') b.hidden = !bdOK;
  });
  if (tab.view === 'breakdown') renderBreakdown(tab);
  else if (tab.view === 'pivot' && pivotOK) renderPivot(tab);
  else renderLong(tab);
  const f = Object.fromEntries(Object.entries(tab.filters).filter(([, v]) => v !== ''));
  writeURL({ run: G.run, c: tab.key, view: tab.view !== 'long' ? tab.view : '', f: Object.keys(f).length ? JSON.stringify(f) : '', cmp: G.cmpRun || '' });
}

function renderTabs() {
  $('sheetTabs').innerHTML = G.tabs.map((t, i) => `
    <span class="sheet-tab ${i === G.active ? 'active' : ''}" data-i="${i}" title="${esc(t.comp ? `${t.comp.tipo} · ${t.comp.doc || t.key}` : t.key)}">
      ${esc(t.key)}<button class="x" data-close="${i}" title="Close">×</button></span>`).join('');
}

function showTableEmpty(msg) {
  $('tblTitle').textContent = 'Table'; $('tblSubtitle').innerHTML = '';
  $('compInfo').innerHTML = ''; $('formula').innerHTML = ''; $('toolbarSub').innerHTML = ''; $('pager').innerHTML = '';
  $('tableBody').innerHTML = `<div class="run-report-empty"><span class="ph-icon">▦</span><p>${msg}</p></div>`;
}

function renderInfo(tab) {
  const c = tab.comp, F = tab.data.fields;
  $('tblTitle').textContent = c.name;
  $('tblSubtitle').innerHTML = `<span class="type-badge type-${c.tipo}">${c.tipo}</span>`;
  const main = F.valor_fisico && G.values !== 'pyomo' ? F.valor_fisico : (c.tipo === 'Constraint' ? F.residuo : F.valor);
  const finite = (main || []).filter(v => typeof v === 'number');
  const nNull = (main || []).length - finite.length;
  let mn = null, mx = null;
  finite.forEach(v => { if (mn === null || v < mn) mn = v; if (mx === null || v > mx) mx = v; });
  const what = c.tipo === 'Constraint' ? 'residual' : (main === F.valor_fisico ? 'physical' : 'value');
  const scale = c.escala_por_elemento ? 'per element' : c.escala_conocida ? (c.escala === null ? '—' : num(c.escala)) : 'unknown';
  const dims = c.dims.length ? c.dims.map(d => `<code>${esc(d.set)}</code>${d.regla && d.fisico ? ` <span style="color:var(--text-muted)">(${esc(d.regla)})</span>` : ''}`).join(' × ') : 'scalar';
  let viol = '';
  if (F.residuo) {
    const nv = F.residuo.filter(v => v !== null && v > VIOL_TOL).length;
    viol = `<span><b>violated</b><span class="${nv ? 'badge-fail' : 'badge-ok'}">${nv} &gt; ${VIOL_TOL}</span></span>`;
  }
  let bounds = '';
  if (c.tipo === 'Var' && F.lb) {
    let nb = 0;
    for (let r = 0; r < c.n; r++) if (atBound(tab, r)) nb++;
    if (nb) bounds = `<span><b>at bounds</b><span class="badge-warn">${nb}</span></span>`;
  }
  $('compInfo').innerHTML = `
    ${c.doc ? `<span class="ci-doc">${esc(c.doc)}</span>` : ''}
    <span><b>stage</b>${esc(c.etapa || 'global')}${c.bloque ? ` · block <code>${esc(c.bloque)}</code>` : ''}</span>
    <span><b>indices</b>${dims}</span>
    <span><b>units</b>${esc(c.unidades || '—')}</span>
    <span><b>scale</b>${scale}</span>
    <span><b>n</b><span class="stat">${c.n}</span></span>
    <span><b>${what}</b><span class="stat">${num(mn)} … ${num(mx)}</span>${nNull ? ` · <span class="stat">${nNull}</span> empty` : ''}</span>
    ${viol}${bounds}
    ${c.sentido ? `<span><b>sense</b>${esc(c.sentido)}</span>` : ''}
    ${compareInfo(tab)}`;
}

function compareInfo(tab) {
  if (!G.cmpRun || !tab.cmp) return '';
  const rA = list.get(G.run), rB = list.get(G.cmpRun);
  const legend = `<span class="cmp-legend-inline"><b class="tagA">A</b>${esc(rA ? fmtDate(rA.date) : G.run)} &nbsp; <b class="tagB">B</b>${esc(rB ? fmtDate(rB.date) : G.cmpRun)}</span>`;
  if (tab.cmp.loading) return `<span class="ci-doc">${legend} · loading B…</span>`;
  if (tab.cmp.error) return `<span class="ci-doc">${legend} · <span class="badge-warn">${esc(tab.cmp.error)}</span></span>`;
  const n = tab.comp.n, cols = compareColumns(tab);
  let nd = 0, maxD = 0, maxP = 0;
  for (let r = 0; r < n; r++) {
    if (differs(tab, r)) nd++;
    const d = cols[2].get(r), p = cols[3].get(r);
    if (typeof d === 'number' && Math.abs(d) > Math.abs(maxD)) maxD = d;
    if (typeof p === 'number' && Math.abs(p) > Math.abs(maxP)) maxP = p;
  }
  const unmatched = n - tab.cmp.matched;
  return `<span class="ci-doc">${legend}
    · <b>matched</b><span class="stat">${tab.cmp.matched}/${n}</span>${tab.cmp.nB !== n || unmatched ? ` <span class="badge-warn">B has ${tab.cmp.nB} elements · grids differ?</span>` : ''}
    · <b>differ</b><span class="${nd ? 'badge-warn' : 'badge-ok'}">${nd}</span>
    · <b>max |Δ|</b><span class="stat">${num(maxD)}</span>
    · <b>max |Δ%|</b><span class="stat">${num(maxP)} %</span></span>`;
}

// ---------------------------------------------------------------------------
// Vista larga
// ---------------------------------------------------------------------------

function renderLong(tab) {
  const c = tab.comp, F = tab.data.fields;
  const vcols = valueColumns(tab);
  let sub = c.dims.map((d, i) => {
    const opts = uniqueVals(tab, i);
    if (opts.length < 2) return '';
    return `<div class="tool-group"><span class="tool-label">${esc(dimLabel(tab, i).main)}</span>
      <select class="panel-select" data-filter="${i}"><option value="">all</option>
      ${opts.map(o => `<option value="${esc(o.key)}" ${tab.filters[i] === o.key ? 'selected' : ''}>${esc(o.label)}</option>`).join('')}</select></div>`;
  }).join('');
  const comparing = cmpActive(tab);
  if (comparing) {
    const cf = cmpField(tab);
    sub += `<div class="tool-group"><span class="tool-label">Compare field</span><select class="panel-select" id="selCmpField">
        ${cmpFields(tab).map(k => `<option value="${k}" ${k === cf ? 'selected' : ''}>${FIELD_LABEL[k]}</option>`).join('')}</select></div>
      <div class="tool-group"><label class="label-check"><input type="checkbox" id="chkDiff" ${tab.onlyDiff ? 'checked' : ''}> only rows that differ</label></div>`;
  } else if (F.residuo) {
    const nv = F.residuo.filter(v => v !== null && v > VIOL_TOL).length;
    sub += `<div class="tool-group"><span class="tool-label">Elements</span><div class="cfg-tabs" id="segRowStat" title="violated = residual > ${VIOL_TOL}">
      ${[['all', 'All', c.n], ['viol', '✗ Violated', nv], ['ok', '✓ Satisfied', c.n - nv]].map(([k, l, n]) =>
        `<button class="cfg-tab ${tab.rowStat === k ? 'active' : ''}" data-rs="${k}">${l} ${n}</button>`).join('')}</div></div>`;
  }
  $('toolbarSub').innerHTML = sub;
  $('toolbarSub').querySelectorAll('[data-filter]').forEach(s => s.onchange = () => { tab.filters[s.dataset.filter] = s.value; tab.page = 0; render(); });
  $('toolbarSub').querySelectorAll('[data-rs]').forEach(b => b.onclick = () => { tab.rowStat = b.dataset.rs; tab.page = 0; render(); });
  if ($('chkDiff')) $('chkDiff').onchange = e => { tab.onlyDiff = e.target.checked; tab.page = 0; render(); };
  if ($('selCmpField')) $('selCmpField').onchange = e => { tab.cmpField = e.target.value; render(); };

  const rows = longRows(tab);
  const pages = Math.max(1, Math.ceil(rows.length / PAGE));
  tab.page = Math.min(tab.page, pages - 1);
  const slice = rows.slice(tab.page * PAGE, (tab.page + 1) * PAGE);

  const th = (key, lab, cls = '') => {
    const s = tab.sort === key ? (tab.dir > 0 ? ' ▲' : ' ▼') : '';
    return `<th class="${cls} ${tab.sort === key ? 'sorted' : ''}" data-sort="${key}" title="Sort">${esc(lab.main)}${s}${lab.sub ? `<span class="th-sub">${esc(lab.sub)}</span>` : ''}</th>`;
  };
  let html = `<table class="data-table"><thead><tr><th class="rownum">#</th>`;
  c.dims.forEach((_, i) => { html += th('@' + i, dimLabel(tab, i), 'idx-col'); });
  vcols.forEach(col => { html += th(col.key, col, col.cmp ? 'cmp-col' : ''); });
  html += `</tr></thead><tbody>`;
  for (const r of slice) {
    const viol = !comparing && F.residuo && F.residuo[r] !== null && F.residuo[r] > VIOL_TOL;
    const bnd = comparing ? null : atBound(tab, r);
    const dif = comparing && differs(tab, r);
    html += `<tr data-r="${r}" class="${viol ? 'row-viol' : ''} ${bnd ? 'row-bound' : ''} ${tab.sel === r ? 'row-sel' : ''}"><td class="rownum">${r + 1}</td>`;
    c.dims.forEach((_, i) => { const v = idxValue(tab, i, r); html += `<td class="idx-col">${typeof v === 'number' ? num(v) : esc(v)}</td>`; });
    vcols.forEach(col => {
      const v = col.get(r);
      if (col.diff && dif) html += cellHTML(v, 'cmp-diff');
      else if (bnd && (col.key === 'valor' || col.key === 'valor_fisico' || col.key === bnd)) html += cellHTML(v, 'at-bound', `at ${bnd}`);
      else html += cellHTML(v);
    });
    html += `</tr>`;
  }
  if (!slice.length) html += `<tr><td class="results-empty" colspan="${1 + c.dims.length + vcols.length}">No rows match the filters.</td></tr>`;
  html += `</tbody></table>`;
  $('tableBody').innerHTML = html;
  $('tableBody').querySelectorAll('th[data-sort]').forEach(h => h.onclick = () => {
    const k = h.dataset.sort;
    if (tab.sort === k) { if (tab.dir > 0) tab.dir = -1; else { tab.sort = null; tab.dir = 1; } }
    else { tab.sort = k; tab.dir = 1; }
    render();
  });
  $('tableBody').querySelector('tr.row-sel')?.scrollIntoView({ block: 'nearest' });

  $('pager').innerHTML = rows.length > PAGE ? `
      <button class="btn-sm" id="pgPrev" ${tab.page === 0 ? 'disabled' : ''}>‹ Prev</button>
      <span>page ${tab.page + 1} / ${pages}</span>
      <button class="btn-sm" id="pgNext" ${tab.page >= pages - 1 ? 'disabled' : ''}>Next ›</button>
      <span class="spacer"></span><span>${rows.length} of ${c.n} rows</span>`
    : `<span class="spacer"></span><span>${rows.length} of ${c.n} rows${c.regla && !comparing ? ' · click a row to point the formula links at it' : ''}</span>`;
  if ($('pgPrev')) $('pgPrev').onclick = () => { tab.page--; render(); };
  if ($('pgNext')) $('pgNext').onclick = () => { tab.page++; render(); };
}

// ---------------------------------------------------------------------------
// Pivote
// ---------------------------------------------------------------------------

/** Campos que el pivote puede mostrar, cada uno con get(r). Comparando: A, B, Δ y Δ% del
 *  campo de comparacion. */
function pivotFields(tab) {
  if (cmpActive(tab)) {
    const f = FIELD_LABEL[cmpField(tab)];
    return compareColumns(tab).map(c => ({ key: c.key, label: `${c.main} · ${f}`, get: c.get }));
  }
  const F = tab.data.fields;
  return Object.keys(FIELD_LABEL).filter(k => F[k]).map(k => ({ key: k, label: FIELD_LABEL[k], get: r => F[k][r] }));
}

function pivotField(tab) {
  const fs = pivotFields(tab);
  return fs.find(f => f.key === tab.pivot.field) || fs.find(f => f.key === 'cmpD') || fs[0];
}

function defaultPivot(tab) {
  const sets = tab.comp.dims.map(d => d.set);
  const row = sets.indexOf('A') >= 0 ? sets.indexOf('A') : Math.max(0, tab.comp.dims.findIndex(d => d.fisico));
  let col = sets.indexOf('T');
  if (col === row || col < 0) col = sets.findIndex((_, i) => i !== row);
  const field = baseField(tab);
  const fixed = {};
  tab.comp.dims.forEach((_, i) => { if (i !== row && i !== col) fixed[i] = uniqueVals(tab, i)[0].key; });
  tab.pivot = { row, col, field, fixed };
}

function pivotGrid(tab) {
  const P = tab.pivot, n = tab.comp.n, field = pivotField(tab);
  const rowVals = uniqueVals(tab, P.row);
  const colVals = P.col >= 0 ? uniqueVals(tab, P.col) : [{ key: '', label: field.label }];
  const ri = new Map(rowVals.map((v, i) => [v.key, i])), ci = new Map(colVals.map((v, i) => [v.key, i]));
  const grid = rowVals.map(() => colVals.map(() => undefined));
  for (let r = 0; r < n; r++) {
    let ok = true;
    for (const [i, val] of Object.entries(P.fixed)) if (String(tab.data.idx[+i][r]) !== val) { ok = false; break; }
    if (!ok) continue;
    const a = ri.get(String(tab.data.idx[P.row][r]));
    const b = P.col >= 0 ? ci.get(String(tab.data.idx[P.col][r])) : 0;
    grid[a][b] = field.get(r);
  }
  return { rowVals, colVals, grid, field };
}

function renderPivot(tab) {
  const c = tab.comp, P = tab.pivot;
  const dimOpt = (sel, allowNone) => (allowNone ? `<option value="-1" ${sel === -1 ? 'selected' : ''}>(none)</option>` : '') +
    c.dims.map((d, i) => `<option value="${i}" ${sel === i ? 'selected' : ''}>${esc(dimLabel(tab, i).main)}</option>`).join('');
  let sub = `
    <div class="tool-group"><span class="tool-label">Rows</span><select class="panel-select" id="pvRow">${dimOpt(P.row)}</select></div>
    <div class="tool-group"><span class="tool-label">Columns</span><select class="panel-select" id="pvCol">${dimOpt(P.col, true)}</select></div>
    <div class="tool-group"><span class="tool-label">Field</span><select class="panel-select" id="pvField">${pivotFields(tab).map(f => `<option value="${f.key}" ${f.key === pivotField(tab).key ? 'selected' : ''}>${esc(f.label)}</option>`).join('')}</select></div>`;
  if (cmpActive(tab)) {
    const cf = cmpField(tab);
    sub += `<div class="tool-group"><span class="tool-label">Compare field</span><select class="panel-select" id="selCmpField">
      ${cmpFields(tab).map(k => `<option value="${k}" ${k === cf ? 'selected' : ''}>${FIELD_LABEL[k]}</option>`).join('')}</select></div>`;
  }
  c.dims.forEach((_, i) => {
    if (i === P.row || i === P.col) return;
    sub += `<div class="tool-group"><span class="tool-label">${esc(dimLabel(tab, i).main)} =</span><select class="panel-select" data-fixed="${i}">
      ${uniqueVals(tab, i).map(o => `<option value="${esc(o.key)}" ${P.fixed[i] === o.key ? 'selected' : ''}>${esc(o.label)}</option>`).join('')}</select></div>`;
  });
  $('toolbarSub').innerHTML = sub;
  const regrid = () => {
    const fixed = {};
    c.dims.forEach((_, i) => { if (i !== P.row && i !== P.col) fixed[i] = P.fixed[i] ?? uniqueVals(tab, i)[0].key; });
    P.fixed = fixed;
    render();
  };
  $('pvRow').onchange = e => { const v = +e.target.value; if (v === P.col) P.col = P.row; P.row = v; regrid(); };
  $('pvCol').onchange = e => { const v = +e.target.value; if (v === P.row) P.row = P.col >= 0 ? P.col : P.row; P.col = v; regrid(); };
  $('pvField').onchange = e => { P.field = e.target.value; render(); };
  if ($('selCmpField')) $('selCmpField').onchange = e => { tab.cmpField = e.target.value; render(); };
  $('toolbarSub').querySelectorAll('[data-fixed]').forEach(s => s.onchange = () => { P.fixed[s.dataset.fixed] = s.value; render(); });

  const { rowVals, colVals, grid, field } = pivotGrid(tab);
  const isDiff = field.key === 'cmpD' || field.key === 'cmpP';
  const rl = dimLabel(tab, P.row), cl = P.col >= 0 ? dimLabel(tab, P.col) : null;
  let html = `<table class="data-table pivot-table"><thead><tr><th class="corner">${esc(rl.main)} \\ ${cl ? esc(cl.main) : ''}</th>`;
  colVals.forEach(cv => { html += `<th>${esc(cv.label)}</th>`; });
  html += `</tr></thead><tbody>`;
  rowVals.forEach((rv, a) => {
    html += `<tr><td class="row-head">${esc(rv.label)}</td>`;
    colVals.forEach((_, b) => {
      const v = grid[a][b];
      const hot = isDiff && typeof v === 'number' && v !== 0;
      html += v === undefined ? '<td class="null"></td>' : cellHTML(v, hot ? 'cmp-diff' : '');
    });
    html += `</tr>`;
  });
  html += `</tbody></table>`;
  $('tableBody').innerHTML = html;
  $('pager').innerHTML = `<span class="spacer"></span><span>${rowVals.length} × ${colVals.length} cells · field: ${esc(field.label)}</span>`;
}

// ---------------------------------------------------------------------------
// Desglose (Breakdown): terminos de la funcion objetivo / de una Expression escalar
// ---------------------------------------------------------------------------

/** slot 'bd' = esta corrida; slot 'bdB' = corrida B, desglosada con la REGLA de esta corrida
 *  y los valores de B (asi funciona aunque B se haya exportado sin formulas). */
async function ensureBd(tab, slot, run) {
  tab[slot] = { run, loading: true };
  try {
    const values = run !== G.run ? `&values=${encodeURIComponent(run)}` : '';
    const d = await api(`/api/model/${encodeURIComponent(G.run)}/breakdown?key=${encodeURIComponent(tab.key)}${values}`);
    if (!tab[slot] || tab[slot].run !== run) return;
    tab[slot] = { run, data: d };
  } catch (e) {
    if (!tab[slot] || tab[slot].run !== run) return;
    tab[slot] = { run, error: e.message };
  }
  if (cur() === tab) render();
}

/** Filas del arbol en orden, con la contribucion de cada nodo a la raiz (signo acumulado).
 *  all=true ignora los nodos colapsados (para exportar y para emparejar con B). */
function bdRows(tab, data, all = false) {
  const rows = [];
  const walk = (n, depth, acc) => {
    const s = acc * n.sign;
    rows.push({ n, depth, contrib: n.kind === 'factor' || n.value === null ? null : s * n.value });
    if (all || !tab.bdClosed.has(n.id)) n.children.forEach(c => walk(c, depth + 1, s));
  };
  walk(data.root, 0, 1);
  return rows;
}

function renderBreakdown(tab) {
  if (!tab.bd || tab.bd.run !== G.run) ensureBd(tab, 'bd', G.run);
  if (G.cmpRun && (!tab.bdB || tab.bdB.run !== G.cmpRun)) ensureBd(tab, 'bdB', G.cmpRun);
  if (!G.cmpRun) tab.bdB = null;
  $('toolbarSub').innerHTML = `<div class="tool-group">
      <button class="btn-sm" id="bdOpen">⊞ Expand all</button><button class="btn-sm" id="bdClose">⊟ Collapse all</button></div>
    <div class="tool-group"><span class="formula-ctx">contribution = value × the signs of its parents · % of |${esc(tab.comp.name)}|</span></div>`;
  $('bdOpen').onclick = () => { tab.bdClosed.clear(); render(); };
  $('bdClose').onclick = () => { if (tab.bd && tab.bd.data) bdRows(tab, tab.bd.data, true).forEach(r => { if (r.depth >= 1 && r.n.children.length) tab.bdClosed.add(r.n.id); }); render(); };

  const S = tab.bd;
  if (!S || S.loading) { $('tableBody').innerHTML = '<div class="run-report-empty"><p>Loading…</p></div>'; $('pager').innerHTML = ''; return; }
  if (S.error) { $('tableBody').innerHTML = `<div class="error-box">${esc(S.error)}</div>`; $('pager').innerHTML = ''; return; }
  const d = S.data, rows = bdRows(tab, d);
  const B = tab.bdB && tab.bdB.data ? new Map(bdRows(tab, tab.bdB.data, true).map(r => [r.n.id, r])) : null;
  const denom = Math.abs(d.value) || 1;
  const maxAbs = Math.max(1e-300, ...rows.filter(r => r.depth > 0 && r.contrib !== null).map(r => Math.abs(r.contrib)));
  const unitOf = n => n.unidades || (n.kind === 'factor' ? '' : d.unidades || '');

  let check = `<b>${esc(tab.comp.name)}</b>${d.sense ? ` (${esc(d.sense)})` : ''} = <span class="stat">${num(d.value)}</span>
    · Σ of the terms = <span class="stat">${num(d.evaluated)}</span>
    ${d.verified ? '<span class="badge-ok">✓ matches Pyomo</span>' : '<span class="badge-warn">⚠ does not match the Pyomo value</span>'}`;
  if (G.cmpRun) {
    const rB = list.get(G.cmpRun);
    check += tab.bdB && tab.bdB.data ? ` · <span class="cmp-legend-inline"><b class="tagB">B</b>${esc(rB ? fmtDate(rB.date) : '')}</span> = <span class="stat">${num(tab.bdB.data.value)}</span>
        ${tab.bdB.data.verified ? '<span class="badge-ok">✓</span>' : '<span class="badge-warn">⚠ B terms do not add up to its value</span>'}`
      : tab.bdB && tab.bdB.error ? ` · <span class="badge-warn">B: ${esc(tab.bdB.error)}</span>` : ' · loading B…';
  }

  let html = `<div class="bd-check">${check}</div><table class="data-table bd-table"><thead><tr>
    <th class="idx-col">Term</th><th class="idx-col">Description</th><th>Value</th><th class="idx-col">Units</th>
    <th>Contribution${B ? ' A' : ''}</th>${B ? '<th class="cmp-col">Contribution B</th><th>Δ = B − A</th>' : ''}<th>% of |${esc(tab.comp.name)}|</th><th class="idx-col"></th></tr></thead><tbody>`;
  for (const { n, depth, contrib } of rows) {
    const open = !tab.bdClosed.has(n.id);
    const tw = n.children.length ? `<span class="tw" data-tw="${esc(n.id)}">${open ? '▾' : '▸'}</span>` : '<span class="tw"></span>';
    const sg = depth === 0 ? '' : `<span class="sg ${n.sign > 0 ? 'pos' : 'neg'}">${n.sign > 0 ? '+' : '−'}</span>`;
    const label = n.key && G.byKey.has(n.key) ? refLink(n.key, n.label, null) : esc(n.label);
    const cls = depth === 0 ? 'bd-root' : n.kind === 'local' ? 'bd-local' : n.kind === 'factor' ? 'bd-factor' : '';
    const share = contrib === null || depth === 0 ? null : 100 * contrib / denom;
    const bar = contrib === null || depth === 0 ? '' :
      `<span class="${contrib >= 0 ? 'pos' : 'neg'}" style="width:${Math.max(1, 100 * Math.abs(contrib) / maxAbs).toFixed(1)}%" title="${contrib >= 0 ? 'adds to' : 'subtracts from'} ${esc(tab.comp.name)}"></span>`;
    const desc = n.kind === 'local' ? 'local variable of the rule' : n.kind === 'term' ? 'product of the factors below' : (n.doc || '');
    html += `<tr class="${cls}"><td class="bd-term" style="padding-left:${10 + depth * 18}px">${tw}${sg}${label}</td>
      <td class="bd-doc" title="${esc(desc)}">${esc(desc)}</td>${cellHTML(n.value)}<td class="idx-col">${esc(unitOf(n))}</td>${cellHTML(contrib)}`;
    if (B) {
      const b = B.get(n.id), cb = b ? b.contrib : null;
      const dlt = typeof cb === 'number' && typeof contrib === 'number' ? cb - contrib : null;
      html += cellHTML(cb) + cellHTML(dlt, dlt !== null && Math.abs(dlt) > DIFF_RTOL * Math.max(1, Math.abs(contrib)) ? 'cmp-diff' : '');
    }
    html += `${share === null ? '<td></td>' : `<td>${num(share)} %</td>`}<td class="bd-bar">${bar}</td></tr>`;
  }
  html += '</tbody></table>';
  $('tableBody').innerHTML = html;
  $('pager').innerHTML = `<span class="spacer"></span><span>${bdRows(tab, d, true).length} terms · click ▸ to expand, a name to open its table</span>`;
}

// ---------------------------------------------------------------------------
// Formula
// ---------------------------------------------------------------------------

const hasFormula = c => c.tipo === 'Expression' || c.tipo === 'Constraint' || c.tipo === 'Objective';

function ruleSignature(src) {
  const m = src.match(/def\s+\w+\s*\(([^)]*)\)/);
  if (!m) return { blk: null, params: [] };
  const ps = m[1].split(',').map(s => s.trim().split('=')[0].trim()).filter(Boolean);
  return { blk: ps[0] || null, params: ps.slice(1).filter(p => !p.startsWith('*')) };
}

/** blk_ox.rho_g -> 'ox.rho_g' ; m.L -> 'L' ; m.ox.T_hat -> 'ox.T_hat'. Devuelve la clave del
 *  componente y cuantas partes de la cadena consumio, o null. */
function resolveChain(prefix, chain, comp, blkVar) {
  if (NOT_MODEL.has(prefix)) return null;
  for (let len = chain.length; len >= 1; len--) {
    const rest = chain.slice(0, len).join('.');
    const blockFirst = prefix === blkVar;
    const cands = comp.bloque ? (blockFirst ? [`${comp.bloque}.${rest}`, rest] : [rest, `${comp.bloque}.${rest}`]) : [rest];
    const key = cands.find(k => G.byKey.has(k));
    if (key) return { key, used: len };
  }
  return null;
}

/** Texto entre corchetes que empieza en src[pos] === '[' (balanceado), o null. */
function bracketAt(src, pos) {
  if (src[pos] !== '[') return null;
  let depth = 0;
  for (let j = pos; j < src.length && src[j] !== '\n'; j++) {
    if (src[j] === '[' || src[j] === '(') depth++;
    else if (src[j] === ']' || src[j] === ')') { depth--; if (depth === 0) return src.slice(pos + 1, j); }
  }
  return null;
}

function splitArgs(s) {
  const out = []; let depth = 0, start = 0;
  for (let j = 0; j < s.length; j++) {
    const ch = s[j];
    if (ch === '(' || ch === '[') depth++;
    else if (ch === ')' || ch === ']') depth--;
    else if (ch === ',' && depth === 0) { out.push(s.slice(start, j).trim()); start = j + 1; }
  }
  out.push(s.slice(start).trim());
  return out;
}

/** Valor de un argumento de la regla en el contexto de la fila seleccionada:
 *  'z' -> el indice z de esa fila; 'blk.R.last()' / '.first()' / '.at(k)' -> punto del Set. */
function argValue(arg, ctx, comp) {
  if (arg in ctx) return ctx[arg];
  const m = arg.match(/^\w+\.(\w+)\.(first|last)\(\)$/) || arg.match(/^\w+\.(\w+)\.at\((\d+)\)$/);
  if (m) {
    const set = G.sets[`${comp.etapa}.${m[1]}`] || G.sets[m[1]];
    if (!set) return undefined;
    const v = set.valor;
    return m[2] === 'first' ? v[0] : m[2] === 'last' ? v[v.length - 1] : v[+m[2] - 1];
  }
  return undefined;
}

function refLink(key, label, filters) {
  const c = G.byKey.get(key);
  const f = filters && Object.keys(filters).length ? JSON.stringify(filters) : '';
  const href = `/tables?run=${encodeURIComponent(G.run)}&c=${encodeURIComponent(key)}${f ? `&f=${encodeURIComponent(f)}` : ''}`;
  const title = `${c.tipo} · ${c.doc || key}${f ? ` · ${Object.values(filters).join(', ')}` : ''}  (Ctrl+click: new browser tab)`;
  return `<a class="ref ref-${c.tipo}" href="${href}" data-key="${esc(key)}" data-f="${esc(f)}" title="${esc(title)}">${esc(label)}</a>`;
}

function plain(s) { return esc(s).replace(KW, '<span class="kw">$1</span>'); }

function renderRuleHTML(tab) {
  const c = tab.comp, src = c.regla;
  const { blk, params } = ruleSignature(src);
  const ctx = {};
  if (tab.sel !== null) params.forEach((p, j) => { if (j < c.dims.length) ctx[p] = tab.data.idx[j][tab.sel]; });
  const re = /\b([A-Za-z_]\w*)((?:\.[A-Za-z_]\w*)+)/g;
  return src.split('\n').map(line => {
    const hash = line.indexOf('#');
    const code = hash >= 0 ? line.slice(0, hash) : line;
    const comment = hash >= 0 ? line.slice(hash) : '';
    let out = '', last = 0, m;
    re.lastIndex = 0;
    while ((m = re.exec(code))) {
      const chain = m[2].slice(1).split('.');
      const hit = resolveChain(m[1], chain, c, blk);
      if (!hit) continue;
      out += plain(code.slice(last, m.index));
      let filters = null;
      const consumed = `${m[1]}.${chain.slice(0, hit.used).join('.')}`;
      const after = m.index + consumed.length;
      if (hit.used === chain.length) {
        const inside = bracketAt(code, after);
        if (inside !== null && tab.sel !== null) {
          filters = {};
          splitArgs(inside).forEach((a, j) => { const v = argValue(a, ctx, c); if (v !== undefined && v !== null) filters[j] = String(v); });
        }
      }
      out += `${esc(m[1])}.${refLink(hit.key, chain.slice(0, hit.used).join('.'), filters)}`;
      out += plain(code.slice(after, m.index + m[0].length));
      last = m.index + m[0].length;
    }
    out += plain(code.slice(last));
    return out + (comment ? `<span class="cm">${esc(comment)}</span>` : '');
  }).join('\n');
}

let katexLoad = null;
function loadKatex() {
  if (window.katex) return Promise.resolve(true);
  if (!katexLoad) {
    katexLoad = new Promise(resolve => {
      const css = document.createElement('link');
      css.rel = 'stylesheet'; css.href = 'https://cdn.jsdelivr.net/npm/katex@0.16.11/dist/katex.min.css';
      document.head.appendChild(css);
      const js = document.createElement('script');
      js.src = 'https://cdn.jsdelivr.net/npm/katex@0.16.11/dist/katex.min.js';
      js.onload = () => resolve(true); js.onerror = () => resolve(false);
      document.head.appendChild(js);
    });
  }
  return katexLoad;
}

function ctxText(tab) {
  if (tab.sel === null) return 'select a row of the table to resolve the indices of the links';
  const parts = tab.comp.dims.map((_, i) => `${dimLabel(tab, i).main.split(' ')[0]} = ${labelOfIdx(tab, i, tab.sel)}`);
  return `row #${tab.sel + 1}${parts.length ? ' · ' + parts.join(', ') : ''}`;
}

function renderFormula(tab) {
  const c = tab.comp, box = $('formula');
  const views = [];
  if (hasFormula(c)) {
    if (c.latex) views.push(['latex', c.latex_fuente === 'pyomo' ? 'LaTeX (Pyomo)' : 'LaTeX']);
    if (c.regla) views.push(['rule', 'Rule (Python)']);
  }
  if (!tab.fview || !views.some(v => v[0] === tab.fview)) tab.fview = views.length ? views[0][0] : null;
  const refs = [...(G.refsIn.get(c.key) || [])].sort();

  let head = `<div class="formula-head"><span class="tool-label">Formula</span>`;
  if (views.length) {
    head += `<div class="cfg-tabs">${views.map(([v, l]) => `<button class="cfg-tab ${tab.fview === v ? 'active' : ''}" data-fv="${v}">${l}</button>`).join('')}</div>`;
    head += `<span class="formula-ctx">${esc(tab.fview === 'latex' ? (c.latex_fuente === 'pyomo' ? 'generated by Pyomo' : 'from the rule · links are in the Rule view') : ctxText(tab))}</span>`;
  } else {
    head += `<span class="formula-ctx">${c.tipo === 'Var' ? 'a variable has no formula: it is determined by the constraints it appears in'
      : c.tipo === 'Param' ? 'a parameter is data, not a formula'
      : 'no readable rule for this component (built without a rule function, or generated by pyomo.dae)'}</span>`;
  }
  head += `<span class="spacer" style="flex:1"></span>
    ${views.length ? '<button class="btn-sm" id="btnCopyFormula" title="Copy the formula text">⧉ Copy</button>' : ''}
    <button class="btn-sm" id="btnFold">${G.formulaCollapsed ? '▸ Show' : '▾ Hide'}</button></div>`;

  let body = '';
  if (views.length) {
    body = tab.fview === 'rule'
      ? `<div class="formula-body">${renderRuleHTML(tab)}</div>`
      : `<div class="formula-body formula-latex" id="latexBox">${esc(c.latex)}</div>`;
  }
  const appears = `<div class="appears"><span class="tool-label">Appears in</span>${refs.length
    ? refs.map(k => refLink(k, k, null)).join('')
    : '<span class="formula-ctx">no exported rule mentions this component</span>'}</div>`;

  box.className = 'formula' + (G.formulaCollapsed ? ' collapsed' : '');
  box.innerHTML = head + body + appears;
  box.querySelectorAll('[data-fv]').forEach(b => b.onclick = () => { tab.fview = b.dataset.fv; renderFormula(tab); });
  $('btnFold').onclick = () => { G.formulaCollapsed = !G.formulaCollapsed; renderFormula(tab); };
  if ($('btnCopyFormula')) $('btnCopyFormula').onclick = async () => {
    const text = tab.fview === 'rule' ? c.regla : c.latex;
    try { await navigator.clipboard.writeText(text); $('btnCopyFormula').textContent = '✓ Copied'; } catch (_) { $('btnCopyFormula').textContent = '✗'; }
  };

  if (tab.fview === 'latex') {
    loadKatex().then(ok => {
      const el = $('latexBox');
      if (!ok || !el || cur() !== tab) return;
      const tex = c.latex.replace(/\\begin\{equation\*?\}|\\end\{equation\*?\}/g, '').trim();
      try { katex.render(tex, el, { displayMode: true, throwOnError: false }); } catch (_) { /* queda el texto */ }
    });
  }
}

/** Quien menciona a quien, a partir del codigo de las reglas. */
function buildRefs() {
  G.refsIn = new Map();
  const re = /\b([A-Za-z_]\w*)((?:\.[A-Za-z_]\w*)+)/g;
  for (const c of G.comps) {
    if (!c.regla) continue;
    const { blk } = ruleSignature(c.regla);
    const code = c.regla.split('\n').map(l => { const h = l.indexOf('#'); return h >= 0 ? l.slice(0, h) : l; }).join('\n');
    let m;
    re.lastIndex = 0;
    while ((m = re.exec(code))) {
      const hit = resolveChain(m[1], m[2].slice(1).split('.'), c, blk);
      if (!hit || hit.key === c.key) continue;
      if (!G.refsIn.has(hit.key)) G.refsIn.set(hit.key, new Set());
      G.refsIn.get(hit.key).add(c.key);
    }
  }
}

// ---------------------------------------------------------------------------
// Copiar / CSV (precision completa, toda la tabla)
// ---------------------------------------------------------------------------

function exportMatrix(tab) {
  const full = v => v === null || v === undefined ? '' : typeof v === 'boolean' ? (v ? 'TRUE' : 'FALSE') : String(v);
  const c = tab.comp;
  if (tab.view === 'breakdown' && tab.bd && tab.bd.data) {
    const B = tab.bdB && tab.bdB.data ? new Map(bdRows(tab, tab.bdB.data, true).map(r => [r.n.id, r])) : null;
    const head = ['Term', 'Sign', 'Description', 'Value', 'Units', B ? 'Contribution A' : 'Contribution', ...(B ? ['Contribution B', 'Delta B-A'] : [])];
    return [head, ...bdRows(tab, tab.bd.data, true).map(({ n, depth, contrib }) => {
      const b = B && B.get(n.id), cb = b ? b.contrib : null;
      return ['  '.repeat(depth) + n.label, depth ? (n.sign > 0 ? '+' : '-') : '', n.kind === 'local' ? 'local variable' : (n.doc || ''),
              full(n.value), n.unidades || '', full(contrib),
              ...(B ? [full(cb), full(typeof cb === 'number' && typeof contrib === 'number' ? cb - contrib : null)] : [])];
    })];
  }
  if (tab.view === 'pivot' && c.dims.length) {
    const P = tab.pivot, { rowVals, colVals, grid } = pivotGrid(tab);
    const raw = (i, key) => { const r = tab.data.idx[i].findIndex(v => String(v) === key); return full(idxValue(tab, i, r)); };
    const head = [`${dimLabel(tab, P.row).main} \\ ${P.col >= 0 ? dimLabel(tab, P.col).main : ''}`, ...colVals.map(cv => P.col >= 0 ? raw(P.col, cv.key) : cv.label)];
    return [head, ...rowVals.map((rv, a) => [raw(P.row, rv.key), ...grid[a].map(full)])];
  }
  const vcols = valueColumns(tab);
  const head = [...c.dims.map((_, i) => dimLabel(tab, i).main), ...vcols.map(col => col.sub ? `${col.main} (${col.sub})` : col.main)];
  return [head, ...longRows(tab).map(r => [...c.dims.map((_, i) => full(idxValue(tab, i, r))), ...vcols.map(col => full(col.get(r)))])];
}

async function copyTable() {
  const tab = cur();
  if (!tab || !tab.data) return;
  const tsv = exportMatrix(tab).map(row => row.join('\t')).join('\n');
  const b = $('btnCopy');
  try { await navigator.clipboard.writeText(tsv); b.textContent = '✓ Copied'; }
  catch (_) { b.textContent = '✗ Could not copy'; }
  setTimeout(() => { b.textContent = '⧉ Copy for Excel'; }, 1500);
}

function downloadCsv() {
  const tab = cur();
  if (!tab || !tab.data) return;
  const q = s => /[",\n]/.test(s) ? `"${s.replace(/"/g, '""')}"` : s;
  const csv = exportMatrix(tab).map(row => row.map(q).join(',')).join('\r\n');
  const a = document.createElement('a');
  a.href = URL.createObjectURL(new Blob(['﻿' + csv], { type: 'text/csv;charset=utf-8' }));
  a.download = `${G.run.replace(/\.h5$/i, '')}__${tab.key.replace(/[^\w.-]+/g, '_')}${tab.view === 'pivot' ? '_pivot' : ''}.csv`;
  a.click();
  setTimeout(() => URL.revokeObjectURL(a.href), 1000);
}

// ---------------------------------------------------------------------------
// Lista de componentes, hojas y carga
// ---------------------------------------------------------------------------

const isViol = c => c.tipo === 'Constraint' && c.n_viol > 0;

function filteredComps() {
  const q = G.query.toLowerCase();
  const rows = G.comps.filter(c =>
    (G.type === 'all' || c.tipo === G.type) &&
    (G.stage === 'all' || (G.stage === 'global' ? !c.etapa : c.etapa === G.stage)) &&
    (G.type !== 'Constraint' || G.cstat === 'all' || (G.cstat === 'viol' ? isViol(c) : !isViol(c))) &&
    (!q || `${c.key} ${c.doc}`.toLowerCase().includes(q)));
  // violadas: la peor primero
  if (G.type === 'Constraint' && G.cstat === 'viol') rows.sort((a, b) => (b.max_res || 0) - (a.max_res || 0));
  return rows;
}

function renderComps() {
  const counts = { all: G.comps.length };
  G.comps.forEach(c => { counts[c.tipo] = (counts[c.tipo] || 0) + 1; });
  let chips = ['all', ...TYPES].filter(t => counts[t]).map(t =>
    `<button class="chip ${t === G.type ? 'active' : ''}" data-type="${t}">${t === 'all' ? 'All' : t} ${counts[t]}</button>`).join('');
  if (G.type === 'Constraint') {   // segunda fila: violadas / satisfechas
    const cons = G.comps.filter(c => c.tipo === 'Constraint' && (G.stage === 'all' || (G.stage === 'global' ? !c.etapa : c.etapa === G.stage)));
    const nv = cons.filter(isViol).length;
    chips += `<span class="chip-break"></span>` + [['all', 'All', cons.length, ''], ['viol', '✗ Violated', nv, 'chip-viol'], ['ok', '✓ Satisfied', cons.length - nv, 'chip-ok']]
      .map(([k, l, n, cls]) => `<button class="chip ${cls} ${G.cstat === k ? 'active' : ''}" data-cstat="${k}" title="violated = some element with residual > ${VIOL_TOL}">${l} ${n}</button>`).join('');
  }
  $('typeChips').innerHTML = chips;
  $('typeChips').querySelectorAll('[data-type]').forEach(b => b.onclick = () => { G.type = b.dataset.type; renderComps(); });
  $('typeChips').querySelectorAll('[data-cstat]').forEach(b => b.onclick = () => { G.cstat = b.dataset.cstat; renderComps(); });
  const rows = filteredComps(), tab = cur();
  $('compCount').textContent = `${rows.length} of ${G.comps.length}`;
  $('compBody').innerHTML = rows.length ? rows.map(c => `
    <tr data-key="${esc(c.key)}" class="${tab && tab.key === c.key ? 'run-active' : ''}" title="${esc(c.doc || c.key)}">
      <td>${esc(c.key)}</td><td><span class="type-badge type-${c.tipo}">${c.tipo}</span>${isViol(c)
        ? ` <span class="badge-fail" title="${c.n_viol} of ${c.n} elements with residual > ${VIOL_TOL} · worst ${num(c.max_res)}">✗ ${c.n_viol}</span>` : ''}</td><td class="num">${c.n}</td>
    </tr>`).join('') : `<tr><td colspan="3" class="results-empty">No components match.</td></tr>`;
}

/** Abre un componente: en la hoja activa (lista) o en una hoja nueva (enlaces). Si ya hay una
 *  hoja con ese componente, va a ella. `filters` = {posicion de dimension: valor} para
 *  enfocar el elemento (los enlaces de la formula). */
async function openComp(key, { inNewTab = false, filters = null } = {}) {
  if (!G.byKey.has(key)) return;
  let i = G.tabs.findIndex(t => t.key === key);
  if (i < 0) {
    if (inNewTab || G.active < 0) {
      G.tabs.push(newTab(key));
      if (G.tabs.length > MAX_TABS) G.tabs.splice(G.tabs.findIndex((_, j) => j !== G.tabs.length - 1), 1);
      i = G.tabs.length - 1;
    } else {
      G.tabs[G.active] = newTab(key);
      i = G.active;
    }
  }
  G.active = i;
  const tab = cur();
  if (filters) tab.pending = filters;
  renderComps();
  $('compBody').querySelector('tr.run-active')?.scrollIntoView({ block: 'nearest' });
  if (!tab.data) {
    render();
    try {
      const d = await api(`/api/model/${encodeURIComponent(G.run)}/component?key=${encodeURIComponent(key)}`);
      tab.data = d;
      if (tab.comp.dims.length) defaultPivot(tab);
    } catch (e) {
      if (cur() === tab) showTableEmpty(`<span class="error-box">${esc(e.message)}</span>`);
      return;
    }
  }
  applyPending(tab);
  if (cur() === tab) render();
}

function applyPending(tab) {
  if (!tab.pending || !tab.data) return;
  const f = {};
  for (const [i, val] of Object.entries(tab.pending)) {
    if (+i >= tab.comp.dims.length) continue;
    const opt = uniqueVals(tab, +i).find(o => {
      const r = tab.data.idx[+i].findIndex(v => String(v) === o.key);
      return idxMatch(tab, +i, r, val);
    });
    if (opt) f[i] = opt.key;
  }
  tab.pending = null;
  tab.filters = f; tab.view = 'long'; tab.page = 0; tab.sort = null; tab.rowStat = 'all';
  const rows = longRows(tab);
  tab.sel = rows.length === 1 ? rows[0] : null;
}

function closeTab(i) {
  G.tabs.splice(i, 1);
  if (G.active >= G.tabs.length) G.active = G.tabs.length - 1;
  else if (i < G.active) G.active--;
  renderComps();
  render();
}

async function loadRun(name, want) {
  const seq = ++G.seq;
  G.run = name; G.tabs = []; G.active = -1;
  list.render();
  fillCmpSelect();  // B no puede ser la misma corrida que A
  setHeaderStatus('running', 'Loading…');
  $('compBody').innerHTML = `<tr><td colspan="3" class="results-empty">Loading…</td></tr>`;
  render();
  try {
    const d = await api(`/api/model/${encodeURIComponent(name)}/components`);
    if (seq !== G.seq) return;
    G.comps = d.components.sort((a, b) => a.key.localeCompare(b.key));
    G.byKey = new Map(G.comps.map(c => [c.key, c]));
    G.sets = d.sets;
    buildRefs();
    const stages = [...new Set(G.comps.map(c => c.etapa).filter(Boolean))];
    $('compStage').innerHTML = `<option value="all">all stages</option><option value="global">global</option>` +
      stages.map(s => `<option value="${s}">${s}</option>`).join('');
    if (G.stage !== 'all' && G.stage !== 'global' && !stages.includes(G.stage)) G.stage = 'all';
    $('compStage').value = G.stage;
    const nForm = G.comps.filter(c => c.regla).length;
    $('compMeta').textContent = `${(d.meta || {}).duales_disponibles ? 'duals ✓' : 'no duals'} · ${nForm ? `formulas ✓ (${nForm})` : 'no formulas'}`;
    const info = list.get(name);
    const [cls, txt] = statusInfo(info ? info.status : '');
    setHeaderStatus(cls, `${info ? info.label : ''} · ${txt}`);
    $('headerSub').textContent = name;
    renderComps();
    render();
    if (want && want.c && G.byKey.has(want.c)) {
      await openComp(want.c, { filters: want.f });
      if (want.view === 'pivot' && cur()) { cur().view = 'pivot'; render(); }
    }
  } catch (e) {
    if (seq !== G.seq) return;
    $('compBody').innerHTML = `<tr><td colspan="3"><div class="error-box">${esc(e.message)}</div></td></tr>`;
    setHeaderStatus('fail', 'Read error');
  }
}

// ---------------------------------------------------------------------------
// Arranque
// ---------------------------------------------------------------------------

function bindSeg(id, apply) {
  document.querySelectorAll(`#${id} .cfg-tab`).forEach(b => b.onclick = () => {
    if (b.disabled) return;
    apply(b.dataset.v);
    document.querySelectorAll(`#${id} .cfg-tab`).forEach(x => x.classList.toggle('active', x === b));
    render();
  });
}

(async function init() {
  bindExpandButtons();
  bindSeg('segView', v => { if (cur()) cur().view = v; });
  bindSeg('segValues', v => {
    G.values = v;
    // en el pivote, cambiar Pyomo/Physical cambia tambien el campo mostrado
    G.tabs.forEach(t => {
      if (t.data && t.comp.tipo !== 'Constraint' && (t.pivot.field === 'valor' || t.pivot.field === 'valor_fisico')) {
        t.pivot.field = t.data.fields.valor_fisico && v !== 'pyomo' ? 'valor_fisico' : 'valor';
      }
    });
  });
  bindSeg('segIdx', v => { G.idxMode = v; });
  $('selDigits').onchange = e => { G.digits = e.target.value; render(); };
  $('selCmp').addEventListener('focus', fillCmpSelect);  // la lista de corridas se refresca sola
  $('selCmp').onchange = e => {
    G.cmpRun = e.target.value || null;
    e.target.classList.toggle('on', !!G.cmpRun);
    G.tabs.forEach(t => { t.onlyDiff = false; if (t.pivot.field.startsWith('cmp')) t.pivot.field = t.data ? baseField(t) : 'valor'; });
    render();
  };
  $('btnCopy').onclick = copyTable;
  $('btnCsv').onclick = downloadCsv;
  $('compSearch').addEventListener('input', e => { G.query = e.target.value; renderComps(); });
  $('compStage').onchange = e => { G.stage = e.target.value; renderComps(); };
  $('compBody').addEventListener('click', e => {
    const tr = e.target.closest('tr[data-key]');
    if (tr) openComp(tr.dataset.key, { inNewTab: e.ctrlKey || e.metaKey });
  });
  $('sheetTabs').addEventListener('click', e => {
    const x = e.target.closest('[data-close]');
    if (x) { e.stopPropagation(); closeTab(+x.dataset.close); return; }
    const t = e.target.closest('[data-i]');
    if (t) { G.active = +t.dataset.i; renderComps(); render(); }
  });
  $('sheetTabs').addEventListener('auxclick', e => {  // clic con la rueda = cerrar, como en el navegador
    const t = e.target.closest('[data-i]');
    if (t && e.button === 1) closeTab(+t.dataset.i);
  });
  // tabla: enlaces (desglose) -> hoja nueva; ▸/▾ -> expandir; fila -> la formula apunta a ese elemento
  $('tableBody').addEventListener('click', e => {
    const tab = cur();
    if (!tab) return;
    const a = e.target.closest('a.ref');
    if (a) {
      if (e.ctrlKey || e.metaKey || e.shiftKey || e.button !== 0) return;
      e.preventDefault();
      openComp(a.dataset.key, { inNewTab: true });
      return;
    }
    const tw = e.target.closest('[data-tw]');
    if (tw) {
      const id = tw.dataset.tw;
      tab.bdClosed.has(id) ? tab.bdClosed.delete(id) : tab.bdClosed.add(id);
      render();
      return;
    }
    const tr = e.target.closest('tr[data-r]');
    if (!tr) return;
    const r = +tr.dataset.r;
    tab.sel = tab.sel === r ? null : r;
    $('tableBody').querySelectorAll('tr.row-sel').forEach(x => x.classList.remove('row-sel'));
    if (tab.sel !== null) tr.classList.add('row-sel');
    renderFormula(tab);
  });
  // enlaces de la formula: clic = hoja nueva; Ctrl/Shift/rueda = el navegador abre su pestaña
  $('formula').addEventListener('click', e => {
    const a = e.target.closest('a.ref');
    if (!a || e.ctrlKey || e.metaKey || e.shiftKey || e.button !== 0) return;
    e.preventDefault();
    openComp(a.dataset.key, { inNewTab: true, filters: a.dataset.f ? JSON.parse(a.dataset.f) : {} });
  });
  document.addEventListener('keydown', e => {
    const tag = (e.target.tagName || '').toLowerCase();
    const typing = tag === 'input' || tag === 'select' || tag === 'textarea';
    if (e.key === '/' && !typing) { e.preventDefault(); $('compSearch').focus(); return; }
    if (typing || !G.comps.length || (e.key !== 'ArrowDown' && e.key !== 'ArrowUp')) return;
    const rows = filteredComps();
    if (!rows.length) return;
    const tab = cur();
    const i = tab ? rows.findIndex(c => c.key === tab.key) : -1;
    const j = Math.min(Math.max(i + (e.key === 'ArrowDown' ? 1 : -1), 0), rows.length - 1);
    e.preventDefault();
    openComp(rows[j].key);
  });
  render();
  const u = readURL();
  let f = null;
  try { f = u.f ? JSON.parse(u.f) : null; } catch (_) { f = null; }
  const runs = await list.load();
  const cand = runs.filter(r => r.has_model);
  const target = (u.run && cand.find(r => r.name === u.run)) || cand[0];
  if (u.cmp && cand.some(r => r.name === u.cmp)) G.cmpRun = u.cmp;
  if (target) loadRun(target.name, { c: u.c, f, view: u.view });
})();
