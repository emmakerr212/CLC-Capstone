/* =========================================================
   Explorador: una corrida, perfiles en z con sliders de t y r.
   Toda la corrida viene en un solo JSON (h5_runs.load_run), los arreglos con el tiempo
   como primer eje: V(k)[t][z] o V(k)[t][r][z].
   ========================================================= */

const $ = id => document.getElementById(id);
const S = { name: null, t: 0, r: 0, view: 'profiles', timer: null, seq: 0, raf: 0 };
let R = null;

const V = k => R.vars[k];
const bulk = () => S.r === R.BULK;
const stageAt = i => R.stages[R.stage_idx[i]];
const fname = what => `${(S.name || 'clc').replace(/\.h5$/i, '').replace(/[^\w.-]+/g, '_')}_${what}`;
const X_MASK = [-1, 2];  // X fuera de su etapa puede ser basura numerica (division ~0): se oculta

const list = new RunList({
  multi: false,
  onPick: name => { if (name !== S.name) loadRun(name); },
  isSelected: n => n === S.name,
  initialKind: readURL().kind,
});

// ---------------------------------------------------------------------------
// Carga
// ---------------------------------------------------------------------------

function prepare(d) {
  d.nt = d.t.length; d.nr = d.r.length; d.nz = d.z.length; d.BULK = d.nr;
  const L = (d.summary.kpis.find(k => k.label === 'L') || {}).value;
  d.zLabel = (d.z[d.nz - 1] <= 1.0001 && L && parseFloat(L) > 1.5) ? 'z (dimensionless)' : 'z (m)';
  d.allGas = [];
  d.gas_by_stage.forEach(g => g.forEach(s => { if (!d.allGas.includes(s)) d.allGas.push(s); }));
  d.cache = {};
  return d;
}

async function loadRun(name, pos = {}) {
  const seq = ++S.seq;
  stopPlay();
  S.name = name;
  list.render();
  setHeaderStatus('running', 'Loading…');
  try {
    const data = await api('/api/run/' + encodeURIComponent(name));
    if (seq !== S.seq) return;
    R = prepare(data);
  } catch (e) {
    if (seq !== S.seq) return;
    R = null;
    showEmpty(`<div class="error-box">Could not read <b>${esc(name)}</b><br>${esc(e.message)}</div>`);
    setHeaderStatus('fail', 'Read error');
    return;
  }
  const clampInt = (v, lo, hi, dflt) => { const n = parseInt(v, 10); return Number.isFinite(n) ? Math.min(Math.max(n, lo), hi) : dflt; };
  S.t = clampInt(pos.t, 0, R.nt - 1, R.stages[0].i1);   // por defecto: final de la primera etapa
  S.r = clampInt(pos.r, 0, R.BULK, R.nr - 1);           // por defecto: superficie de la particula

  $('tSlider').max = R.nt - 1; $('tSlider').value = S.t;
  $('rSlider').max = R.BULK;   $('rSlider').value = S.r;
  buildStageStrip();
  buildMapOptions();
  renderSummary();

  const [cls, txt] = statusInfo(R.info.status);
  setHeaderStatus(cls, `${R.info.label} · ${txt}`);
  $('headerSub').textContent = R.name;
  showEmpty(null);
  renderAll();
  list.scrollToActive();
}

function showEmpty(html) {
  const on = html !== null;
  $('mainEmpty').hidden = !on;
  if (on) $('mainEmpty').innerHTML = html;
  $('controls').hidden = on;
  document.querySelectorAll('.view').forEach(v => { v.hidden = on; });
  if (on) $('summaryBody').innerHTML = html;
}

// ---------------------------------------------------------------------------
// Controles
// ---------------------------------------------------------------------------

function buildStageStrip() {
  const n = R.nt;
  $('stageStrip').innerHTML = R.stages.map(s =>
    `<span style="flex:${s.i1 - s.i0 + 1};background:${s.color}"></span>`).join('');
  $('stageLabels').innerHTML = R.stages.map(s =>
    `<span style="flex:${s.i1 - s.i0 + 1}" title="${esc(s.name)}: ${fmt(s.t0)}–${fmt(s.t1)} s"><b>${esc(s.name)}</b> ${fmt(s.t0)}–${fmt(s.t1)} s</span>`).join('');
  $('stageStrip').title = `${n} time points`;
}

function updateLabels() {
  const st = stageAt(S.t);
  $('tStage').textContent = st.name;
  $('tStage').style.background = st.color;
  $('tValue').textContent = `t = ${fmt(R.t_local[S.t])} s`;
  $('tSub').textContent = (R.stages.length > 1 ? `global ${fmt(R.t[S.t])} s · ` : '') + `step ${S.t + 1}/${R.nt}`;
  $('rValue').textContent = bulk() ? 'Bulk (gas phase)' : `r = ${R.r[S.r].toExponential(2)} · node ${S.r + 1}/${R.nr}`;
  $('tSlider').value = S.t;
  $('rSlider').value = S.r;
  writeURL({ run: S.name, t: S.t, r: S.r, view: S.view !== 'profiles' ? S.view : '', kind: list.kind !== 'all' ? list.kind : '' });
}

function setT(i) { if (!R) return; S.t = Math.min(Math.max(i, 0), R.nt - 1); scheduleRender(); }
function setR(i) { if (!R) return; S.r = Math.min(Math.max(i, 0), R.BULK); scheduleRender(); }

function scheduleRender() {
  if (S.raf) return;
  S.raf = requestAnimationFrame(() => { S.raf = 0; renderAll(); });
}

function togglePlay() { S.timer ? stopPlay() : startPlay(); }
function startPlay() {
  if (!R) return;
  if (S.t >= R.nt - 1) S.t = 0;
  $('btnPlay').textContent = '❚❚'; $('btnPlay').classList.add('active');
  S.timer = setInterval(() => {
    if (S.t >= R.nt - 1) { stopPlay(); return; }
    setT(S.t + 1);
  }, 160);
}
function stopPlay() {
  clearInterval(S.timer); S.timer = null;
  $('btnPlay').textContent = '▶'; $('btnPlay').classList.remove('active');
}

function setView(v) {
  S.view = v;
  document.querySelectorAll('#viewTabs .cfg-tab').forEach(b => b.classList.toggle('active', b.dataset.view === v));
  document.querySelectorAll('.view').forEach(el => el.classList.toggle('active', el.dataset.view === v));
  if (R) renderAll();
}

function nearestIndex(x) {
  let best = 0, bd = Infinity;
  R.t.forEach((t, i) => { const d = Math.abs(t - x); if (d < bd) { bd = d; best = i; } });
  return best;
}

/** Clic en una grafica contra el tiempo -> mover el slider de t a ese instante. */
function bindTimeClick(div) {
  if (div._timeClick) return;
  div._timeClick = true;
  div.on('plotly_click', ev => { if (ev.points && ev.points.length) { stopPlay(); setT(nearestIndex(ev.points[0].x)); } });
}

// ---------------------------------------------------------------------------
// Graficas — perfiles en z
// ---------------------------------------------------------------------------

function renderAll() {
  if (!R) return;
  updateLabels();
  if (S.view === 'profiles') { renderGas(); renderSolid(); renderTemp(); renderConv(); renderFlow(); }
  else if (S.view === 'map') renderMap();
  else renderOutlet();
}

function memo(key, fn) { return key in R.cache ? R.cache[key] : (R.cache[key] = fn()); }

function rangeOver(keys, i0, i1, cacheKey) {
  return memo(cacheKey, () => extent(keys.filter(k => V(k)).map(k => V(k).slice(i0, i1 + 1))));
}

function note(id, text) { $(id).hidden = !text; $(id).textContent = text || ''; }

function renderGas() {
  const id = 'plotGas', div = $(id);
  const si = R.stage_idx[S.t], st = R.stages[si], gas = R.gas_by_stage[si];
  const traces = gas.map((sp, i) => {
    const arr = bulk() ? V('C_' + sp) : V('Cci_' + sp);
    const y = arr ? (bulk() ? arr[S.t] : arr[S.t][S.r]) : [];
    return line(sp, pretty(sp), R.z, y, colorFor(GAS_COLORS, sp, i), id);
  });
  const shown = gas.filter(sp => !hiddenSet(id).has(sp));
  const keys = (shown.length ? shown : gas).flatMap(sp => ['C_' + sp, 'Cci_' + sp]);
  // rango de la ETAPA: estable mientras el slider se mueve dentro de ella
  const range = rangeOver(keys, st.i0, st.i1, `gas|${si}|${keys}`);
  $('subGas').textContent = `${st.name} · ${bulk() ? 'bulk' : 'inside the particle, r = ' + R.r[S.r].toExponential(2)}`;
  Plotly.react(div, traces, baseLayout(R.zLabel, 'C (mol/m³)', { yaxis: axis('C (mol/m³)', { range }) }), plotConfig(fname('gas')));
  bindLegend(div, renderGas);
}

function renderSolid() {
  const id = 'plotSolid', div = $(id);
  let traces = [], range = [0, 1];
  if (!R.solids.length) note('noteSolid', 'This H5 has no solid concentrations (Ccs_*).');
  else if (bulk()) note('noteSolid', 'Not applicable in bulk: the solid lives inside the particle. Move the radius slider.');
  else {
    note('noteSolid', '');
    traces = R.solids.map((s, i) => line(s, pretty(s), R.z, V('Ccs_' + s)[S.t][S.r], colorFor(SOLID_COLORS, s, i), id));
    const shown = R.solids.filter(s => !hiddenSet(id).has(s));
    const keys = (shown.length ? shown : R.solids).map(s => 'Ccs_' + s);
    range = rangeOver(keys, 0, R.nt - 1, `solid|${keys}`);
  }
  $('subSolid').textContent = bulk() ? '' : `r = ${R.r[S.r].toExponential(2)}`;
  Plotly.react(div, traces, baseLayout(R.zLabel, 'C solid (mol/m³)', { yaxis: axis('C solid (mol/m³)', { range }) }), plotConfig(fname('solid')));
  bindLegend(div, renderSolid);
}

function renderTemp() {
  const id = 'plotTemp', div = $(id);
  const traces = [line('T', 'T bulk', R.z, V('T')[S.t], CAT[0], id)];
  if (!bulk() && V('Tc')) {
    traces.push(line('Tc', `T particle (r = ${R.r[S.r].toExponential(2)})`, R.z, V('Tc')[S.t][S.r], CAT[1], id, { line: { width: 2, color: CAT[1], dash: 'dash' } }));
  }
  const keys = ['T', 'Tc'].filter(k => V(k) && !hiddenSet(id).has(k));
  const range = rangeOver(keys.length ? keys : ['T', 'Tc'], 0, R.nt - 1, `temp|${keys}`);
  $('subTemp').textContent = stageAt(S.t).name;
  Plotly.react(div, traces, baseLayout(R.zLabel, 'T (K)', { yaxis: axis('T (K)', { range }) }), plotConfig(fname('temperature')));
  bindLegend(div, renderTemp);
}

function renderConv() {
  const id = 'plotConv', div = $(id);
  let traces = [], masked = false;
  if (!R.x_keys.length) note('noteConv', 'This H5 has no conversion (X).');
  else if (bulk()) note('noteConv', 'Not applicable in bulk: conversion depends on the particle radius.');
  else {
    note('noteConv', '');
    traces = R.x_keys.map((k, i) => {
      const y = V(k)[S.t][S.r].map(v => {
        if (v === null || v < X_MASK[0] || v > X_MASK[1]) { if (v !== null) masked = true; return null; }
        return v;
      });
      const label = k === 'X' ? 'X' : `X ${pretty(k.slice(2))}`;
      return line(k, label, R.z, y, colorFor(SOLID_COLORS, k.slice(2), i), id);
    });
  }
  $('subConv').textContent = `${stageAt(S.t).name}${masked ? ' · values outside [−1, 2] hidden' : ''}`;
  Plotly.react(div, traces, baseLayout(R.zLabel, 'X (–)', { yaxis: axis('X (–)', { range: [-0.05, 1.05] }) }), plotConfig(fname('conversion')));
  bindLegend(div, renderConv);
}

function renderFlow() {
  const id = 'plotFlow', div = $(id);
  const traces = R.feed_species.map((sp, i) =>
    line(sp, pretty(sp), R.t, V('Fin_' + sp), colorFor(GAS_COLORS, sp, i), id, { line: { width: 2, color: colorFor(GAS_COLORS, sp, i), shape: 'vh' } }));
  if (V('m_total') && R.feed_species.length > 1) {
    traces.push(line('total', 'Total', R.t, V('m_total'), INK.muted, id, { line: { width: 2, color: INK.muted, dash: 'dash', shape: 'vh' } }));
  }
  note('noteFlow', traces.length ? '' : 'This H5 has no inlet mass flow (m_total / w_in_*).');
  const { shapes, annotations } = stageShapes(R.stages, R.t[S.t]);
  Plotly.react(div, traces, baseLayout('global t (s)', 'Flow (kg/s)', {
    yaxis: axis('Flow (kg/s)', { rangemode: 'tozero' }), shapes, annotations,
  }), plotConfig(fname('inlet_flow')));
  bindLegend(div, renderFlow);
  bindTimeClick(div);
}

// ---------------------------------------------------------------------------
// Mapa z–t
// ---------------------------------------------------------------------------

function buildMapOptions() {
  const sel = $('mapVar'), prev = sel.value;
  const groups = [
    ['Temperature', [['T', 'T bulk'], ['Tc', 'T particle']]],
    ['Gas — bulk', R.allGas.map(s => ['C_' + s, pretty(s)])],
    ['Gas — in particle', R.allGas.map(s => ['Cci_' + s, pretty(s) + ' (particle)'])],
    ['Solid', R.solids.map(s => ['Ccs_' + s, pretty(s)])],
    ['Conversion', R.x_keys.map(k => [k, k === 'X' ? 'X' : 'X ' + pretty(k.slice(2))])],
  ];
  sel.innerHTML = groups.map(([g, opts]) => {
    const ok = opts.filter(([k]) => V(k));
    return ok.length ? `<optgroup label="${g}">${ok.map(([k, l]) => `<option value="${k}">${l}</option>`).join('')}</optgroup>` : '';
  }).join('');
  sel.value = [...sel.options].some(o => o.value === prev) ? prev : 'T';
}

const UNITS = k => k.startsWith('X') ? 'X (–)' : (k === 'T' || k === 'Tc') ? 'T (K)' : 'C (mol/m³)';

function renderMap() {
  const div = $('plotMap'), key = $('mapVar').value, arr = V(key);
  if (!arr) return;
  const is3D = Array.isArray(arr[0][0]) || (arr.find(a => a && Array.isArray(a[0])) !== undefined);
  const rIdx = is3D ? (bulk() ? R.nr - 1 : S.r) : null;
  const isX = key.startsWith('X');
  const z = R.z.map((_, iz) => R.t.map((_, it) => {
    const row = is3D ? (arr[it] && arr[it][rIdx]) : arr[it];
    let v = row ? row[iz] : null;
    if (isX && v !== null && (v < X_MASK[0] || v > X_MASK[1])) v = null;
    return v;
  }));
  const trace = {
    type: 'heatmap', x: R.t, y: R.z, z, colorscale: SEQ_SCALE, zsmooth: false,
    zmin: isX ? 0 : undefined, zmax: isX ? 1 : undefined,
    colorbar: { title: { text: UNITS(key), side: 'right', font: { size: 11, color: INK.muted } }, thickness: 12, outlinewidth: 0, tickfont: { size: 10 } },
    hovertemplate: 't = %{x:.4g} s<br>z = %{y:.3g}<br>%{z:.4g}<extra></extra>',
  };
  const { shapes, annotations } = stageShapes(R.stages, R.t[S.t]);
  shapes.forEach(s => { if (s.type === 'rect') { s.opacity = 0; } });
  R.stages.slice(0, -1).forEach(s => shapes.push({ type: 'line', xref: 'x', yref: 'paper', x0: s.t1, x1: s.t1, y0: 0, y1: 1, line: { color: '#ffffff', width: 2 } }));
  const where = is3D ? (bulk() ? 'r = surface (not applicable in bulk, using the last node)' : `r = ${R.r[rIdx].toExponential(2)}`) : 'bulk';
  $('subMap').textContent = where + (R.stages.length > 1 ? ' · blank where the variable does not exist in that stage' : '');
  Plotly.react(div, [trace], baseLayout('global t (s)', R.zLabel, {
    hovermode: 'closest', showlegend: false, shapes, annotations, margin: { l: 56, r: 14, t: 22, b: 40 },
  }), plotConfig(fname('map_' + key)));
  bindTimeClick(div);
}

// ---------------------------------------------------------------------------
// Salida del lecho
// ---------------------------------------------------------------------------

function renderOutlet() {
  const iL = R.nz - 1, log = $('outScale').value === 'log';
  const { shapes, annotations } = stageShapes(R.stages, R.t[S.t]);

  const idG = 'plotOutGas', divG = $(idG);
  const tracesG = R.allGas.filter(sp => V('C_' + sp)).map((sp, i) =>
    line(sp, pretty(sp), R.t, V('C_' + sp).map(row => {
      const v = row ? row[iL] : null;
      return log && v !== null && v <= 0 ? null : v;
    }), colorFor(GAS_COLORS, sp, i), idG, { connectgaps: false }));
  $('subOutGas').textContent = `bulk at z = L (${fmt(R.z[iL])})` + (log ? ' · values ≤ 0 omitted' : '');
  Plotly.react(divG, tracesG, baseLayout('global t (s)', 'C (mol/m³)', {
    yaxis: axis('C (mol/m³)', { type: log ? 'log' : 'linear', exponentformat: 'power' }), shapes, annotations,
  }), plotConfig(fname('outlet_gas')));
  bindLegend(divG, renderOutlet);
  bindTimeClick(divG);

  const idT = 'plotOutT', divT = $(idT);
  const tracesT = [
    line('Tout', 'T outlet (z = L)', R.t, V('T').map(row => row ? row[iL] : null), CAT[0], idT),
    line('Tin', 'T inlet (z = 0)', R.t, V('T').map(row => row ? row[0] : null), INK.muted, idT, { line: { width: 2, color: INK.muted, dash: 'dash' } }),
  ];
  $('subOutT').textContent = 'bulk';
  Plotly.react(divT, tracesT, baseLayout('global t (s)', 'T (K)', { shapes, annotations }), plotConfig(fname('T_outlet')));
  bindLegend(divT, renderOutlet);
  bindTimeClick(divT);
}

// ---------------------------------------------------------------------------
// Resumen
// ---------------------------------------------------------------------------

function renderSummary() {
  const sm = R.summary;
  const tone = { good: 'is-good', warn: 'is-warn', danger: 'is-danger' };
  let html = `<div class="file-name">${esc(R.name)}</div>`;
  html += `<div class="kpi-grid">${sm.kpis.map(k => `
    <div class="kpi-tile ${tone[k.tone] || ''}" title="${esc(k.label)}: ${esc(k.value)} ${esc(k.sub)}">
      <div class="kpi-label">${esc(k.label)}</div>
      <div class="kpi-value">${esc(k.value)}</div>
      <div class="kpi-sub">${esc(k.sub) || '&nbsp;'}</div>
    </div>`).join('')}</div>`;

  if (sm.decisions.length) {
    html += `<div class="section-head">Decisions</div>
      <table class="results-table"><thead><tr><th>Variable</th><th style="text-align:right">Value</th><th>Status</th></tr></thead><tbody>
      ${sm.decisions.map(d => `<tr><td class="run-date">${esc(d.name)}</td><td class="num">${fmt(d.value, 6)}</td>
        <td>${d.free === null ? '' : d.free ? '<span class="badge-ok">free</span>' : '<span class="badge-neutral">fixed</span>'}</td></tr>`).join('')}
      </tbody></table>`;
  }

  html += `<div class="section-head">Stages</div>
    <table class="results-table"><thead><tr><th>Stage</th><th style="text-align:right">global t (s)</th><th style="text-align:right">t nodes</th><th>Gas</th></tr></thead><tbody>
    ${R.stages.map((s, i) => `<tr><td><span class="swatch" style="background:${s.color}"></span> ${esc(s.name)}</td>
      <td class="num">${fmt(s.t0)}–${fmt(s.t1)}</td><td class="num">${s.i1 - s.i0 + 1}</td>
      <td>${R.gas_by_stage[i].map(pretty).join(', ')}</td></tr>`).join('')}
    </tbody></table>`;

  html += `<div class="section-head">H5 attributes</div>`;
  html += sm.groups.map(g => `<details class="form-group" ${g.title === 'Decisions' ? 'open' : ''}>
      <summary>${esc(g.title)} <span style="font-weight:400">(${g.items.length})</span></summary>
      <table class="attr-table"><tbody>${g.items.map(([k, v]) => `<tr><td>${esc(k)}</td><td>${esc(fmt(v, 6))}</td></tr>`).join('')}</tbody></table>
    </details>`).join('');

  html += `<div class="summary-actions">
      <button class="btn-sm btn-load" id="btnReveal">📁 Show in folder</button>
      <button class="btn-sm" id="btnCopy">⧉ Copy name</button>
      <a class="btn-sm" href="/compare?runs=${encodeURIComponent(R.name)}" style="text-decoration:none">⚖ Compare with…</a>
    </div>`;
  $('summaryBody').innerHTML = html;
  $('btnReveal').onclick = () => api('/api/reveal', { method: 'POST', headers: { 'Content-Type': 'application/json' }, body: JSON.stringify({ name: R.name }) }).catch(e => alert(e.message));
  $('btnCopy').onclick = async () => {
    try { await navigator.clipboard.writeText(R.name); $('btnCopy').textContent = '✓ Copied'; }
    catch (_) { $('btnCopy').textContent = '✗ Could not copy'; }
    setTimeout(() => { const b = $('btnCopy'); if (b) b.textContent = '⧉ Copy name'; }, 1500);
  };
}

// ---------------------------------------------------------------------------
// Arranque
// ---------------------------------------------------------------------------

function bindControls() {
  $('tSlider').addEventListener('input', e => { stopPlay(); setT(+e.target.value); });
  $('rSlider').addEventListener('input', e => setR(+e.target.value));
  $('btnPrev').onclick = () => { stopPlay(); setT(S.t - 1); };
  $('btnNext').onclick = () => { stopPlay(); setT(S.t + 1); };
  $('btnPlay').onclick = togglePlay;
  $('mapVar').onchange = () => renderMap();
  $('outScale').onchange = () => renderOutlet();
  document.querySelectorAll('#viewTabs .cfg-tab').forEach(b => b.onclick = () => setView(b.dataset.view));

  document.addEventListener('keydown', e => {
    const tag = (e.target.tagName || '').toLowerCase();
    if (!R || tag === 'select' || tag === 'textarea' || (tag === 'input' && e.target.type !== 'range')) return;
    const step = e.shiftKey ? 5 : 1;
    if (e.key === 'ArrowRight') { stopPlay(); setT(S.t + step); }
    else if (e.key === 'ArrowLeft') { stopPlay(); setT(S.t - step); }
    else if (e.key === 'ArrowUp') setR(S.r + 1);
    else if (e.key === 'ArrowDown') setR(S.r - 1);
    else if (e.key === ' ') togglePlay();
    else if (e.key === 'Home') setT(0);
    else if (e.key === 'End') setT(R.nt - 1);
    else return;
    e.preventDefault();
  });
}

(async function init() {
  bindExpandButtons();
  bindControls();
  showEmpty('<span class="ph-icon">📈</span><p>Pick a run from the list to see its profiles.</p>');
  $('summaryBody').innerHTML = '<div class="run-report-empty"><span class="ph-icon">📋</span><p>Pick a run from the list.</p></div>';
  const u = readURL();
  if (u.view) setView(u.view);
  const runs = await list.load();
  const target = (u.run && runs.find(r => r.name === u.run)) || runs.find(r => !r.error);
  if (target) loadRun(target.name, { t: u.t, r: u.r });
})();
