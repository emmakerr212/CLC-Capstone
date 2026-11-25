/* =========================================================
   Comparar: hasta 4 corridas superpuestas contra el tiempo global.
   Cada corrida conserva su color mientras siga seleccionada (el color sigue a la corrida,
   no a su posicion): quitar una no repinta las demas.
   ========================================================= */

const $ = id => document.getElementById(id);
const MAX_RUNS = 4;
const LETTERS = ['A', 'B', 'C', 'D'];
const C = { slots: [null, null, null, null], data: {}, seq: 0 };

const selected = () => C.slots.filter(Boolean);
const slotOf = name => C.slots.indexOf(name);

const list = new RunList({
  multi: true,
  onPick: name => toggle(name),
  isSelected: n => slotOf(n) >= 0,
  colorOf: n => RUN_COLORS[slotOf(n)],
});

function toggle(name) {
  const run = list.get(name);
  if (!run || run.error) return;
  const i = slotOf(name);
  if (i >= 0) C.slots[i] = null;
  else {
    const free = C.slots.indexOf(null);
    if (free < 0) { alert(`At most ${MAX_RUNS} runs: remove one before adding another.`); return; }
    C.slots[free] = name;
  }
  list.render();
  refresh();
}

function trapzWeights(x) {
  const n = x.length, w = new Array(n).fill(0);
  for (let i = 0; i < n - 1; i++) { const h = (x[i + 1] - x[i]) / 2; w[i] += h; w[i + 1] += h; }
  return w;
}

function prepareRun(d) {
  d.nt = d.t.length; d.nr = d.r.length; d.nz = d.z.length;
  d.allGas = [];
  d.gas_by_stage.forEach(g => g.forEach(s => { if (!d.allGas.includes(s)) d.allGas.push(s); }));
  // promedio volumetrico en la particula (peso r^2) y promedio axial (trapecio en z)
  const wr = trapzWeights(d.r).map((w, i) => w * d.r[i] * d.r[i]);
  const sr = wr.reduce((a, b) => a + b, 0) || 1;
  const wz = trapzWeights(d.z), sz = wz.reduce((a, b) => a + b, 0) || 1;
  d.solidAvg = {};
  d.solids.forEach(s => {
    const A = d.vars['Ccs_' + s];
    d.solidAvg[s] = A.map(tr => {
      if (!tr) return null;
      let acc = 0;
      for (let iz = 0; iz < d.nz; iz++) {
        let rad = 0;
        for (let ir = 0; ir < d.nr; ir++) {
          const v = tr[ir] ? tr[ir][iz] : null;
          if (v === null) return null;
          rad += wr[ir] * v;
        }
        acc += wz[iz] * rad / sr;
      }
      return acc / sz;
    });
  });
  return d;
}

async function refresh() {
  const seq = ++C.seq;
  const names = selected();
  writeURL({ runs: names.join('|') });
  if (!names.length) { show(false, 'Check 2 to 4 runs in the list.'); setHeaderStatus('neutral', 'No runs'); return; }
  setHeaderStatus('running', 'Loading…');
  const missing = names.filter(n => !C.data[n]);
  try {
    const loaded = await Promise.all(missing.map(n => api('/api/run/' + encodeURIComponent(n))));
    loaded.forEach((d, i) => { C.data[missing[i]] = prepareRun(d); });
  } catch (e) {
    if (seq !== C.seq) return;
    show(false, `<div class="error-box">${esc(e.message)}</div>`);
    setHeaderStatus('fail', 'Read error');
    return;
  }
  if (seq !== C.seq) return;
  setHeaderStatus('ok', `${names.length} run${names.length > 1 ? 's' : ''}`);
  show(true);
  buildSelectors();
  renderAll();
}

function show(on, msg) {
  $('cmpContent').hidden = !on;
  $('mainEmpty').hidden = on;
  if (!on && msg) $('mainEmpty').innerHTML = msg.startsWith('<') ? msg : `<span class="ph-icon">⚖</span><p>${msg}</p>`;
}

function runs() { return C.slots.map((n, i) => n && { name: n, i, d: C.data[n], color: RUN_COLORS[i] }).filter(x => x && x.d); }
function runLabel(r) { return `${LETTERS[r.i]} · ${fmtDate(r.d.info.date)} · ${r.d.info.label}`; }

function fillSelect(sel, values, labelFn, preferred) {
  const prev = sel.value;
  sel.innerHTML = values.map(v => `<option value="${v}">${labelFn(v)}</option>`).join('');
  sel.value = values.includes(prev) ? prev : (values.find(v => v === preferred) || values[0] || '');
}

function buildSelectors() {
  const gas = [], sol = [];
  runs().forEach(r => {
    r.d.allGas.forEach(s => { if (!gas.includes(s)) gas.push(s); });
    r.d.solids.forEach(s => { if (!sol.includes(s)) sol.push(s); });
  });
  fillSelect($('selSpecies'), gas, pretty, 'O2');
  fillSelect($('selSolid'), sol, pretty, 'Fe2O3');
}

function renderAll() {
  const rs = runs();
  $('cmpLegend').innerHTML = rs.map(r =>
    `<span><span class="swatch" style="background:${r.color}"></span><b>${LETTERS[r.i]}</b> ${esc(r.d.name)}</span>`).join('');
  renderSeries('cmpTout', rs, r => r.d.vars.T.map(row => row ? row[r.d.nz - 1] : null), 'T (K)');
  const sp = $('selSpecies').value;
  $('subCout').textContent = sp ? `${pretty(sp)} bulk at z = L` : '';
  renderSeries('cmpCout', rs, r => {
    const A = r.d.vars['C_' + sp];
    return A ? A.map(row => row ? row[r.d.nz - 1] : null) : null;
  }, 'C (mol/m³)');
  renderSeries('cmpFlow', rs, r => r.d.vars.m_total || null, 'Flow (kg/s)', { line: { shape: 'vh' } }, { rangemode: 'tozero' });
  const s = $('selSolid').value;
  $('subSolidAvg').textContent = s ? `${pretty(s)}, averaged over the particle (r² weight) and over z` : '';
  renderSeries('cmpSolid', rs, r => r.d.solidAvg[s] || null, 'C solid (mol/m³)');
  renderTable(rs);
}

function renderSeries(id, rs, yOf, yTitle, traceExtra = {}, yExtra = {}) {
  const div = $(id);
  const traces = [];
  rs.forEach(r => {
    const y = yOf(r);
    if (!y) return;
    const t = line(r.name, runLabel(r), r.d.t, y, r.color, id, { connectgaps: false, hovertemplate: '%{y:.4g}' });
    if (traceExtra.line) Object.assign(t.line, traceExtra.line);
    traces.push(t);
  });
  const layout = baseLayout('global t (s)', yTitle, { yaxis: axis(yTitle, yExtra), hovermode: 'x' });
  if (!traces.length) layout.annotations = [{ text: 'None of the runs has this variable', xref: 'paper', yref: 'paper', x: 0.5, y: 0.5, showarrow: false, font: { color: INK.muted } }];
  Plotly.react(div, traces, layout, plotConfig('comparison_' + id));
  bindLegend(div, renderAll);
}

function renderTable(rs) {
  const rows = [];  // [etiqueta, [valor por corrida]]
  const add = (label, vals) => rows.push([label, vals]);
  add('File', rs.map(r => r.d.name));
  add('Type', rs.map(r => r.d.info.label + (r.d.info.tipo ? ` (${r.d.info.tipo})` : '')));
  add('Status', rs.map(r => statusInfo(r.d.info.status)[1]));
  const kpiLabels = [];
  rs.forEach(r => r.d.summary.kpis.forEach(k => { if (k.label !== 'Solver' && !kpiLabels.includes(k.label)) kpiLabels.push(k.label); }));
  kpiLabels.forEach(l => add(l, rs.map(r => { const k = r.d.summary.kpis.find(x => x.label === l); return k ? `${k.value} ${k.sub}`.trim() : '—'; })));
  const decKeys = [];
  rs.forEach(r => (r.d.summary.groups.find(g => g.title === 'Decisions') || { items: [] }).items.forEach(([k]) => { if (!decKeys.includes(k)) decKeys.push(k); }));
  decKeys.forEach(k => add(k, rs.map(r => {
    const g = r.d.summary.groups.find(x => x.title === 'Decisions');
    const it = g && g.items.find(([kk]) => kk === k);
    return it ? fmt(it[1], 6) : '—';
  })));

  $('cmpTable').innerHTML = `<thead><tr><th></th>${rs.map(r => `<th><span class="swatch" style="background:${r.color}"></span>${LETTERS[r.i]}</th>`).join('')}</tr></thead>
    <tbody>${rows.map(([label, vals]) => {
      const diff = rs.length > 1 && label !== 'File' && new Set(vals).size > 1;
      return `<tr><td>${esc(label)}</td>${vals.map(v => `<td class="${diff ? 'diff' : ''} ${label === 'File' ? 'run-date' : ''}">${esc(v)}</td>`).join('')}</tr>`;
    }).join('')}</tbody>`;
}

(async function init() {
  bindExpandButtons();
  $('selSpecies').onchange = renderAll;
  $('selSolid').onchange = renderAll;
  $('btnClear').onclick = () => { C.slots = [null, null, null, null]; list.render(); refresh(); };
  const all = await list.load();
  const want = (readURL().runs || '').split('|').filter(n => all.some(r => r.name === n && !r.error)).slice(0, MAX_RUNS);
  want.forEach((n, i) => { C.slots[i] = n; });
  list.render();
  refresh();
})();
