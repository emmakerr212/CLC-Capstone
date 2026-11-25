/* =========================================================
   Comun a Explorador y Comparar: colores, tema de Plotly, lista de corridas,
   leyendas con estado y paneles expandibles.
   ========================================================= */

// Paleta categorica (orden fijo, validado para pares adyacentes). El color sigue a la
// entidad: una especie tiene SIEMPRE el mismo color, este en la etapa que este.
const CAT = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300', '#6250d6', '#e34948'];
const GAS_COLORS = { CH4: CAT[0], CO: CAT[1], H2: CAT[2], N2: CAT[3], CO2: CAT[4], H2O: CAT[5], O2: CAT[7], gas: CAT[0] };
const SOLID_COLORS = { Fe3O4: CAT[0], Fe2O3: CAT[1] };
const RUN_COLORS = [CAT[0], CAT[1], CAT[2], CAT[6]];
const INK = { primary: '#1a1f2e', mid: '#374151', muted: '#6b7280', grid: '#eef0f3', axis: '#d1d5db' };
const SEQ_SCALE = [[0, '#cde2fb'], [0.17, '#9ec5f4'], [0.33, '#6da7ec'], [0.5, '#3987e5'],
                   [0.67, '#256abf'], [0.83, '#184f95'], [1, '#0d366b']];

const SUB = { O2: 'O₂', N2: 'N₂', CH4: 'CH₄', CO2: 'CO₂', H2: 'H₂', H2O: 'H₂O', CO: 'CO',
              Fe3O4: 'Fe₃O₄', Fe2O3: 'Fe₂O₃', gas: 'gas' };
const pretty = s => SUB[s] || s;

function colorFor(table, key, i) { return table[key] || CAT[i % CAT.length]; }

// ---------------------------------------------------------------------------
// Utilidades
// ---------------------------------------------------------------------------

async function api(url, opts) {
  const res = await fetch(url, opts);
  let body = null;
  try { body = await res.json(); } catch (_) { /* sin cuerpo */ }
  if (!res.ok) throw new Error((body && body.error) || `${res.status} ${res.statusText}`);
  return body;
}

function esc(s) {
  return String(s ?? '').replace(/[&<>"']/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c]));
}

function fmt(v, digits = 4) {
  if (v === null || v === undefined || v === '') return '—';
  if (typeof v === 'boolean') return v ? 'yes' : 'no';
  if (typeof v !== 'number') return Array.isArray(v) ? v.map(x => fmt(x, digits)).join(', ') : String(v);
  if (!isFinite(v)) return '—';
  if (v !== 0 && (Math.abs(v) >= 1e6 || Math.abs(v) < 1e-3)) return v.toExponential(Math.max(digits - 2, 1));
  return Number(v.toPrecision(digits)).toLocaleString('en-US', { maximumFractionDigits: 6 });
}

function fmtDate(iso) { return iso ? iso.replace('T', ' ').slice(0, 16) : '—'; }

/** min/max de arreglos anidados ignorando null, con margen; null si no hay datos. */
function extent(arrays, padFrac = 0.05) {
  let lo = Infinity, hi = -Infinity;
  const walk = a => {
    if (a === null || a === undefined) return;
    if (Array.isArray(a)) { for (const x of a) walk(x); return; }
    if (a < lo) lo = a;
    if (a > hi) hi = a;
  };
  arrays.forEach(walk);
  if (!isFinite(lo)) return null;
  const pad = hi > lo ? (hi - lo) * padFrac : (Math.abs(hi) * 0.1 || 1);
  return [lo - pad, hi + pad];
}

const STATUS_MAP = {
  optimal: ['ok', 'optimal'], ok: ['ok', 'ok'], locallyOptimal: ['ok', 'locally optimal'],
  infeasible: ['fail', 'infeasible'], error: ['fail', 'unreadable'],
  maxIterations: ['warn', 'max. iter.'], warning: ['warn', 'warning'],
};
function statusInfo(s) { return STATUS_MAP[s] || (s ? ['warn', s] : ['neutral', '—']); }
function statusBadge(s) { const [cls, txt] = statusInfo(s); return `<span class="badge-${cls}">${esc(txt)}</span>`; }

function setHeaderStatus(status, text) {
  const badge = document.getElementById('statusBadge');
  if (!badge) return;
  const cls = { ok: 'done', fail: 'error', warn: 'warn', neutral: '', running: 'running' }[status] ?? '';
  badge.className = 'status-badge ' + cls;
  document.getElementById('statusText').textContent = text;
}

function readURL() { return Object.fromEntries(new URLSearchParams(location.search)); }
function writeURL(params) {
  const q = new URLSearchParams();
  for (const [k, v] of Object.entries(params)) if (v !== undefined && v !== null && v !== '') q.set(k, v);
  history.replaceState(null, '', `${location.pathname}?${q}`);
}

// ---------------------------------------------------------------------------
// Plotly: tema, leyendas con estado, sombreado de etapas
// ---------------------------------------------------------------------------

function axis(title, extra = {}) {
  return Object.assign({
    title: { text: title, font: { size: 11, color: INK.muted } },
    gridcolor: INK.grid, zeroline: false, linecolor: INK.axis, ticks: 'outside', tickcolor: INK.axis,
    ticklen: 4, automargin: true, tickfont: { size: 10.5, color: INK.mid },
  }, extra);
}

function baseLayout(xTitle, yTitle, extra = {}) {
  return Object.assign({
    font: { family: "'Segoe UI', system-ui, sans-serif", size: 11, color: INK.mid },
    paper_bgcolor: '#ffffff', plot_bgcolor: '#ffffff',
    margin: { l: 56, r: 14, t: 26, b: 40 },
    xaxis: axis(xTitle), yaxis: axis(yTitle),
    legend: { orientation: 'h', x: 0, y: 1.02, yanchor: 'bottom', font: { size: 11, color: INK.mid } },
    hovermode: 'x unified',
    hoverlabel: { bgcolor: '#ffffff', bordercolor: INK.axis, font: { size: 11, color: INK.primary } },
    showlegend: true,
  }, extra);
}

function plotConfig(filename) {
  return {
    displaylogo: false, responsive: true,
    modeBarButtonsToRemove: ['lasso2d', 'select2d', 'autoScale2d', 'toggleSpikelines'],
    toImageButtonOptions: { format: 'png', scale: 3, filename: filename || 'clc_plot' },
  };
}

/** Curvas ocultas por grafica: clic en la leyenda = mostrar/ocultar, doble clic = aislar.
 *  Se guarda aqui (y no en Plotly) para que sobreviva a cada Plotly.react y para poder
 *  ajustar el rango del eje solo a lo visible. */
const HIDDEN = {};
function hiddenSet(id) { return HIDDEN[id] || (HIDDEN[id] = new Set()); }
function vis(id, key) { return hiddenSet(id).has(key) ? 'legendonly' : true; }

function bindLegend(div, rerender) {
  if (div._legendBound) return;
  div._legendBound = true;
  const id = div.id;
  let clickTimer = null;
  div.on('plotly_legendclick', ev => {
    const key = ev.data[ev.curveNumber].meta;
    if (key === undefined) return false;
    clearTimeout(clickTimer);
    clickTimer = setTimeout(() => {  // esperar por si es doble clic
      const h = hiddenSet(id);
      h.has(key) ? h.delete(key) : h.add(key);
      rerender();
    }, 320);  // > doubleClickDelay de Plotly (300 ms)
    return false;
  });
  div.on('plotly_legenddoubleclick', ev => {
    clearTimeout(clickTimer);
    const key = ev.data[ev.curveNumber].meta;
    const keys = ev.data.map(t => t.meta).filter(k => k !== undefined);
    const h = hiddenSet(id);
    const isolated = keys.every(k => k === key ? !h.has(k) : h.has(k));
    h.clear();
    if (!isolated) keys.forEach(k => { if (k !== key) h.add(k); });
    rerender();
    return false;
  });
}

function stageShapes(stages, tCursor) {
  const shapes = [], annotations = [];
  if (stages.length > 1) {
    stages.forEach(s => {
      shapes.push({ type: 'rect', xref: 'x', yref: 'paper', x0: s.t0, x1: s.t1, y0: 0, y1: 1,
                    fillcolor: s.color, opacity: 0.06, line: { width: 0 }, layer: 'below' });
      annotations.push({ x: (s.t0 + s.t1) / 2, y: 1, xref: 'x', yref: 'paper', yanchor: 'top', showarrow: false,
                         text: s.name, font: { size: 10.5, color: INK.muted } });
    });
  }
  if (tCursor !== undefined && tCursor !== null) {
    shapes.push({ type: 'line', xref: 'x', yref: 'paper', x0: tCursor, x1: tCursor, y0: 0, y1: 1,
                  line: { color: INK.mid, width: 1.5, dash: 'dot' } });
  }
  return { shapes, annotations };
}

function line(key, name, x, y, color, plotId, extra = {}) {
  const nPts = x.length;
  return Object.assign({
    type: 'scatter', mode: nPts <= 30 ? 'lines+markers' : 'lines',
    x, y, name, meta: key, visible: vis(plotId, key),
    line: { width: 2, color }, marker: { size: 7, color, line: { width: 1.5, color: '#ffffff' } },
    hovertemplate: '%{y:.4g}',
  }, extra);
}

// ---------------------------------------------------------------------------
// Paneles expandibles
// ---------------------------------------------------------------------------

function bindExpandButtons() {
  document.querySelectorAll('[data-expand]').forEach(btn => {
    btn.addEventListener('click', () => {
      const panel = document.getElementById(btn.dataset.expand);
      const on = panel.classList.toggle('panel-expanded');
      btn.classList.toggle('active', on);
      btn.textContent = on ? '⤡' : '⤢';
      btn.title = on ? 'Restore' : 'Expand';
      panel.querySelectorAll('.js-plotly-plot').forEach(p => Plotly.Plots.resize(p));
    });
  });
  document.addEventListener('keydown', e => {
    if (e.key !== 'Escape') return;
    document.querySelectorAll('.panel-expanded').forEach(p => {
      p.querySelector('[data-expand]')?.click();
    });
  });
}

// ---------------------------------------------------------------------------
// Lista de corridas (filtros por tipo, busqueda, refresco automatico)
// ---------------------------------------------------------------------------

const KINDS = [
  ['all', 'All'], ['cycle', 'Full cycle'], ['oxpurge', 'Ox + Purge'], ['ox', 'Oxidation'],
  ['purge', 'Purge'], ['red', 'Reduction'], ['error', 'Unreadable'],
];

class RunList {
  /** opts: {multi, onPick(name), isSelected(name), colorOf(name), initialKind} */
  constructor(opts) {
    this.opts = opts;
    this.runs = [];
    this.kind = opts.initialKind || 'all';
    this.query = '';
    this.known = null;   // nombres vistos en la primera carga; lo que aparezca despues = "nuevo"
    this.fresh = new Set();
    this.body = document.getElementById('runsBody');
    this.chips = document.getElementById('kindChips');
    document.getElementById('runSearch').addEventListener('input', e => { this.query = e.target.value.toLowerCase(); this.render(); });
    document.getElementById('btnRefresh').addEventListener('click', () => this.load());
    this.body.addEventListener('click', e => {
      const tr = e.target.closest('tr[data-name]');
      if (tr) this.opts.onPick(tr.dataset.name, e);
    });
    setInterval(() => { if (!document.hidden) this.load(true); }, 15000);
  }

  async load(quiet) {
    try {
      const data = await api('/api/runs');
      const names = new Set(data.runs.map(r => r.name));
      if (this.known) names.forEach(n => { if (!this.known.has(n)) this.fresh.add(n); });
      else this.known = names;
      const changed = JSON.stringify(data.runs.map(r => [r.name, r.mtime])) !== JSON.stringify(this.runs.map(r => [r.name, r.mtime]));
      this.runs = data.runs;
      const up = document.getElementById('runsUpdated');
      if (up) up.textContent = 'read at ' + new Date().toLocaleTimeString('en-GB', { hour: '2-digit', minute: '2-digit', second: '2-digit' });
      if (changed || !quiet) this.render();
      return this.runs;
    } catch (e) {
      if (!quiet) this.body.innerHTML = `<tr><td colspan="4"><div class="error-box">${esc(e.message)}</div></td></tr>`;
      return this.runs;
    }
  }

  get(name) { return this.runs.find(r => r.name === name); }

  /** Corridas candidatas: todas, o solo las que pasan opts.filter (p.ej. con /modelo). */
  base() { return this.opts.filter ? this.runs.filter(this.opts.filter) : this.runs; }

  filtered() {
    return this.base().filter(r =>
      (this.kind === 'all' || r.kind === this.kind) &&
      (!this.query || `${r.name} ${r.label} ${r.status} ${r.tipo} ${r.date}`.toLowerCase().includes(this.query)));
  }

  render() {
    const counts = {};
    const base = this.base();
    base.forEach(r => { counts[r.kind] = (counts[r.kind] || 0) + 1; });
    counts.all = base.length;
    this.chips.innerHTML = KINDS.filter(([k]) => counts[k]).map(([k, label]) =>
      `<button class="chip ${k === this.kind ? 'active' : ''}" data-kind="${k}">${label} ${counts[k]}</button>`).join('');
    this.chips.querySelectorAll('.chip').forEach(c => c.onclick = () => { this.kind = c.dataset.kind; this.render(); });

    const rows = this.filtered();
    document.getElementById('runsCount').textContent = `${rows.length} of ${base.length}`;
    const multi = this.opts.multi;
    const cols = multi ? 4 : 3;
    if (!rows.length) {
      const msg = base.length ? 'No runs match.' : (this.opts.emptyMessage || 'No runs.');
      this.body.innerHTML = `<tr><td colspan="${cols}" class="results-empty">${msg}</td></tr>`;
      return;
    }
    this.body.innerHTML = rows.map(r => {
      const sel = this.opts.isSelected(r.name);
      const color = multi && sel ? this.opts.colorOf(r.name) : null;
      const tipo = r.tipo ? `<span class="tipo"> · ${r.tipo === 'optimization' ? 'opt' : r.tipo === 'simulation' ? 'sim' : esc(r.tipo)}</span>` : '';
      return `<tr data-name="${esc(r.name)}" class="${sel ? 'run-active' : ''} ${r.error ? 'run-error' : ''}" title="${esc(r.error || r.name)}">
        ${multi ? `<td class="check-cell"><input type="checkbox" ${sel ? 'checked' : ''} ${r.error ? 'disabled' : ''} tabindex="-1">${color ? ` <span class="swatch" style="background:${color}"></span>` : ''}</td>` : ''}
        <td><span class="run-date">${fmtDate(r.date)}</span>${this.fresh.has(r.name) ? '<span class="tag-new">new</span>' : ''}<span class="run-file">${esc(r.name)}</span></td>
        <td><span class="kind-tag">${esc(r.label)}${tipo}</span></td>
        <td>${statusBadge(r.status)}</td>
      </tr>`;
    }).join('');
  }

  scrollToActive() {
    this.body.querySelector('tr.run-active')?.scrollIntoView({ block: 'nearest' });
  }
}
