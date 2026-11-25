"""
Lectura y normalizacion de los H5 de resultados de CLC para el dashboard.

En Runs/ conviven tres familias de archivos y las tres se llevan a la MISMA forma, asi el
front-end no necesita saber de formatos:

  formato 2   attrs 'formato'=2 y 'etapas' en la raiz, un grupo por etapa
              (<etapa>/parametros/ejes, <etapa>/variables) + grupo 'decisiones'.
              OX_PURGE, RED, COMPLETE_CYCLE (desde 2026-09-23).
  formato 1   parametros/ejes + variables con especies (C_O2, Cci_O2, ...), una sola etapa
              (CLC_OX_results de septiembre). La etapa se deduce del nombre del archivo.
  legacy      junio-julio: un solo gas sin nombre de especie ('C', 'Cc'), ejes en 'axes/' o en
              'parametros/ejes', y a veces grupos 'funcion_objetivo' / 'restricciones'. Ese gas se
              muestra como la especie 'gas'.

Forma normalizada (load_run):
  stages        [{key, name, color, i0, i1, t0, t1}]  i0..i1 = indices globales de la etapa
  t, t_local    tiempo global (etapas encadenadas) y tiempo local de cada etapa
  stage_idx     indice de etapa de cada instante
  z, r
  vars          {nombre: arreglo} con el TIEMPO COMO PRIMER EJE: (nt,), (nt, nz) o (nt, nr, nz).
                Una variable que no existe en alguna etapa queda en NaN ahi (p.ej. C_CH4 en ox).
  gas_by_stage, solids, x_keys, feed_species
  summary       tarjetas (kpis), tabla de decisiones y todos los attrs agrupados

Igual que en CYCLE_dashboard.py, el t=0 de cada etapa > 0 se salta: es el mismo estado que el t
final de la anterior y repetirlo duplicaria un instante en el eje global.
"""

import os
import re
from datetime import datetime
from pathlib import Path

import h5py
import numpy as np

STAGE_NAMES = {'ox': 'Oxidation', 'purge': 'Purge', 'red': 'Reduction'}
TIPO_EN = {'optimizacion': 'optimization', 'simulacion': 'simulation'}
GROUP_TITLES = {'parametros': 'Parameters', 'funcion_objetivo': 'Objective function', 'restricciones': 'Constraints'}
STAGE_COLORS = {'ox': '#2a78d6', 'purge': '#1baf7a', 'red': '#eb6834'}
PREFERRED_GAS = ['O2', 'CH4', 'CO', 'H2', 'N2', 'CO2', 'H2O', 'gas']
PREFERRED_SOLID = ['Fe3O4', 'Fe2O3']

# legacy: un solo gas sin especie
LEGACY_RENAME = {'C': 'C_gas', 'Cc': 'Cci_gas'}


# ---------------------------------------------------------------------------
# Utilidades
# ---------------------------------------------------------------------------

def _decode(v):
    if isinstance(v, bytes):
        return v.decode(errors='replace')
    if isinstance(v, np.ndarray):
        if v.dtype.kind in 'OSU':
            return [_decode(x) for x in v.tolist()]
        return v.tolist() if v.size <= 64 else f'<arreglo {v.shape}>'
    if isinstance(v, (np.bool_,)):
        return bool(v)
    if isinstance(v, np.integer):
        return int(v)
    if isinstance(v, np.floating):
        return None if not np.isfinite(v) else float(v)
    if isinstance(v, (list, tuple)):
        return [_decode(x) for x in v]
    return v


def _attrs(obj):
    return {k: _decode(v) for k, v in obj.attrs.items()}


def _ordered(found, preferred):
    return [s for s in preferred if s in found] + sorted(s for s in found if s not in preferred)


def to_json_array(a, sig=6):
    """Arreglo -> listas anidadas para JSON: 6 cifras significativas y NaN/inf -> null."""
    a = np.asarray(a, dtype=float)
    finite = np.isfinite(a)
    out = np.zeros_like(a)
    with np.errstate(divide='ignore', invalid='ignore'):
        mag = np.floor(np.log10(np.abs(np.where(finite & (a != 0), a, 1.0))))
        scale = 10.0 ** (sig - 1 - mag)
        out[finite] = np.round(a[finite] * scale[finite]) / scale[finite]
    obj = out.astype(object)
    obj[~finite] = None
    return obj.tolist()


def _stage_from_name(name):
    low = name.lower()
    if 'purge' in low or 'purga' in low:
        return 'purge'
    if re.search(r'(^|[_\s-])red([_\s.-]|$)', low) or 'reduc' in low:
        return 'red'
    return 'ox'


def _date_from_name(name, path):
    m = re.search(r'(20\d{6})_(\d{6})', name)
    if m:
        try:
            return datetime.strptime(m.group(1) + m.group(2), '%Y%m%d%H%M%S')
        except ValueError:
            pass
    return datetime.fromtimestamp(os.path.getmtime(path))


# ---------------------------------------------------------------------------
# Deteccion de formato + lectura cruda
# ---------------------------------------------------------------------------

def _detect(f):
    if int(f.attrs.get('formato', 0) or 0) == 2 and 'etapas' in f.attrs:
        return 2
    if 'variables' in f and 'parametros' in f and any(k.startswith('C_') for k in f['variables']):
        return 1
    return 0  # legacy


_AXES_LAYOUTS = [  # (grupo, nombre de t, de z, de r) — en orden de preferencia
    ('parametros/ejes', 't', 'z', 'r'),
    ('axes', 't', 'z', 'r'),
    ('parametros/discretizacion', 't_pts', 'y_pts', 'r_pts'),
    ('/', 't_pts', 'y_pts', 'r_pts'),
]
_AXIS_NAMES = {'t', 'z', 'r', 't_pts', 'y_pts', 'r_pts'}


def _axes(f):
    """-> (t, z, r, grupo de ejes o None). El grupo se devuelve para recoger perfiles guardados
    junto a los ejes (p.ej. m_N2 en la purga legacy)."""
    for grp, kt, kz, kr in _AXES_LAYOUTS:
        if grp in f and all(k in f[grp] for k in (kt, kz, kr)):
            g = f[grp]
            return g[kt][:], g[kz][:], g[kr][:], g
    raise ValueError('The file has no t/z/r axes (empty or interrupted export?).')


def _read_raw(f, fmt, name):
    """-> (stages_raw, groups) ; stages_raw = [(key, {t, z, r, data})]."""
    groups = []  # [(titulo, dict de attrs)] para el panel de resumen
    root = _attrs(f)
    if root:
        groups.append(('General', root))

    if fmt == 2:
        stages = [str(s) for s in _decode(f.attrs['etapas'])]
        if 'decisiones' in f:
            groups.append(('Decisions', _attrs(f['decisiones'])))
        raw = []
        for st in stages:
            g = f[st]
            par = g['parametros']
            groups.append((f'Parameters — {STAGE_NAMES.get(st, st)}', _attrs(par)))
            ax = par['ejes']
            raw.append((st, dict(t=ax['t'][:], z=ax['z'][:], r=ax['r'][:],
                                 data={k: g['variables'][k][:] for k in g['variables']})))
        return raw, groups

    st = _stage_from_name(name)
    t, z, r, ax = _axes(f)
    if fmt == 1:
        data = {k: f['variables'][k][:] for k in f['variables']}
    else:
        src = f['variables'] if 'variables' in f else f
        data = {LEGACY_RENAME.get(k, k): src[k][:] for k in src
                if isinstance(src[k], h5py.Dataset) and k not in _AXIS_NAMES}
        # perfiles de entrada guardados junto a los ejes (p.ej. m_N2 en la purga legacy)
        for k in ax:
            if k not in _AXIS_NAMES and isinstance(ax[k], h5py.Dataset) and ax[k].shape == (len(t),):
                data.setdefault(k, ax[k][:])
    for gname in ('parametros', 'funcion_objetivo', 'restricciones'):
        if gname in f and isinstance(f[gname], h5py.Group):
            a = _attrs(f[gname])
            if a:
                groups.append((GROUP_TITLES[gname], a))
    return [(st, dict(t=t, z=z, r=r, data=data))], groups


# ---------------------------------------------------------------------------
# Listado (solo attrs; rapido)
# ---------------------------------------------------------------------------

_INFO_CACHE = {}


def _kind_label(stages, name, tipo):
    s = tuple(stages)
    label = {('ox', 'purge', 'red'): 'Full cycle', ('ox', 'purge'): 'Ox + Purge'}.get(
        s, ' + '.join(STAGE_NAMES.get(x, x) for x in s))
    kind = {('ox', 'purge', 'red'): 'cycle', ('ox', 'purge'): 'oxpurge'}.get(s, s[0] if len(s) == 1 else 'other')
    return kind, label


def run_info(path):
    path = Path(path)
    st_ = path.stat()
    key = (str(path), st_.st_mtime, st_.st_size)
    if key in _INFO_CACHE:
        return _INFO_CACHE[key]
    name = path.name
    info = dict(name=name, size=st_.st_size, mtime=st_.st_mtime,
                date=_date_from_name(name, path).isoformat(timespec='seconds'))
    try:
        with h5py.File(path, 'r') as f:
            fmt = _detect(f)
            info['has_model'] = 'modelo' in f  # exportacion completa del modelo (visor de tablas)
            if fmt == 2:
                stages = [str(s) for s in _decode(f.attrs['etapas'])]
                status = _decode(f.attrs.get('solver_status', ''))
                tipo = _decode(f.attrs.get('tipo_corrida', ''))
            else:
                _axes(f)  # un archivo sin ejes (exportacion interrumpida) se marca ya en el listado
                stages = [_stage_from_name(name)]
                par = _attrs(f['parametros']) if 'parametros' in f else {}
                root = _attrs(f)
                status = next((d[k] for k in ('solver_status', 'termination', 'status')
                               for d in (par, root) if d.get(k)), '')
                tipo = ''
        if not tipo:
            low = name.lower()
            tipo = 'optimizacion' if 'opt' in low else ('simulacion' if 'sim' in low else '')
        kind, label = _kind_label(stages, name, tipo)
        tipo = TIPO_EN.get(str(tipo), str(tipo))
        info.update(format=fmt, stages=stages, status=str(status), tipo=tipo, kind=kind, label=label)
    except Exception as e:  # archivo vacio/corrupto: se lista igual, marcado con error
        info.update(format=None, stages=[], status='error', tipo='', kind='error', label='Unreadable',
                    error=f'{type(e).__name__}: {e}')
    _INFO_CACHE[key] = info
    return info


def list_runs(runs_dir):
    runs = [run_info(p) for p in Path(runs_dir).glob('*.h5')]
    runs.sort(key=lambda d: (d['date'], d['name']), reverse=True)
    return runs


# ---------------------------------------------------------------------------
# Carga completa + normalizacion
# ---------------------------------------------------------------------------

def _to_time_first(a):
    a = np.asarray(a, dtype=float)
    if a.ndim == 1:
        return a
    if a.ndim == 2:
        return a.T                      # (nz, nt) -> (nt, nz)
    if a.ndim == 3:
        return a.transpose(2, 1, 0)     # (nz, nr, nt) -> (nt, nr, nz)
    raise ValueError(f'Variables with {a.ndim} dimensions are not supported.')


def _num(v):
    try:
        x = float(v)
        return x if np.isfinite(x) else None
    except (TypeError, ValueError):
        return None


def _build_summary(info, groups, stages_out, z, r, nt, decisions):
    flat = {}
    for _, g in groups:
        for k, v in g.items():
            flat.setdefault(k, v)

    kpis = []

    def add(label, value, sub='', fmt='{:.3g}', tone=''):
        if value is None or value == '':
            return
        if isinstance(value, (int, float)):
            value = fmt.format(value)
        kpis.append(dict(label=label, value=str(value), sub=sub, tone=tone))

    status = info.get('status', '')
    tone = 'good' if status in ('optimal', 'ok', 'locallyOptimal') else ('danger' if status in ('infeasible', 'error') else ('warn' if status else ''))
    add('Solver', status or '—', info.get('tipo', ''), tone=tone)
    add('L', _num(flat.get('L')), 'm')
    add('Di', _num(flat.get('Di')), 'm')
    for s in stages_out:
        add(f'tf {s["name"].lower()}', s['t1'] - s['t0'], 's', fmt='{:.6g}')
    add('Grid', f'{len(z)}×{len(r)}×{nt}', 'z × r × t')
    add('DOF', _num(flat.get('dof', flat.get('degrees_of_freedom'))), '', fmt='{:.0f}')
    add('Solve time', _num(flat.get('solve_time_s')), 's', fmt='{:.1f}')
    obj = _num(flat.get('obj_usd_yr'))
    add('Objective', obj, 'USD/yr', fmt='{:,.0f}')
    return dict(kpis=kpis, groups=[dict(title=t, items=[[k, v] for k, v in g.items()]) for t, g in groups],
                decisions=decisions)


def _read_decisions(f):
    if 'decisiones' not in f:
        return []
    g = f['decisiones']
    if 'decision_nombres' not in g or 'decision_valor' not in g:
        return []
    names = _decode(g['decision_nombres'][:])
    vals = g['decision_valor'][:]
    free = g['decision_libre'][:] if 'decision_libre' in g else [None] * len(names)
    return [dict(name=str(n), value=_num(v), free=None if l is None else bool(l))
            for n, v, l in zip(names, vals, free)]


def load_run(path):
    path = Path(path)
    info = run_info(path)
    if info.get('error'):
        raise ValueError(info['error'])
    with h5py.File(path, 'r') as f:
        fmt = _detect(f)
        raw, groups = _read_raw(f, fmt, path.name)
        decisions = _read_decisions(f)

    z = np.asarray(raw[0][1]['z'], dtype=float)
    r = np.asarray(raw[0][1]['r'], dtype=float)
    for st, d in raw[1:]:
        if d['z'].shape != z.shape or d['r'].shape != r.shape or not (
                np.allclose(d['z'], z) and np.allclose(d['r'], r)):
            raise ValueError(f"The z/r grids of '{st}' differ from those of '{raw[0][0]}'.")

    # tiempo global
    t_glob, t_loc, stage_idx, stages_out = [], [], [], []
    offset = 0.0
    for n, (st, d) in enumerate(raw):
        t = np.asarray(d['t'], dtype=float)
        sl = slice(None) if n == 0 else slice(1, None)
        i0 = len(t_glob)
        t_glob.extend((t[sl] + offset).tolist())
        t_loc.extend(t[sl].tolist())
        stage_idx.extend([n] * len(t[sl]))
        stages_out.append(dict(key=st, name=STAGE_NAMES.get(st, st), color=STAGE_COLORS.get(st, '#898781'),
                               i0=i0, i1=len(t_glob) - 1, t0=offset if n == 0 else t_glob[i0 - 1],
                               t1=t_glob[-1]))
        offset += float(t[-1])
    nt = len(t_glob)

    # variables: union de todas las etapas, NaN donde una etapa no la tiene
    shapes = {}
    for st, d in raw:
        for k, v in d['data'].items():
            if np.asarray(v).dtype.kind not in 'fiub':
                continue
            shapes.setdefault(k, _to_time_first(v).shape[1:])
    vars_ = {}
    for k, shp in shapes.items():
        parts = []
        for n, (st, d) in enumerate(raw):
            sl = slice(None) if n == 0 else slice(1, None)
            n_loc = len(d['t'][sl])
            v = d['data'].get(k)
            if v is not None and _to_time_first(v).shape[1:] == shp and _to_time_first(v).shape[0] == len(d['t']):
                parts.append(_to_time_first(v)[sl])
            else:
                parts.append(np.full((n_loc,) + shp, np.nan))
        vars_[k] = np.concatenate(parts, axis=0)

    gas_by_stage = []
    for st, d in raw:
        found = {k[2:] for k in d['data'] if k.startswith('C_')}
        gas_by_stage.append(_ordered(found, PREFERRED_GAS))
    solids = _ordered({k[4:] for k in vars_ if k.startswith('Ccs_')}, PREFERRED_SOLID)
    x_keys = [f'X_{s}' for s in _ordered({k[2:] for k in vars_ if k.startswith('X_')}, PREFERRED_SOLID)]
    if not x_keys and 'X' in vars_:
        x_keys = ['X']

    # flujo masico de entrada [kg/s] por especie alimentada = m_total * w_in_<i>
    feed_species = []
    if fmt == 0 and 'm_total' not in vars_ and 'm_N2' in vars_:
        # purga legacy: solo se alimenta N2, sin fracciones w_in_*
        vars_['m_total'] = vars_['m_N2']
        vars_['Fin_N2'] = vars_['m_N2']
        feed_species = ['N2']
    elif 'm_total' in vars_:
        found = {k[5:] for k in vars_ if k.startswith('w_in_') and np.nanmax(np.abs(vars_[k]), initial=0) > 1e-9}
        feed_species = _ordered(found, PREFERRED_GAS)
        for sp in feed_species:
            vars_[f'Fin_{sp}'] = vars_['m_total'] * np.nan_to_num(vars_[f'w_in_{sp}'])

    summary = _build_summary(info, groups, stages_out, z, r, nt, decisions)
    return dict(
        name=path.name, info=info, stages=stages_out,
        t=to_json_array(t_glob), t_local=to_json_array(t_loc), stage_idx=stage_idx,
        z=to_json_array(z), r=to_json_array(r),
        gas_by_stage=gas_by_stage, solids=solids, x_keys=x_keys, feed_species=feed_species,
        vars={k: to_json_array(v) for k, v in vars_.items()},
        summary=summary,
    )


_RUN_CACHE = {}


def load_run_cached(path, max_items=16):
    path = Path(path)
    key = (str(path), path.stat().st_mtime)
    if key not in _RUN_CACHE:
        if len(_RUN_CACHE) >= max_items:
            _RUN_CACHE.pop(next(iter(_RUN_CACHE)))
        _RUN_CACHE[key] = load_run(path)
    return _RUN_CACHE[key]
