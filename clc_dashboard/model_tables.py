"""
Lectura del grupo /modelo (exportacion completa del modelo de Pyomo) para el visor de tablas.

Formato (version_export = 1), escrito por h5_export_full.py desde los notebooks:

  /modelo/<componente>        un subgrupo por Var/Expression/Param/Constraint/Objective
      attrs: nombre, tipo, bloque, etapa, doc, indices, n, unidades, escala,
             escala_conocida, escala_por_elemento
      datasets (largo n): idx_0.., valor, valor_fisico, y segun el tipo
             lb/ub/fijo (Var) o cuerpo/inferior/superior/residuo/activa/dual (Constraint)
  /modelo/_sets/<etapa>.<Set> los puntos de cada Set, con valor (adimensional) y
             valor_fisico + unidades_fisicas + regla para los continuos (A, R, T)

Aqui los valores NO se redondean: el visor promete el valor exacto de cada elemento (el
redondeo es solo de presentacion, en el navegador). El CSV exporta lo que llega en el JSON,
que para float64 es la representacion exacta mas corta.
"""

from pathlib import Path

import h5py
import numpy as np

import rule_latex

EXTRA_FIELDS = {
    'Var': ['lb', 'ub', 'fijo'],
    'Constraint': ['cuerpo', 'inferior', 'superior', 'residuo', 'dual', 'activa'],
}


def _s(v):
    return v.decode(errors='replace') if isinstance(v, bytes) else v


def _attr(g, k, default=None):
    v = g.attrs.get(k, default)
    if isinstance(v, bytes):
        return v.decode(errors='replace')
    if isinstance(v, np.ndarray):
        return [_s(x) for x in v.tolist()]
    if isinstance(v, np.bool_):
        return bool(v)
    if isinstance(v, np.integer):
        return int(v)
    if isinstance(v, np.floating):
        return None if not np.isfinite(v) else float(v)
    return v


def _jsonable(a):
    """Arreglo -> lista para JSON sin redondear; NaN/inf -> None, bytes -> str."""
    a = np.asarray(a)
    if a.dtype.kind in 'OSU':
        return [_s(x) for x in a.tolist()]
    if a.dtype.kind == 'b':
        return a.astype(bool).tolist()
    a = a.astype(float)
    out = a.astype(object)
    out[~np.isfinite(a)] = None
    return out.tolist()


def has_model(f):
    return 'modelo' in f and isinstance(f['modelo'], h5py.Group)


def _sets(M):
    sets = {}
    if '_sets' not in M:
        return sets
    for key, g in M['_sets'].items():
        valor = g['valor'][:]
        sets[key] = dict(
            key=key, nombre=_attr(g, 'nombre', key.split('.')[-1]), etapa=_attr(g, 'etapa', ''),
            regla=_attr(g, 'regla', ''), unidades=_attr(g, 'unidades_fisicas', ''),
            categorico=valor.dtype.kind in 'OSU',
            valor=_jsonable(valor),
            valor_fisico=_jsonable(g['valor_fisico'][:]) if 'valor_fisico' in g else None,
        )
    return sets


def _find_set(sets, etapa, name):
    """Set con ese nombre para la etapa del componente; si el componente es global se usa el
    Set sin prefijo, o el de cualquier etapa siempre que todas coincidan (A y R son iguales en
    todas las etapas; T no, porque t = T·tf depende de la etapa)."""
    if etapa and f'{etapa}.{name}' in sets:
        return sets[f'{etapa}.{name}']
    if name in sets:
        return sets[name]
    cands = [s for k, s in sets.items() if k.split('.')[-1] == name]
    if cands and all(c['valor_fisico'] == cands[0]['valor_fisico'] for c in cands):
        return cands[0]
    return None


def _formula(g, key, idx):
    """Formula de Expression/Constraint/Objective: el codigo de la regla y su LaTeX.

    El LaTeX sale de la REGLA (rule_latex: con los nombres del modelo). El de Pyomo
    (attr 'latex') sustituye los parametros por numeros, asi que solo se usa cuando el
    componente no tiene regla (p.ej. Expressions construidas con expr=, como D_ft = 3.28 Di).
    La expresion expandida por elemento ('expresion') no se usa."""
    regla = _attr(g, 'regla_fuente', '') if _attr(g, 'regla_fuente_ok', False) else ''
    latex, fuente = '', ''
    if regla:
        latex = rule_latex.to_latex(regla, _attr(g, 'tipo', ''), key.split('.')[-1], [str(s) for s in idx],
                                    sense=_attr(g, 'sentido')) or ''
        fuente = 'rule' if latex else ''
    elif _attr(g, 'latex_ok', False):
        latex, fuente = _attr(g, 'latex', ''), 'pyomo'
    return dict(regla=regla, latex=latex, latex_fuente=fuente)


VIOL_TOL = 1e-6  # residuo por encima de esto = elemento violado (el mismo umbral del navegador)


def _violation(g):
    """Restricciones: cuantos elementos violan (residuo > VIOL_TOL) y el peor residuo, para
    poder filtrar la lista por violadas / satisfechas sin abrir cada una."""
    if _attr(g, 'tipo', '') != 'Constraint' or 'residuo' not in g:
        return {}
    r = g['residuo'][:].astype(float)
    finite = r[np.isfinite(r)]
    return dict(n_viol=int((finite > VIOL_TOL).sum()),
                max_res=float(finite.max()) if finite.size else None)


def list_components(path):
    """Metadatos de todos los componentes (sin datos) + Sets."""
    with h5py.File(path, 'r') as f:
        if not has_model(f):
            return dict(has_model=False, components=[], sets={})
        M = f['modelo']
        sets = _sets(M)
        comps = []
        for name, g in M.items():
            if name == '_sets' or not isinstance(g, h5py.Group):
                continue
            idx = _attr(g, 'indices', []) or []
            if not isinstance(idx, list):
                idx = [idx]
            etapa = _attr(g, 'etapa', '')
            dims = []
            for s in idx:
                st = _find_set(sets, etapa, str(s))
                dims.append(dict(set=str(s), unidades=st['unidades'] if st else '', regla=st['regla'] if st else '',
                                 fisico=bool(st and st['valor_fisico'] is not None and not st['categorico'])))
            comps.append(dict(
                name=_attr(g, 'nombre', name), key=name, tipo=_attr(g, 'tipo', ''),
                bloque=_attr(g, 'bloque', ''), etapa=etapa, doc=_attr(g, 'doc', ''),
                dims=dims, n=_attr(g, 'n', 0), unidades=_attr(g, 'unidades', ''),
                escala=_attr(g, 'escala'), escala_conocida=_attr(g, 'escala_conocida', False),
                escala_por_elemento=_attr(g, 'escala_por_elemento', False),
                tiene_fisico='valor_fisico' in g,
                sentido=_attr(g, 'sentido'), es_derivada=_attr(g, 'es_derivada'),
                **_formula(g, name, idx),
                **_violation(g),
            ))
        meta = {k: _attr(M, k) for k in M.attrs}
    return dict(has_model=True, components=comps, sets=sets, meta=meta)


def read_component(path, key):
    with h5py.File(path, 'r') as f:
        if not has_model(f) or key == '_sets' or key not in f['modelo']:
            raise KeyError(f"Component '{key}' not found in /modelo.")
        M = f['modelo']
        g = M[key]
        sets = _sets(M)
        tipo = _attr(g, 'tipo', '')
        etapa = _attr(g, 'etapa', '')
        idx_names = _attr(g, 'indices', []) or []
        if not isinstance(idx_names, list):
            idx_names = [idx_names]
        idx, idx_fis = [], []
        for i, s in enumerate(idx_names):
            raw = g[f'idx_{i}'][:]
            idx.append(_jsonable(raw))
            st = _find_set(sets, etapa, str(s))
            if st and not st['categorico'] and st['valor_fisico'] is not None and raw.dtype.kind not in 'OSU':
                sv = np.asarray(st['valor'], dtype=float)
                sf = np.asarray(st['valor_fisico'], dtype=float)
                # cada indice -> punto del Set mas cercano (los indices SON puntos del Set)
                pos = np.abs(raw[:, None].astype(float) - sv[None, :]).argmin(axis=1)
                ok = np.isclose(raw.astype(float), sv[pos], rtol=1e-9, atol=1e-12)
                phys = np.where(ok, sf[pos], np.nan)
                idx_fis.append(_jsonable(phys))
            else:
                idx_fis.append(None)
        fields = {'valor': _jsonable(g['valor'][:])}
        if 'valor_fisico' in g:
            fields['valor_fisico'] = _jsonable(g['valor_fisico'][:])
        for k in EXTRA_FIELDS.get(tipo, []):
            if k in g:
                fields[k] = _jsonable(g[k][:])
        return dict(key=key, tipo=tipo, idx=idx, idx_fis=idx_fis, fields=fields)


def model_flag(path):
    try:
        with h5py.File(Path(path), 'r') as f:
            return has_model(f)
    except Exception:
        return False
