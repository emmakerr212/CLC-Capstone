"""
Desglose de un componente ESCALAR con regla (tipicamente la funcion objetivo) en sus terminos.

La regla del objetivo es, por ejemplo:
    CAPEX  = m.capital_annualized
    OPEX   = m.stream_cost_year + m.OC_cost + m.Compressor_cost + m.Fan_cost + m.Heating_cost
    INCOME = m.total_outlet_mass_CO2_year*m.co2_incentive + m.Turbine_income
    PROFIT = INCOME - (CAPEX + OPEX)
    return PROFIT

y el desglose es un arbol:
    PROFIT
      + INCOME           -> + total_outlet_mass_CO2_year·co2_incentive, + Turbine_income
      - CAPEX            -> + capital_annualized
      - OPEX             -> + stream_cost_year, + OC_cost, + Compressor_cost, ...

  - Los terminos se separan por sumas y restas (el signo se acumula hacia abajo).
  - Las variables locales (CAPEX, OPEX...) se expanden con su definicion y se EVALUAN con los
    valores exportados de los componentes; la raiz se verifica contra el valor de Pyomo.
  - Un componente sin regla (construido con expr=, p.ej. Fan_cost) se desglosa por etapa SOLO
    si existen componentes ox./purge./red. con el mismo nombre que suman exactamente su valor
    (verificado numericamente, nunca supuesto).

Los valores son los de Pyomo (attr 'valor'): para el objetivo y sus terminos de costo no hay
escalamiento, asi que estan en sus unidades fisicas.
"""

import ast
import textwrap

import h5py
import numpy as np

import rule_latex

STAGES = ('ox', 'purge', 'red')
RTOL = 1e-6


class _Unsup(Exception):
    pass


def _s(v):
    return v.decode(errors='replace') if isinstance(v, bytes) else v


def _scalars(M):
    """{clave: (valor, tipo, unidades, doc)} de todos los componentes escalares."""
    out = {}
    for k, g in M.items():
        if k == '_sets' or int(g.attrs.get('n', 0)) != 1:
            continue
        v = float(g['valor'][0])
        out[k] = (v if np.isfinite(v) else None, _s(g.attrs.get('tipo', '')),
                  _s(g.attrs.get('unidades', '')), _s(g.attrs.get('doc', '')))
    return out


class Breakdown:
    def __init__(self, M, key, M_values=None):
        """M: /modelo de donde sale la REGLA. M_values: /modelo de donde salen los VALORES (por
        defecto el mismo). Separarlos permite desglosar una corrida B exportada sin formulas con
        la regla de A: el arbol es identico (mismos ids) y la raiz se verifica contra la FO de B."""
        self.M = M
        self.Mv = M_values if M_values is not None else M
        self.key = key
        self.block = key.rsplit('.', 1)[0] if '.' in key else ''
        self.vals = _scalars(self.Mv)
        g = M[key]
        if not g.attrs.get('regla_fuente_ok', False):
            raise ValueError('this component has no rule source to break down')
        src = textwrap.dedent(_s(g.attrs['regla_fuente']))
        fn = next(n for n in ast.parse(src).body if isinstance(n, ast.FunctionDef))
        self.blkvar = fn.args.args[0].arg if fn.args.args else None
        if len(fn.args.args) > 1 or fn.args.vararg:
            raise ValueError('only scalar (non-indexed) rules can be broken down')
        self.defs, self.ret = {}, None
        for st in fn.body:
            if isinstance(st, ast.Assign) and len(st.targets) == 1 and isinstance(st.targets[0], ast.Name):
                self.defs[st.targets[0].id] = st.value
            elif isinstance(st, ast.Return):
                self.ret = st.value
                break
            elif isinstance(st, ast.Expr) and isinstance(st.value, ast.Constant):
                continue
            else:
                raise ValueError(f'unsupported statement in rule: {type(st).__name__}')
        if self.ret is None:
            raise ValueError('rule has no return')
        self.tr = rule_latex.Translator([], STAGES)

    # --- componentes -------------------------------------------------------
    def comp_key(self, node):
        """blk.X / m.X / m.ox.X -> clave del componente, o None."""
        if not isinstance(node, ast.Attribute):
            return None
        try:
            names = rule_latex._chain(node)
        except rule_latex.Unsupported:
            return None
        rest = '.'.join(names[1:])
        cands = ([f'{self.block}.{rest}', rest] if names[0] == self.blkvar else [rest, f'{self.block}.{rest}']) \
            if self.block else [rest]
        return next((c for c in cands if c in self.vals), None)

    # --- evaluacion --------------------------------------------------------
    def ev(self, n):
        if isinstance(n, ast.Constant) and isinstance(n.value, (int, float)):
            return float(n.value)
        if isinstance(n, ast.Name) and n.id in self.defs:
            return self.ev(self.defs[n.id])
        if isinstance(n, ast.Attribute):
            k = self.comp_key(n)
            if k is None or self.vals[k][0] is None:
                raise _Unsup('component')
            return self.vals[k][0]
        if isinstance(n, ast.UnaryOp) and isinstance(n.op, (ast.USub, ast.UAdd)):
            v = self.ev(n.operand)
            return -v if isinstance(n.op, ast.USub) else v
        if isinstance(n, ast.BinOp):
            a, b = self.ev(n.left), self.ev(n.right)
            op = n.op
            if isinstance(op, ast.Add): return a + b
            if isinstance(op, ast.Sub): return a - b
            if isinstance(op, ast.Mult): return a * b
            if isinstance(op, ast.Div): return a / b
            if isinstance(op, ast.Pow): return a ** b
        raise _Unsup(type(n).__name__)

    def try_ev(self, n):
        try:
            return self.ev(n)
        except (_Unsup, ZeroDivisionError, OverflowError):
            return None

    def latex(self, n):
        try:
            return self.tr.ex(n)[0]
        except Exception:
            return ''

    # --- arbol ---------------------------------------------------------------
    def terms(self, n, sign):
        """Separa n en terminos de suma con su signo: [(nodo_ast, signo)]."""
        if isinstance(n, ast.BinOp) and isinstance(n.op, (ast.Add, ast.Sub)):
            right = sign if isinstance(n.op, ast.Add) else -sign
            return self.terms(n.left, sign) + self.terms(n.right, right)
        if isinstance(n, ast.UnaryOp) and isinstance(n.op, ast.USub):
            return self.terms(n.operand, -sign)
        return [(n, sign)]

    def node(self, n, sign, path):
        """Nodo del arbol. `path` es el id del nodo: posicion + nombre desde la raiz, estable entre
        corridas con la misma regla (asi se emparejan A y B al comparar)."""
        if isinstance(n, ast.Name) and n.id in self.defs:      # variable local: se expande
            d = self.defs[n.id]
            return dict(id=path, label=n.id, latex=rule_latex.sym(n.id), kind='local',
                        sign=sign, value=self.try_ev(n), formula=self.latex(d),
                        children=[self.node(t, s, f'{path}/{i}') for i, (t, s) in enumerate(self.terms(d, 1))])
        k = self.comp_key(n)
        if k is not None:                                        # componente
            v, tipo, unit, doc = self.vals[k]
            return dict(id=path, label=k, latex=self.latex(n), kind='component', key=k, tipo=tipo,
                        unidades=unit, doc=doc, sign=sign, value=v, children=self.stage_split(k, path))
        # termino compuesto (producto, cociente...): se evalua; sus factores se muestran como detalle
        factors = []
        if isinstance(n, ast.BinOp) and isinstance(n.op, (ast.Mult, ast.Div)):
            for i, f in enumerate(self.factors(n)):
                fk = self.comp_key(f)
                factors.append(dict(id=f'{path}/{i}', label=fk or self.label(f), latex=self.latex(f),
                                    kind='factor', key=fk, tipo=self.vals[fk][1] if fk else '',
                                    unidades=self.vals[fk][2] if fk else '', doc=self.vals[fk][3] if fk else '',
                                    sign=1, value=self.try_ev(f), children=[]))
        return dict(id=path, label=self.label(n), latex=self.latex(n), kind='term',
                    sign=sign, value=self.try_ev(n), children=factors)

    def label(self, n):
        """Texto legible de un termino: sin 'm.' ni la variable del bloque, × en vez de *."""
        s = ast.unparse(n)
        for pre in {'m', self.blkvar} - {None}:
            s = s.replace(f'{pre}.', '')
        return s.replace(' * ', ' × ').replace('**', '^')

    def factors(self, n):
        if isinstance(n, ast.BinOp) and isinstance(n.op, ast.Mult):
            return self.factors(n.left) + self.factors(n.right)
        return [n]

    def stage_split(self, k, path):
        """Componente global sin regla -> sus partes por etapa, solo si suman exactamente."""
        if '.' in k:
            return []
        v = self.vals[k][0]
        low = k.lower()
        parts = [f'{st}.{c.split(".", 1)[1]}' for st in STAGES for c in self.vals
                 if c.startswith(st + '.') and c.split('.', 1)[1].lower() == low]
        parts = list(dict.fromkeys(parts))
        if len(parts) < 2 or v is None or any(self.vals[p][0] is None for p in parts):
            return []
        total = sum(self.vals[p][0] for p in parts)
        if abs(total - v) > RTOL * max(1.0, abs(v)):
            return []
        return [dict(id=f'{path}/{p.split(".")[0]}', label=p, latex='', kind='component', key=p, tipo=self.vals[p][1],
                     unidades=self.vals[p][2], doc=self.vals[p][3], sign=1, value=self.vals[p][0],
                     children=[], stage=p.split('.')[0]) for p in parts]

    def build(self):
        g = self.M[self.key]
        if self.key not in self.Mv:
            raise KeyError(f"'{self.key}' does not exist in the run that provides the values")
        pyomo_val = float(self.Mv[self.key]['valor'][0])
        kids = self.terms(self.ret, 1)
        if len(kids) == 1 and isinstance(kids[0][0], ast.Name) and kids[0][0].id in self.defs:
            root = self.node(kids[0][0], 1, 'r')                 # return PROFIT -> la raiz es PROFIT
        else:
            root = dict(id='r', label=self.key, latex='', kind='local', sign=1, value=self.try_ev(self.ret),
                        formula=self.latex(self.ret),
                        children=[self.node(t, s, f'r/{i}') for i, (t, s) in enumerate(kids)])
        check = root['value']
        ok = check is not None and abs(check - pyomo_val) <= RTOL * max(1.0, abs(pyomo_val))
        return dict(key=self.key, tipo=_s(g.attrs.get('tipo', '')), sense=_s(g.attrs.get('sentido', '')),
                    unidades=_s(g.attrs.get('unidades', '')), value=pyomo_val, evaluated=check, verified=ok,
                    root=root)


def breakdown(path, key, values_path=None):
    """Desglose de `key` con la regla de `path` y los valores de `values_path` (o de `path`)."""
    with h5py.File(path, 'r') as f:
        M = f['modelo']
        if key not in M:
            raise KeyError(f"Component '{key}' not found in /modelo.")
        if values_path is None:
            return Breakdown(M, key).build()
        with h5py.File(values_path, 'r') as fv:
            if 'modelo' not in fv:
                raise KeyError('the run that provides the values has no /modelo')
            return Breakdown(M, key, fv['modelo']).build()
