"""
Codigo de una regla de Pyomo  ->  LaTeX legible, con los NOMBRES del modelo.

El LaTeX de pyomo.contrib.latex_printer se construye desde el modelo ya armado: falla con
reglas con if/for y, cuando funciona, sustituye los parametros por numeros. Aqui se parte del
codigo fuente de la regla (attr 'regla_fuente' del H5) y se traduce con el arbol sintactico
de Python (ast):

    def _deltaP_(blk_ox, tao):
        return m.L*sum(blk_ox.hz*blk_ox.dPdz[z,tao] for z in blk_ox.A)
  ->
    \\mathrm{dP}(\\tau) = L \\, \\sum_{z \\in A} \\mathrm{hz} \\, \\mathrm{dPdz}(z,\\tau)

Convenciones:
  - blk.X / m.X / m.ox.X  ->  simbolo de X (el bloque no se escribe)
  - X[a,b]                ->  X(a,b)   (parentesis: evita dobles subindices como rho_g_{z})
  - X_hat                 ->  \\hat{X};  rho_g -> \\rho_{g};  tao -> \\tau
  - dX_hat_dA, dX_hat_dA2 ->  derivada parcial (primera / segunda) de \\hat{X} respecto de A
  - S.first() / S.last()  ->  \\min S / \\max S   (los Sets continuos estan ordenados)
  - if ...: return Constraint.Skip  ->  "except ..." en la linea del dominio
  - variables locales (Accum = ...) -> la ecuacion final primero y luego "with" + definiciones

Si una regla usa algo que no se sabe traducir, to_latex devuelve None y el dashboard no
muestra LaTeX para ese componente (nunca un LaTeX inventado).
"""

import ast
import re
import textwrap

GREEK = {
    'alpha': r'\alpha', 'beta': r'\beta', 'gamma': r'\gamma', 'delta': r'\delta',
    'epsilon': r'\varepsilon', 'eps': r'\varepsilon', 'zeta': r'\zeta', 'eta': r'\eta',
    'theta': r'\theta', 'kappa': r'\kappa', 'lambda': r'\lambda', 'lam': r'\lambda', 'mu': r'\mu',
    'nu': r'\nu', 'xi': r'\xi', 'pi': r'\pi', 'rho': r'\rho', 'sigma': r'\sigma', 'tau': r'\tau',
    'tao': r'\tau', 'phi': r'\phi', 'chi': r'\chi', 'psi': r'\psi', 'omega': r'\omega',
    'Gamma': r'\Gamma', 'Delta': r'\Delta', 'Theta': r'\Theta', 'Lambda': r'\Lambda',
    'Sigma': r'\Sigma', 'Phi': r'\Phi', 'Psi': r'\Psi', 'Omega': r'\Omega',
}
MODULES = {'pyo', 'np', 'math', 'pe', 'pyomo', 'numpy'}
FUNCS = {'exp': r'\exp', 'log': r'\ln', 'log10': r'\log_{10}', 'sin': r'\sin', 'cos': r'\cos',
         'tan': r'\tan', 'tanh': r'\tanh', 'sinh': r'\sinh', 'cosh': r'\cosh', 'atan': r'\arctan'}
CMP = {ast.Eq: '=', ast.LtE: r'\le', ast.GtE: r'\ge', ast.Lt: '<', ast.Gt: '>', ast.NotEq: r'\ne',
       ast.In: r'\in', ast.NotIn: r'\notin'}

ATOM, POW, UNARY, MUL, ADD, CMPP = 5, 4, 3, 2, 1, 0


class Unsupported(Exception):
    pass


def _base(s):
    if s in GREEK:
        return GREEK[s]
    return s if len(s) == 1 else rf'\mathrm{{{s}}}'


def sym(name):
    """Nombre de Python -> simbolo LaTeX."""
    m = re.match(r'^d([A-Z]\w*?)_d([A-Z])(2?)$', name)
    if m:
        top, var, two = sym(m.group(1)), _base(m.group(2)), m.group(3)
        return (rf'\frac{{\partial^{{2}} {top}}}{{\partial {var}^{{2}}}}' if two
                else rf'\frac{{\partial {top}}}{{\partial {var}}}')
    parts = name.split('_')
    hat = len(parts) > 1 and parts[-1] == 'hat'
    if hat:
        parts = parts[:-1]
    base, subs = parts[0], [p for p in parts[1:] if p]
    if not base:
        return rf'\mathrm{{{name}}}'
    # nombres "de palabras" (capital_annualized, Fan_cost, co2_incentive): una sola palabra,
    # no base + subindice; rho_g, MW_mix, C_feed siguen como simbolo con subindice
    if base not in GREEK and len(base) > 2 and any(len(p) > 3 for p in subs):
        word = '\\_'.join([base] + subs)
        return rf'\hat{{\mathrm{{{word}}}}}' if hat else rf'\mathrm{{{word}}}'
    out = _base(base)
    if hat:
        out = rf'\hat{{{out}}}'
    if subs:
        out += '_{' + ','.join(GREEK.get(p, p if len(p) == 1 else rf'\mathrm{{{p}}}') for p in subs) + '}'
    return out


def _num(v):
    if isinstance(v, bool) or v is None:
        return rf'\text{{{v}}}'
    if isinstance(v, int):
        return str(v)
    if isinstance(v, float):
        if v != 0 and (abs(v) < 1e-3 or abs(v) >= 1e5):
            mant, exp = f'{v:.6g}'.split('e') if 'e' in f'{v:.6g}' else (f'{v:.6g}', '0')
            return rf'{mant} \times 10^{{{int(exp)}}}'
        return f'{v:.10g}'
    if isinstance(v, str):
        return rf'\text{{{v}}}'
    raise Unsupported(f'constant {v!r}')


def _chain(node):
    """Attribute anidado -> lista de nombres ['blk_ox', 'R', 'last']."""
    out = []
    while isinstance(node, ast.Attribute):
        out.append(node.attr)
        node = node.value
    if not isinstance(node, ast.Name):
        raise Unsupported('attribute on a non-name')
    out.append(node.id)
    return out[::-1]


class Translator:
    def __init__(self, params, blocks):
        self.params = set(params)       # indices de la regla (z, tao, i, ...)
        self.blocks = set(blocks)       # nombres de bloque que se omiten (ox, purge, red)

    def wrap(self, node, minp):
        s, p = self.ex(node)
        return rf'\left({s}\right)' if p < minp else s

    def ex(self, n):
        if isinstance(n, ast.Constant):
            return _num(n.value), ATOM
        if isinstance(n, ast.Name):
            return (GREEK.get(n.id, n.id) if n.id in self.params else sym(n.id)), ATOM
        if isinstance(n, ast.Attribute):
            names = _chain(n)
            rest = [x for x in names[1:] if x not in self.blocks]
            if names[0] in MODULES:
                if names[-1] == 'pi':
                    return r'\pi', ATOM
                raise Unsupported(f'module attribute {".".join(names)}')
            if not rest:
                raise Unsupported('bare block')
            return sym(rest[-1]), ATOM
        if isinstance(n, ast.Subscript):
            base, _ = self.ex(n.value)
            sl = n.slice
            elts = sl.elts if isinstance(sl, ast.Tuple) else [sl]
            return f'{base}({",".join(self.ex(e)[0] for e in elts)})', ATOM
        if isinstance(n, ast.BinOp):
            return self.binop(n)
        if isinstance(n, ast.UnaryOp):
            if isinstance(n.op, ast.USub):
                return '-' + self.wrap(n.operand, MUL), UNARY
            if isinstance(n.op, ast.UAdd):
                return self.ex(n.operand)
            if isinstance(n.op, ast.Not):
                return r'\lnot ' + self.wrap(n.operand, UNARY), UNARY
            raise Unsupported('unary op')
        if isinstance(n, ast.Compare):
            parts = [self.wrap(n.left, ADD)]
            for op, right in zip(n.ops, n.comparators):
                if type(op) not in CMP:
                    raise Unsupported('comparison')
                parts += [CMP[type(op)], self.wrap(right, ADD)]
            return ' '.join(parts), CMPP
        if isinstance(n, ast.BoolOp):
            word = r'\text{ and }' if isinstance(n.op, ast.And) else r'\text{ or }'
            return word.join(self.wrap(v, CMPP) for v in n.values), CMPP
        if isinstance(n, (ast.Tuple, ast.List, ast.Set)):        # z not in (z0, z1)
            return r'\{' + ', '.join(self.ex(e)[0] for e in n.elts) + r'\}', ATOM
        if isinstance(n, ast.IfExp):
            return (rf'\begin{{cases}} {self.ex(n.body)[0]} & \text{{if }} {self.ex(n.test)[0]} \\ '
                    rf'{self.ex(n.orelse)[0]} & \text{{otherwise}} \end{{cases}}'), ATOM
        if isinstance(n, ast.Call):
            return self.call(n)
        raise Unsupported(type(n).__name__)

    def binop(self, n):
        op = n.op
        if isinstance(op, ast.Add):
            return f'{self.wrap(n.left, ADD)} + {self.wrap(n.right, ADD)}', ADD
        if isinstance(op, ast.Sub):
            return f'{self.wrap(n.left, ADD)} - {self.wrap(n.right, MUL)}', ADD
        if isinstance(op, ast.Mult):
            # una suma a la DERECHA de un producto se lee sin parentesis (L \sum_z ...);
            # a la izquierda si los necesita, o el factor quedaria "dentro" de la suma
            right = self.ex(n.right)[0] if _is_sum(n.right) else self.wrap(n.right, MUL)
            return rf'{self.wrap(n.left, MUL)} \, {right}', MUL
        if isinstance(op, ast.Div):
            return rf'\frac{{{self.ex(n.left)[0]}}}{{{self.ex(n.right)[0]}}}', ATOM
        if isinstance(op, ast.Pow):
            return f'{{{self.wrap(n.left, ATOM)}}}^{{{self.ex(n.right)[0]}}}', POW
        raise Unsupported(type(op).__name__)

    def call(self, n):
        f = n.func
        args = n.args
        if isinstance(f, ast.Name) and f.id == 'sum' and len(args) == 1 and isinstance(args[0], ast.GeneratorExp):
            g = args[0]
            under = []
            for comp in g.generators:
                tgt = comp.target
                names = [e.id for e in tgt.elts] if isinstance(tgt, ast.Tuple) else [tgt.id]
                self.params.update(names)
                under.append(rf'{",".join(GREEK.get(x, x) for x in names)} \in {self.ex(comp.iter)[0]}')
                under += [self.ex(c)[0] for c in comp.ifs]
            body = self.wrap(g.elt, MUL)
            sub = under[0] if len(under) == 1 else r'\substack{' + r' \\ '.join(under) + '}'
            return rf'\sum_{{{sub}}} {body}', ADD
        if isinstance(f, ast.Attribute):
            names = _chain(f)
            meth = names[-1]
            # blk.R.last() / blk.A.first() / blk.A.at(k)
            if meth in ('first', 'last') and not args:
                return rf'\{"min" if meth == "first" else "max"} {self.ex(f.value)[0]}', ATOM
            if meth == 'at' and len(args) == 1:
                return rf'{self.ex(f.value)[0]}_{{({self.ex(args[0])[0]})}}', ATOM
            if meth in ('prev', 'next') and len(args) == 1:      # T.prev(tao): punto anterior del Set
                return rf'\operatorname{{{meth}}}_{{{self.ex(f.value)[0]}}}\left({self.ex(args[0])[0]}\right)', ATOM
            if names[0] in MODULES:
                return self.func(meth, args)
            raise Unsupported(f'method {meth}')
        if isinstance(f, ast.Name):
            return self.func(f.id, args)
        raise Unsupported('call')

    def func(self, name, args):
        a = [self.ex(x)[0] for x in args]
        if name in FUNCS and len(a) == 1:
            return rf'{FUNCS[name]}\left({a[0]}\right)', ATOM
        if name == 'sqrt' and len(a) == 1:
            return rf'\sqrt{{{a[0]}}}', ATOM
        if name in ('abs', 'fabs') and len(a) == 1:
            return rf'\left|{a[0]}\right|', ATOM
        if name in ('min', 'max'):
            return rf'\{name}\left({", ".join(a)}\right)', ATOM
        if name == 'value' and len(a) == 1:
            return a[0], ATOM
        return rf'\operatorname{{{name}}}\left({", ".join(a)}\right)', ATOM


def _is_sum(node):
    return (isinstance(node, ast.Call) and isinstance(node.func, ast.Name) and node.func.id == 'sum'
            and len(node.args) == 1 and isinstance(node.args[0], ast.GeneratorExp))


def _is_skip(node):
    return isinstance(node, ast.Attribute) and node.attr in ('Skip', 'Infeasible', 'Feasible')


def to_latex(src, tipo, local_name, dim_sets, blocks=('ox', 'purge', 'red'), sense=None):
    """LaTeX (KaTeX) de la regla, o None si no se puede traducir con certeza."""
    try:
        tree = ast.parse(textwrap.dedent(src))
    except SyntaxError:
        return None
    fn = next((n for n in tree.body if isinstance(n, (ast.FunctionDef, ast.Lambda))), None)
    if not isinstance(fn, ast.FunctionDef):
        return None
    a = fn.args
    if a.vararg is not None or a.kwarg is not None:
        return None                       # reglas genericas (p.ej. discretizacion de pyomo.dae)
    params = [x.arg for x in a.args][1:]  # el primero es el bloque
    tr = Translator(params, blocks)
    try:
        defs, skips, early, final = [], [], [], None
        for st in fn.body:
            if isinstance(st, ast.Expr) and isinstance(st.value, ast.Constant):
                continue                                         # docstring
            if isinstance(st, ast.If) and len(st.body) == 1 and isinstance(st.body[0], ast.Return) and not st.orelse:
                if _is_skip(st.body[0].value):
                    skips.append(tr.ex(st.test)[0])
                elif tipo == 'Expression':                       # if i not in I: return 0  -> por casos
                    early.append((tr.ex(st.test)[0], tr.ex(st.body[0].value)[0]))
                else:
                    raise Unsupported('conditional return')
                continue
            if isinstance(st, ast.Assign) and len(st.targets) == 1 and isinstance(st.targets[0], (ast.Name, ast.Tuple)):
                t = st.targets[0]
                if isinstance(t, ast.Tuple):                   # z0, z1, z2 = A.at(1), A.at(2), A.at(3)
                    if not (isinstance(st.value, ast.Tuple) and len(st.value.elts) == len(t.elts)):
                        raise Unsupported('tuple assignment')
                    for tn, tv in zip(t.elts, st.value.elts):
                        defs.append((tr.ex(tn)[0], tr.ex(tv)[0]))
                else:
                    defs.append((sym(t.id), tr.ex(st.value)[0]))
                continue
            if isinstance(st, ast.Return) and final is None:
                final = st.value
                continue
            raise Unsupported(type(st).__name__)
        if final is None or _is_skip(final):
            return None

        idx = ','.join(GREEK.get(p, p) for p in params)
        lhs_name = sym(local_name) + (f'({idx})' if params else '')
        if tipo == 'Constraint':
            if isinstance(final, ast.Tuple):                     # (lb, expr, ub)
                parts = [tr.ex(e)[0] for e in final.elts if not (isinstance(e, ast.Constant) and e.value is None)]
                first = r' \le '.join(parts[:1]) + r' &\le ' + r' \le '.join(parts[1:])
            elif isinstance(final, ast.Compare) and len(final.ops) == 1 and type(final.ops[0]) in CMP:
                first = f'{tr.wrap(final.left, ADD)} &{CMP[type(final.ops[0])]} {tr.wrap(final.comparators[0], ADD)}'
            else:
                return None
        elif tipo == 'Objective':
            first = rf'\{"max" if sense == "max" else "min"} \quad & {tr.ex(final)[0]}'
        else:
            rhs = tr.ex(final)[0]
            if early:
                cases = r' \\ '.join(rf'{v} & \text{{if }} {c}' for c, v in early)
                rhs = rf'\begin{{cases}} {cases} \\ {rhs} & \text{{otherwise}} \end{{cases}}'
            first = f'{lhs_name} &= {rhs}'

        lines = [first]
        for i, (lhs, rhs) in enumerate(defs):
            lines.append((r'\text{with}\quad ' if i == 0 else '') + f'{lhs} &= {rhs}')
        if params and len(dim_sets) == len(params):
            dom = r',\ '.join(f'{GREEK.get(p, p)} \\in {_base(s)}' for p, s in zip(params, dim_sets))
            dom = rf'& \forall\, {dom}'
            if skips:
                dom += r',\ \text{except } ' + r'\text{ or }'.join(skips)
            lines.append(dom)
        return r'\begin{aligned} ' + r' \\ '.join(lines) + r' \end{aligned}'
    except (Unsupported, AttributeError, TypeError):
        return None
