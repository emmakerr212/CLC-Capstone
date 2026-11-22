"""
h5_export_full.py

Exporta TODOS los componentes de un modelo Pyomo (Var, Expression, Param, Constraint,
Objective) a un grupo /modelo dentro de un archivo h5py ya abierto, sin tocar nada de
lo que la celda de exportacion de cada notebook ya escribe hoy (grupos <etapa>/*,
decisiones, attrs de la raiz). Tambien exporta /modelo/_sets con los Sets de cada
etapa, adimensionales y (cuando se puede deducir) en unidades fisicas.

Uso:
    import pyomo.environ as pyo
    asegurar_sufijo_dual(modelo)     # antes de resolver
    # ... resolver modelo ...
    with h5py.File(fpath, 'w') as f:
        ... exportar_etapa(...) de siempre, decisiones, etc. (sin cambios) ...
        exportar_modelo_completo(modelo, f, etapa_de_bloque={"ox": "ox", "purge": "purge", "red": "red"})
        exportar_sets_fisicos(modelo, f, etapa_de_bloque={"ox": "ox", "purge": "purge", "red": "red"})

No modifica el modelo ni sus ecuaciones. La unica escritura permitida sobre el modelo
es, si no existe ya, el sufijo de duales (`asegurar_sufijo_dual`).
"""

import re
import inspect
import numpy as np
import h5py
import pyomo.environ as pyo
import pyomo.dae as pdae

try:
    from pyomo.contrib.latex_printer.latex_printer import latex_printer as _latex_printer
    _LATEX_DISPONIBLE = True
    # El impresor reporta sus fallas (muy frecuentes: cualquier regla con if/for,
    # que es la mayoria de nuestras Constraint) via logging, no con una excepcion
    # silenciosa -- sin esto, cada export llena la salida de cientos de "ERROR:".
    # El try/except de _info_regla_y_latex ya captura el fallo real (latex_ok=False);
    # esto solo apaga el ruido del log, no cambia el resultado.
    import logging as _logging
    # El mensaje "templatizing the rule..." sale de pyomo.core.expr.template_expr
    # (la maquinaria general de templatizing que usa el impresor), no del propio
    # modulo latex_printer -- hay que apagar ESE logger especificamente.
    _logging.getLogger("pyomo.core.expr.template_expr").setLevel(_logging.CRITICAL)
    _logging.getLogger("pyomo.contrib.latex_printer").setLevel(_logging.CRITICAL)
except Exception:
    _LATEX_DISPONIBLE = False

# Si el latex generado pasa de esto, lo trato como "no util" (latex_ok=False) en vez
# de guardar un attr de HDF5 de cientos de KB -- le pasa a componentes como el
# objetivo, donde el impresor expande TODA la cadena de sub-expresiones.
_LATEX_MAX_CHARS = 20000


def _info_regla_y_latex(comp):
    """(regla_fuente, regla_fuente_ok, latex, latex_ok) para Expression/Constraint/Objective.
    regla_fuente: codigo fuente de la funcion 'rule' original, via comp._rule._fcn
    (solo existe si el componente se construyo con rule=; con expr= no hay funcion que
    mostrar). latex: forma generica (indexada) via pyomo.contrib.latex_printer; falla
    (y se deja '') para reglas con control de flujo (if/for) que el impresor no puede
    templatizar -- eso es una limitacion conocida del impresor, no un bug nuestro."""
    regla_fuente, regla_fuente_ok = "", False
    rule = getattr(comp, "_rule", None)
    fcn = getattr(rule, "_fcn", None) if rule is not None else None
    if fcn is not None:
        try:
            regla_fuente = inspect.getsource(fcn)
            regla_fuente_ok = True
        except Exception:
            pass

    latex_str, latex_ok = "", False
    if _LATEX_DISPONIBLE:
        try:
            s = _latex_printer(comp)
            if len(s) <= _LATEX_MAX_CHARS:
                latex_str, latex_ok = s, True
        except Exception:
            pass

    return regla_fuente, regla_fuente_ok, latex_str, latex_ok


# ---------------------------------------------------------------------------
# Reglas de escalamiento -- UNA sola fuente de verdad, derivada de como el
# propio modelo define Cnd_i, Cnd_s y Tnd (ver CLC COMPLETE CYCLE.ipynb,
# bloques ox/purge/red). Si un nombre de componente no esta aqui ni en
# UNIDADES_FIJAS mas abajo, su escala queda NaN / escala_conocida=False --
# nunca se inventa un factor.
# ---------------------------------------------------------------------------

def _especie_de(idx_tuple):
    """La primera posicion del indice es casi siempre la especie/solido/corriente
    en este modelo (C_hat[i,z,tao], Ccs_hat[s,z,r,tao], etc.)."""
    if isinstance(idx_tuple, tuple) and len(idx_tuple) > 0:
        return idx_tuple[0]
    return idx_tuple


def _regla_Cnd_i(blk, idx_tuple):
    sp = _especie_de(idx_tuple)
    try:
        return float(pyo.value(blk.Cnd_i[sp])), "mol/m3"
    except Exception:
        return float("nan"), ""


def _regla_Cnd_s(blk, idx_tuple):
    sp = _especie_de(idx_tuple)
    try:
        return float(pyo.value(blk.Cnd_s[sp])), "mol/m3"
    except Exception:
        return float("nan"), ""


def _regla_Tnd(blk, idx_tuple):
    try:
        return float(pyo.value(blk.Tnd)), "K"
    except Exception:
        return float("nan"), ""


def _regla_Cnd_i_sobre_tf(blk, idx_tuple):
    sp = _especie_de(idx_tuple)
    try:
        return float(pyo.value(blk.Cnd_i[sp]) / pyo.value(blk.tf)), "mol/m3/s"
    except Exception:
        return float("nan"), ""


def _regla_Cnd_i_sobre_L(blk, idx_tuple):
    sp = _especie_de(idx_tuple)
    try:
        return float(pyo.value(blk.Cnd_i[sp]) / pyo.value(blk.model().L)), "mol/m4"
    except Exception:
        return float("nan"), ""


def _regla_Cnd_i_sobre_L2(blk, idx_tuple):
    sp = _especie_de(idx_tuple)
    try:
        return float(pyo.value(blk.Cnd_i[sp]) / pyo.value(blk.model().L) ** 2), "mol/m5"
    except Exception:
        return float("nan"), ""


def _regla_Cnd_i_sobre_Rp(blk, idx_tuple):
    sp = _especie_de(idx_tuple)
    try:
        return float(pyo.value(blk.Cnd_i[sp]) / pyo.value(blk.Rp)), "mol/m4"
    except Exception:
        return float("nan"), ""


def _regla_Cnd_i_sobre_Rp2(blk, idx_tuple):
    sp = _especie_de(idx_tuple)
    try:
        return float(pyo.value(blk.Cnd_i[sp]) / pyo.value(blk.Rp) ** 2), "mol/m5"
    except Exception:
        return float("nan"), ""


def _regla_Tnd_sobre_tf(blk, idx_tuple):
    try:
        return float(pyo.value(blk.Tnd) / pyo.value(blk.tf)), "K/s"
    except Exception:
        return float("nan"), ""


def _regla_Tnd_sobre_L(blk, idx_tuple):
    try:
        return float(pyo.value(blk.Tnd) / pyo.value(blk.model().L)), "K/m"
    except Exception:
        return float("nan"), ""


def _regla_Tnd_sobre_L2(blk, idx_tuple):
    try:
        return float(pyo.value(blk.Tnd) / pyo.value(blk.model().L) ** 2), "K/m2"
    except Exception:
        return float("nan"), ""


def _regla_Tnd_sobre_Rp(blk, idx_tuple):
    try:
        return float(pyo.value(blk.Tnd) / pyo.value(blk.Rp)), "K/m"
    except Exception:
        return float("nan"), ""


def _regla_Tnd_sobre_Rp2(blk, idx_tuple):
    try:
        return float(pyo.value(blk.Tnd) / pyo.value(blk.Rp) ** 2), "K/m2"
    except Exception:
        return float("nan"), ""


def _regla_Cnd_s_sobre_tf(blk, idx_tuple):
    sp = _especie_de(idx_tuple)
    try:
        return float(pyo.value(blk.Cnd_s[sp]) / pyo.value(blk.tf)), "mol/m3/s"
    except Exception:
        return float("nan"), ""


# nombre_local_exacto -> (funcion_regla, varia_por_elemento)
# varia_por_elemento=True cuando la escala depende de la especie del indice
# (Cnd_i/Cnd_s distinto por especie); False cuando es un solo valor para todo
# el componente (Tnd es igual sin importar el indice).
REGLAS_ESCALA = {
    "C_hat":        (_regla_Cnd_i, True),
    "C_feed_hat":   (_regla_Cnd_i, True),
    "Cci_hat":      (_regla_Cnd_i, True),
    "Ccs_hat":      (_regla_Cnd_s, True),
    "T_hat":        (_regla_Tnd, False),
    "Tc_hat":       (_regla_Tnd, False),
    "T_in_hat":     (_regla_Tnd, False),
    "dC_hat_dT":    (_regla_Cnd_i_sobre_tf, True),
    "dCci_hat_dT":  (_regla_Cnd_i_sobre_tf, True),
    "dC_hat_dA":    (_regla_Cnd_i_sobre_L, True),
    "dC_hat_dA2":   (_regla_Cnd_i_sobre_L2, True),
    "dCci_hat_dR":  (_regla_Cnd_i_sobre_Rp, True),
    "dCci_hat_dR2": (_regla_Cnd_i_sobre_Rp2, True),
    "dT_hat_dT":    (_regla_Tnd_sobre_tf, False),
    "dT_hat_dA":    (_regla_Tnd_sobre_L, False),
    "dT_hat_dA2":   (_regla_Tnd_sobre_L2, False),
    "dTc_hat_dT":   (_regla_Tnd_sobre_tf, False),
    "dTc_hat_dR":   (_regla_Tnd_sobre_Rp, False),
    "dTc_hat_dR2":  (_regla_Tnd_sobre_Rp2, False),
    "dCcs_hat_dT":  (_regla_Cnd_s_sobre_tf, True),
}

# ---------------------------------------------------------------------------
# Componentes SIN adimensionalizar (ya son fisicos desde que se definieron),
# cuyo doc no trae la unidad entre corchetes -- se las asigno a mano, leyendo
# la formula/doc real de cada una (no son una suposicion: son exactamente lo
# que dice su definicion en el notebook). escala=1 para todas estas.
#
# IMPORTANTE (nota de consistencia, ya reportada antes): Fan_cost, Heating_cost,
# Compressor_cost, Turbine_income, stream_cost_year y sus sub-terminos NO
# aplican CAD_conversion en su formula (a diferencia de capital/cRM_OC que si
# la aplican) -- por eso su unidad real queda como "$/año" generico (moneda
# ambigua, tal cual esta en el codigo), no "$CAD/año". No es un numero
# inventado, es fiel a la formula; la inconsistencia de moneda entre terminos
# es un problema del modelo, ya señalado aparte.
# ---------------------------------------------------------------------------
UNIDADES_FIJAS = {
    # factores y constantes adimensionales
    "CAD_conversion": "-", "CEPCI25": "-", "CEPCI94": "-", "MS": "-",
    "Fc": "-", "Fm": "-", "Fp": "-", "USD_yr_conv": "-", "tn_kg_conv": "kg/ton",
    "i": "-", "AF": "1/año",
    "eta_comp": "-", "eta_fan": "-", "eta_turbina": "-",
    # tiempo / lifetime / ciclos
    "n": "años", "t_op": "h", "t_cycle": "h", "n_cycles": "1/año", "n_OC_changes": "1/año",
    # geometria derivada
    "D_ft": "ft", "L_ft": "ft",
    # precios y condiciones de referencia
    "P_atm": "bar", "T_ambiente": "K", "Pe": "$/J", "c_heat": "$/J",
    "cRM": "$USD(2014)/ton", "cRM_OC": "$CAD/kg",
    "air_cost": "$/kg", "n2_cost": "$/kg", "ch4_cost": "$/kg", "co_cost": "$/kg", "h2_cost": "$/kg",
    "co2_incentive": "$/kg",
    # CAPEX (si aplica CAD_conversion en su formula)
    "capital": "$CAD", "capital_annualized": "$CAD/año",
    # OPEX / INCOME por etapa (moneda ambigua, ver nota arriba)
    "fan_cost": "$/año", "compressor_cost": "$/año", "heating_cost": "$/año", "turbine_income": "$/año",
    "Fan_cost": "$/año", "Compressor_cost": "$/año", "Heating_cost": "$/año", "Turbine_income": "$/año",
    "OC_cost": "$/año", "stream_cost_year": "$/año",
    "total_air_cost_year": "$/año", "total_n2_cost_year": "$/año", "total_ch4_cost_year": "$/año",
    "total_co_cost_year": "$/año", "total_h2_cost_year": "$/año",
    # masas/moles de alimentacion y salida, por ciclo y por año
    "total_feed_mass_air": "kg", "total_feed_mass_n2": "kg", "total_feed_mass_ch4": "kg",
    "total_feed_mass_co": "kg", "total_feed_mass_h2": "kg",
    "total_feed_mass_air_year": "kg/año", "total_feed_mass_n2_year": "kg/año",
    "total_feed_mass_ch4_year": "kg/año", "total_feed_mass_co_year": "kg/año", "total_feed_mass_h2_year": "kg/año",
    "total_outlet_mass_O2_year": "kg/año", "total_outlet_mass_N2_year": "kg/año",
    "total_outlet_mass_CH4_year": "kg/año", "total_outlet_mass_CO_year": "kg/año",
    "total_outlet_mass_H2_year": "kg/año", "total_outlet_mass_CO2_year": "kg/año",
    "total_outlet_mass_H2O_year": "kg/año",
    "total_outlet_moles_O2_year": "mol/año", "total_outlet_moles_N2_year": "mol/año",
    "total_outlet_moles_CH4_year": "mol/año", "total_outlet_moles_CO_year": "mol/año",
    "total_outlet_moles_H2_year": "mol/año", "total_outlet_moles_CO2_year": "mol/año",
    "total_outlet_moles_H2O_year": "mol/año",
}

_RE_UNIDAD_EN_DOC = re.compile(r"\[([^\[\]]+)\]\s*$")


def _unidad_desde_doc(doc):
    if not doc:
        return ""
    m = _RE_UNIDAD_EN_DOC.search(doc.strip())
    return m.group(1) if m else ""


def _info_escala(comp, local_name, claves, blk, tipo):
    """Devuelve (escala_attr, unidades, escala_conocida, escala_por_elemento, vf_calculator).
    vf_calculator(j, k, valor_j) -> valor_fisico en el elemento j (indice k), o NaN."""
    if tipo in ("Constraint", "Objective"):
        # El "valor" de una restriccion/objetivo es un residuo o una suma de
        # terminos con escalas mezcladas -- no hay un factor unico honesto.
        return float("nan"), "", False, False, None

    regla_info = REGLAS_ESCALA.get(local_name)
    if regla_info is not None:
        regla_fn, varia = regla_info
        blk_real = blk if blk is not None else comp.parent_block()
        if varia:
            def vf_calc(j, k, valor_j, _fn=regla_fn, _blk=blk_real):
                e, _u = _fn(_blk, k)
                return valor_j * e if not np.isnan(e) else float("nan")
            # unidades: tomo la de cualquier elemento valido (son todas la misma unidad fisica)
            _, unidades = regla_fn(blk_real, claves[0] if claves else None)
            return float("nan"), unidades, True, True, vf_calc
        else:
            escala, unidades = regla_fn(blk_real, claves[0] if claves else None)
            if np.isnan(escala):
                return float("nan"), unidades, False, False, None
            return escala, unidades, True, False, (lambda j, k, v, e=escala: v * e)

    if local_name in UNIDADES_FIJAS:
        unidad = UNIDADES_FIJAS[local_name]
        return 1.0, unidad, True, False, (lambda j, k, v: v * 1.0)

    doc = comp.doc if hasattr(comp, "doc") else ""
    unidad = _unidad_desde_doc(doc)
    if unidad:
        return 1.0, unidad, True, False, (lambda j, k, v: v * 1.0)

    return float("nan"), "", False, False, None


# ---------------------------------------------------------------------------
# Utilidades genericas de introspeccion Pyomo
# ---------------------------------------------------------------------------

def asegurar_sufijo_dual(modelo):
    """Agrega m.dual (IMPORT) si no existe ya. Es la UNICA escritura permitida
    sobre el modelo antes de resolver."""
    if not hasattr(modelo, "dual"):
        modelo.dual = pyo.Suffix(direction=pyo.Suffix.IMPORT)
    return modelo


def _nombre_set(s, pos):
    nombre = s.getname() if hasattr(s, "getname") else str(s)
    if not nombre or nombre.startswith("Any") or nombre.startswith("_"):
        return f"dim{pos}"
    return nombre


def _indices_y_sets(comp):
    if not comp.is_indexed():
        return []
    idx_set = comp.index_set()
    try:
        subsets = list(idx_set.subsets())
    except Exception:
        subsets = [idx_set]
    return [_nombre_set(s, i) for i, s in enumerate(subsets)]


def _es_derivada(comp):
    return isinstance(comp, pdae.DerivativeVar)


def _valor_seguro(x):
    try:
        v = pyo.value(x, exception=False)
    except Exception:
        v = None
    return float("nan") if v is None else float(v)


def _escribir_dataset(grupo, nombre, datos, dtype=None):
    datos = np.asarray(datos, dtype=dtype)
    kwargs = {}
    if datos.size > 1000:
        kwargs["compression"] = "gzip"
    grupo.create_dataset(nombre, data=datos, **kwargs)


def _bloque_y_etapa(comp, etapa_de_bloque):
    parent = comp.parent_block()
    modelo_raiz = comp.model()
    if parent is modelo_raiz:
        return "", ""
    nombre_bloque = parent.local_name
    etapa = etapa_de_bloque.get(nombre_bloque, "")
    return nombre_bloque, etapa


def _tipo_de(comp):
    if isinstance(comp, pyo.Objective):
        return "Objective"
    if isinstance(comp, pyo.Constraint):
        return "Constraint"
    if isinstance(comp, pyo.Var):
        return "Var"
    if isinstance(comp, pyo.Expression):
        return "Expression"
    if isinstance(comp, pyo.Param):
        return "Param"
    return None


# ---------------------------------------------------------------------------
# Exportador principal: /modelo/<nombre_completo>
# ---------------------------------------------------------------------------

def exportar_modelo_completo(modelo, archivo_h5, etapa_de_bloque=None):
    etapa_de_bloque = etapa_de_bloque or {}
    hay_duales = hasattr(modelo, "dual") and len(modelo.dual) > 0

    g_modelo = archivo_h5.require_group("modelo")

    tipos_pyomo = (pyo.Var, pyo.Expression, pyo.Param, pyo.Constraint, pyo.Objective)
    n_componentes = 0

    for tipo_pyomo in tipos_pyomo:
        for comp in modelo.component_objects(tipo_pyomo, active=None, descend_into=True):
            tipo = _tipo_de(comp)
            if tipo is None:
                continue

            nombre_completo = comp.getname(fully_qualified=True).replace("/", "|")
            local_name = comp.local_name
            bloque, etapa = _bloque_y_etapa(comp, etapa_de_bloque)

            claves = list(comp.keys()) if comp.is_indexed() else [None]
            n = len(claves)
            nombres_dims = _indices_y_sets(comp)

            g = g_modelo.create_group(nombre_completo)
            g.attrs["nombre"] = nombre_completo
            g.attrs["tipo"] = tipo
            g.attrs["bloque"] = bloque
            g.attrs["etapa"] = etapa
            g.attrs["doc"] = comp.doc if comp.doc else ""
            g.attrs["indices"] = nombres_dims
            g.attrs["n"] = n

            if claves[0] is not None:
                ndims = len(claves[0]) if isinstance(claves[0], tuple) else 1
                for d in range(ndims):
                    valores_dim = [k[d] if isinstance(k, tuple) else k for k in claves]
                    try:
                        arr = np.asarray(valores_dim, dtype=float)
                        _escribir_dataset(g, f"idx_{d}", arr, dtype=np.float64)
                    except (ValueError, TypeError):
                        arr = np.array([str(v) for v in valores_dim], dtype=object)
                        g.create_dataset(f"idx_{d}", data=arr, dtype=h5py.string_dtype())

            valores = np.full(n, np.nan, dtype=np.float64)
            for j, k in enumerate(claves):
                data = comp[k] if comp.is_indexed() else comp
                valores[j] = _valor_seguro(data)
            _escribir_dataset(g, "valor", valores, dtype=np.float64)

            blk_padre = comp.parent_block()
            escala_attr, unidades, escala_conocida, escala_por_elem, vf_calc = _info_escala(
                comp, local_name, claves, blk_padre, tipo
            )

            if vf_calc is not None:
                vf = np.array([vf_calc(j, k, valores[j]) for j, k in enumerate(claves)], dtype=np.float64)
                _escribir_dataset(g, "valor_fisico", vf, dtype=np.float64)
                if escala_por_elem and np.all(np.isnan(vf)):
                    escala_conocida = False

            g.attrs["unidades"] = unidades
            g.attrs["escala"] = escala_attr
            g.attrs["escala_conocida"] = bool(escala_conocida)
            g.attrs["escala_por_elemento"] = bool(escala_por_elem)

            if tipo in ("Expression", "Constraint", "Objective"):
                regla_fuente, regla_fuente_ok, latex_str, latex_ok = _info_regla_y_latex(comp)
                g.attrs["regla_fuente"] = regla_fuente
                g.attrs["regla_fuente_ok"] = regla_fuente_ok
                g.attrs["latex"] = latex_str
                g.attrs["latex_ok"] = latex_ok

                expresiones = []
                for k in claves:
                    data = comp[k] if comp.is_indexed() else comp
                    try:
                        expresiones.append(str(data.expr))
                    except Exception:
                        expresiones.append("")
                arr_expr = np.array(expresiones, dtype=object)
                g.create_dataset("expresion", data=arr_expr, dtype=h5py.string_dtype(), compression="gzip")

            if tipo == "Var":
                lb = np.full(n, np.nan); ub = np.full(n, np.nan)
                fijo = np.zeros(n, dtype=bool)
                for j, k in enumerate(claves):
                    data = comp[k] if comp.is_indexed() else comp
                    lb[j] = data.lb if data.lb is not None else float("nan")
                    ub[j] = data.ub if data.ub is not None else float("nan")
                    fijo[j] = bool(data.fixed)
                _escribir_dataset(g, "lb", lb, dtype=np.float64)
                _escribir_dataset(g, "ub", ub, dtype=np.float64)
                _escribir_dataset(g, "fijo", fijo, dtype=bool)
                g.attrs["es_derivada"] = _es_derivada(comp)

            elif tipo == "Constraint":
                cuerpo = np.full(n, np.nan); inferior = np.full(n, np.nan)
                superior = np.full(n, np.nan); residuo = np.full(n, np.nan)
                activa = np.zeros(n, dtype=bool)
                dual = np.full(n, np.nan)
                for j, k in enumerate(claves):
                    data = comp[k] if comp.is_indexed() else comp
                    b = _valor_seguro(data.body)
                    lo = _valor_seguro(data.lower) if data.lower is not None else float("nan")
                    up = _valor_seguro(data.upper) if data.upper is not None else float("nan")
                    cuerpo[j] = b; inferior[j] = lo; superior[j] = up
                    viol_lo = (lo - b) if not np.isnan(lo) else -np.inf
                    viol_up = (b - up) if not np.isnan(up) else -np.inf
                    residuo[j] = max(viol_lo, viol_up, 0.0)
                    activa[j] = bool(data.active)
                    if hay_duales:
                        try:
                            dual[j] = modelo.dual.get(data, float("nan"))
                        except Exception:
                            dual[j] = float("nan")
                _escribir_dataset(g, "cuerpo", cuerpo, dtype=np.float64)
                _escribir_dataset(g, "inferior", inferior, dtype=np.float64)
                _escribir_dataset(g, "superior", superior, dtype=np.float64)
                _escribir_dataset(g, "residuo", residuo, dtype=np.float64)
                _escribir_dataset(g, "activa", activa, dtype=bool)
                _escribir_dataset(g, "dual", dual, dtype=np.float64)

            elif tipo == "Objective":
                g.attrs["sentido"] = "max" if comp.sense == pyo.maximize else "min"
                g.attrs["activa"] = bool(comp.active)

            n_componentes += 1

    g_modelo.attrs["version_export"] = 1
    g_modelo.attrs["n_componentes"] = n_componentes
    g_modelo.attrs["duales_disponibles"] = bool(hay_duales)

    return n_componentes


# ---------------------------------------------------------------------------
# /modelo/_sets -- Sets de cada etapa, adimensionales y (si se puede) fisicos
# ---------------------------------------------------------------------------

# nombre_local_del_set -> (unidad_fisica, funcion(blk) -> factor_de_escala)
_SETS_NUMERICOS = {
    "A": ("m", lambda blk: pyo.value(blk.model().L)),       # z = A_norm * L
    "R": ("m", lambda blk: pyo.value(blk.Rp)),               # r = R_norm * Rp
    "T": ("s", lambda blk: pyo.value(blk.tf)),               # t = T_norm * tf
}
_SETS_TEXTO = {"I", "S", "F", "GI", "PS", "PG", "X"}  # X = reacciones (p.ej. "R4"), no confundir con Xconv


def exportar_sets_fisicos(modelo, archivo_h5, etapa_de_bloque=None):
    """Escribe /modelo/_sets/<etapa>.<nombre_set> con el valor adimensional y,
    cuando se puede deducir de la definicion del modelo (L, tf, Rp de cada
    etapa), el valor fisico. Devuelve la lista de sets sin conversion fisica
    conocida."""
    etapa_de_bloque = etapa_de_bloque or {}
    g_sets = archivo_h5.require_group("modelo/_sets")
    sin_fisica = []

    vistos = set()
    for blk_name, etapa in etapa_de_bloque.items():
        blk = getattr(modelo, blk_name, None)
        if blk is None:
            continue
        sets_del_bloque = list(blk.component_objects(pyo.Set, descend_into=False)) + \
                           list(blk.component_objects(pdae.ContinuousSet, descend_into=False))
        for s in sets_del_bloque:
            nombre_set = s.local_name
            clave = f"{blk_name}.{nombre_set}"
            if clave in vistos:
                continue
            vistos.add(clave)

            valores = list(s)
            g = g_sets.create_group(clave)
            g.attrs["nombre"] = nombre_set
            g.attrs["etapa"] = etapa

            if nombre_set in _SETS_NUMERICOS:
                unidad, factor_fn = _SETS_NUMERICOS[nombre_set]
                arr = np.asarray(valores, dtype=np.float64)
                _escribir_dataset(g, "valor", arr, dtype=np.float64)
                try:
                    factor = factor_fn(blk)
                    g.attrs["unidades_fisicas"] = unidad
                    regla_txt = {"A": "z = A·L", "R": "r = R·Rp", "T": "t = T·tf"}[nombre_set]
                    g.attrs["regla"] = regla_txt
                    _escribir_dataset(g, "valor_fisico", arr * factor, dtype=np.float64)
                except Exception:
                    g.attrs["unidades_fisicas"] = unidad
                    g.attrs["regla"] = "no se pudo evaluar"
                    sin_fisica.append(clave)
            elif nombre_set in _SETS_TEXTO:
                arr = np.array([str(v) for v in valores], dtype=object)
                g.create_dataset("valor", data=arr, dtype=h5py.string_dtype())
                g.attrs["unidades_fisicas"] = ""
                g.attrs["regla"] = "categorico, no aplica conversion"
            else:
                # set numerico o de texto no previsto -- lo guardo igual, sin fisica
                try:
                    arr = np.asarray(valores, dtype=np.float64)
                    _escribir_dataset(g, "valor", arr, dtype=np.float64)
                except (ValueError, TypeError):
                    arr = np.array([str(v) for v in valores], dtype=object)
                    g.create_dataset("valor", data=arr, dtype=h5py.string_dtype())
                g.attrs["unidades_fisicas"] = ""
                g.attrs["regla"] = "set no reconocido, sin regla de conversion"
                sin_fisica.append(clave)

    return sin_fisica
