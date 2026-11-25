"""
Dashboard de CLC (Flask) — explorador y comparador de los H5 de resultados.

Paginas:
  /           portada con los destinos
  /explorer   una corrida: perfiles en z (gas, solido, T, X) con sliders de t y r, mapas z-t,
              series a la salida, flujo de entrada y el resumen de parametros/decisiones
  /compare    hasta 4 corridas superpuestas contra el tiempo global + tabla de parametros

Toda la lectura/normalizacion de los H5 vive en h5_runs.py; este archivo solo conecta rutas.
Los datos de una corrida se mandan completos al navegador (el H5 mas grande pesa ~0.5 MB en JSON),
asi que mover los sliders no vuelve a pedir nada al servidor.

Uso:  python app.py   (o run_dashboard.bat)  ->  http://localhost:5050
Carpeta de corridas: ../Runs, o la variable de entorno CLC_RUNS_DIR.
"""

import os
import subprocess
import sys
from functools import lru_cache
from pathlib import Path

from flask import Flask, Response, abort, jsonify, render_template, request

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import breakdown  # noqa: E402
import h5_runs  # noqa: E402
import model_tables  # noqa: E402

RUNS_DIR = Path(os.environ.get('CLC_RUNS_DIR', HERE.parent / 'Runs')).resolve()
PORT = int(os.environ.get('CLC_DASHBOARD_PORT', 5050))

app = Flask(__name__)
app.config['TEMPLATES_AUTO_RELOAD'] = True
# que el navegador revalide JS/CSS en cada carga: si no, tras actualizar el dashboard puede
# seguir mostrando la version vieja hasta un Ctrl+F5
app.config['SEND_FILE_MAX_AGE_DEFAULT'] = 0
app.json.sort_keys = False
app.json.ensure_ascii = False


def _run_path(name):
    """Solo archivos .h5 que esten directamente dentro de RUNS_DIR (nada de rutas relativas)."""
    if not name or Path(name).name != name or not name.lower().endswith('.h5'):
        abort(404)
    path = RUNS_DIR / name
    if not path.is_file():
        abort(404)
    return path


# ---------------------------------------------------------------------------
# Paginas
# ---------------------------------------------------------------------------

@app.route('/')
def landing():
    runs = h5_runs.list_runs(RUNS_DIR)
    by_kind = {}
    for r in runs:
        by_kind[r['kind']] = by_kind.get(r['kind'], 0) + 1
    return render_template('landing.html', page='home', n_runs=len(runs), by_kind=by_kind,
                           latest=runs[0] if runs else None, runs_dir=str(RUNS_DIR))


@app.route('/explorer')
def explorer():
    return render_template('explorer.html', page='explorer')


@app.route('/compare')
def compare():
    return render_template('compare.html', page='compare')


@app.route('/tables')
def tables():
    return render_template('tables.html', page='tables')


# ---------------------------------------------------------------------------
# API
# ---------------------------------------------------------------------------

@app.route('/api/runs')
def api_runs():
    return jsonify(runs=h5_runs.list_runs(RUNS_DIR), runs_dir=str(RUNS_DIR))


@app.route('/api/run/<path:name>')
def api_run(name):
    path = _run_path(name)
    try:
        return jsonify(h5_runs.load_run_cached(path))
    except Exception as e:  # el front-end muestra el mensaje en el panel
        return jsonify(error=f'{type(e).__name__}: {e}'), 422


@app.route('/api/model/<path:name>/components')
def api_model_components(name):
    path = _run_path(name)
    try:
        return jsonify(model_tables.list_components(path))
    except Exception as e:
        return jsonify(error=f'{type(e).__name__}: {e}'), 422


@app.route('/api/model/<path:name>/component')
def api_model_component(name):
    path = _run_path(name)
    try:
        return jsonify(model_tables.read_component(path, request.args.get('key', '')))
    except Exception as e:
        return jsonify(error=f'{type(e).__name__}: {e}'), 422


@app.route('/api/model/<path:name>/breakdown')
def api_model_breakdown(name):
    """Terminos de un componente escalar con regla (la funcion objetivo), evaluados y verificados."""
    path = _run_path(name)
    values = request.args.get('values')  # otra corrida: su FO desglosada con la regla de esta
    values_path = _run_path(values) if values else None
    try:
        return jsonify(breakdown.breakdown(path, request.args.get('key', ''), values_path))
    except Exception as e:
        return jsonify(error=f'{type(e).__name__}: {e}'), 422


@app.route('/api/reveal', methods=['POST'])
def api_reveal():
    """Abre el Explorador de Windows con el archivo seleccionado."""
    path = _run_path((request.get_json(silent=True) or {}).get('name', ''))
    if sys.platform == 'win32':
        subprocess.Popen(['explorer', '/select,', str(path)])
        return jsonify(ok=True)
    return jsonify(ok=False, error='Only available on Windows.'), 400


@lru_cache(maxsize=1)
def _plotly_js():
    from plotly.offline import get_plotlyjs  # el mismo plotly.js que trae el paquete de Python
    return get_plotlyjs()


@app.route('/vendor/plotly.min.js')
def plotly_js():
    return Response(_plotly_js(), mimetype='application/javascript',
                    headers={'Cache-Control': 'public, max-age=86400'})


def _already_running():
    """El puerto ya responde: hay un dashboard abierto (p.ej. otro doble clic al acceso directo)."""
    import socket
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        s.settimeout(0.3)
        return s.connect_ex(('127.0.0.1', PORT)) == 0


if __name__ == '__main__':
    import threading
    import webbrowser
    url = f'http://localhost:{PORT}'
    open_browser = '--open' in sys.argv
    if _already_running():
        print(f'[CLC dashboard] Already running at {url}')
        if open_browser:
            webbrowser.open(url)
        sys.exit(0)
    print(f'[CLC dashboard] Runs folder: {RUNS_DIR}')
    print(f'[CLC dashboard] Open: {url}')
    if open_browser:  # el navegador se abre cuando el servidor ya esta escuchando
        threading.Timer(1.0, webbrowser.open, args=(url,)).start()
    app.run(host='127.0.0.1', port=PORT, debug=False, threaded=True)
