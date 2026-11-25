@echo off
cd /d "%~dp0"

rem --- Python: conda env "UoW" (the one with h5py + flask + plotly). ---
set "PY=%USERPROFILE%\miniconda3\envs\UoW\python.exe"
if not exist "%PY%" (
    echo ERROR: conda env "UoW" not found at %PY%
    pause
    exit /b 1
)

echo Starting the CLC dashboard...
echo It opens in the browser by itself when ready: http://localhost:5050
echo Keep this window open while you use it. Ctrl+C or close the window to stop it.
echo.
"%PY%" app.py --open
if errorlevel 1 pause
