@echo off
REM ============================================================================
REM  Protein Backbones as Piecewise Cornu Spirals on T^2 (v3)
REM  Pure-numpy DSSP (no external binary needed) + dihedral fallback
REM ============================================================================

echo.
echo ========================================================================
echo   PROTEIN BACKBONES AS PIECEWISE CORNU SPIRALS ON T^2  (v3)
echo   Pure-numpy DSSP ^| No external binary needed
echo ========================================================================
echo.

REM ── Check for Python ──
python --version >nul 2>&1
if errorlevel 1 (
    echo ERROR: Python not found in PATH.
    echo Please install Python 3.8+ from https://www.python.org/
    pause
    exit /b 1
)

echo [1/4] Installing Python dependencies...
pip install numpy scipy matplotlib requests --quiet --upgrade
if errorlevel 1 (
    echo WARNING: Some packages may have failed. Script will continue.
)

echo.
echo [2/4] Creating directory structure...
if not exist "data" mkdir data
if not exist "data\pdb_files" mkdir data\pdb_files
if not exist "results" mkdir results
if not exist "results\figures" mkdir results\figures

echo.
echo [3/4] Running Cornu spiral analysis pipeline (v3)...
echo   This will:
echo     - Download 23 PDB structures from RCSB
echo     - Run embedded pure-numpy DSSP for SS assignment (no binary needed)
echo     - Compute differential geometry on the Ramachandran torus T^2
echo     - Test winding-rate constancy within SS elements
echo     - Classify segments as logarithmic/Fermat/Cornu spirals
echo     - Compute per-protein torus knot (p,q) descriptors
echo     - Test correlations with folding rate
echo     - Generate 6 publication-quality figures
echo.

python cornu_spirals_on_T2.py

echo.
echo [4/4] Analysis complete!
echo.
echo   Results:
echo     Report:     results\cornu_analysis_report.md
echo     Stats:      results\stats_summary.json
echo     Proteins:   results\protein_results.json
echo     Figures:    results\figures\
echo     Log:        cornu_analysis.log
echo.

if exist "results\cornu_analysis_report.md" (
    echo Opening report...
    start "" "results\cornu_analysis_report.md"
)

pause
