@echo off
REM ═══════════════════════════════════════════════════════════════
REM  BACKBONE STRAIN ANALYZER — SELF TEST
REM  Runs 4 scenarios to verify the tool works correctly.
REM ═══════════════════════════════════════════════════════════════
REM
REM  Place this .bat in the same folder as strain_analyzer.py
REM  and your .cif files. Requires Python 3 with numpy and scipy.
REM
REM  Test 1: ER-alpha agonist vs antagonist (just run through)
REM  Test 2: HIV-PR apo vs indinavir (select option B at each step)
REM  Test 3: ERK phosphorylation (select last option at each step)
REM  Test 4: Thermolysin with chain E (non-default chain)
REM ═══════════════════════════════════════════════════════════════

echo.
echo ================================================================
echo   SELF TEST — Backbone Strain Analyzer
echo   %date% %time%
echo ================================================================
echo.

REM ── Check dependencies ──
python --version >nul 2>&1
if %errorlevel% neq 0 (
    echo ERROR: Python not found. Install Python 3.8+ and add to PATH.
    pause
    exit /b 1
)

python -c "import numpy; import scipy" >nul 2>&1
if %errorlevel% neq 0 (
    echo ERROR: numpy/scipy not found. Run: pip install numpy scipy
    pause
    exit /b 1
)

if not exist strain_analyzer.py (
    echo ERROR: strain_analyzer.py not found in current directory.
    pause
    exit /b 1
)

REM ── Check CIF files ──
set MISSING=0
for %%F in (1ERE 3ERT 1HHP 1HSG 1ERK 2ERK 4TMN 8TLN) do (
    if not exist %%F.cif (
        echo MISSING: %%F.cif
        set MISSING=1
    )
)
if %MISSING% equ 1 (
    echo.
    echo Download missing files from https://files.rcsb.org/download/XXXX.cif
    pause
    exit /b 1
)

echo   All dependencies OK. All CIF files found.
echo.

REM ═══════════════════════════════════════════════════════════════
REM  TEST 1: ER-alpha — Just press Enter through everything
REM  Estradiol (agonist) vs Tamoxifen (antagonist), Chain A
REM  Expected: dGini positive (stiffens), ~75%% distal
REM ═══════════════════════════════════════════════════════════════
echo ================================================================
echo   TEST 1: ER-alpha — Tamoxifen (skip all options)
echo   Expected: Stiffened, substantial effect, ~75%% distal
echo ================================================================
echo.

(
echo.
echo.
echo.
echo.
echo.
echo.
echo.
echo.
) | python strain_analyzer.py 1ERE.cif 3ERT.cif A

echo.
echo ────────────────────────────────────────────────────────────────
echo   TEST 1 COMPLETE. Check: Did it say STIFFENED?
echo ────────────────────────────────────────────────────────────────
echo.
echo.

REM ═══════════════════════════════════════════════════════════════
REM  TEST 2: HIV-PR — Select option B at each step
REM  Apo vs Indinavir, Chain A
REM  Expected: dGini ~+0.034, flap residues in hotspots
REM ═══════════════════════════════════════════════════════════════
echo ================================================================
echo   TEST 2: HIV-PR — Indinavir (option B at each step)
echo   Expected: Stiffened, flap residues 44-47 in hotspots
echo ================================================================
echo.

(
echo B
echo B
echo B
echo B
echo B
echo B
echo B
echo B
) | python strain_analyzer.py 1HHP.cif 1HSG.cif A

echo.
echo ────────────────────────────────────────────────────────────────
echo   TEST 2 COMPLETE. Check: No NaN? Wrapped angles shown?
echo   Check: Residues 44-47 in top 20 hotspots?
echo ────────────────────────────────────────────────────────────────
echo.
echo.

REM ═══════════════════════════════════════════════════════════════
REM  TEST 3: ERK — Select last option (C where available, else B)
REM  Unphosphorylated vs Phosphorylated, Chain A
REM  Expected: dGini positive, 95%% distal, phosphorylation demo
REM ═══════════════════════════════════════════════════════════════
echo ================================================================
echo   TEST 3: ERK — Phosphorylation (last option at each step)
echo   Expected: Stiffened, ~95%% distal
echo ================================================================
echo.

(
echo B
echo B
echo C
echo B
echo B
echo B
echo C
echo B
) | python strain_analyzer.py 1ERK.cif 2ERK.cif A

echo.
echo ────────────────────────────────────────────────────────────────
echo   TEST 3 COMPLETE. Check: Frenet explanation displayed?
echo   Check: Did it say STIFFENED?
echo ────────────────────────────────────────────────────────────────
echo.
echo.

REM ═══════════════════════════════════════════════════════════════
REM  TEST 4: Thermolysin — Chain E (non-default)
REM  Apo vs Phosphoramidon, Chain E
REM  Expected: dGini positive, 100%% distal, small sphere
REM ═══════════════════════════════════════════════════════════════
echo ================================================================
echo   TEST 4: Thermolysin — Chain E (non-default chain)
echo   Expected: Stiffened, 100%% distal, sphere ~8%%
echo ================================================================
echo.

(
echo A
echo.
echo.
echo.
echo.
echo A
echo A
echo A
) | python strain_analyzer.py 4TMN.cif 8TLN.cif E

echo.
echo ────────────────────────────────────────────────────────────────
echo   TEST 4 COMPLETE. Check: Chain E parsed correctly?
echo   Check: Ligand identified? Mode classified?
echo ────────────────────────────────────────────────────────────────
echo.
echo.

REM ═══════════════════════════════════════════════════════════════
REM  SUMMARY
REM ═══════════════════════════════════════════════════════════════
echo ================================================================
echo   ALL 4 TESTS COMPLETE
echo ================================================================
echo.
echo   Test 1 (ER-alpha):     Tamoxifen, skip all options
echo   Test 2 (HIV-PR):       Indinavir, option B everywhere
echo   Test 3 (ERK):          Phosphorylation, last options
echo   Test 4 (Thermolysin):  Chain E, option A at start/end
echo.
echo   Review output above for:
echo     - No NaN values in any ranking
echo     - All phi/psi values 0-180 degrees (circular-wrapped)
echo     - Terminal/invalid residues excluded from common count
echo     - Correct chain parsing (Test 4)
echo     - Consistent direction with atlas v15
echo.
echo   Expected results:
echo     Test 1: dGini ~ +0.11   STIFFENED  75%% distal
echo     Test 2: dGini ~ +0.034  STIFFENED  75%% distal
echo     Test 3: dGini ~ +0.028  STIFFENED  95%% distal
echo     Test 4: dGini ~ +0.020  STIFFENED  100%% distal
echo.
pause
