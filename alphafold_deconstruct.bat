@echo off
setlocal enabledelayedexpansion

:: ============================================================================
:: ALPHAFOLD GEOMETRIC DECONSTRUCTION CONTROLLER
:: Torus curvature analysis on AlphaFold predicted structures
:: 
:: Usage:
::   alphafold_deconstruct.bat                    - Show this help
::   alphafold_deconstruct.bat single P00533      - Analyze one protein
::   alphafold_deconstruct.bat batch ids.txt      - Analyze proteins from file
::   alphafold_deconstruct.bat search EGFR        - Search UniProt, then analyze
::   alphafold_deconstruct.bat human              - Full human proteome (~11 hrs)
::   alphafold_deconstruct.bat demo               - Run demo (insulin, p53, EGFR)
::   alphafold_deconstruct.bat status             - Show barcode database stats
::   alphafold_deconstruct.bat view P00533        - View barcode for a protein
::   alphafold_deconstruct.bat compare P00533 P04637  - Compare two proteins
::   alphafold_deconstruct.bat defects P00533     - Show curvature defects only
:: ============================================================================

set SCRIPT_DIR=%~dp0
set PYTHON=python
set PIPELINE=%SCRIPT_DIR%cornu_spirals_on_T2.py
set RESULTS_DIR=%SCRIPT_DIR%results
set BARCODE_DIR=%RESULTS_DIR%\barcodes
set DB_FILE=%RESULTS_DIR%\alphafold_barcode_db.json

:: Check python
where %PYTHON% >nul 2>&1
if errorlevel 1 (
    echo [ERROR] Python not found. Install Python 3.8+ and add to PATH.
    exit /b 1
)

:: Check pipeline exists
if not exist "%PIPELINE%" (
    echo [ERROR] Pipeline not found: %PIPELINE%
    echo         Place cornu_spirals_on_T2.py in the same directory as this .bat
    exit /b 1
)

:: Route commands
if "%1"=="" goto :help
if /i "%1"=="help" goto :help
if /i "%1"=="single" goto :single
if /i "%1"=="batch" goto :batch
if /i "%1"=="search" goto :search
if /i "%1"=="human" goto :human
if /i "%1"=="demo" goto :demo
if /i "%1"=="status" goto :status
if /i "%1"=="view" goto :view
if /i "%1"=="compare" goto :compare
if /i "%1"=="defects" goto :defects
echo [ERROR] Unknown command: %1
goto :help


:: ============================================================================
:help
:: ============================================================================
echo.
echo  ================================================================
echo   ALPHAFOLD GEOMETRIC DECONSTRUCTION
echo   Torus curvature analysis pipeline
echo  ================================================================
echo.
echo  COMMANDS:
echo.
echo    single ^<UniProt_ID^>         Analyze a single protein
echo    batch ^<file.txt^>            Analyze proteins from ID list
echo    search ^<gene_name^>          Search UniProt by gene name
echo    human                        Full human proteome (~20k, ~11 hrs)
echo    demo                         Quick demo (3 well-known proteins)
echo.
echo    status                       Show barcode database statistics
echo    view ^<UniProt_ID^>           Display barcode for one protein
echo    compare ^<ID1^> ^<ID2^>         Compare two protein barcodes
echo    defects ^<UniProt_ID^>        Show curvature defects (Gaussian peaks)
echo.
echo  EXAMPLES:
echo.
echo    alphafold_deconstruct single P00533
echo    alphafold_deconstruct batch my_targets.txt
echo    alphafold_deconstruct demo
echo    alphafold_deconstruct view P00533
echo    alphafold_deconstruct compare P00533 P04637
echo.
echo  ID LIST FORMAT (one per line, # for comments):
echo    # Kinases
echo    P00533
echo    P04637
echo    P01308
echo.
goto :eof


:: ============================================================================
:single
:: ============================================================================
if "%2"=="" (
    echo [ERROR] Provide a UniProt ID.  Example: alphafold_deconstruct single P00533
    exit /b 1
)
echo.
echo  Analyzing: %2
echo  ----------------------------------------
%PYTHON% "%PIPELINE%" alphafold %2
if errorlevel 1 (
    echo [ERROR] Analysis failed for %2
    exit /b 1
)
echo.
echo  Done. Barcode saved to: %BARCODE_DIR%\%2_barcode.json
goto :eof


:: ============================================================================
:batch
:: ============================================================================
if "%2"=="" (
    echo [ERROR] Provide an ID list file.  Example: alphafold_deconstruct batch targets.txt
    exit /b 1
)
if not exist "%2" (
    echo [ERROR] File not found: %2
    exit /b 1
)
:: Count non-comment lines
set /a COUNT=0
for /f "usebackq eol=# tokens=*" %%a in ("%2") do set /a COUNT+=1
echo.
echo  Batch processing: %2
echo  Proteins: %COUNT%
echo  Estimated time: ~%COUNT% x 2 seconds
echo  ----------------------------------------
echo.
set START_TIME=%TIME%
%PYTHON% "%PIPELINE%" alphafold @%2
echo.
echo  Started:  %START_TIME%
echo  Finished: %TIME%
echo  Barcodes: %BARCODE_DIR%\
echo  Database: %DB_FILE%
goto :eof


:: ============================================================================
:search
:: ============================================================================
if "%2"=="" (
    echo [ERROR] Provide a gene name.  Example: alphafold_deconstruct search EGFR
    exit /b 1
)
echo.
echo  Searching UniProt for: %2
echo  ----------------------------------------
%PYTHON% -c "import requests,json; r=requests.get('https://rest.uniprot.org/uniprotkb/search',params={'query':'gene:%2 AND organism_id:9606 AND reviewed:true','format':'json','size':'5'},timeout=15); data=r.json(); entries=data.get('results',[]); [print(f'  {e[\"primaryAccession\"]:10s}  {e.get(\"proteinDescription\",{}).get(\"recommendedName\",{}).get(\"fullName\",{}).get(\"value\",\"?\")}  ({e.get(\"sequence\",{}).get(\"length\",\"?\")} aa)') for e in entries] if entries else print('  No results found.')"
if errorlevel 1 (
    echo  [Note] Search requires internet. You can manually find IDs at uniprot.org
    exit /b 0
)
echo.
echo  To analyze, run:  alphafold_deconstruct single ^<UniProt_ID^>
goto :eof


:: ============================================================================
:human
:: ============================================================================
echo.
echo  ================================================================
echo   FULL HUMAN PROTEOME ANALYSIS
echo  ================================================================
echo.
echo  This will:
echo    1. Download ~20,000 AlphaFold structures (~50 GB)
echo    2. Compute torus curvature for each protein
echo    3. Generate geometric barcodes
echo    4. Build searchable barcode database
echo.
echo  Estimated time: ~11 hours
echo  Estimated disk:  ~55 GB (structures + barcodes)
echo.
set /p CONFIRM="  Proceed? (y/n): "
if /i not "%CONFIRM%"=="y" (
    echo  Aborted.
    goto :eof
)
echo.
echo  Step 1: Downloading human proteome ID list...
%PYTHON% -c "import requests; r=requests.get('https://rest.uniprot.org/uniprotkb/stream',params={'query':'organism_id:9606 AND reviewed:true','format':'list'},timeout=60,stream=True); f=open('human_proteome_ids.txt','w'); [f.write(line.decode()+'\n') if isinstance(line,bytes) else f.write(line+'\n') for line in r.iter_lines()]; f.close(); print(f'  Downloaded {sum(1 for _ in open(\"human_proteome_ids.txt\"))} IDs')"
if errorlevel 1 (
    echo  [ERROR] Failed to download proteome list. Check internet connection.
    exit /b 1
)
echo.
echo  Step 2: Running analysis...
set START_TIME=%TIME%
%PYTHON% "%PIPELINE%" alphafold @human_proteome_ids.txt
echo.
echo  ================================================================
echo   HUMAN PROTEOME ANALYSIS COMPLETE
echo   Started:  %START_TIME%
echo   Finished: %TIME%
echo   Database: %DB_FILE%
echo  ================================================================
goto :eof


:: ============================================================================
:demo
:: ============================================================================
echo.
echo  DEMO: Analyzing 3 well-known proteins
echo  ----------------------------------------
echo    P01308 - Insulin
echo    P04637 - Tumor protein p53
echo    P00533 - Epidermal growth factor receptor (EGFR)
echo  ----------------------------------------
echo.
%PYTHON% "%PIPELINE%" alphafold P01308 P04637 P00533
goto :eof


:: ============================================================================
:status
:: ============================================================================
echo.
echo  BARCODE DATABASE STATUS
echo  ----------------------------------------
if not exist "%BARCODE_DIR%" (
    echo  No barcodes found. Run 'alphafold_deconstruct demo' to get started.
    goto :eof
)
:: Count barcode files
set /a BC_COUNT=0
for %%f in ("%BARCODE_DIR%\*_barcode.json") do set /a BC_COUNT+=1
echo  Barcodes:  %BC_COUNT%
echo  Location:  %BARCODE_DIR%
if exist "%DB_FILE%" (
    echo  Database:  %DB_FILE%
)
echo.
if %BC_COUNT% gtr 0 (
    echo  Recent barcodes:
    for /f "delims=" %%f in ('dir /b /o-d "%BARCODE_DIR%\*_barcode.json" 2^>nul ^| head 2^>nul') do (
        echo    %%f
    )
)
goto :eof


:: ============================================================================
:view
:: ============================================================================
if "%2"=="" (
    echo [ERROR] Provide a UniProt ID.  Example: alphafold_deconstruct view P00533
    exit /b 1
)
set BC_FILE=%BARCODE_DIR%\%2_barcode.json
if not exist "%BC_FILE%" (
    echo  [ERROR] No barcode found for %2
    echo  Run:  alphafold_deconstruct single %2
    exit /b 1
)
echo.
echo  GEOMETRIC BARCODE: %2
echo  ========================================
%PYTHON% -c "import json; b=json.load(open(r'%BC_FILE%')); print(f'  UniProt:     {b.get(\"uniprot_id\",\"?\")}'); print(f'  Length:      {b.get(\"length\",0)} residues'); print(f'  SS:          H={b[\"ss_composition\"][\"H\"]} E={b[\"ss_composition\"][\"E\"]} C={b[\"ss_composition\"][\"C\"]}'); print(f'  p-winding:   {b[\"p\"]:+.3f}'); print(f'  q-winding:   {b[\"q\"]:+.3f}'); print(f'  |Q|:         {b[\"Q_magnitude\"]:.3f}'); print(f'  Segments:    {b[\"n_segments\"]}'); print(f'  Geodesic:    {b[\"n_geodesic\"]} ({b[\"frac_geodesic\"]:.0%%})'); print(f'  Defects:     {b[\"n_defect\"]} ({b[\"frac_defect\"]:.0%%})'); print(); print('  SEGMENT BARCODES:'); print('  {:>4s} {:>3s} {:>15s} {:>5s} {:>7s} {:>5s}'.format('Idx','SS','Model','Len','Kappa','R2')); print('  '+'-'*50); [print(f'  {i+1:4d} {s[\"ss_type\"]:>3s} {s[\"model_class\"]:>15s} {s[\"length\"]:5d} {s[\"mean_kappa\"]:+7.3f} {s[\"R2\"]:5.3f}') for i,s in enumerate(b.get('segments',[]))]"
goto :eof


:: ============================================================================
:compare
:: ============================================================================
if "%2"=="" (
    echo [ERROR] Provide two UniProt IDs.  Example: alphafold_deconstruct compare P00533 P04637
    exit /b 1
)
if "%3"=="" (
    echo [ERROR] Provide two UniProt IDs.  Example: alphafold_deconstruct compare P00533 P04637
    exit /b 1
)
set BC1=%BARCODE_DIR%\%2_barcode.json
set BC2=%BARCODE_DIR%\%3_barcode.json
if not exist "%BC1%" (
    echo  [ERROR] No barcode for %2. Run: alphafold_deconstruct single %2
    exit /b 1
)
if not exist "%BC2%" (
    echo  [ERROR] No barcode for %3. Run: alphafold_deconstruct single %3
    exit /b 1
)
echo.
echo  GEOMETRIC COMPARISON: %2 vs %3
echo  ========================================
%PYTHON% -c "import json,math; a=json.load(open(r'%BC1%')); b=json.load(open(r'%BC2%')); dp=a['p']-b['p']; dq=a['q']-b['q']; dQ=math.sqrt(dp**2+dq**2); print(f'  {'':20s} {'%2':>12s} {'%3':>12s} {'Delta':>10s}'); print('  '+'-'*58); print(f'  {\"Length\":20s} {a[\"length\"]:12d} {b[\"length\"]:12d} {a[\"length\"]-b[\"length\"]:+10d}'); print(f'  {\"p-winding\":20s} {a[\"p\"]:+12.3f} {b[\"p\"]:+12.3f} {dp:+10.3f}'); print(f'  {\"q-winding\":20s} {a[\"q\"]:+12.3f} {b[\"q\"]:+12.3f} {dq:+10.3f}'); print(f'  {\"|Q|\":20s} {a[\"Q_magnitude\"]:12.3f} {b[\"Q_magnitude\"]:12.3f} {a[\"Q_magnitude\"]-b[\"Q_magnitude\"]:+10.3f}'); print(f'  {\"Segments\":20s} {a[\"n_segments\"]:12d} {b[\"n_segments\"]:12d} {a[\"n_segments\"]-b[\"n_segments\"]:+10d}'); print(f'  {\"Geodesic frac\":20s} {a[\"frac_geodesic\"]:11.0%%} {b[\"frac_geodesic\"]:11.0%%}'); print(f'  {\"Defect frac\":20s} {a[\"frac_defect\"]:11.0%%} {b[\"frac_defect\"]:11.0%%}'); print(); print(f'  Winding distance (Euclidean in p,q): {dQ:.3f}'); fg_diff=abs(a['frac_geodesic']-b['frac_geodesic']); print(f'  Geodesic similarity: {1-fg_diff:.0%%}')"
goto :eof


:: ============================================================================
:defects
:: ============================================================================
if "%2"=="" (
    echo [ERROR] Provide a UniProt ID.  Example: alphafold_deconstruct defects P00533
    exit /b 1
)
set BC_FILE=%BARCODE_DIR%\%2_barcode.json
if not exist "%BC_FILE%" (
    echo  [ERROR] No barcode for %2. Run: alphafold_deconstruct single %2
    exit /b 1
)
echo.
echo  CURVATURE DEFECTS: %2
echo  ========================================
%PYTHON% -c "import json; b=json.load(open(r'%BC_FILE%')); defects=[s for s in b.get('segments',[]) if s['model_class'] in ('gauss_peak','quadratic','damped_oscillation')]; print(f'  Total segments: {b[\"n_segments\"]}'); print(f'  Defects found:  {len(defects)}'); print(); [print(f'  [{i+1}] {s[\"ss_type\"]} seg, residues {s[\"start_idx\"]}-{s[\"end_idx\"]} ({s[\"length\"]} res)') or print(f'      Model: {s[\"model_class\"]}  kappa={s[\"mean_kappa\"]:+.3f}  R2={s[\"R2\"]:.3f}') or (print(f'      Peak: A={s[\"params\"][\"A\"]:.2f}  sigma={s[\"params\"][\"sigma\"]:.3f}  peak_res={s[\"params\"][\"peak_residue\"]}') if s['model_class']=='gauss_peak' and 'A' in s.get('params',{}) else None) or print() for i,s in enumerate(defects)] if defects else print('  No curvature defects detected.'); print(f'  Geodesic fraction: {b[\"frac_geodesic\"]:.0%%} (higher = smoother backbone on T2)')"
goto :eof
