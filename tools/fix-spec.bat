@echo off
setlocal

REM ============================================================
REM  Krisite - fix docs\SPEC-phase0.md editorial leftovers
REM
REM  Usage (run from anywhere):
REM    tools\fix-spec.bat            show diff and apply (keeps .bak)
REM    tools\fix-spec.bat --check    show diff only, do not write
REM    tools\fix-spec.bat --no-backup
REM
REM  NOTE: Keep this file ASCII-only. CMD tracks the batch file by
REM        byte offset, so multi-byte characters can corrupt parsing.
REM        All Japanese text lives in fix-spec.py instead.
REM ============================================================

set "SCRIPT_DIR=%~dp0"
set "PYSCRIPT=%SCRIPT_DIR%fix-spec.py"

if not exist "%PYSCRIPT%" (
    echo [ERROR] Not found: %PYSCRIPT%
    exit /b 1
)

REM --- Find a Python interpreter -----------------------------
set "PY="
for %%P in (py python python3) do (
    if not defined PY (
        where %%P >nul 2>nul && set "PY=%%P"
    )
)

if not defined PY (
    echo [ERROR] No Python interpreter found on PATH ^(tried: py, python, python3^).
    echo.
    echo         Option A - install Python from https://www.python.org/downloads/
    echo                    and re-run this script.
    echo.
    echo         Option B - run it inside the sandbox, which already has python3:
    echo                    claude_sandbox
    echo                    then, at the Claude Code prompt, type:
    echo                        ! python3 tools/fix-spec.py
    echo.
    echo         Option C - edit docs\SPEC-phase0.md by hand. The four spots are
    echo                    listed in the docstring at the top of tools\fix-spec.py.
    exit /b 1
)

REM --- Make the console UTF-8 so Japanese output is readable --
set "OLDCP="
for /f "tokens=2 delims=:" %%C in ('chcp') do set "OLDCP=%%C"
chcp 65001 >nul 2>nul
set "PYTHONIOENCODING=utf-8"
set "PYTHONUTF8=1"

"%PY%" "%PYSCRIPT%" %*
set "RC=%ERRORLEVEL%"

if defined OLDCP chcp %OLDCP% >nul 2>nul

if not "%RC%"=="0" (
    echo.
    echo [FAILED] fix-spec.py exited with code %RC%. Nothing was written.
)

exit /b %RC%
