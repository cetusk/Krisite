@echo off
setlocal EnableDelayedExpansion

REM ============================================================
REM  Claude Code / Docker Sandbox launcher
REM
REM  Usage:
REM    claude_sandbox            start (create if not exists)
REM    claude_sandbox continue   continue the last session  (alias: c)
REM    claude_sandbox resume     pick a session to resume   (alias: r)
REM    claude_sandbox bypass     start with bypass permissions
REM    claude_sandbox ls         list sandboxes
REM    claude_sandbox stop       stop this sandbox
REM    claude_sandbox rm         remove this sandbox
REM
REM  NOTE: Keep this file ASCII-only. CMD tracks the batch file by
REM        byte offset, so multi-byte characters can corrupt parsing.
REM ============================================================

set "WORKSPACE=C:\Users\owner\Desktop\Apps\DockerSandboxWs\Krisite"
set "SANDBOX_NAME=claude-krisite"

REM --- Detect CLI: prefer new "sbx", fall back to legacy "docker sandbox" ---
set "NEWCLI=0"
where sbx >nul 2>&1
if %errorlevel%==0 set "NEWCLI=1"

if not exist "%WORKSPACE%" (
    echo [ERROR] Workspace not found: %WORKSPACE%
    exit /b 1
)

set "MODE=%~1"

REM --- Management subcommands ---
if /i "%MODE%"=="ls"    goto :do_ls
if /i "%MODE%"=="stop"  goto :do_stop
if /i "%MODE%"=="rm"    goto :do_rm

REM --- Build agent args (passed to Claude Code after "--") ---
set "AGENT_ARGS="
if /i "%MODE%"==""         goto :args_ok
if /i "%MODE%"=="continue" set "AGENT_ARGS=-- -c" & goto :args_ok
if /i "%MODE%"=="c"        set "AGENT_ARGS=-- -c" & goto :args_ok
if /i "%MODE%"=="resume"   set "AGENT_ARGS=-- --resume" & goto :args_ok
if /i "%MODE%"=="r"        set "AGENT_ARGS=-- --resume" & goto :args_ok
if /i "%MODE%"=="bypass"   set "AGENT_ARGS=-- --dangerously-skip-permissions" & goto :args_ok

echo [ERROR] Unknown argument: %MODE%
echo         Valid: continue ^| resume ^| bypass ^| ls ^| stop ^| rm
exit /b 1

:args_ok
echo Sandbox : %SANDBOX_NAME%
echo Workdir : %WORKSPACE%
if not "!AGENT_ARGS!"=="" echo Args    : !AGENT_ARGS!
echo.

if "%NEWCLI%"=="1" goto :run_new
goto :run_legacy


REM ============================================================
REM  New CLI (sbx): run creates or re-attaches by name
REM ============================================================
:run_new
cd /d "%WORKSPACE%"
sbx run claude --name %SANDBOX_NAME% !AGENT_ARGS!
goto :end


REM ============================================================
REM  Legacy CLI (docker sandbox): create and re-attach differ
REM ============================================================
:run_legacy
docker sandbox ls 2>nul | findstr /C:"%SANDBOX_NAME%" >nul 2>&1
if %errorlevel%==0 goto :legacy_attach

echo Creating a new sandbox...
docker sandbox run --name %SANDBOX_NAME% claude "%WORKSPACE%" !AGENT_ARGS!
goto :end

:legacy_attach
echo Re-attaching to the existing sandbox...
docker sandbox run %SANDBOX_NAME% !AGENT_ARGS!
goto :end


REM ============================================================
REM  Management
REM ============================================================
:do_ls
if "%NEWCLI%"=="1" (sbx ls) else (docker sandbox ls)
goto :end

:do_stop
echo Stopping %SANDBOX_NAME% ...
if "%NEWCLI%"=="1" (sbx stop %SANDBOX_NAME%) else (docker sandbox stop %SANDBOX_NAME%)
goto :end

:do_rm
echo [WARNING] This removes the sandbox "%SANDBOX_NAME%".
echo           Files in the workspace stay on the host, but anything
echo           installed only inside the container will be lost.
choice /c YN /m "Continue"
if errorlevel 2 goto :end
if "%NEWCLI%"=="1" (sbx rm %SANDBOX_NAME%) else (docker sandbox rm %SANDBOX_NAME%)
goto :end

:end
endlocal