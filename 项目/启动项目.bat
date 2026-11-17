@echo off
chcp 65001 >nul
title BioNova - Bioinformatics Analysis Platform
cd /d "%~dp0"
echo ============================================
echo   BioNova - Bioinformatics Analysis Platform
echo ============================================
where python >nul 2>nul
if errorlevel 1 (
    echo [INFO] python not found in PATH, trying "py" launcher...
    py -X utf8 server.py
) else (
    python -X utf8 server.py
)
if errorlevel 1 pause
