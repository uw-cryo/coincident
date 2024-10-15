@echo off

REM Check the shell being used and set paths accordingly
IF "%SHELL%"=="" (
    echo "Running activation for Windows Command Prompt..."
    set PATH=%PATH%;C:\Program Files\Git\bin  REM Adjust to your Git Bash path
) ELSE (
    echo "Running activation for other shells..."
    set PATH=%PATH%;C:\Path\To\Other\Shells\bin  REM Adjust for other shells if necessary
)

REM Common commands for all shells
echo Activating the project environment...
REM Any additional setup can go here, e.g., setting virtualenv paths, etc.

REM Example of activating a VM:
REM call C:\path\to\your\venv\Scripts\activate