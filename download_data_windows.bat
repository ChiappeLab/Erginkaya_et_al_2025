@echo off
setlocal enabledelayedexpansion


:: Dropbox folder URL (fixed to use dl=1 for direct download)
set "FOLDER_1=https://zenodo.org/records/14967806/files/data.zip?download=1"

:: Define the ZIP filename (in script directory)
set "ZIP_FILE=data.zip"

:: Download folder ZIP
echo Downloading folder ZIP...
curl -L -o "%ZIP_FILE%" "%FOLDER_1%"

:: Verify the ZIP was downloaded
if not exist "%ZIP_FILE%" (
    echo Error: The ZIP file was not downloaded. Please check the link or internet connection.
    pause
    exit /b 1
)


:: Extract ZIP file to script directory (ensure process fully completes)
echo Extracting ZIP file to script directory ...
start /wait powershell -Command "Expand-Archive -Path '%ZIP_FILE%' -DestinationPath '.' -Force"

:: Wait a bit longer to ensure extraction is fully complete and file lock is released
timeout /t 5 /nobreak >nul


:: Verify extraction was successful (check for at least one known folder)
if exist "data" (
    echo Extraction successful.
    del /f /q "%ZIP_FILE%" 2>nul
    if exist "%ZIP_FILE%" (
        echo Warning: ZIP file could not be removed. Please delete it manually.
        pause
        exit /b 1
    )
    echo Download and extraction complete. Files are saved in the script directory.
    pause
    exit /b 0
) else (
    echo Error: Extraction failed or expected folders not found.
    pause
    exit /b 1
)
