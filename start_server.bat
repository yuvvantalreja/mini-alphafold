@echo off

REM Protein Structure Viewer Startup Script for Windows
echo 🧬 Starting Protein Structure Viewer...

REM Check if virtual environment exists
if not exist "venv" (
    echo ❌ Virtual environment not found. Please run setup first:
    echo    python -m venv venv
    echo    venv\Scripts\activate
    echo    pip install -r requirements.txt
    pause
    exit /b 1
)

REM Activate virtual environment
echo 🔧 Activating virtual environment...
call venv\Scripts\activate

REM Check if dependencies are installed
python -c "import flask" 2>nul
if errorlevel 1 (
    echo ❌ Dependencies not installed. Installing now...
    pip install -r requirements.txt
)

REM Create uploads directory if it doesn't exist
if not exist "backend\uploads" mkdir backend\uploads

REM Start the server
echo 🚀 Starting Flask server...
echo 📱 Open http://localhost:5000 in your browser
echo 💡 Press Ctrl+C to stop the server
echo.

cd backend && python app.py

pause
