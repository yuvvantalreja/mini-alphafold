#!/bin/bash

# Protein Structure Viewer Startup Script
echo "🧬 Starting Protein Structure Viewer..."

# Check if virtual environment exists
if [ ! -d "venv" ]; then
    echo "❌ Virtual environment not found. Please run setup first:"
    echo "   python -m venv venv"
    echo "   source venv/bin/activate"
    echo "   pip install -r requirements.txt"
    exit 1
fi

# Activate virtual environment
echo "🔧 Activating virtual environment..."
source venv/bin/activate

# Check if dependencies are installed
if ! python -c "import flask" 2>/dev/null; then
    echo "❌ Dependencies not installed. Installing now..."
    pip install -r requirements.txt
fi

# Create uploads directory if it doesn't exist
mkdir -p backend/uploads

# Start the server
echo "🚀 Starting Flask server..."
echo "📱 Open http://localhost:5000 in your browser"
echo "💡 Press Ctrl+C to stop the server"
echo ""

cd backend && python app.py
