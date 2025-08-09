# Protein Structure Viewer

A simplified AlphaFold-style system for predicting and visualizing protein structures from amino acid sequences with interactive 3D visualization.

## Features

- **Interactive 3D Protein Visualization** - Powered by 3Dmol.js for smooth, responsive 3D rendering
- **Multiple Visualization Styles** - Cartoon, stick, spacefill, and backbone representations
- **Drag & Drop File Upload** - Easy PDB file upload with validation
- **Sample Proteins** - Pre-loaded sample structures for quick testing
- **Protein Information Panel** - Detailed structure statistics and metadata
- **Responsive Design** - Modern UI that works on desktop and mobile
- **Real-time Controls** - Interactive rotation, zoom, and color scheme controls

## Tech Stack

- **Backend**: Python, Flask, BioPython
- **Frontend**: HTML/CSS/JavaScript, React (CDN), 3Dmol.js
- **Visualization**: 3Dmol.js for molecular rendering
- **Styling**: Modern CSS with gradients and blur effects

## Installation

### Prerequisites

- Python 3.8 or higher
- pip package manager

### Setup

1. **Clone or download the project**
   ```bash
   cd protein-structure-viewer
   ```

2. **Create and activate virtual environment**
   ```bash
   python -m venv venv
   
   # On macOS/Linux:
   source venv/bin/activate
   
   # On Windows:
   venv\Scripts\activate
   ```

3. **Install dependencies**
   ```bash
   pip install -r requirements.txt
   ```

4. **Run the application**
   ```bash
   python backend/app.py
   ```

5. **Open your browser**
   Navigate to `http://localhost:5000`

## Usage

### Upload PDB Files

1. **Drag & Drop**: Simply drag a PDB file onto the upload area
2. **Browse**: Click the upload area to browse and select a file
3. **Supported Formats**: `.pdb` and `.ent` files

### Sample Proteins

The application comes with several sample proteins for testing:

- **Hemoglobin** - Oxygen-carrying protein from red blood cells
- **Insulin** - Hormone that regulates blood sugar
- **Lysozyme** - Antimicrobial enzyme found in tears and saliva

Click on any sample to load it immediately.

### Visualization Controls

#### Representation Styles

- **Cartoon** - Smooth ribbon representation showing secondary structure
- **Stick** - Shows individual bonds between atoms
- **Spacefill** - Van der Waals sphere representation
- **Backbone** - Shows only the protein backbone

#### Color Schemes

- **Spectrum** - Rainbow coloring from N to C terminus
- **By Chain** - Different colors for each protein chain
- **By Residue** - Colors based on residue type
- **By Element** - Standard element-based coloring (CPK)

#### Interactive Controls

- **Rotate**: Click and drag to rotate the structure
- **Zoom**: Mouse wheel or pinch gestures to zoom in/out
- **Reset View**: Button to return to the original view
- **Fullscreen**: Expand the viewer to fullscreen mode

## API Endpoints

### Upload PDB File
```
POST /api/upload
Content-Type: multipart/form-data
Body: file (PDB file)
```

### Load Sample Protein
```
GET /api/sample/<sample_name>
```

### List Available Samples
```
GET /api/samples
```

## File Structure

```
protein-structure-viewer/
├── backend/
│   ├── app.py              # Flask server
│   └── pdb_parser.py       # PDB file parsing logic
├── frontend/               # (Reserved for future React build)
├── templates/
│   └── index.html          # Main HTML template
├── static/
│   ├── css/
│   │   └── styles.css      # Application styling
│   └── js/
│       └── app.js          # React application
├── sample_data/            # Sample PDB files
│   ├── hemoglobin.pdb
│   ├── insulin.pdb
│   └── lysozyme.pdb
├── venv/                   # Python virtual environment
├── requirements.txt        # Python dependencies
└── README.md              # This file
```

## Key Components

### Backend (Python/Flask)

- **PDB Parser**: Comprehensive parser that extracts atoms, bonds, chains, and metadata
- **File Upload**: Secure file handling with validation
- **Sample Management**: API for loading pre-configured sample structures
- **CORS Support**: Cross-origin resource sharing for frontend integration

### Frontend (React/3Dmol.js)

- **Protein Viewer**: 3D visualization component using 3Dmol.js
- **File Upload**: Drag & drop interface with progress feedback
- **Controls Panel**: Interactive controls for visualization customization
- **Info Panel**: Displays protein statistics and metadata
- **Responsive Layout**: Mobile-friendly design with sidebar and main viewer

## Development

### Adding New Sample Proteins

1. Place PDB files in the `sample_data/` directory
2. The application automatically discovers and lists them
3. File names become the sample identifiers

### Customizing Visualization

The visualization system supports extensive customization through 3Dmol.js:

- Add new representation styles in the `applyVisualizationStyle` function
- Modify color schemes in the visualization controls
- Extend the PDB parser for additional metadata extraction

### Extending the Backend

The Flask backend is modular and can be extended with:

- Additional file format support (mmCIF, MOL2, etc.)
- Protein structure prediction capabilities
- Database integration for PDB ID lookup
- Advanced analysis features

## Browser Compatibility

- Chrome 60+
- Firefox 55+
- Safari 12+
- Edge 79+

WebGL support is required for 3D visualization.

## Troubleshooting

### Common Issues

1. **File Upload Fails**
   - Ensure file is a valid PDB format
   - Check file size (max 16MB)
   - Verify file permissions

2. **3D Viewer Not Loading**
   - Check browser WebGL support
   - Ensure JavaScript is enabled
   - Try refreshing the page

3. **Sample Proteins Not Showing**
   - Verify `sample_data/` directory exists
   - Check file permissions
   - Restart the Flask server

### Performance Tips

- Large protein structures (>10,000 atoms) may load slowly
- Use cartoon representation for better performance with large structures
- Close other browser tabs for optimal performance

## Contributing

This is a demonstration project showcasing modern web-based protein visualization. Feel free to fork and extend for your own use cases.

## License

This project is for educational and research purposes. Protein structure data should be properly attributed to their original sources.

## Acknowledgments

- **3Dmol.js** - WebGL-accelerated molecular visualization
- **BioPython** - Python tools for computational biology
- **PDB Database** - Protein structure data repository
- **Flask** - Python web framework
- **React** - Frontend JavaScript library
