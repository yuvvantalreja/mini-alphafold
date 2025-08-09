from flask import Flask, request, jsonify, render_template
from flask_cors import CORS
import os
import json
from werkzeug.utils import secure_filename
from pdb_parser import PDBParser

app = Flask(__name__, template_folder='../templates', static_folder='../static')
CORS(app)

# Configuration
UPLOAD_FOLDER = 'uploads'
ALLOWED_EXTENSIONS = {'pdb', 'ent'}
MAX_FILE_SIZE = 16 * 1024 * 1024  # 16MB

app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
app.config['MAX_CONTENT_LENGTH'] = MAX_FILE_SIZE

# Ensure upload directory exists
os.makedirs(UPLOAD_FOLDER, exist_ok=True)

def allowed_file(filename):
    return '.' in filename and \
           filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS

@app.route('/')
def index():
    """Serve the main application page"""
    return render_template('index.html')

@app.route('/api/upload', methods=['POST'])
def upload_file():
    """Handle PDB file uploads"""
    try:
        if 'file' not in request.files:
            return jsonify({'error': 'No file provided'}), 400
        
        file = request.files['file']
        if file.filename == '':
            return jsonify({'error': 'No file selected'}), 400
        
        if file and allowed_file(file.filename):
            filename = secure_filename(file.filename)
            filepath = os.path.join(app.config['UPLOAD_FOLDER'], filename)
            file.save(filepath)
            
            # Parse the PDB file
            parser = PDBParser()
            try:
                structure_data = parser.parse_pdb_file(filepath)
                return jsonify({
                    'success': True,
                    'filename': filename,
                    'structure': structure_data
                })
            except Exception as e:
                return jsonify({'error': f'Failed to parse PDB file: {str(e)}'}), 400
            finally:
                # Clean up uploaded file
                if os.path.exists(filepath):
                    os.remove(filepath)
        
        return jsonify({'error': 'Invalid file type. Please upload a PDB file.'}), 400
    
    except Exception as e:
        return jsonify({'error': f'Upload failed: {str(e)}'}), 500

@app.route('/api/sample/<sample_name>')
def load_sample(sample_name):
    """Load a sample protein structure"""
    try:
        sample_path = os.path.join('../sample_data', f'{sample_name}.pdb')
        if not os.path.exists(sample_path):
            return jsonify({'error': 'Sample not found'}), 404
        
        parser = PDBParser()
        structure_data = parser.parse_pdb_file(sample_path)
        
        return jsonify({
            'success': True,
            'filename': f'{sample_name}.pdb',
            'structure': structure_data
        })
    
    except Exception as e:
        return jsonify({'error': f'Failed to load sample: {str(e)}'}), 500

@app.route('/api/samples')
def list_samples():
    """List available sample protein structures"""
    try:
        sample_dir = '../sample_data'
        if not os.path.exists(sample_dir):
            return jsonify({'samples': []})
        
        samples = []
        for filename in os.listdir(sample_dir):
            if filename.endswith('.pdb'):
                sample_name = filename[:-4]  # Remove .pdb extension
                samples.append({
                    'name': sample_name,
                    'filename': filename,
                    'title': sample_name.replace('_', ' ').title()
                })
        
        return jsonify({'samples': samples})
    
    except Exception as e:
        return jsonify({'error': f'Failed to list samples: {str(e)}'}), 500

if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0', port=8000)
