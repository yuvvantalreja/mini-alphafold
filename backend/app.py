

import matplotlib
matplotlib.use('Agg')

from flask import Flask, request, jsonify, render_template, redirect, url_for
from flask_cors import CORS
import os
import json
from werkzeug.utils import secure_filename
from pdb_parser import PDBParser
import chatbot
import pdb_script


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

@app.route('/api/sample/<sample_path>')
def load_sample(sample_path):
    """Load a sample protein structure"""
    try:
        # sample_path = os.path.join('../sample_data', f'{sample_name}.pdb')
        if not os.path.exists(sample_path):
            return jsonify({'error': 'Sample not found'}), 404
        
        parser = PDBParser()
        structure_data = parser.parse_pdb_file(sample_path)
        
        return jsonify({
            'success': True,
            'filename': f'{sample_path}',
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

@app.route('/api/chat', methods=['POST'])
def chat():
    data = request.get_json()
    user_message = data.get('message', '')
    return jsonify({'response': chatbot.ChatBot(user_message)})
    #return jsonify({'response': "Wasgood"})

@app.route('/api/sequence', methods=['POST'])
def receive_sequence():
    data = request.get_json()
    sequence = data.get('sequence', '')
    print(f"Received sequence: {sequence}")
    jobname = pdb_script.GeneratePDB(sequence)

    if jobname == -1:
        return jsonify({'error': 'Invalid sequence provided'}), 400

    # Search for pdb files in the jobname directory
    pdb_files = [f for f in os.listdir(jobname) if f.endswith('.pdb')]
    if not pdb_files:
        return jsonify({'error': 'No PDB files generated for the provided sequence'}), 404

    finalpath = os.path.join(jobname, pdb_files[0])

    # Parse the PDB file and return structure data directly
    parser = PDBParser()
    structure_data = parser.parse_pdb_file(finalpath)
    return jsonify({
        'success': True,
        'filename': finalpath,
        'structure': structure_data
    })
    

if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0', port=8000)




