from flask import Flask, request, jsonify, render_template
from flask_cors import CORS
import os
import json
from werkzeug.utils import secure_filename
from pdb_parser import PDBParser

from interaction import ComplexPredictor

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

# Initialize complex predictor and a simple in-memory cache for sample PDBs
predictor = ComplexPredictor(cache_dir=os.path.join('uploads', '.cache'))
_sample_cache = {}

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

@app.route('/api/interaction_available')
def interaction_available():
    """Check if interaction prediction is available"""
    return jsonify({'available': True})

@app.route('/api/predict_interaction', methods=['POST'])
def predict_interaction():
    """Predict protein-protein complex poses between two input PDBs or sample names.

    Request JSON body options:
    - pdb_a: PDB text for protein A
    - pdb_b: PDB text for protein B
    - sample_a: optional sample name to load from sample_data instead of pdb_a
    - sample_b: optional sample name to load from sample_data instead of pdb_b
    - num_poses: int, default 5
    - refine: bool, default False (requires OpenMM if True)
    """
    
    try:
        payload = request.get_json(force=True, silent=True) or {}
        sample_a = payload.get('sample_a')
        sample_b = payload.get('sample_b')
        pdb_a = payload.get('pdb_a')
        pdb_b = payload.get('pdb_b')
        num_poses = int(payload.get('num_poses', 5))
        refine = bool(payload.get('refine', False))

        if sample_a and not pdb_a:
            pdb_a = _load_sample_pdb_text(sample_a)
        if sample_b and not pdb_b:
            pdb_b = _load_sample_pdb_text(sample_b)

        if not pdb_a or not pdb_b:
            return jsonify({'error': 'Provide both pdb_a and pdb_b (or sample_a/sample_b).'}), 400

        poses = predictor.predict(pdb_a, pdb_b, num_poses=num_poses, refine=refine)
        response_poses = []
        for pose in poses:
            response_poses.append({
                'pose_id': pose.pose_id,
                'pdb': pose.pdb,
                'scores': {
                    'confidence': pose.scores.confidence,
                    'binding_energy': pose.scores.binding_energy,
                    'buried_surface_area': pose.scores.buried_surface_area,
                    'contact_count': pose.scores.contact_count,
                    'hydrogen_bonds': pose.scores.hydrogen_bonds,
                    'salt_bridges': pose.scores.salt_bridges,
                    'shape_complementarity': pose.scores.shape_complementarity,
                    'final_score': pose.scores.final_score,
                }
            })

        return jsonify({
            'success': True,
            'num_poses': len(response_poses),
            'poses': response_poses
        })
    except Exception as e:
        return jsonify({'error': f'Prediction failed: {str(e)}'}), 500


def _load_sample_pdb_text(sample_name: str) -> str:
    if sample_name in _sample_cache:
        return _sample_cache[sample_name]
    sample_path = os.path.join('../sample_data', f'{sample_name}.pdb')
    if not os.path.exists(sample_path):
        raise FileNotFoundError(f'Sample {sample_name} not found')
    with open(sample_path, 'r') as f:
        text = f.read()
    _sample_cache[sample_name] = text
    return text

if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0', port=8000)
