
try:
    import matplotlib
    matplotlib.use('Agg')
except Exception:
    # Optional: plotting not required for core API
    pass

from flask import Flask, request, jsonify, render_template, redirect, url_for
from flask_cors import CORS
import os
import re
import json
import requests
import time
import uuid
import hashlib
import logging
from werkzeug.utils import secure_filename
from pdb_parser import PDBParser
import chatbot
try:
    import pdb_script
except Exception:
    pdb_script = None
from similarity import (
    compute_properties,
    is_valid_protein_sequence,
    smith_waterman,
    percent_identity,
)
# Import docking functionality from parent directory
import sys
import os
parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, parent_dir)
try:
    from predict_and_rank import run_pipeline as run_docking_pipeline
    print("Successfully imported docking functionality")
except Exception as e:
    print(f"Warning: Could not import docking functionality: {e}")
    run_docking_pipeline = None


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

# --------------------------- Logging setup ---------------------------
logger = logging.getLogger("similarity")
if not logger.handlers:
    logger.setLevel(logging.INFO)
    _h = logging.StreamHandler()
    _h.setFormatter(logging.Formatter('%(asctime)s %(levelname)s %(message)s'))
    logger.addHandler(_h)

def _seq_preview(seq: str) -> str:
    s = ''.join([c for c in (seq or '') if c.isalpha()])
    if len(s) <= 14:
        return s
    return f"{s[:7]}...{s[-4:]}"

def _seq_sha1(seq: str) -> str:
    try:
        return hashlib.sha1((seq or '').encode()).hexdigest()[:12]
    except Exception:
        return ""

# No local DB initialization; RCSB is the sole similarity source now.


# --------------------------- RCSB helpers ---------------------------
RCSB_SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query?json"
# RCSB Data API expects split path: /polymer_entity/{entry_id}/{entity_id}
# and /polymer_entity_instance/{entry_id}/{asym_id}
RCSB_ENTITY_URL = "https://data.rcsb.org/rest/v1/core/polymer_entity/{entry_id}/{entity_id}"
RCSB_ENTITY_INSTANCE_URL = "https://data.rcsb.org/rest/v1/core/polymer_entity_instance/{entry_id}/{asym_id}"
_ENTITY_CACHE: dict[str, dict] = {}


def _rcsb_fetch_entity(identifier: str, request_id: str | None = None) -> dict | None:
    """Fetch polymer_entity JSON by resolving common identifier forms.
    Accepted forms:
      - polymer entity composite id: 1ABC_1  -> GET /polymer_entity/1ABC/1
      - instance id: 1ABC.A                  -> GET /polymer_entity_instance/1ABC/A, then map to entity and fetch /polymer_entity/1ABC/{id}
    Casing of entry id is case-insensitive per RCSB.
    """
    if identifier in _ENTITY_CACHE:
        return _ENTITY_CACHE[identifier]

    def _get_entity_split(entry_id: str, ent_id: str) -> dict | None:
        for e in (entry_id, entry_id.lower(), entry_id.upper()):
            try:
                t0 = time.time()
                url = RCSB_ENTITY_URL.format(entry_id=e, entity_id=ent_id)
                resp = requests.get(url, timeout=15)
                dt = round((time.time() - t0) * 1000)
                logger.info(json.dumps({'event': 'rcsb_fetch_entity', 'rid': request_id, 'entry_id': e, 'entity_id': ent_id, 'status': resp.status_code, 'dt_ms': dt}))
                if resp.status_code == 200:
                    return resp.json()
            except Exception as ex:
                logger.warning(json.dumps({'event': 'rcsb_fetch_entity.exception', 'rid': request_id, 'entry_id': e, 'entity_id': ent_id, 'error': str(ex)}))
        return None

    def _get_instance_split(entry_id: str, asym_id: str) -> dict | None:
        for e in (entry_id, entry_id.lower(), entry_id.upper()):
            for a in (asym_id, asym_id.upper()):
                try:
                    t0 = time.time()
                    url = RCSB_ENTITY_INSTANCE_URL.format(entry_id=e, asym_id=a)
                    resp = requests.get(url, timeout=15)
                    dt = round((time.time() - t0) * 1000)
                    logger.info(json.dumps({'event': 'rcsb_fetch_instance', 'rid': request_id, 'entry_id': e, 'asym_id': a, 'status': resp.status_code, 'dt_ms': dt}))
                    if resp.status_code == 200:
                        return resp.json()
                except Exception as ex:
                    logger.warning(json.dumps({'event': 'rcsb_fetch_instance.exception', 'rid': request_id, 'entry_id': e, 'asym_id': a, 'error': str(ex)}))
        return None

    # Case 1: composite polymer entity id like 1AKI_1
    if '_' in identifier and '.' not in identifier:
        entry, ent = identifier.split('_', 1)
        data = _get_entity_split(entry, ent)
        if data:
            _ENTITY_CACHE[identifier] = data
            return data

    # Case 2: instance id like 1AKI.A -> resolve to entity
    if '.' in identifier:
        entry, asym = identifier.split('.', 1)
        inst = _get_instance_split(entry, asym)
        if inst:
            cid = inst.get('rcsb_polymer_entity_instance_container_identifiers', {})
            peid = cid.get('polymer_entity_id')  # typically numeric string, e.g., "1"
            if peid:
                data2 = _get_entity_split(entry, peid)
                if data2:
                    _ENTITY_CACHE[identifier] = data2
                    return data2

    # As a last resort, try interpreting identifier as already split by accidental forms
    # e.g., "1AKI/1" passed in; handle gracefully
    if '/' in identifier:
        parts = identifier.split('/')
        if len(parts) == 2:
            data3 = _get_entity_split(parts[0], parts[1])
            if data3:
                _ENTITY_CACHE[identifier] = data3
                return data3

    logger.warning(json.dumps({'event': 'rcsb_fetch_entity.not_found', 'rid': request_id, 'id': identifier}))
    return None


def _extract_entity_metadata(j: dict) -> tuple[str, str | None, str | None, str | None]:
    # name/description
    name = (
        j.get('entity_poly', {}).get('pdbx_description')
        or j.get('rcsb_polymer_entity', {}).get('polymer_entity_name')
        or j.get('rcsb_polymer_entity_container_identifiers', {}).get('entity_name')
        or 'Protein entity'
    )
    # organism
    org = None
    srcs = j.get('rcsb_polymer_entity', {}).get('source_organism', []) or []
    if isinstance(srcs, list) and srcs:
        sci = srcs[0].get('scientific_name') or srcs[0].get('name')
        if sci:
            org = sci
    # function/description
    func = j.get('rcsb_polymer_entity', {}).get('description') or None
    # sequence
    seq = (
        j.get('rcsb_polymer_entity', {}).get('sequence')
        or j.get('entity_poly', {}).get('pdbx_seq_one_letter_code_can')
        or j.get('entity_poly', {}).get('pdbx_seq_one_letter_code')
    )
    if seq:
        seq = ''.join([c for c in seq if c.isalpha()])
    return name, org, func, seq


def _sanitize_protein_seq(seq: str) -> str:
    """Map non-standard/ambiguous residues to reasonable substitutes for scoring.
    B->N, Z->Q, J->L, U->C, O->K, X/others->A. Uppercase and keep A-Z only.
    """
    if not seq:
        return ''
    s = ''.join([c for c in seq.upper() if c.isalpha()])
    mapping = {
        'B': 'N',  # Asx
        'Z': 'Q',  # Glx
        'J': 'L',  # Leu/Ile ambiguous
        'U': 'C',  # Selenocysteine -> Cys
        'O': 'K',  # Pyrrolysine -> Lys
    }
    out = []
    for c in s:
        if c in 'ACDEFGHIKLMNPQRSTVWY':
            out.append(c)
        elif c in mapping:
            out.append(mapping[c])
        else:
            out.append('A')  # unknown -> Ala
    return ''.join(out)


def rcsb_sequence_search(query_seq: str, rows: int = 10, request_id: str | None = None) -> list[dict]:
    # Try multiple payload variants for maximum compatibility with RCSB schema
    qseq = ''.join([c for c in query_seq.upper() if c.isalpha()])
    logger.info(
        json.dumps({
            'event': 'rcsb_sequence_search.start',
            'rid': request_id,
            'len': len(qseq),
            'sha1': _seq_sha1(qseq),
            'preview': _seq_preview(qseq),
            'rows': rows
        })
    )
    targets = [
        ("polymer_entity", ["polymer_entity"]),
        ("polymer_entity_instance", ["polymer_entity_instance"]),
        # Prefer instance return_type first for pdb_protein_sequence so IDs like 1AKI.A are returned
        ("pdb_protein_sequence", ["polymer_entity_instance", "polymer_entity"]),
    ]
    variants = []
    for target, return_types in targets:
        for ret in return_types:
            # Group form
            base_group = {
            "query": {
                "type": "group",
                "logical_operator": "and",
                "nodes": [{
                    "type": "terminal",
                    "service": "sequence",
                    "parameters": {
                        "target": target,
                        "value": qseq
                    }
                }]
            },
            "return_type": ret
            }
            # Terminal form
            base_terminal = {
            "query": {
                "type": "terminal",
                "service": "sequence",
                "parameters": {
                    "target": target,
                    "value": qseq
                }
            },
            "return_type": ret
            }
            # Variants per target
            # A: group with evalue + identity + sequence_type
            a = json.loads(json.dumps(base_group))
            a["query"]["nodes"][0]["parameters"].update({"evalue_cutoff": 10.0, "identity_cutoff": 0.3, "sequence_type": "protein"})
            variants.append(a)
            # B: group with evalue only
            b = json.loads(json.dumps(base_group))
            b["query"]["nodes"][0]["parameters"].update({"evalue_cutoff": 10.0})
            variants.append(b)
            # C: terminal with evalue + identity + sequence_type
            c = json.loads(json.dumps(base_terminal))
            c["query"]["parameters"].update({"evalue_cutoff": 10.0, "identity_cutoff": 0.3, "sequence_type": "protein"})
            variants.append(c)
            # D: terminal minimal
            d = json.loads(json.dumps(base_terminal))
            variants.append(d)

    errors = []
    saw_204 = False
    for payload in variants:
        meta = {}
        try:
            # Add pager via request_options only if schema allows; start with simplest
            meta = {
                'target': payload.get('query', {}).get('parameters', {}).get('target') or 
                          (payload.get('query', {}).get('nodes', [{}])[0].get('parameters', {}).get('target') if payload.get('query', {}).get('nodes') else None),
                'return_type': payload.get('return_type')
            }
            t0 = time.time()
            resp = requests.post(RCSB_SEARCH_URL, json=payload, timeout=20)
            dt = round((time.time() - t0) * 1000)
            logger.info(json.dumps({'event': 'rcsb_sequence_search.try', 'rid': request_id, 'status': resp.status_code, 'dt_ms': dt, **meta}))
            if resp.status_code == 204:
                # No content -> no hits
                saw_204 = True
                logger.info(json.dumps({'event': 'rcsb_sequence_search.no_content', 'rid': request_id, **meta}))
                continue
            if resp.status_code != 200:
                errors.append(f"HTTP {resp.status_code}: {resp.text[:200]}")
                continue
            data = resp.json() or {}
            result_set = data.get('result_set', []) or []
            logger.info(json.dumps({'event': 'rcsb_sequence_search.ok', 'rid': request_id, 'count': len(result_set), **meta}))
            if result_set:
                sample_ids = [r.get('identifier') for r in result_set[:5]]
                logger.info(json.dumps({'event': 'rcsb_sequence_search.sample', 'rid': request_id, 'ids': sample_ids}))
            # If there are many results and rows was requested, try pagination if supported
            if result_set and rows and len(result_set) < rows:
                pass  # keep simple; RCSB returns top-ranked by default
            return result_set
        except Exception as e:
            logger.warning(json.dumps({'event': 'rcsb_sequence_search.exception', 'rid': request_id, 'error': str(e), **meta}))
            errors.append(str(e))
            continue
    # No variant worked
    if saw_204:
        logger.info(json.dumps({'event': 'rcsb_sequence_search.empty', 'rid': request_id}))
        return []
    logger.error(json.dumps({'event': 'rcsb_sequence_search.fail', 'rid': request_id, 'errors': errors[:3]}))
    raise requests.HTTPError("RCSB search failed: " + " | ".join(errors))


def rcsb_similarity_top3(query_seq: str, request_id: str | None = None) -> tuple[list[dict], dict]:
    # Fetch candidates from RCSB, then compute alignment-based identity locally, return top 3
    t_start = time.time()
    candidates = rcsb_sequence_search(query_seq, rows=20, request_id=request_id)
    logger.info(json.dumps({'event': 'similarity_top3.candidates', 'rid': request_id, 'count': len(candidates)}))
    results = []
    skipped = {'no_entity': 0, 'no_seq': 0}
    for c in candidates:
        ent_id = c.get('identifier')
        score = c.get('score', 0)
        if not ent_id:
            continue
        j = _rcsb_fetch_entity(ent_id, request_id=request_id)
        if not j:
            skipped['no_entity'] += 1
            continue
        name, org, func, seq = _extract_entity_metadata(j)
        if not seq:
            skipped['no_seq'] += 1
            continue
        subj_seq = _sanitize_protein_seq(seq)
        if not subj_seq:
            skipped['no_seq'] += 1
            continue
        # Align and compute identity/coverage
        t0 = time.time()
        s, qa, sb, _ = smith_waterman(query_seq, subj_seq)
        dt = round((time.time() - t0) * 1000)
        aln_len = sum(1 for x, y in zip(qa, sb) if x != '-' or y != '-')
        cov = aln_len / max(1, len(query_seq))
        pid = percent_identity(qa, sb)
        logger.info(json.dumps({
            'event': 'similarity_top3.align', 'rid': request_id,
            'entity': ent_id, 'score': int(score) if isinstance(score, (int, float)) else 0,
            'seq_len': len(seq), 'aln_len': aln_len, 'coverage': round(cov, 3), 'identity_percent': round(pid, 2), 'dt_ms': dt
        }))
        results.append({
            'id': ent_id,
            'name': name,
            'organism': org,
            'function': func,
            'length': len(subj_seq),
            'score': int(score) if isinstance(score, (int, float)) else 0,
            'identity_percent': round(pid, 2),
            'alignment_length': aln_len,
            'coverage': round(cov, 3),
            '_seq': subj_seq,
        })
    # Rank: score -> coverage -> identity
    results.sort(key=lambda h: (h['score'], h['coverage'], h['identity_percent']), reverse=True)
    top3 = results[:3]
    top_props = compute_properties(top3[0]['_seq']) if top3 else {}
    # Remove internal _seq from hits
    for h in top3:
        h.pop('_seq', None)
    logger.info(json.dumps({'event': 'similarity_top3.done', 'rid': request_id, 'top_count': len(top3), 'dt_ms': round((time.time() - t_start) * 1000), 'skipped': skipped}))
    return top3, top_props

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

@app.route('/api/sample/<path:sample_path>')
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
    # If the message looks like or contains a protein sequence, augment context with similarity results
    context = None
    # Find longest AA-only substring (ACDEFGHIKLMNPQRSTVWY) length >= 10
    aa_regex = re.compile(r'[ACDEFGHIKLMNPQRSTVWY]{10,}', re.IGNORECASE)
    matches = aa_regex.findall(user_message)
    seq_candidate = max(matches, key=len).upper() if matches else None
    if seq_candidate and is_valid_protein_sequence(seq_candidate) and 10 <= len(seq_candidate) <= 100:
        try:
            # Use RCSB remote search only
            rid = uuid.uuid4().hex[:8]
            logger.info(json.dumps({'event': 'chat.similarity.start', 'rid': rid, 'len': len(seq_candidate), 'sha1': _seq_sha1(seq_candidate), 'preview': _seq_preview(seq_candidate)}))
            hits_json, props = rcsb_similarity_top3(seq_candidate, request_id=rid)
            logger.info(json.dumps({'event': 'chat.similarity.ok', 'rid': rid, 'hits': len(hits_json)}))
            context = (
                f"Similarity results for query (length {len(seq_candidate)}):\n"
                f"Top hits: {hits_json}\n"
                f"Top-hit properties: {props}\n"
                f"Use this as background; respond concisely and do not over-claim."
            )
        except Exception as e:
            logger.error(json.dumps({'event': 'chat.similarity.error', 'error': str(e)}))
            context = f"Similarity search error: {e}"
    return jsonify({'response': chatbot.ChatBot(user_message, context=context)})
    #return jsonify({'response': "Wasgood"})

@app.route('/api/similarity', methods=['POST'])
def similarity_api():
    try:
        data = request.get_json(force=True)
        query = data.get('sequence', '')
        rid = uuid.uuid4().hex[:8]
        logger.info(json.dumps({'event': 'similarity_api.request', 'rid': rid, 'len': len(query or ''), 'sha1': _seq_sha1(query or ''), 'preview': _seq_preview(query or '')}))
        if not query:
            return jsonify({'error': 'Missing sequence'}), 400
        if not is_valid_protein_sequence(query):
            logger.warning(json.dumps({'event': 'similarity_api.invalid_seq', 'rid': rid}))
            return jsonify({'error': 'Invalid protein sequence (use standard amino acids).'}), 400
        # if len(query) > 100:
        #     logger.warning(json.dumps({'event': 'similarity_api.too_long', 'rid': rid, 'len': len(query)}))
        #     return jsonify({'error': 'Sequence too long for this search (max 100 aa).'}), 400
        if len(query) < 10:
            logger.warning(json.dumps({'event': 'similarity_api.too_short', 'rid': rid, 'len': len(query)}))
            return jsonify({'error': 'Sequence too short for similarity search (min 10 aa).'}), 400
        # Use RCSB remote similarity search only
        t0 = time.time()
        hits_json, props = rcsb_similarity_top3(query, request_id=rid)
        dt = round((time.time() - t0) * 1000)
        logger.info(json.dumps({'event': 'similarity_api.response', 'rid': rid, 'hits': len(hits_json), 'dt_ms': dt}))
        return jsonify({'success': True, 'hits': hits_json, 'top_properties': props})
    except requests.HTTPError as he:
        logger.error(json.dumps({'event': 'similarity_api.rcsb_http_error', 'error': str(he)}))
        return jsonify({'error': f'RCSB HTTP error: {he}'}), 502
    except Exception as e:
        logger.error(json.dumps({'event': 'similarity_api.error', 'error': str(e)}))
        return jsonify({'error': f'Similarity search failed: {e}'}), 500

@app.route('/api/sequence', methods=['POST'])
def receive_sequence():
    data = request.get_json()
    sequence = data.get('sequence', '')
    print(f"Received sequence: {sequence}")
    if pdb_script is None:
        return jsonify({'error': 'Structure generation backend not available (missing colabfold/pdb_script).'}), 503
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

@app.route('/api/dock', methods=['POST'])
def protein_dock():
    """Handle protein-protein docking requests"""
    try:
        if run_docking_pipeline is None:
            return jsonify({'error': 'Docking functionality not available'}), 503
            
        # Check if we have file uploads
        if 'target_file' in request.files and 'ligand_file' in request.files:
            target_file = request.files['target_file']
            ligand_file = request.files['ligand_file']
            
            if target_file.filename == '' or ligand_file.filename == '':
                return jsonify({'error': 'Both target and ligand files must be provided'}), 400
            
            if not (allowed_file(target_file.filename) and allowed_file(ligand_file.filename)):
                return jsonify({'error': 'Both files must be PDB format'}), 400
                
            # Save uploaded files
            target_filename = secure_filename(target_file.filename)
            ligand_filename = secure_filename(ligand_file.filename)
            target_path = os.path.join(app.config['UPLOAD_FOLDER'], f"dock_target_{uuid.uuid4().hex[:8]}_{target_filename}")
            ligand_path = os.path.join(app.config['UPLOAD_FOLDER'], f"dock_ligand_{uuid.uuid4().hex[:8]}_{ligand_filename}")
            
            target_file.save(target_path)
            ligand_file.save(ligand_path)
            
        else:
            # Check for file paths in JSON data
            data = request.get_json()
            if not data:
                return jsonify({'error': 'No files or data provided'}), 400
                
            target_path = data.get('target_path')
            ligand_path = data.get('ligand_path')
            
            if not target_path or not ligand_path:
                return jsonify({'error': 'Target and ligand file paths must be provided'}), 400
                
            # Check if files exist
            if not os.path.exists(target_path) or not os.path.exists(ligand_path):
                return jsonify({'error': 'One or both specified files do not exist'}), 404
        
        # Get docking parameters
        data = data if 'data' in locals() else request.form
        num_poses = int(data.get('num_poses', 5))
        seed = int(data.get('seed', 42))
        
        # Create output directory for this docking run
        dock_id = uuid.uuid4().hex[:8]
        out_dir = os.path.join(app.config['UPLOAD_FOLDER'], f'docking_{dock_id}')
        
        # Run docking pipeline
        try:
            df = run_docking_pipeline(target_path, ligand_path, out_dir, num_poses=num_poses, seed=seed)
            
            # Parse the best pose for structure data
            best_pose_path = os.path.join(out_dir, 'best_pose.pdb')
            if os.path.exists(best_pose_path):
                parser = PDBParser()
                best_structure = parser.parse_pdb_file(best_pose_path)
            else:
                best_structure = None
            
            # Convert DataFrame to dict for JSON response
            results = df.to_dict('records')
            
            # Add file paths for download
            for i, result in enumerate(results):
                pose_id = result['Pose ID']
                pose_filename = f'pose_{pose_id}.pdb'
                result['pose_download_path'] = f'/api/download/{dock_id}/{pose_filename}'
            
            response = {
                'success': True,
                'dock_id': dock_id,
                'results': results,
                'best_pose_structure': best_structure,
                'best_pose_download': f'/api/download/{dock_id}/best_pose.pdb',
                'summary_download': f'/api/download/{dock_id}/pose_ranking.csv'
            }
            
            return jsonify(response)
            
        except Exception as e:
            return jsonify({'error': f'Docking failed: {str(e)}'}), 500
            
        finally:
            # Clean up uploaded files if they were temporary
            if 'target_file' in request.files:
                try:
                    if os.path.exists(target_path):
                        os.remove(target_path)
                    if os.path.exists(ligand_path):
                        os.remove(ligand_path)
                except Exception:
                    pass
    
    except Exception as e:
        return jsonify({'error': f'Docking request failed: {str(e)}'}), 500

@app.route('/api/download/<dock_id>/<filename>')
def download_docking_file(dock_id, filename):
    """Download docking result files"""
    try:
        file_path = os.path.join(app.config['UPLOAD_FOLDER'], f'docking_{dock_id}', filename)
        if not os.path.exists(file_path):
            return jsonify({'error': 'File not found'}), 404
            
        from flask import send_file
        return send_file(file_path, as_attachment=True)
        
    except Exception as e:
        return jsonify({'error': f'Download failed: {str(e)}'}), 500

@app.route('/api/ligand/upload', methods=['POST'])
def upload_ligand():
    """Handle ligand file uploads (SDF, PDB, MOL, MOL2)"""
    try:
        if 'ligand_file' not in request.files:
            return jsonify({'error': 'No ligand file provided'}), 400
        
        file = request.files['ligand_file']
        if file.filename == '':
            return jsonify({'error': 'No file selected'}), 400
        
        # Check file extension
        filename = secure_filename(file.filename)
        file_ext = filename.rsplit('.', 1)[1].lower() if '.' in filename else ''
        
        if file_ext not in ['sdf', 'pdb', 'mol', 'mol2']:
            return jsonify({'error': 'Invalid file type. Please upload SDF, PDB, MOL, or MOL2 files.'}), 400
        
        # Save file temporarily
        ligand_id = uuid.uuid4().hex[:8]
        ligand_filename = f"ligand_{ligand_id}_{filename}"
        ligand_path = os.path.join(app.config['UPLOAD_FOLDER'], ligand_filename)
        file.save(ligand_path)
        
        try:
            # Parse ligand file based on format
            if file_ext == 'pdb':
                # Use existing PDB parser
                parser = PDBParser()
                structure_data = parser.parse_pdb_file(ligand_path)
            else:
                # For SDF/MOL files, create a basic structure representation
                # This is a simplified parser - in production you'd want RDKit
                structure_data = parse_ligand_file(ligand_path, file_ext)
            
            # Store the file path for later use in docking
            structure_data['_ligand_file_path'] = ligand_path
            
            return jsonify({
                'success': True,
                'ligand_id': ligand_id,
                'filename': filename,
                'structure': structure_data
            })
            
        except Exception as e:
            # Clean up file on error
            if os.path.exists(ligand_path):
                os.remove(ligand_path)
            return jsonify({'error': f'Failed to parse ligand file: {str(e)}'}), 400
    
    except Exception as e:
        return jsonify({'error': f'Ligand upload failed: {str(e)}'}), 500

def parse_ligand_file(file_path: str, file_ext: str) -> dict:
    """Basic ligand file parser for SDF/MOL formats"""
    atoms = []
    bonds = []
    
    try:
        with open(file_path, 'r') as f:
            lines = f.readlines()
        
        if file_ext in ['sdf', 'mol']:
            # Simple MOL/SDF parser
            if len(lines) < 4:
                raise ValueError("Invalid MOL/SDF file format")
            
            # Parse counts line (line 4)
            counts_line = lines[3].strip()
            if len(counts_line) >= 6:
                try:
                    atom_count = int(counts_line[:3].strip())
                    bond_count = int(counts_line[3:6].strip())
                except ValueError:
                    raise ValueError("Invalid counts line in MOL file - non-numeric values")
            else:
                raise ValueError("Invalid counts line in MOL file")
            
            # Parse atoms
            for i in range(4, 4 + atom_count):
                if i >= len(lines):
                    break
                line = lines[i]
                if len(line) >= 31:
                    x = float(line[0:10].strip())
                    y = float(line[10:20].strip()) 
                    z = float(line[20:30].strip())
                    element = line[31:34].strip()
                    
                    atoms.append({
                        'serial': i - 3,
                        'name': element + str(i - 3),
                        'res_name': 'LIG',
                        'chain': 'L',
                        'res_seq': 1,
                        'x': x,
                        'y': y,
                        'z': z,
                        'element': element,
                        'occupancy': 1.0,
                        'temp_factor': 20.0
                    })
        
        elif file_ext == 'mol2':
            # Basic MOL2 parser
            in_atoms = False
            atom_count = 0
            
            for line in lines:
                line = line.strip()
                if line == '@<TRIPOS>ATOM':
                    in_atoms = True
                    continue
                elif line.startswith('@<TRIPOS>') and line != '@<TRIPOS>ATOM':
                    in_atoms = False
                    continue
                
                if in_atoms and line:
                    parts = line.split()
                    if len(parts) >= 6:
                        atom_count += 1
                        atoms.append({
                            'serial': atom_count,
                            'name': parts[1],
                            'res_name': 'LIG',
                            'chain': 'L', 
                            'res_seq': 1,
                            'x': float(parts[2]),
                            'y': float(parts[3]),
                            'z': float(parts[4]),
                            'element': parts[5].split('.')[0],
                            'occupancy': 1.0,
                            'temp_factor': 20.0
                        })
        
        if not atoms:
            raise ValueError("No atoms found in ligand file")
        
        return {
            'header': {
                'title': 'LIGAND',
                'resolution': None,
                'experiment_type': 'LIGAND',
                'organism': None
            },
            'atoms': atoms,
            'stats': {
                'atom_count': len(atoms),
                'chain_count': 1,
                'residue_count': 1,
                'bond_count': len(bonds)
            }
        }
        
    except Exception as e:
        raise ValueError(f"Error parsing ligand file: {str(e)}")
    

if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0', port=8000)




