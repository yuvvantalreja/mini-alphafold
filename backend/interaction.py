import os
import io
import json
import hashlib
import math
import random
from dataclasses import dataclass, asdict
from typing import List, Dict, Optional, Tuple, Any

from Bio.PDB import PDBParser as BioPDBParser, PDBIO, Select, NeighborSearch, ShrakeRupley
from Bio.PDB.StructureBuilder import StructureBuilder
from Bio.PDB.vectors import rotaxis2m, Vector


def _sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode('utf-8')).hexdigest()


def _ensure_dir(path: str) -> None:
    os.makedirs(path, exist_ok=True)


@dataclass
class PoseScores:
    confidence: float
    binding_energy: float
    buried_surface_area: float
    contact_count: int
    hydrogen_bonds: int
    salt_bridges: int
    shape_complementarity: float
    final_score: float


@dataclass
class PoseResult:
    pose_id: int
    pdb: str
    scores: PoseScores


class ComplexPredictor:
    """Protein-protein interaction pose prediction with optional refinement and analysis.

    Primary path tries to call an external EquiDock/EquiBind-like CLI if available.
    Fallback path uses simple rigid-body sampling and heuristic scoring.
    """

    def __init__(self, cache_dir: str = 'backend_cache'):
        self.cache_dir = os.path.abspath(cache_dir)
        _ensure_dir(self.cache_dir)

    # --------------------------- Public API ---------------------------
    def predict(
        self,
        pdb_a_text: str,
        pdb_b_text: str,
        num_poses: int = 5,
        refine: bool = False,
        random_seed: Optional[int] = 42,
    ) -> List[PoseResult]:
        cache_key = self._make_cache_key(pdb_a_text, pdb_b_text, num_poses, refine)
        cached = self._read_cache(cache_key)
        if cached is not None:
            return cached

        # Try external predictor first
        poses = self._try_external_predictor(pdb_a_text, pdb_b_text, num_poses=num_poses)
        if not poses:
            poses = self._fallback_sampling(pdb_a_text, pdb_b_text, num_poses=num_poses, random_seed=random_seed)

        # Optional refinement
        if refine:
            poses = [self._refine_pose(pose) for pose in poses]

        # Score and rank
        poses = self._score_and_rank_poses(poses)

        self._write_cache(cache_key, poses)
        return poses

    # --------------------------- Cache ---------------------------
    def _make_cache_key(self, a: str, b: str, n: int, refine: bool) -> str:
        digest = _sha256_text(a)[:12] + '_' + _sha256_text(b)[:12]
        return f"ppip_{digest}_n{n}_r{int(refine)}"

    def _cache_path(self, key: str) -> str:
        return os.path.join(self.cache_dir, f"{key}.json")

    def _read_cache(self, key: str) -> Optional[List[PoseResult]]:
        path = self._cache_path(key)
        if not os.path.exists(path):
            return None
        try:
            with open(path, 'r') as f:
                data = json.load(f)
            poses = [PoseResult(p['pose_id'], p['pdb'], PoseScores(**p['scores'])) for p in data]
            return poses
        except Exception:
            return None

    def _write_cache(self, key: str, poses: List[PoseResult]) -> None:
        path = self._cache_path(key)
        try:
            with open(path, 'w') as f:
                json.dump([{'pose_id': p.pose_id, 'pdb': p.pdb, 'scores': asdict(p.scores)} for p in poses], f)
        except Exception:
            pass

    # --------------------------- External predictor hook ---------------------------
    def _try_external_predictor(self, pdb_a_text: str, pdb_b_text: str, num_poses: int) -> List[PoseResult]:
        """Hook for EquiDock/EquiBind CLI if user installed it. Return [] if unavailable.

        Protocol: If environment variable EQUIDOCK_CLI is set to an executable that accepts:
            equidock_cli --receptor a.pdb --ligand b.pdb --num_poses N --out out_dir
        and writes pose_{i}.pdb and a JSON scores file, we will parse and return.
        """
        import shutil
        import subprocess
        cli = os.environ.get('EQUIDOCK_CLI') or os.environ.get('EQUIBIND_CLI')
        if not cli or not shutil.which(cli):
            return []

        tmpdir = os.path.join(self.cache_dir, f"tmp_{random.randint(1, 1_000_000)}")
        _ensure_dir(tmpdir)
        a_path = os.path.join(tmpdir, 'a.pdb')
        b_path = os.path.join(tmpdir, 'b.pdb')
        with open(a_path, 'w') as f:
            f.write(pdb_a_text)
        with open(b_path, 'w') as f:
            f.write(pdb_b_text)

        out_dir = os.path.join(tmpdir, 'out')
        _ensure_dir(out_dir)

        try:
            cmd = [cli, '--receptor', a_path, '--ligand', b_path, '--num_poses', str(num_poses), '--out', out_dir]
            subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        except Exception:
            # Cleanup and bail
            try:
                import shutil as _sh
                _sh.rmtree(tmpdir)
            except Exception:
                pass
            return []

        poses: List[PoseResult] = []
        # Read PDB poses
        for i in range(1, num_poses + 1):
            pose_path = os.path.join(out_dir, f'pose_{i}.pdb')
            if not os.path.exists(pose_path):
                continue
            with open(pose_path, 'r') as f:
                pdb_text = f.read()
            # Scores JSON optional
            scores_json = os.path.join(out_dir, f'pose_{i}.json')
            if os.path.exists(scores_json):
                with open(scores_json, 'r') as f:
                    s = json.load(f)
                scores = PoseScores(
                    confidence=float(s.get('confidence', 0.5)),
                    binding_energy=float(s.get('binding_energy', 0.0)),
                    buried_surface_area=float(s.get('bsa', 0.0)),
                    contact_count=int(s.get('contact_count', 0)),
                    hydrogen_bonds=int(s.get('hydrogen_bonds', 0)),
                    salt_bridges=int(s.get('salt_bridges', 0)),
                    shape_complementarity=float(s.get('shape_complementarity', 0.5)),
                    final_score=float(s.get('final_score', 0.0)),
                )
            else:
                scores = PoseScores(0.5, 0.0, 0.0, 0, 0, 0, 0.5, 0.0)
            poses.append(PoseResult(i, pdb_text, scores))

        # Cleanup temp dir
        try:
            import shutil as _sh
            _sh.rmtree(tmpdir)
        except Exception:
            pass

        return poses

    # --------------------------- Fallback sampling ---------------------------
    def _fallback_sampling(self, pdb_a_text: str, pdb_b_text: str, num_poses: int, random_seed: Optional[int]) -> List[PoseResult]:
        random.seed(random_seed or 42)
        struct_a = self._read_structure_from_text(pdb_a_text, struct_id='A')
        struct_b_original = self._read_structure_from_text(pdb_b_text, struct_id='B')
        center_a = self._calc_structure_center(struct_a)
        center_b_orig = self._calc_structure_center(struct_b_original)

        poses: List[PoseResult] = []
        base_distance = 8.0  # Å separation between surfaces before random adjustment
        for i in range(num_poses):
            # Fresh copy of B each iteration
            struct_b = self._read_structure_from_text(pdb_b_text, struct_id=f'B{i}')
            # Recentre B to origin relative to A
            translated_b = self._transform_structure(struct_b, translation=-(self._calc_structure_center(struct_b) - center_a))
            # Random small rotation
            axis = Vector(random.random(), random.random(), random.random())
            angle = random.uniform(0, math.pi)
            rot = rotaxis2m(angle, axis)
            rotated_b = self._rotate_structure(translated_b, rot)

            # Translate along vector from A center to ensure contact
            direction = Vector(0.0, 0.0, 1.0)  # arbitrary; viewer-independent
            jitter = Vector(random.uniform(-2, 2), random.uniform(-2, 2), random.uniform(-2, 2))
            # Manually scale direction vector and add jitter
            scaled_direction = Vector(direction[0] * base_distance, direction[1] * base_distance, direction[2] * base_distance)
            translation = Vector(scaled_direction[0] + jitter[0], scaled_direction[1] + jitter[1], scaled_direction[2] + jitter[2])
            pose_struct = self._translate_structure(rotated_b, translation)

            # Combine A and B into one structure, with chains labeled distinctly
            combined = self._combine_structures(struct_a, pose_struct)
            pdb_text = self._write_structure_to_pdb_text(combined)

            # Placeholder scores; will be recomputed in _score_and_rank_poses
            placeholder_scores = PoseScores(0.5, 0.0, 0.0, 0, 0, 0, 0.5, 0.0)
            poses.append(PoseResult(i + 1, pdb_text, placeholder_scores))

        return poses

    # --------------------------- Refinement (OpenMM optional) ---------------------------
    def _refine_pose(self, pose: PoseResult) -> PoseResult:
        try:
            from openmm import LangevinMiddleIntegrator, Platform
            from openmm.app import PDBFile, Modeller, ForceField, Simulation
            from openmm.unit import kelvin, picoseconds, femtoseconds, nanometer
        except Exception:
            # OpenMM not available; return pose unchanged
            return pose

        pdb = PDBFile(io.StringIO(pose.pdb))
        modeller = Modeller(pdb.topology, pdb.positions)

        # Use Amber14 with GB implicit solvent
        forcefield = ForceField('amber14-all.xml', 'amber14/gbn2.xml')
        system = forcefield.createSystem(
            modeller.topology,
            nonbondedMethod=None,
            constraints=None,
            removeCMMotion=True,
        )

        integrator = LangevinMiddleIntegrator(300*kelvin, 1.0/picoseconds, 2.0*femtoseconds)
        platform = Platform.getPlatformByName('Reference')
        simulation = Simulation(modeller.topology, system, integrator, platform)
        simulation.context.setPositions(modeller.positions)

        # Energy minimization followed by a very short MD to relax clashes
        try:
            simulation.minimizeEnergy(maxIterations=200)
            simulation.step(1500)  # ~3 ps
        except Exception:
            pass

        # Get refined positions back to PDB text
        refined_io = PDBIO()
        refined_io.set_structure(self._bio_structure_from_simulation(simulation))
        buf = io.StringIO()
        refined_io.save(buf)
        refined_pdb = buf.getvalue()
        return PoseResult(pose.pose_id, refined_pdb, pose.scores)

    def _bio_structure_from_simulation(self, simulation) -> Any:
        """Convert OpenMM simulation state to Bio.PDB Structure for saving."""
        state = simulation.context.getState(getPositions=True)
        positions = state.getPositions(asNumpy=True)
        # Re-parse original topology via PDBIO approach using Modeller topology
        # Simplest path: dump current positions to PDB text via OpenMM app and re-read with Bio.PDB
        from openmm.app import PDBFile
        buf = io.StringIO()
        PDBFile.writeFile(simulation.topology, positions, buf)
        buf.seek(0)
        return self._read_structure_from_text(buf.getvalue(), struct_id='REF')

    # --------------------------- Scoring and analysis ---------------------------
    def _score_and_rank_poses(self, poses: List[PoseResult]) -> List[PoseResult]:
        scored: List[Tuple[PoseResult, PoseScores]] = []
        for pose in poses:
            scores = self._analyze_pose(pose.pdb)
            # Combine with binding energy heuristic into final score
            final = 0.3 * scores.buried_surface_area \
                    + 0.2 * scores.contact_count \
                    + 0.2 * scores.hydrogen_bonds \
                    + 0.2 * scores.shape_complementarity \
                    + 0.1 * scores.binding_energy
            scores.final_score = float(final)
            scored.append((pose, scores))

        # Normalize confidence by rank (higher score => higher confidence)
        scored.sort(key=lambda ps: ps[1].final_score, reverse=True)
        max_score = scored[0][1].final_score if scored else 1.0
        min_score = scored[-1][1].final_score if scored else 0.0
        score_span = max(max_score - min_score, 1e-6)

        ranked: List[PoseResult] = []
        for rank, (pose, scores) in enumerate(scored, start=1):
            confidence = (scores.final_score - min_score) / score_span
            scores.confidence = float(max(0.0, min(1.0, confidence)))
            pose.scores = scores
            ranked.append(pose)

        return ranked

    def _analyze_pose(self, pdb_text: str) -> PoseScores:
        structure = self._read_structure_from_text(pdb_text, struct_id='POSE')
        # Identify two chains (assume first two chains belong to partners)
        model = next(structure.get_models())
        chains = list(model.get_chains())
        if len(chains) < 2:
            # Not enough chains; treat as single with duplicates
            return PoseScores(0.5, 0.0, 0.0, 0, 0, 0, 0.5, 0.0)

        chain_a, chain_b = chains[0], chains[1]

        # Contacts
        atoms_a = [a for a in chain_a.get_atoms() if a.element != 'H']
        atoms_b = [a for a in chain_b.get_atoms() if a.element != 'H']
        ns = NeighborSearch(atoms_a + atoms_b)
        contact_pairs = set()
        for atom in atoms_a:
            neighbors = ns.search(atom.coord, 4.0, level='A')
            for nb in neighbors:
                if nb.get_parent().get_parent().id == chain_b.id:
                    # residue-residue level contact
                    key = (atom.get_parent().id, nb.get_parent().id)
                    contact_pairs.add(key)
        contact_count = len(contact_pairs)

        # Hydrogen bonds (simple distance criterion N-O <= 3.5 Å across chains)
        n_atoms_a = [a for a in atoms_a if a.element == 'N']
        o_atoms_b = [a for a in atoms_b if a.element == 'O']
        n_atoms_b = [a for a in atoms_b if a.element == 'N']
        o_atoms_a = [a for a in atoms_a if a.element == 'O']
        hbonds = 0
        for na in n_atoms_a:
            for ob in o_atoms_b:
                distance_vec = Vector(na.coord) - Vector(ob.coord)
                if distance_vec.norm() <= 3.5:
                    hbonds += 1
        for nb in n_atoms_b:
            for oa in o_atoms_a:
                distance_vec = Vector(nb.coord) - Vector(oa.coord)
                if distance_vec.norm() <= 3.5:
                    hbonds += 1

        # Salt bridges: Lys/Arg/His (positives) with Asp/Glu (negatives) heavy atoms within 4 Å
        pos_res = {'LYS', 'ARG', 'HIS'}
        neg_res = {'ASP', 'GLU'}
        pos_atoms = [a for a in atoms_a if a.get_parent().get_resname() in pos_res] + \
                    [a for a in atoms_b if a.get_parent().get_resname() in pos_res]
        neg_atoms = [a for a in atoms_a if a.get_parent().get_resname() in neg_res] + \
                    [a for a in atoms_b if a.get_parent().get_resname() in neg_res]
        salt_bridges = 0
        for pa in pos_atoms:
            for na in neg_atoms:
                if pa.get_parent().get_parent().id == na.get_parent().get_parent().id:
                    continue  # same chain
                distance_vec = Vector(pa.coord) - Vector(na.coord)
                if distance_vec.norm() <= 4.0:
                    salt_bridges += 1

        # Buried Surface Area using Shrake-Rupley
        sr = ShrakeRupley()
        # ASA for isolated chains
        structure_isolated = self._read_structure_from_text(self._write_structure_to_pdb_text(structure), 'ISO')
        model_iso = next(structure_isolated.get_models())
        chain_iso_a, chain_iso_b = list(model_iso.get_chains())[0:2]
        sr.compute(chain_iso_a, level='A')
        sr.compute(chain_iso_b, level='A')
        asa_a = sum(getattr(atom, 'sasa', 0.0) for atom in chain_iso_a.get_atoms())
        asa_b = sum(getattr(atom, 'sasa', 0.0) for atom in chain_iso_b.get_atoms())

        # ASA for complex
        sr.compute(model, level='A')
        asa_complex = sum(getattr(atom, 'sasa', 0.0) for atom in model.get_atoms())
        bsa = max(0.0, asa_a + asa_b - asa_complex)

        # Shape complementarity heuristic: normalized contact density around interface
        shape_comp = 0.0
        if contact_count > 0:
            # Normalize by approximate size (sqrt of atom counts) and cap to 1.0
            norm = math.sqrt(max(1, len(atoms_a))) + math.sqrt(max(1, len(atoms_b)))
            shape_comp = min(1.0, (contact_count / max(1.0, norm)))

        # Pseudo binding energy: more negative is better
        binding_energy = (
            -0.05 * bsa
            -0.8 * hbonds
            -0.6 * salt_bridges
            -0.2 * contact_count
        )

        return PoseScores(
            confidence=0.5,
            binding_energy=float(binding_energy),
            buried_surface_area=float(bsa),
            contact_count=int(contact_count),
            hydrogen_bonds=int(hbonds),
            salt_bridges=int(salt_bridges),
            shape_complementarity=float(shape_comp),
            final_score=0.0,
        )

    # --------------------------- Bio.PDB helpers ---------------------------
    def _read_structure_from_text(self, pdb_text: str, struct_id: str) -> Any:
        parser = BioPDBParser(QUIET=True)
        handle = io.StringIO(pdb_text)
        structure = parser.get_structure(struct_id, handle)
        return structure

    def _write_structure_to_pdb_text(self, structure: Any) -> str:
        io_obj = PDBIO()
        io_obj.set_structure(structure)
        buf = io.StringIO()
        io_obj.save(buf)
        return buf.getvalue()

    def _calc_structure_center(self, structure: Any) -> Vector:
        coords = [Vector(a.coord) for a in structure.get_atoms() if a.element != 'H']
        if not coords:
            return Vector(0.0, 0.0, 0.0)
        cx = sum(v[0] for v in coords) / len(coords)
        cy = sum(v[1] for v in coords) / len(coords)
        cz = sum(v[2] for v in coords) / len(coords)
        return Vector(cx, cy, cz)

    def _transform_structure(self, structure: Any, translation: Vector) -> Any:
        return self._translate_structure(structure, translation)

    def _translate_structure(self, structure: Any, translation: Vector) -> Any:
        import numpy as np
        # Convert Vector to numpy array for atom.transform
        translation_array = np.array([translation[0], translation[1], translation[2]])
        for atom in structure.get_atoms():
            atom.transform(np.eye(3), translation_array)
        return structure

    def _rotate_structure(self, structure: Any, rotation_matrix) -> Any:
        import numpy as np
        # Use numpy array for translation (zero translation for rotation)
        zero_translation = np.array([0.0, 0.0, 0.0])
        for atom in structure.get_atoms():
            atom.transform(rotation_matrix, zero_translation)
        return structure

    def _combine_structures(self, struct_a: Any, struct_b: Any) -> Any:
        # Simple approach: write both to PDB and re-parse combined text
        pdb_a = self._write_structure_to_pdb_text(struct_a)
        pdb_b = self._write_structure_to_pdb_text(struct_b)
        combined = []
        for line in pdb_a.splitlines():
            if line.startswith(('ATOM', 'HETATM', 'TER', 'MODEL', 'ENDMDL', 'TITLE', 'REMARK')):
                combined.append(line)
        # Ensure B chain identifiers for second structure
        for line in pdb_b.splitlines():
            if line.startswith(('ATOM', 'HETATM')):
                # Set chain to B at column 22 (index 21)
                line = line[:21] + 'B' + line[22:]
                combined.append(line)
        combined.append('END')
        return self._read_structure_from_text('\n'.join(combined), 'COMB')


