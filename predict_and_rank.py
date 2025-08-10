#!/usr/bin/env python3
"""
Predict and rank protein–protein binding poses using rigid-body sampling 
and simple CPU-friendly scoring metrics.

- Inputs: target.pdb (receptor, fixed), ligand.pdb (moving protein)
- Outputs: outputs/pose_*.pdb, outputs/best_pose.pdb, and a printed ranking table

Implementation notes:
- Generates diverse rigid-body poses (random rotations/translations) as a quick CPU baseline
- Preprocesses PDBs via pdbfixer: removes waters, fills missing residues/atoms, adds hydrogens.
- Binding energy computed via OpenMM GBSA (OBC2) potentials: ΔE = E_complex − (E_target + E_ligand)
  using the pose coordinates without minimization for speed.
- Hydrogen bonds computed via a simple Biopython-based geometric heuristic for inter-chain pairs.
- Interface area computed from SASA via Biopython Shrake-Rupley: IA = ASA_A + ASA_B − ASA_complex.
- Normalization maps each metric to [0,1], inverting those where lower is better.

This script is designed for small proteins (<100 aa each) and runs quickly on CPU.
"""

import argparse
import os
import sys
import shutil
import math
from dataclasses import dataclass
from typing import List, Tuple, Optional

import numpy as np
import pandas as pd

# Biopython for structure IO, transformations, SASA
from Bio.PDB import PDBParser, PDBIO, StructureBuilder, Superimposer
from Bio.PDB.Polypeptide import is_aa
from Bio.PDB.SASA import ShrakeRupley

# pdbfixer for PDB cleanup
try:
    from pdbfixer import PDBFixer
except Exception as e:
    PDBFixer = None

# OpenMM imports: try openmm first (>=7.7/8), then fallback to simtk
try:
    import openmm
    import openmm.app as app
    from openmm import unit
except Exception:
    try:
        import simtk.openmm as openmm  # type: ignore
        import simtk.openmm.app as app  # type: ignore
        from simtk import unit  # type: ignore
    except Exception:
        openmm = None
        app = None
        unit = None


@dataclass
class PoseMetrics:
    pose_id: int
    binding_energy: float  # lower better
    hydrogen_bonds: int    # higher better
    interface_area: float  # higher better
    final_score: Optional[float] = None


# -----------------------------
# Utility: filesystem helpers
# -----------------------------

def ensure_dir(path: str) -> None:
    os.makedirs(path, exist_ok=True)


def safe_copy(src: str, dst: str) -> None:
    ensure_dir(os.path.dirname(dst))
    shutil.copy2(src, dst)


# -----------------------------
# PDB Preprocessing with pdbfixer
# -----------------------------

def preprocess_pdb(input_pdb: str, output_pdb: str, ph: float = 7.0) -> str:
    """Clean a PDB: remove waters/heterogens, add missing residues/atoms, add hydrogens.
    Returns the cleaned PDB path.
    """
    if PDBFixer is None or app is None:
        # If pdbfixer or openmm.app isn't available, just pass-through the file.
        # This keeps pipeline runnable, but quality may degrade.
        safe_copy(input_pdb, output_pdb)
        return output_pdb

    try:
        fixer = PDBFixer(filename=input_pdb)
        # Remove waters (HOH) and heterogens
        fixer.removeHeterogens(keepWater=False)

        # Add missing residues and atoms
        fixer.findMissingResidues()
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()

        # Add hydrogens
        fixer.addMissingHydrogens(ph)

        # Save
        with open(output_pdb, 'w') as f:
            app.PDBFile.writeFile(fixer.topology, fixer.positions, f, keepIds=True)
        return output_pdb
    except Exception as e:
        # If PDBFixer fails (e.g., on ligand files), just copy the original
        print(f"[WARN] PDBFixer failed on {input_pdb}: {e}. Using original file.")
        safe_copy(input_pdb, output_pdb)
        return output_pdb


# -----------------------------
# Structure helpers (Biopython)
# -----------------------------

def load_structure(pdb_path: str, struct_id: str) -> "StructureBuilder.Structure":
    parser = PDBParser(QUIET=True)
    return parser.get_structure(struct_id, pdb_path)


def write_structure(structure, out_path: str) -> None:
    io = PDBIO()
    io.set_structure(structure)
    ensure_dir(os.path.dirname(out_path))
    io.save(out_path)


def get_atoms_by_chain(structure, chain_id: str):
    for model in structure:
        for chain in model:
            if chain.id == chain_id:
                return list(chain.get_atoms())
    return []


def transform_chain(structure, chain_id: str, rotation: np.ndarray, translation: np.ndarray) -> None:
    """Apply rigid transform to a single chain in-place."""
    assert rotation.shape == (3, 3)
    assert translation.shape == (3,)
    for model in structure:
        for chain in model:
            if chain.id != chain_id:
                continue
            for atom in chain.get_atoms():
                coord = atom.get_coord()
                new_coord = rotation.dot(coord) + translation
                atom.set_coord(new_coord)


def merge_structures_as_chains(target_struct, ligand_struct, target_chain_id: str = 'A', ligand_chain_id: str = 'B'):
    """Create a new Biopython structure containing two chains: target(A) and ligand(B)."""
    builder = StructureBuilder.StructureBuilder()
    builder.init_structure('complex')
    builder.init_model(0)
    builder.init_seg(' ')  # Initialize segment with default ID

    # Helper to copy a chain with new ID
    def copy_chain(src_chain, new_chain_id):
        builder.init_chain(new_chain_id)
        # residue serials
        res_id = 1
        for residue in src_chain:
            # Only keep standard amino acids to avoid waters/ligands
            if not is_aa(residue, standard=True):
                continue
            hetfield, resseq, icode = residue.id
            builder.init_residue(residue.resname, hetfield, res_id, icode)
            res_id += 1
            for atom in residue:
                # atom properties
                builder.init_atom(atom.name, atom.get_coord(), atom.bfactor, atom.occupancy, atom.altloc, atom.fullname, atom.serial_number, element=atom.element)

    # Find first chains
    target_chain = next((c for c in target_struct.get_chains()), None)
    ligand_chain = next((c for c in ligand_struct.get_chains()), None)
    if target_chain is None or ligand_chain is None:
        raise ValueError("Could not find chains in input structures.")

    copy_chain(target_chain, target_chain_id)
    copy_chain(ligand_chain, ligand_chain_id)
    return builder.get_structure()


# -----------------------------
# Pose prediction (rigid body sampling only)
# -----------------------------

def random_rotation_matrix(rng: np.random.Generator) -> np.ndarray:
    """Generate a random 3x3 rotation matrix using the Haar measure."""
    u1, u2, u3 = rng.random(3)
    q1 = math.sqrt(1 - u1) * math.sin(2 * math.pi * u2)
    q2 = math.sqrt(1 - u1) * math.cos(2 * math.pi * u2)
    q3 = math.sqrt(u1) * math.sin(2 * math.pi * u3)
    q4 = math.sqrt(u1) * math.cos(2 * math.pi * u3)
    # Quaternion to rotation matrix
    R = np.array([
        [1 - 2*(q3**2 + q4**2), 2*(q2*q3 - q1*q4),     2*(q2*q4 + q1*q3)],
        [2*(q2*q3 + q1*q4),     1 - 2*(q2**2 + q4**2), 2*(q3*q4 - q1*q2)],
        [2*(q2*q4 - q1*q3),     2*(q3*q4 + q1*q2),     1 - 2*(q2**2 + q3**2)]
    ])
    return R


def center_structure(structure) -> np.ndarray:
    """Center the whole structure at the origin; return original centroid."""
    coords = []
    for atom in structure.get_atoms():
        coords.append(atom.get_coord())
    coords = np.array(coords)
    centroid = coords.mean(axis=0)
    for atom in structure.get_atoms():
        atom.set_coord(atom.get_coord() - centroid)
    return centroid


def generate_poses(
    target_pdb: str,
    ligand_pdb: str,
    out_dir: str,
    num_poses: int = 5,
    seed: int = 42,
) -> List[str]:
    """Generate diverse rigid-body poses of ligand around target.

    Returns list of output PDB file paths (each containing both proteins).
    """
    ensure_dir(out_dir)

    # Pre-center both structures for stable transforms
    target_struct = load_structure(target_pdb, 'target')
    ligand_struct = load_structure(ligand_pdb, 'ligand')
    center_structure(target_struct)
    center_structure(ligand_struct)

    # Random rigid-body placements of ligand around target with a brief contact adjustment
    rng = np.random.default_rng(seed)

    # Compute target bounding radius for translation scaling
    t_coords = np.array([a.get_coord() for a in target_struct.get_atoms()])
    t_radius = np.linalg.norm(t_coords, axis=1).max() if len(t_coords) > 0 else 20.0

    out_paths: List[str] = []
    for i in range(num_poses):
        # Make fresh copies to avoid cumulative transforms
        t = load_structure(target_pdb, f'target_{i}')
        l = load_structure(ligand_pdb, f'ligand_{i}')

        # Random rotation + moderate translation of ligand
        R = random_rotation_matrix(rng)
        # Translation scaled by target radius (place near the target surface)
        direction = rng.normal(size=3)
        direction = direction / (np.linalg.norm(direction) + 1e-8)
        distance = t_radius * (0.6 + 0.6 * rng.random())  # between 0.6R and 1.2R
        T = direction * distance
        transform_chain(l, chain_id=next(l.get_chains()).id, rotation=R, translation=T)

        # Brief contact adjustment: nudge ligand closer to target by minimizing average distance
        try:
            t_coords = np.array([a.get_coord() for a in t.get_atoms()], dtype=float)
            l_atoms = list(l.get_atoms())
            l_coords = np.array([a.get_coord() for a in l_atoms], dtype=float)
            if len(t_coords) > 0 and len(l_coords) > 0:
                # Compute vector from ligand centroid to nearest target point approximation
                from scipy.spatial import cKDTree
                kdt = cKDTree(t_coords)
                dists, idxs = kdt.query(l_coords, k=1)
                # Move ligand a small step toward target to encourage interface (avoid full overlap)
                step = min(max(float(np.median(dists)) * 0.5, 1.0), t_radius * 0.5)
                # Direction from ligand centroid to target centroid
                l_cent = l_coords.mean(axis=0)
                t_cent = t_coords[idxs].mean(axis=0) if isinstance(idxs, np.ndarray) else t_coords[int(idxs)]
                dir_vec = t_cent - l_cent
                nrm = np.linalg.norm(dir_vec) + 1e-8
                delta = (dir_vec / nrm) * step
                for atom in l_atoms:
                    atom.set_coord(atom.get_coord() + delta)
        except Exception:
            # If scipy not available or anything fails, skip adjustment
            pass

        # Merge as chains A (target) and B (ligand)
        complex_struct = merge_structures_as_chains(t, l, target_chain_id='A', ligand_chain_id='B')
        out_path = os.path.join(out_dir, f"pose_{i+1}.pdb")
        write_structure(complex_struct, out_path)
        out_paths.append(out_path)

    return out_paths


# -----------------------------
# Scoring: Binding energy (OpenMM)
# -----------------------------

def _energy_of_structure_with_openmm(structure, forcefield) -> float:
    """Compute potential energy (kJ/mol) for a Biopython structure using OpenMM GBSA.
    Uses the current atomic coordinates without minimization for speed.
    """
    # Convert to OpenMM Topology/Positions via PDBIO -> temp PDB -> read via app.PDBFile
    import tempfile

    with tempfile.TemporaryDirectory() as tmpd:
        tmp_pdb = os.path.join(tmpd, 'tmp.pdb')
        write_structure(structure, tmp_pdb)
        pdb = app.PDBFile(tmp_pdb)

    modeller = app.Modeller(pdb.topology, pdb.positions)
    # Add hydrogens at pH 7
    modeller.addHydrogens(forcefield, pH=7.0)

    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=app.NoCutoff,
        constraints=None,
        removeCMMotion=False,
    )

    # Add GBSA implicit solvent (OBC2)
    try:
        from openmm.app.internal.customgbforces import GBSAOBC2  # type: ignore
        obc = GBSAOBC2(modeller.topology, system)
        system.addForce(obc)
    except Exception:
        # If unavailable, continue without GBSA
        pass

    integrator = openmm.LangevinIntegrator(300*unit.kelvin, 1.0/unit.picoseconds, 0.002*unit.picoseconds)
    platform = openmm.Platform.getPlatformByName('Reference') if 'Reference' in [openmm.Platform.getPlatform(i).getName() for i in range(openmm.Platform.getNumPlatforms())] else openmm.Platform.getPlatformByName('CPU')

    context = openmm.Context(system, integrator, platform)
    context.setPositions(modeller.positions)
    state = context.getState(getEnergy=True)
    energy = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)

    del context, integrator, system, modeller
    return float(energy)


def compute_binding_energy_deltaE(complex_struct, chainA_id: str = 'A', chainB_id: str = 'B') -> float:
    """Compute ΔE = E_complex − (E_A + E_B) using OpenMM. Returns kJ/mol.
    """
    # Fast LJ-based proxy fallback if OpenMM is unavailable or parameterization fails
    def _proxy_interaction_energy(struct) -> float:
        """Approximate inter-chain interaction energy with a simple 12-6 LJ sum.
        Units are arbitrary but consistent across poses; good enough for ranking.
        """
        # Basic per-element LJ params (sigma in Å, epsilon in kJ/mol-like units)
        lj_params = {
            'H': (1.20, 0.02),
            'C': (1.70, 0.12),
            'N': (1.55, 0.20),
            'O': (1.52, 0.20),
            'S': (1.80, 0.25),
            'P': (1.80, 0.20),
        }

        def get_params(elem: Optional[str]):
            if not elem:
                return (1.70, 0.10)
            e = elem.upper()
            return lj_params.get(e, (1.70, 0.10))

        atomsA = list(struct[0][chainA_id].get_atoms())
        atomsB = list(struct[0][chainB_id].get_atoms())
        if not atomsA or not atomsB:
            return 0.0

        coordsA = np.array([a.get_coord() for a in atomsA], dtype=float)
        coordsB = np.array([a.get_coord() for a in atomsB], dtype=float)
        elemsA = [getattr(a.element, 'upper', lambda: str(a.element).upper())() if a.element is not None else None for a in atomsA]
        elemsB = [getattr(a.element, 'upper', lambda: str(a.element).upper())() if a.element is not None else None for a in atomsB]

        # Cutoff for interactions to keep it fast
        cutoff = 8.0
        cutoff2 = cutoff * cutoff

        E = 0.0
        # Simple O(N*M) loop is ok for small proteins; if large, could spatially bin
        for i, ra in enumerate(coordsA):
            sa, ea = get_params(elemsA[i])
            for j, rb in enumerate(coordsB):
                dr = ra - rb
                r2 = float(dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2])
                if r2 > cutoff2 or r2 < 1e-6:
                    continue
                sb, eb = get_params(elemsB[j])
                sigma = 0.5 * (sa + sb)
                epsilon = math.sqrt(ea * eb)
                inv_r = 1.0 / math.sqrt(r2)
                sr = sigma * inv_r
                sr6 = sr**6
                sr12 = sr6 * sr6
                E += 4.0 * epsilon * (sr12 - sr6)
        return float(E)

    if openmm is None or app is None or unit is None:
        # No OpenMM available; return a placeholder heuristic based on clash penalty
        return _proxy_interaction_energy(complex_struct)

    # Build structures for A, B, and complex limited to standard residues
    def structure_of_chain(struct, chain_id):
        builder = StructureBuilder.StructureBuilder()
        builder.init_structure('s')
        builder.init_model(0)
        builder.init_seg(' ')  # Initialize segment with default ID
        builder.init_chain('X')
        for residue in struct[0][chain_id]:
            hetfield, resseq, icode = residue.id
            builder.init_residue(residue.resname, hetfield, resseq, icode)
            for atom in residue:
                builder.init_atom(atom.name, atom.get_coord(), atom.bfactor, atom.occupancy, atom.altloc, atom.fullname, atom.serial_number, element=atom.element)
        return builder.get_structure()

    try:
        forcefield = app.ForceField('amber14-all.xml', 'amber14/tip3p.xml')

        # Energies
        E_complex = _energy_of_structure_with_openmm(complex_struct, forcefield)
        E_A = _energy_of_structure_with_openmm(structure_of_chain(complex_struct, chainA_id), forcefield)
        E_B = _energy_of_structure_with_openmm(structure_of_chain(complex_struct, chainB_id), forcefield)
        return float(E_complex - (E_A + E_B))
    except Exception:
        # Parameterization likely failed (e.g., nonstandard residue LIG). Use proxy.
        return _proxy_interaction_energy(complex_struct)


# -----------------------------
# Scoring: Hydrogen bonds (Biopython heuristic)
# -----------------------------

def count_hydrogen_bonds(structure, chainA_id: str = 'A', chainB_id: str = 'B', cutoff: float = 3.5) -> int:
    """Count inter-chain hydrogen bonds using a simple heavy-atom distance criterion.
    Donors/acceptors: N or O atoms; pairs across A/B within cutoff.
    Note: Angle not enforced for speed; hydrogens may be missing. This is a heuristic.
    """
    donors_elems = {'N', 'O'}
    acceptors_elems = {'N', 'O'}

    atomsA = [a for a in structure[0][chainA_id].get_atoms() if (a.element is not None and a.element.upper() in donors_elems)]
    atomsB = [a for a in structure[0][chainB_id].get_atoms() if (a.element is not None and a.element.upper() in acceptors_elems)]

    coordsA = np.array([a.get_coord() for a in atomsA], dtype=float)
    coordsB = np.array([a.get_coord() for a in atomsB], dtype=float)

    if len(coordsA) == 0 or len(coordsB) == 0:
        return 0

    # Pairwise distances; to keep fast for <100 aa, direct computation is fine
    dists = np.linalg.norm(coordsA[:, None, :] - coordsB[None, :, :], axis=-1)
    # Count pairs within cutoff; each unique donor-acceptor pair counted once
    hb_pairs = np.argwhere(dists <= cutoff)
    return int(hb_pairs.shape[0])


# -----------------------------
# Scoring: Interface area via SASA (Biopython Shrake-Rupley)
# -----------------------------

def compute_interface_area(structure, chainA_id: str = 'A', chainB_id: str = 'B') -> float:
    """Compute interface area as IA = ASA(A) + ASA(B) − ASA(A∪B)."""
    # Build isolated A and B structures
    def structure_of_chain(struct, chain_id):
        builder = StructureBuilder.StructureBuilder()
        builder.init_structure('s')
        builder.init_model(0)
        builder.init_seg(' ')  # Initialize segment with default ID
        builder.init_chain('X')
        for residue in struct[0][chain_id]:
            hetfield, resseq, icode = residue.id
            builder.init_residue(residue.resname, hetfield, resseq, icode)
            for atom in residue:
                builder.init_atom(atom.name, atom.get_coord(), atom.bfactor, atom.occupancy, atom.altloc, atom.fullname, atom.serial_number, element=atom.element)
        return builder.get_structure()

    # Compute SASA for A, B, and complex
    sr = ShrakeRupley(n_points=960)  # finer grid for stability

    A = structure_of_chain(structure, chainA_id)
    B = structure_of_chain(structure, chainB_id)

    # Complex structure with both chains already present: reuse
    # Ensure SASA (.sasa) attributes get computed
    sr.compute(A, level='A')
    asa_A = sum(atom.sasa for atom in A.get_atoms() if hasattr(atom, 'sasa'))

    sr.compute(B, level='A')
    asa_B = sum(atom.sasa for atom in B.get_atoms() if hasattr(atom, 'sasa'))

    sr.compute(structure, level='A')
    asa_complex = sum(atom.sasa for atom in structure.get_atoms() if hasattr(atom, 'sasa'))

    interface_area = float(asa_A + asa_B - asa_complex)
    return max(interface_area, 0.0)


# -----------------------------
# Normalization and scoring
# -----------------------------

def normalize(values: List[float], invert: bool = False) -> List[float]:
    arr = np.array(values, dtype=float)
    vmin = float(np.min(arr))
    vmax = float(np.max(arr))
    if math.isclose(vmin, vmax):
        # All equal -> neutral 0.5
        return [0.5] * len(values)
    norm = (arr - vmin) / (vmax - vmin)
    if invert:
        norm = 1.0 - norm
    return norm.tolist()


# -----------------------------
# Main orchestration
# -----------------------------

def run_pipeline(target_pdb: str, ligand_pdb: str, out_dir: str, num_poses: int = 5, seed: int = 42) -> pd.DataFrame:
    ensure_dir(out_dir)

    # 1) Preprocess PDBs
    cleaned_target = os.path.join(out_dir, 'target_cleaned.pdb')
    cleaned_ligand = os.path.join(out_dir, 'ligand_cleaned.pdb')
    print("[INFO] Preprocessing input PDBs with pdbfixer (if available)…")
    preprocess_pdb(target_pdb, cleaned_target)
    preprocess_pdb(ligand_pdb, cleaned_ligand)

    # 2) Pose prediction using rigid-body sampling
    print("[INFO] Generating poses using rigid-body sampling…")
    pose_paths = generate_poses(cleaned_target, cleaned_ligand, out_dir=out_dir, num_poses=num_poses, seed=seed)

    # 3) Metrics for each pose
    metrics: List[PoseMetrics] = []
    print("[INFO] Scoring poses (binding energy, hydrogen bonds, interface area)…")
    for i, pose_path in enumerate(pose_paths, start=1):
        struct = load_structure(pose_path, f'pose_{i}')
        try:
            dE = compute_binding_energy_deltaE(struct, chainA_id='A', chainB_id='B')
        except Exception as e:
            print(f"[WARN] Binding energy failed for pose {i}: {e}. Using 0.0.")
            dE = 0.0
        try:
            hb = count_hydrogen_bonds(struct, chainA_id='A', chainB_id='B', cutoff=3.5)
        except Exception as e:
            print(f"[WARN] H-bond count failed for pose {i}: {e}. Using 0.")
            hb = 0
        try:
            ia = compute_interface_area(struct, chainA_id='A', chainB_id='B')
        except Exception as e:
            print(f"[WARN] Interface area failed for pose {i}: {e}. Using 0.0.")
            ia = 0.0
        metrics.append(PoseMetrics(pose_id=i, binding_energy=dE, hydrogen_bonds=hb, interface_area=ia))

    # 4) Normalize metrics and compute final score
    energies = [m.binding_energy for m in metrics]
    hbonds = [m.hydrogen_bonds for m in metrics]
    iface = [m.interface_area for m in metrics]

    nrg_norm = normalize(energies, invert=True)   # lower is better -> invert
    hb_norm = normalize(hbonds, invert=False)
    ia_norm = normalize(iface, invert=False)

    for m, a, b, c in zip(metrics, nrg_norm, hb_norm, ia_norm):
        m.final_score = float((a + b + c) / 3.0)

    # 5) Sort and print table
    metrics_sorted = sorted(metrics, key=lambda x: x.final_score if x.final_score is not None else -1.0, reverse=True)

    df = pd.DataFrame([
        {
            'Pose ID': m.pose_id,
            'Pose Path': pose_paths[m.pose_id - 1],
            'Binding Energy (kJ/mol)': m.binding_energy,
            'Hydrogen Bonds': m.hydrogen_bonds,
            'Interface Area (Å^2)': m.interface_area,
            'Final Score': m.final_score,
        }
        for m in metrics_sorted
    ])

    # 6) Save best pose
    best_pose_src = df.iloc[0]['Pose Path']
    best_pose_dst = os.path.join(out_dir, 'best_pose.pdb')
    safe_copy(best_pose_src, best_pose_dst)

    # Print ranking table (without paths)
    print("\nPose ranking (best first):")
    print(df.drop(columns=['Pose Path']).to_string(index=False, justify='center', float_format=lambda x: f"{x: .3f}"))

    return df


def main():
    parser = argparse.ArgumentParser(description="Predict and rank protein–protein binding poses using rigid-body sampling and normalized scoring.")
    parser.add_argument('--target', required=True, help='Path to target (receptor) PDB file')
    parser.add_argument('--ligand', required=True, help='Path to ligand (moving protein) PDB file')
    parser.add_argument('--outdir', default='outputs', help='Output directory for poses and results')
    parser.add_argument('--num_poses', type=int, default=5, help='Number of poses to generate')
    parser.add_argument('--seed', type=int, default=42, help='Random seed for rigid-body sampler')
    args = parser.parse_args()

    df = run_pipeline(args.target, args.ligand, args.outdir, num_poses=args.num_poses, seed=args.seed)

    # Save CSV summary for reproducibility
    csv_path = os.path.join(args.outdir, 'pose_ranking.csv')
    df.to_csv(csv_path, index=False)
    print(f"\n[INFO] Saved summary to: {csv_path}")
    print(f"[INFO] Best pose saved to: {os.path.join(args.outdir, 'best_pose.pdb')}")


if __name__ == '__main__':
    main()
