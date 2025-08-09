import re
import json
from collections import defaultdict
from typing import Dict, List, Tuple, Optional

class PDBParser:
    """Parser for PDB (Protein Data Bank) files"""
    
    def __init__(self):
        self.amino_acids = {
            'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
            'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'
        }
        
        self.nucleotides = {'A', 'T', 'G', 'C', 'U', 'DA', 'DT', 'DG', 'DC'}
        
        # Atom colors based on CPK coloring scheme
        self.atom_colors = {
            'C': 0x909090,   # Carbon - Gray
            'N': 0x3050F8,   # Nitrogen - Blue
            'O': 0xFF0D0D,   # Oxygen - Red
            'S': 0xFFFF30,   # Sulfur - Yellow
            'P': 0xFF8000,   # Phosphorus - Orange
            'H': 0xFFFFFF,   # Hydrogen - White
            'CA': 0x40E0D0,  # Calcium - Turquoise
            'FE': 0xE06633,  # Iron - Orange-red
            'MG': 0x8AFF00,  # Magnesium - Green
            'ZN': 0x7D80B0,  # Zinc - Blue-gray
        }
        
        # Van der Waals radii (in Angstroms)
        self.vdw_radii = {
            'C': 1.70, 'N': 1.55, 'O': 1.52, 'S': 1.80, 'P': 1.80,
            'H': 1.20, 'CA': 2.31, 'FE': 2.00, 'MG': 1.73, 'ZN': 1.39
        }

    def parse_pdb_file(self, filepath: str) -> Dict:
        """Parse a PDB file and return structured data for visualization"""
        
        atoms = []
        bonds = []
        chains = defaultdict(list)
        residues = defaultdict(list)
        header_info = {}
        
        try:
            with open(filepath, 'r') as file:
                lines = file.readlines()
            
            # Parse header information
            header_info = self._parse_header(lines)
            
            # Parse atoms
            atom_dict = {}
            for line in lines:
                if line.startswith('ATOM') or line.startswith('HETATM'):
                    atom = self._parse_atom_line(line)
                    if atom:
                        atoms.append(atom)
                        atom_dict[atom['serial']] = atom
                        chains[atom['chain']].append(atom)
                        residues[f"{atom['chain']}_{atom['res_seq']}"].append(atom)
            
            # Generate bonds
            bonds = self._generate_bonds(atoms, atom_dict)
            
            # Calculate secondary structure (basic implementation)
            secondary_structure = self._calculate_secondary_structure(residues)
            
            # Calculate center and bounding box
            center, bbox = self._calculate_bounds(atoms)
            
            return {
                'header': header_info,
                'atoms': atoms,
                'bonds': bonds,
                'chains': dict(chains),
                'residues': dict(residues),
                'secondary_structure': secondary_structure,
                'center': center,
                'bounding_box': bbox,
                'stats': {
                    'atom_count': len(atoms),
                    'chain_count': len(chains),
                    'residue_count': len(residues),
                    'bond_count': len(bonds)
                }
            }
            
        except Exception as e:
            raise Exception(f"Error parsing PDB file: {str(e)}")
    
    def _parse_header(self, lines: List[str]) -> Dict:
        """Parse header information from PDB file"""
        header = {
            'title': '',
            'classification': '',
            'deposition_date': '',
            'keywords': '',
            'organism': '',
            'authors': '',
            'resolution': None,
            'experiment_type': ''
        }
        
        for line in lines:
            if line.startswith('TITLE'):
                header['title'] += line[10:].strip() + ' '
            elif line.startswith('HEADER'):
                header['classification'] = line[10:50].strip()
                header['deposition_date'] = line[50:59].strip()
            elif line.startswith('KEYWDS'):
                header['keywords'] += line[10:].strip() + ' '
            elif line.startswith('SOURCE'):
                if 'ORGANISM_SCIENTIFIC:' in line:
                    header['organism'] = line.split('ORGANISM_SCIENTIFIC:')[1].strip().rstrip(';')
            elif line.startswith('AUTHOR'):
                header['authors'] += line[10:].strip() + ' '
            elif line.startswith('REMARK   2 RESOLUTION'):
                try:
                    res_match = re.search(r'(\d+\.\d+)', line)
                    if res_match:
                        header['resolution'] = float(res_match.group(1))
                except ValueError:
                    pass
            elif line.startswith('EXPDTA'):
                header['experiment_type'] = line[10:].strip()
        
        # Clean up strings
        for key in ['title', 'keywords', 'authors']:
            header[key] = header[key].strip()
        
        return header
    
    def _parse_atom_line(self, line: str) -> Optional[Dict]:
        """Parse an ATOM or HETATM line from PDB file"""
        try:
            if len(line) < 54:
                return None
            
            atom = {
                'serial': int(line[6:11].strip()),
                'name': line[12:16].strip(),
                'alt_loc': line[16].strip(),
                'res_name': line[17:20].strip(),
                'chain': line[21].strip(),
                'res_seq': int(line[22:26].strip()),
                'x': float(line[30:38].strip()),
                'y': float(line[38:46].strip()),
                'z': float(line[46:54].strip()),
                'occupancy': float(line[54:60].strip()) if len(line) > 60 and line[54:60].strip() else 1.0,
                'temp_factor': float(line[60:66].strip()) if len(line) > 66 and line[60:66].strip() else 0.0,
                'element': line[76:78].strip() if len(line) > 78 else line[12:16].strip()[:2].strip('0123456789'),
                'is_hetatm': line.startswith('HETATM')
            }
            
            # Clean up element symbol
            atom['element'] = atom['element'].upper()
            if len(atom['element']) > 1:
                atom['element'] = atom['element'][0] + atom['element'][1:].lower()
            
            # Add visualization properties
            atom['color'] = self.atom_colors.get(atom['element'], 0x808080)
            atom['radius'] = self.vdw_radii.get(atom['element'], 1.5)
            atom['is_backbone'] = atom['name'] in ['N', 'CA', 'C', 'O']
            atom['is_amino_acid'] = atom['res_name'] in self.amino_acids
            atom['is_nucleotide'] = atom['res_name'] in self.nucleotides
            
            return atom
            
        except (ValueError, IndexError) as e:
            print(f"Error parsing line: {line.strip()}")
            return None
    
    def _generate_bonds(self, atoms: List[Dict], atom_dict: Dict) -> List[Dict]:
        """Generate bonds based on distance and chemical rules"""
        bonds = []
        bond_cutoffs = {
            ('C', 'C'): 1.8,
            ('C', 'N'): 1.7,
            ('C', 'O'): 1.6,
            ('C', 'S'): 2.1,
            ('N', 'N'): 1.6,
            ('N', 'O'): 1.5,
            ('O', 'O'): 1.4,
            ('S', 'S'): 2.2,
            ('P', 'O'): 1.8,
        }
        
        # Group atoms by residue for efficiency
        residue_atoms = defaultdict(list)
        for atom in atoms:
            residue_key = f"{atom['chain']}_{atom['res_seq']}"
            residue_atoms[residue_key].append(atom)
        
        # Generate bonds within residues and between adjacent residues
        processed_pairs = set()
        
        for residue_key, res_atoms in residue_atoms.items():
            # Bonds within residue
            for i, atom1 in enumerate(res_atoms):
                for j, atom2 in enumerate(res_atoms[i+1:], i+1):
                    pair_key = tuple(sorted([atom1['serial'], atom2['serial']]))
                    if pair_key in processed_pairs:
                        continue
                    
                    if self._should_bond(atom1, atom2, bond_cutoffs):
                        bonds.append({
                            'atom1': atom1['serial'],
                            'atom2': atom2['serial'],
                            'type': 'covalent',
                            'order': 1
                        })
                        processed_pairs.add(pair_key)
            
            # Bonds to next residue (backbone)
            chain, res_seq = residue_key.split('_')
            next_residue_key = f"{chain}_{int(res_seq) + 1}"
            
            if next_residue_key in residue_atoms:
                current_c = next((a for a in res_atoms if a['name'] == 'C'), None)
                next_n = next((a for a in residue_atoms[next_residue_key] if a['name'] == 'N'), None)
                
                if current_c and next_n:
                    pair_key = tuple(sorted([current_c['serial'], next_n['serial']]))
                    if pair_key not in processed_pairs:
                        if self._should_bond(current_c, next_n, bond_cutoffs):
                            bonds.append({
                                'atom1': current_c['serial'],
                                'atom2': next_n['serial'],
                                'type': 'peptide',
                                'order': 1
                            })
                            processed_pairs.add(pair_key)
        
        return bonds
    
    def _should_bond(self, atom1: Dict, atom2: Dict, bond_cutoffs: Dict) -> bool:
        """Determine if two atoms should be bonded"""
        # Skip hydrogen bonds for now
        if atom1['element'] == 'H' or atom2['element'] == 'H':
            return False
        
        # Calculate distance
        dx = atom1['x'] - atom2['x']
        dy = atom1['y'] - atom2['y']
        dz = atom1['z'] - atom2['z']
        distance = (dx*dx + dy*dy + dz*dz)**0.5
        
        # Get bond cutoff
        elements = tuple(sorted([atom1['element'], atom2['element']]))
        cutoff = bond_cutoffs.get(elements, 2.0)
        
        return distance <= cutoff
    
    def _calculate_secondary_structure(self, residues: Dict) -> Dict:
        """Basic secondary structure assignment (placeholder)"""
        # This is a simplified implementation
        # Real secondary structure assignment would use DSSP algorithm
        structure = {}
        
        for residue_key, atoms in residues.items():
            # Simple heuristic based on phi/psi angles (not implemented here)
            # For now, assign random structure types for demonstration
            structure[residue_key] = 'coil'  # Could be 'helix', 'sheet', or 'coil'
        
        return structure
    
    def _calculate_bounds(self, atoms: List[Dict]) -> Tuple[List[float], Dict]:
        """Calculate center point and bounding box of the molecule"""
        if not atoms:
            return [0, 0, 0], {'min': [0, 0, 0], 'max': [0, 0, 0]}
        
        xs = [atom['x'] for atom in atoms]
        ys = [atom['y'] for atom in atoms]
        zs = [atom['z'] for atom in atoms]
        
        center = [
            sum(xs) / len(xs),
            sum(ys) / len(ys),
            sum(zs) / len(zs)
        ]
        
        bbox = {
            'min': [min(xs), min(ys), min(zs)],
            'max': [max(xs), max(ys), max(zs)]
        }
        
        return center, bbox
