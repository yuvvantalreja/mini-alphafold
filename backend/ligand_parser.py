import os
import tempfile
import json
from typing import Dict, List, Tuple, Optional, Union
from pathlib import Path

try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors
    from rdkit.Chem import AllChem
    from meeko import MoleculePreparation, PDBQTMolecule
    import pubchempy as pcp
except ImportError as e:
    print(f"Warning: Some molecular libraries not available: {e}")
    print("Please install: pip install rdkit meeko pubchempy")

class LigandParser:
    """Parser and processor for ligand molecules in various formats"""
    
    def __init__(self):
        self.supported_formats = {
            'sdf': 'SD File',
            'mol': 'MDL Molfile',
            'mol2': 'Mol2 File',
            'smi': 'SMILES',
            'smiles': 'SMILES'
        }
        
        # Drug-like property thresholds (Lipinski's Rule of Five)
        self.lipinski_thresholds = {
            'mw': 500.0,          # Molecular weight
            'logp': 5.0,          # Partition coefficient
            'hbd': 5,             # Hydrogen bond donors
            'hba': 10,            # Hydrogen bond acceptors
            'rotatable_bonds': 10  # Rotatable bonds
        }

    def parse_ligand_file(self, filepath: str) -> Dict:
        """Parse a ligand file and return structured data"""
        file_ext = Path(filepath).suffix.lower().lstrip('.')
        
        if file_ext not in self.supported_formats:
            raise ValueError(f"Unsupported file format: {file_ext}")
        
        try:
            if file_ext in ['sdf']:
                return self._parse_sdf(filepath)
            elif file_ext in ['mol']:
                return self._parse_mol(filepath)
            elif file_ext in ['mol2']:
                return self._parse_mol2(filepath)
            elif file_ext in ['smi', 'smiles']:
                return self._parse_smiles(filepath)
            else:
                raise ValueError(f"Parser not implemented for {file_ext}")
                
        except Exception as e:
            raise Exception(f"Error parsing ligand file: {str(e)}")

    def parse_smiles_string(self, smiles: str, name: str = "ligand") -> Dict:
        """Parse a SMILES string directly"""
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                raise ValueError("Invalid SMILES string")
            
            return self._process_molecule(mol, name, smiles)
            
        except Exception as e:
            raise Exception(f"Error parsing SMILES: {str(e)}")

    def fetch_pubchem_compound(self, identifier: str, identifier_type: str = "name") -> Dict:
        """Fetch compound from PubChem database"""
        try:
            if identifier_type == "name":
                compounds = pcp.get_compounds(identifier, 'name')
            elif identifier_type == "cid":
                compounds = pcp.get_compounds(int(identifier), 'cid')
            elif identifier_type == "smiles":
                compounds = pcp.get_compounds(identifier, 'smiles')
            else:
                raise ValueError("Identifier type must be 'name', 'cid', or 'smiles'")
            
            if not compounds:
                raise ValueError(f"No compounds found for {identifier}")
            
            compound = compounds[0]  # Take the first match
            smiles = compound.canonical_smiles
            
            if not smiles:
                raise ValueError("No SMILES available for this compound")
            
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                raise ValueError("Invalid SMILES from PubChem")
            
            return self._process_molecule(mol, compound.iupac_name or identifier, smiles, compound)
            
        except Exception as e:
            raise Exception(f"Error fetching from PubChem: {str(e)}")

    def _parse_sdf(self, filepath: str) -> Dict:
        """Parse SDF file"""
        supplier = Chem.SDMolSupplier(filepath)
        molecules = []
        
        for i, mol in enumerate(supplier):
            if mol is not None:
                name = mol.GetProp('_Name') if mol.HasProp('_Name') else f"ligand_{i+1}"
                mol_data = self._process_molecule(mol, name)
                molecules.append(mol_data)
        
        if not molecules:
            raise ValueError("No valid molecules found in SDF file")
        
        return {
            'format': 'sdf',
            'molecule_count': len(molecules),
            'molecules': molecules,
            'primary_molecule': molecules[0]  # Use first molecule as primary
        }

    def _parse_mol(self, filepath: str) -> Dict:
        """Parse MOL file"""
        mol = Chem.MolFromMolFile(filepath)
        if mol is None:
            raise ValueError("Invalid MOL file")
        
        name = os.path.splitext(os.path.basename(filepath))[0]
        mol_data = self._process_molecule(mol, name)
        
        return {
            'format': 'mol',
            'molecule_count': 1,
            'molecules': [mol_data],
            'primary_molecule': mol_data
        }

    def _parse_mol2(self, filepath: str) -> Dict:
        """Parse MOL2 file"""
        # RDKit has limited MOL2 support, try to read it
        mol = Chem.MolFromMol2File(filepath)
        if mol is None:
            raise ValueError("Invalid MOL2 file or format not supported")
        
        name = os.path.splitext(os.path.basename(filepath))[0]
        mol_data = self._process_molecule(mol, name)
        
        return {
            'format': 'mol2',
            'molecule_count': 1,
            'molecules': [mol_data],
            'primary_molecule': mol_data
        }

    def _parse_smiles(self, filepath: str) -> Dict:
        """Parse SMILES file"""
        molecules = []
        
        with open(filepath, 'r') as f:
            for i, line in enumerate(f):
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                
                parts = line.split()
                smiles = parts[0]
                name = parts[1] if len(parts) > 1 else f"ligand_{i+1}"
                
                mol = Chem.MolFromSmiles(smiles)
                if mol is not None:
                    mol_data = self._process_molecule(mol, name, smiles)
                    molecules.append(mol_data)
        
        if not molecules:
            raise ValueError("No valid SMILES found in file")
        
        return {
            'format': 'smiles',
            'molecule_count': len(molecules),
            'molecules': molecules,
            'primary_molecule': molecules[0]
        }

    def _process_molecule(self, mol, name: str, smiles: str = None, pubchem_data=None) -> Dict:
        """Process RDKit molecule and extract comprehensive data"""
        
        # Generate 3D coordinates if not present
        if mol.GetNumConformers() == 0:
            mol = Chem.AddHs(mol)
            AllChem.EmbedMolecule(mol, randomSeed=42)
            AllChem.UFFOptimizeMolecule(mol)
        
        # Get SMILES if not provided
        if smiles is None:
            smiles = Chem.MolToSmiles(mol)
        
        # Calculate molecular properties
        properties = self._calculate_properties(mol)
        
        # Get atom information
        atoms = self._extract_atoms(mol)
        
        # Get bond information
        bonds = self._extract_bonds(mol)
        
        # Drug-likeness assessment
        drug_like = self._assess_drug_likeness(properties)
        
        # Prepare for docking (convert to PDBQT)
        pdbqt_data = self._prepare_for_docking(mol)
        
        molecule_data = {
            'name': name,
            'smiles': smiles,
            'formula': Chem.rdMolDescriptors.CalcMolFormula(mol),
            'properties': properties,
            'drug_likeness': drug_like,
            'atoms': atoms,
            'bonds': bonds,
            'conformers': self._extract_conformers(mol),
            'pdbqt': pdbqt_data,
            'center': self._calculate_center(atoms),
            'stats': {
                'atom_count': mol.GetNumAtoms(),
                'bond_count': mol.GetNumBonds(),
                'ring_count': rdMolDescriptors.CalcNumRings(mol),
                'aromatic_ring_count': rdMolDescriptors.CalcNumAromaticRings(mol)
            }
        }
        
        # Add PubChem data if available
        if pubchem_data:
            molecule_data['pubchem'] = {
                'cid': pubchem_data.cid,
                'iupac_name': pubchem_data.iupac_name,
                'molecular_formula': pubchem_data.molecular_formula,
                'molecular_weight': pubchem_data.molecular_weight,
                'synonyms': pubchem_data.synonyms[:5] if pubchem_data.synonyms else []
            }
        
        return molecule_data

    def _calculate_properties(self, mol) -> Dict:
        """Calculate molecular properties"""
        return {
            'molecular_weight': Descriptors.MolWt(mol),
            'logp': Crippen.MolLogP(mol),
            'hbd': rdMolDescriptors.CalcNumHBD(mol),  # Hydrogen bond donors
            'hba': rdMolDescriptors.CalcNumHBA(mol),  # Hydrogen bond acceptors
            'tpsa': rdMolDescriptors.CalcTPSA(mol),   # Topological polar surface area
            'rotatable_bonds': rdMolDescriptors.CalcNumRotatableBonds(mol),
            'formal_charge': Chem.rdmolops.GetFormalCharge(mol),
            'num_heavy_atoms': mol.GetNumHeavyAtoms(),
            'num_heteroatoms': rdMolDescriptors.CalcNumHeteroatoms(mol),
            'fraction_csp3': rdMolDescriptors.CalcFractionCsp3(mol),
            'num_aliphatic_rings': rdMolDescriptors.CalcNumAliphaticRings(mol),
            'num_saturated_rings': rdMolDescriptors.CalcNumSaturatedRings(mol)
        }

    def _extract_atoms(self, mol) -> List[Dict]:
        """Extract atom information"""
        atoms = []
        conf = mol.GetConformer(0)
        
        for i, atom in enumerate(mol.GetAtoms()):
            pos = conf.GetAtomPosition(i)
            atoms.append({
                'index': i,
                'symbol': atom.GetSymbol(),
                'atomic_num': atom.GetAtomicNum(),
                'formal_charge': atom.GetFormalCharge(),
                'hybridization': str(atom.GetHybridization()),
                'is_aromatic': atom.GetIsAromatic(),
                'x': pos.x,
                'y': pos.y,
                'z': pos.z,
                'radius': self._get_vdw_radius(atom.GetSymbol()),
                'color': self._get_atom_color(atom.GetSymbol())
            })
        
        return atoms

    def _extract_bonds(self, mol) -> List[Dict]:
        """Extract bond information"""
        bonds = []
        for bond in mol.GetBonds():
            bonds.append({
                'atom1': bond.GetBeginAtomIdx(),
                'atom2': bond.GetEndAtomIdx(),
                'bond_type': str(bond.GetBondType()),
                'is_aromatic': bond.GetIsAromatic(),
                'order': int(bond.GetBondTypeAsDouble())
            })
        return bonds

    def _extract_conformers(self, mol) -> List[Dict]:
        """Extract conformer information"""
        conformers = []
        for i in range(mol.GetNumConformers()):
            conf = mol.GetConformer(i)
            coords = []
            for j in range(mol.GetNumAtoms()):
                pos = conf.GetAtomPosition(j)
                coords.append([pos.x, pos.y, pos.z])
            conformers.append({
                'id': i,
                'coordinates': coords
            })
        return conformers

    def _assess_drug_likeness(self, properties: Dict) -> Dict:
        """Assess drug-likeness using Lipinski's Rule of Five"""
        violations = 0
        checks = {}
        
        # Check each Lipinski rule
        if properties['molecular_weight'] > self.lipinski_thresholds['mw']:
            violations += 1
            checks['molecular_weight'] = False
        else:
            checks['molecular_weight'] = True
            
        if properties['logp'] > self.lipinski_thresholds['logp']:
            violations += 1
            checks['logp'] = False
        else:
            checks['logp'] = True
            
        if properties['hbd'] > self.lipinski_thresholds['hbd']:
            violations += 1
            checks['hbd'] = False
        else:
            checks['hbd'] = True
            
        if properties['hba'] > self.lipinski_thresholds['hba']:
            violations += 1
            checks['hba'] = False
        else:
            checks['hba'] = True
        
        return {
            'lipinski_violations': violations,
            'lipinski_compliant': violations <= 1,  # Allow 1 violation
            'checks': checks,
            'drug_like_score': max(0, 1 - (violations / 4))  # Simple scoring
        }

    def _prepare_for_docking(self, mol) -> str:
        """Prepare molecule for docking by converting to PDBQT format"""
        try:
            # Use meeko to prepare molecule for AutoDock
            preparator = MoleculePreparation()
            mol_prepared = preparator.prepare(mol)
            
            # Convert to PDBQT string
            pdbqt_string = mol_prepared.write_pdbqt_string()
            return pdbqt_string
            
        except Exception as e:
            print(f"Warning: Could not prepare molecule for docking: {e}")
            return ""

    def _calculate_center(self, atoms: List[Dict]) -> List[float]:
        """Calculate geometric center of molecule"""
        if not atoms:
            return [0, 0, 0]
        
        x_coords = [atom['x'] for atom in atoms]
        y_coords = [atom['y'] for atom in atoms]
        z_coords = [atom['z'] for atom in atoms]
        
        return [
            sum(x_coords) / len(x_coords),
            sum(y_coords) / len(y_coords),
            sum(z_coords) / len(z_coords)
        ]

    def _get_vdw_radius(self, symbol: str) -> float:
        """Get van der Waals radius for atom"""
        radii = {
            'H': 1.20, 'C': 1.70, 'N': 1.55, 'O': 1.52, 'F': 1.47,
            'P': 1.80, 'S': 1.80, 'Cl': 1.75, 'Br': 1.85, 'I': 1.98
        }
        return radii.get(symbol, 1.5)

    def _get_atom_color(self, symbol: str) -> int:
        """Get CPK color for atom"""
        colors = {
            'C': 0x909090, 'N': 0x3050F8, 'O': 0xFF0D0D, 'H': 0xFFFFFF,
            'S': 0xFFFF30, 'P': 0xFF8000, 'F': 0x90E050, 'Cl': 0x1FF01F,
            'Br': 0xA62929, 'I': 0x940094
        }
        return colors.get(symbol, 0x808080)

    def save_as_pdbqt(self, ligand_data: Dict, output_path: str) -> str:
        """Save ligand as PDBQT file for docking"""
        pdbqt_content = ligand_data.get('pdbqt', '')
        if not pdbqt_content:
            raise ValueError("No PDBQT data available for this ligand")
        
        with open(output_path, 'w') as f:
            f.write(pdbqt_content)
        
        return output_path

    def generate_conformers(self, mol, num_conformers: int = 10) -> List[Dict]:
        """Generate multiple conformers for a molecule"""
        mol = Chem.AddHs(mol)
        
        # Generate conformers
        confIds = AllChem.EmbedMultipleConfs(
            mol, 
            numConfs=num_conformers, 
            randomSeed=42,
            clearConfs=True
        )
        
        # Optimize conformers
        for confId in confIds:
            AllChem.UFFOptimizeMolecule(mol, confId=confId)
        
        return self._extract_conformers(mol)
