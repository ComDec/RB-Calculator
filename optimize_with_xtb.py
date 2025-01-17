import os
import subprocess
from rdkit import Chem
from calculate_sascore import select_top_molecules
from read_xyz_files import read_xyz_files

# Function to optimize molecules with xTB
def optimize_with_xtb(molecules, base_dir='optimized_molecules'):
    # Create base directory if it doesn't exist
    if not os.path.exists(base_dir):
        os.makedirs(base_dir)

    for i, mol in enumerate(molecules):
        # Create a directory for each molecule
        mol_dir = os.path.join(base_dir, f'molecule_{i+1}')
        if not os.path.exists(mol_dir):
            os.makedirs(mol_dir)

        # Write molecule to XYZ file
        xyz_filename = os.path.join(mol_dir, 'input.xyz')
        with open(xyz_filename, 'w') as xyz_file:
            xyz_file.write(Chem.MolToXYZBlock(mol))

        # Run xTB optimization
        subprocess.run(['xtb', xyz_filename, '--opt'], cwd=mol_dir)

# Example usage
molecules = read_xyz_files('scaffold_1_30_xyz')
top_molecules = select_top_molecules(molecules)
optimize_with_xtb(top_molecules) 