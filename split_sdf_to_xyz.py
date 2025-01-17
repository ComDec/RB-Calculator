from rdkit import Chem
from rdkit.Chem import AllChem
import os

# Function to split SDF into XYZ files
def split_sdf_to_xyz(sdf_file, output_dir):
    # Create output directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Read molecules from SDF file
    suppl = Chem.SDMolSupplier(sdf_file)
    
    for i, mol in enumerate(suppl):
        if mol is None:
            continue
        
        # Generate 3D coordinates
        # AllChem.EmbedMolecule(mol)
        # AllChem.UFFOptimizeMolecule(mol)

        # Write to XYZ file
        xyz_filename = os.path.join(output_dir, f'molecule_{i+1}.xyz')
        with open(xyz_filename, 'w') as xyz_file:
            xyz_file.write(Chem.MolToXYZBlock(mol))

# Example usage
split_sdf_to_xyz('/home/wangxi/project/ligands_library/scaffold_1_30.sdf', 'scaffold_1_30_xyz') 