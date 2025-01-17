from rdkit import Chem
import os

# Function to read XYZ files and convert to RDKit molecules
def read_xyz_files(xyz_dir):
    molecules = []
    for filename in os.listdir(xyz_dir):
        if filename.endswith('.xyz'):
            filepath = os.path.join(xyz_dir, filename)
            with open(filepath, 'r') as file:
                xyz_data = file.read()
                mol = Chem.MolFromXYZBlock(xyz_data)
                if mol is not None:
                    molecules.append(mol)
    return molecules

# Example usage
molecules = read_xyz_files('scaffold_1_30_xyz')
print(f"Read {len(molecules)} molecules from XYZ files.") 