from rdkit import Chem
from rdkit.Chem import Descriptors
import os

# Function to calculate synthetic accessibility score
# Placeholder function for SAScore calculation
# Replace with actual SAScore calculation if available
def calculate_sascore(mol):
    return Descriptors.MolWt(mol)  # Using molecular weight as a placeholder

# Function to select top N molecules based on SAScore
def select_top_molecules(sdf_file, top_n=20):
    suppl = Chem.SDMolSupplier(sdf_file)
    
    # Calculate SAScore for each molecule
    scores = []
    for mol in suppl:
        if mol is None:
            continue
        score = calculate_sascore(mol)
        scores.append((score, mol))

    # Sort molecules by SAScore
    scores.sort(reverse=True, key=lambda x: x[0])

    # Select top N molecules
    top_molecules = scores[:top_n]
    return [mol for _, mol in top_molecules]

# Import the function to read XYZ files
from read_xyz_files import read_xyz_files

# Example usage
molecules = read_xyz_files('scaffold_1_30_xyz')
top_molecules = select_top_molecules(molecules)
print(f"Selected {len(top_molecules)} top molecules based on SAScore.") 