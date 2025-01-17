import os
import sys
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from pyscf import gto, dft
from run_dft_refinement import cal_mol_energy_via_dft  # Assuming this function is available or defined within the same script.

def read_smiles_csv(file_path):
    """
    Read a CSV file containing SMILES strings and return a list of RDKit molecule objects.
    """
    df = pd.read_csv(file_path)
    smiles_list = df['SMILES']
    mols = [Chem.MolFromSmiles(smiles) for smiles in smiles_list]
    # Add Hydrogens and generate 3D coordinates
    mols = [Chem.AddHs(mol) for mol in mols if mol is not None]
    [AllChem.EmbedMolecule(mol, randomSeed=42) for mol in mols if mol is not None]
    [AllChem.UFFOptimizeMolecule(mol, maxIters=200) for mol in mols if mol is not None]
    return mols

def optimize_with_xtb(mol, save_folder, mol_name):
    """
    Optimize the structure using xtb by saving an XYZ file and running xtb.
    """
    xyz_file = os.path.join(save_folder, f"{mol_name}.xyz")
    with open(xyz_file, "w") as f:
        f.write(Chem.MolToXYZBlock(mol))

    # Setup for xtb optimization
    os.makedirs(save_folder, exist_ok=True)
    xtb_xyz_file = os.path.join(save_folder, f"xtb_{mol_name}.xyz")
    save_folder_path = os.path.abspath(save_folder)
    xtb_cmd = f" cd {save_folder_path} && xtb {mol_name}.xyz --opt --input {mol_name}.xyz --ohess --gfn 2 --chrg 0 --uhf 0 --speed moderate --gbsa none > xtb_output.log"
    os.system(xtb_cmd)

    # Rename optimized file if it exists
    if os.path.exists(os.path.join(save_folder, "xtbopt.xyz")):
        os.rename(os.path.join(save_folder, "xtbopt.xyz"), xtb_xyz_file)

def refine_with_pyscf(mol, save_folder, mol_name):
    """
    Refine the energies using pyscf.
    """
    mol_atoms = [(atom.GetSymbol(), mol.GetConformer().GetAtomPosition(i)) for i, atom in enumerate(mol.GetAtoms())]
    pyscf_mol = gto.M(atom=mol_atoms, basis='def2-SVP')
    pyscf_mf = dft.RKS(pyscf_mol)
    pyscf_mf.xc = 'M06-2X'
    energy = pyscf_mf.kernel()

    # Write the optimized XYZ file
    pyscf_xyz_file = os.path.join(save_folder, f"pyscf_{mol_name}.xyz")
    with open(pyscf_xyz_file, "w") as xyz_file:
        xyz_file.write(Chem.MolToXYZBlock(mol))

    with open(os.path.join(save_folder, f"{mol_name}_pyscf_energy.txt"), "w") as f:
        f.write(f"{energy}\n")

def main(smiles_file, output_dir):
    mols = read_smiles_csv(smiles_file)
    for i, mol in enumerate(mols):
        mol_name = f"molecule_{i}"
        mol_dir = os.path.join(output_dir, mol_name)
        os.makedirs(mol_dir, exist_ok=True)

        optimize_with_xtb(mol, mol_dir, mol_name)
        refine_with_pyscf(mol, mol_dir, mol_name)

if __name__ == "__main__":
    smiles_csv = sys.argv[1]
    output_directory = sys.argv[2]
    main(smiles_csv, output_directory)
