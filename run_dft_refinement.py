import os, re, sys
import pickle
import numpy as np
from pyscf import gto
from pyscf import dft
from rdkit import Chem
from rdkit.Chem import rdMolTransforms
from pyscf.geomopt.geometric_solver import optimize

BASIS = 'def2-SVP'
XC = 'M06-2X'

def cal_mol_energy_via_dft_with_geomopt(mol, dihedrals: list, angle_base: float, mol_name: str, chk_path=None):
    """
        Mol: rdkit molecule object
        dihedrals: list of 4 atom indices
        angle_base: float, the value of dihedral angle
        mol_name: str, the name of the molecule
    """
    atom_coords = mol.GetConformer().GetPositions()
    atom_nums = [atom.GetAtomicNum() for atom in mol.GetAtoms()]

    # create pyscf molecule object
    pyscf_mol = gto.Mole()
    pyscf_mol.atom = [(atom_nums[i], atom_coords[i]) for i in range(mol.GetNumAtoms())]
    
    # set the basis set
    pyscf_mol.build(basis=BASIS)
    
    # geometry optimization
    mf_pre_eq = dft.RKS(pyscf_mol)
    mf_pre_eq.xc = XC
    
    # add constraint for geometry optimization, scan the dihedral angle with a fluctuation of 2 degrees
    constraint_path = f"constraint_{mol_name}"
    with open(constraint_path, "w") as f:
        f.write("$freeze\ndihedral {} {} {} {} {} {} 2".format(*dihedrals, angle_base-2, angle_base+2))
        
    # set max steps to 3, save time
    mf_eq = optimize(mf_pre_eq, maxsteps=3, constraints=constraint_path)
    
    os.system(f"rm {constraint_path}")
    
    # B3LYP functional
    mf_dft = gto.Mole(atom=mf_eq.tostring(), basis=BASIS)
    mf_dft.xc = XC 
    
    if chk_path:
        mf_dft.chkfile = chk_path
        
    e = mf_dft.kernel()
    
    return e, mf_dft

def cal_mol_energy_via_dft(mol, chk_path=None):
    """
        Mol: rdkit molecule object
        dihedrals: list of 4 atom indices
        angle_base: float, the value of dihedral angle
        mol_name: str, the name of the molecule
    """
    atom_coords = mol.GetConformer().GetPositions()
    atom_nums = [atom.GetAtomicNum() for atom in mol.GetAtoms()]

    # create pyscf molecule object
    pyscf_mol = gto.Mole()
    pyscf_mol.atom = [(atom_nums[i], atom_coords[i]) for i in range(mol.GetNumAtoms())]
    
    # set the basis set
    pyscf_mol.build(basis=BASIS)
    mf_dft = dft.RKS(pyscf_mol)
    mf_dft.xc = XC
    
    if chk_path:
        mf_dft.chkfile = chk_path
        
    e = mf_dft.kernel()
    
    return e, mf_dft

def extract_index_from_xtb_input(xtb_input):
    """
        Extract the dihedral indices from xtb input file
    """
    with open(xtb_input, 'r') as f:
        lines = f.readlines()
        
    match = re.search(r'dihedral: (\d+),(\d+),(\d+),(\d+)', lines[3])
    if match:
        numbers = [int(match.group(i)) for i in range(1, 5)]
    return numbers

def calculate_dihedral_angle(mol, atom_indices: list):
    """
        Calculate the dihedral angle of the molecule
    """
    conf = mol.GetConformer()
    angle = rdMolTransforms.GetDihedralDeg(conf, atom_indices[0], atom_indices[1], atom_indices[2], atom_indices[3])
    return angle

def extract_molecules_with_dihedral(dir_path):
    """
        Extract the molecules and dihedral angles from the xtb scan log file, find the max and min energy conformations
    """
    with open(os.path.join(dir_path, "xtbscan.log"), 'r') as file:
        lines = file.readlines()
    index = extract_index_from_xtb_input(os.path.join(dir_path, "xtb_input"))
    
    molecules = []
    molecules_block = []
    angle_list = []
    energy_list = []
    current_molecule = []
    for line in lines:
        if line.split()[0].isdigit():
            if current_molecule:
                molecule_str = ''.join(current_molecule)
                molecule = Chem.MolFromXYZBlock(molecule_str)
                molecules_block.append(molecule_str)
                molecules.append(molecule)
            current_molecule = [line]
        else:
            current_molecule.append(line)
        if line.strip().startswith("energy"):
            energy_list.append(float(line.split()[1]))
        
    for mol in molecules:
        angle_list.append(calculate_dihedral_angle(mol, index))

    return {"alist": np.array(angle_list), "elist": np.array(energy_list), "mols": molecules, "mols_block": molecules_block, "index": index}

def extract_min_max_record(record: dict):
    """
        Extract the max and min energy information from the record
    """
    max_energy_indices = np.argpartition(record['elist'], -2)[-1:]
    min_energy_indices = np.argpartition(record['elist'], 2)[:1]
    
    return {
        "max_e": record['elist'][max_energy_indices],
        "min_e": record['elist'][min_energy_indices],
        'max_angles': record['alist'][max_energy_indices],
        'min_angles': record['alist'][min_energy_indices],
        'max_mols_block': [record['mols_block'][i] for i in max_energy_indices],
        'min_mols_block': [record['mols_block'][i] for i in min_energy_indices],
        "max_mols": [record['mols'][i] for i in max_energy_indices],
        "min_mols": [record['mols'][i] for i in min_energy_indices]
    }
    
def generate_dft_energy(root_path, save_path, mol_dir):
    """
        parsing the results dir, re-calculating the energy of the max and min conformations using DFT
        
        parameters:
        `root_path`: str, the root path of the results
        `save_path`: str, the path to save the results
        `mol_dir`: str, the name of each molecule dir
    """
    file_path = os.path.join(root_path, mol_dir)
    assert os.path.exists(file_path), f"Mol {mol_dir} not exists"
    
    save_dir = os.path.join(save_path, mol_dir)
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
    
    e_list = []
    gto_mol_list = []
    angle_list = []
    all_rec = extract_molecules_with_dihedral(file_path)
    min_max_rec = extract_min_max_record(all_rec)
    
    for i in range(len(min_max_rec["max_mols_block"])):
        max_mol = min_max_rec["max_mols"][i]
        min_mol = min_max_rec["min_mols"][i]
        max_mol_block = min_max_rec["max_mols_block"][i]
        min_mol_block = min_max_rec["min_mols_block"][i]
        
        # cal the energy of the max conformations
        dihedral = min_max_rec["max_angles"][i]
        chk_path = os.path.join(save_dir, f"dft_opt_max_{i}_{dihedral}.chk")
        max_energy, dft_mol = cal_mol_energy_via_dft(max_mol, chk_path)
        e_list.append(max_energy)
        gto_mol_list.append(max_mol)
        angle_list.append(dihedral)
        
        with open(os.path.join(save_dir, f"dft_opt_max_{i}_{dihedral}.xyz"), "w+") as f:
            f.write(max_mol_block)
        
        # cal the energy of the min conformations
        dihedral = min_max_rec["min_angles"][i]
        chk_path = os.path.join(save_dir, f"dft_opt_min_{i}_{dihedral}.chk")
        min_energy, dft_mol = cal_mol_energy_via_dft(min_mol, chk_path)
        e_list.append(min_energy)
        gto_mol_list.append(min_mol)
        angle_list.append(dihedral)
        
        with open(os.path.join(save_dir, f"dft_opt_min_{i}_{dihedral}.xyz"), "w+") as f:
            f.write(min_mol_block)
        
    with open(os.path.join(save_dir, "dft_energy.pkl"), "wb+") as f:
        pickle.dump({"energy": e_list, "mols": gto_mol_list, "angles": angle_list}, f)
        
if __name__ == "__main__":

    generate_dft_energy(sys.argv[1], sys.argv[2], sys.argv[3])