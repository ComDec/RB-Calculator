import os
import sys
import openbabel
from rdkit import Chem
import tqdm
from rdkit.Chem import AllChem

def xyz_to_sdf_to_mol(xyz_file, sdf_file="temp.sdf"):
    

    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("xyz", "sdf")
    mol = openbabel.OBMol()
    if not obConversion.ReadFile(mol, xyz_file):
        raise FileNotFoundError(f"Unable read xyz: {xyz_file}")

    obConversion.WriteFile(mol, sdf_file)
    rdkit_mol = Chem.SDMolSupplier(sdf_file)[0]

    return rdkit_mol

def get_adjacent_atoms(mol, atom_index):
    atom = mol.GetAtomWithIdx(atom_index)  # 获取指定索引的原子
    adjacent_atoms = [nbr.GetIdx() for nbr in atom.GetNeighbors()]  # 获取所有邻接原子的索引
    return adjacent_atoms

def xyz_file_sanity_check(xyz_input_dir, error_case, fail_case):
    f = open(error_case, "w+")
    f1 = open(fail_case, "w+")
    
    for file_name in tqdm.tqdm(os.listdir(xyz_input_dir)):
        xyz_input = os.path.join(xyz_input_dir, file_name)
        xyz_id = os.path.basename(xyz_input).split(".")[0]
        mol = xyz_to_sdf_to_mol(xyz_input)
        if mol is None:
            f1.write(os.path.basename(xyz_input) + "\n")
            continue
        index_a, index_b = int(xyz_id.split("_")[2]), int(xyz_id.split("_")[3])
        index_a_adjs = get_adjacent_atoms(mol, index_a)
        index_b_adjs = get_adjacent_atoms(mol, index_b)
        index_a_adj = index_a + 1
        index_b_adj = index_b - 1
        
        if index_a_adj in index_a_adjs and index_a_adj != index_b:
            if index_b_adj in index_b_adjs and index_b_adj != index_a:
                continue
        else:
            f.write(os.path.basename(xyz_input) + "\n")
    
    f.close()
    f1.close()

def cal_torsion_profile_via_xtb(xyz_input, save_folder):
    xyz_input_base = os.path.basename(xyz_input)
    xyz_id = os.path.basename(xyz_input).split(".")[0]
    
    index_a, index_b = int(xyz_id.split("_")[2]), int(xyz_id.split("_")[3])
    index_a_adj = index_a + 1
    index_b_adj = index_b - 1
    
    dihedral_indices = f"{index_a_adj + 1},{index_a + 1},{index_b + 1},{index_b_adj + 1}"
    
    workdir_dir = os.path.join(save_folder, xyz_id)
    workdir_dir = os.path.abspath(workdir_dir)
    os.makedirs(workdir_dir, exist_ok=True)
    
    with open(f"{workdir_dir}/xtb_input", "w") as f:
        f.write(f"$constrain\nforce constant=0.05\n$scan\ndihedral: {dihedral_indices},0.0; 0.0,360.0,144\n$end")
    
    os.system(f"cp {xyz_input} {workdir_dir}")
    os.system(f"cd {workdir_dir} && xtb {xyz_input_base} --input xtb_input --opt crude")
    
if __name__ == "__main__":
    # cal_torsion_profile_via_xtb(sys.argv[1], sys.argv[2])
    xyz_file_sanity_check(sys.argv[1], sys.argv[2], sys.argv[3])