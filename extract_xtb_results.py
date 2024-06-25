import os
from tqdm import tqdm
import numpy as np
import pickle as pkl
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
from IPython.display import display, Image
    
def read_energy(filename):
    energies = []

    with open(filename, 'r') as file:
        lines = file.readlines()
        step = 0
        for line in lines:
            if line.strip().startswith("energy"):
                energy = float(line.split()[1])
                energies.append(energy)
                step += 1
                
    return np.array(energies), step

def parser_xtb_results(dir_path):
    dir_list = os.listdir(dir_path)
    tqdmer = tqdm(dir_list)
    save_list = []
    for dir_name in tqdmer:
        save_dir = os.path.join(dir_path, dir_name)
        if os.path.exists(os.path.join(save_dir, "xtbopt.xyz")) and os.path.exists(os.path.join(save_dir, "xtbscan.log")):
            mol = Chem.rdmolfiles.MolFromXYZFile(os.path.join(save_dir, "xtbopt.xyz"))
            elist, steps = read_energy(os.path.join(save_dir, "xtbscan.log"))
            ebar = elist.max() - elist.min()
            save_list.append({
                "id": dir_name,
                "mol": mol,
                "energy": ebar
            })
        else:
            continue
    with open("xtb_results.pkl", "wb+") as f:
        pkl.dump(save_list, f)
        
parser_xtb_results("./xtb_cal_mols")