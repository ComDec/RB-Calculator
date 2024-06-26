import pubchempy as pcp
from tqdm import tqdm
import itertools

qudaiji_list=['C', 'N','CCl',
 'CCC#N',
 'CC(=O)O',
 'NC',
 'C(CS)C',
 'CC#N',
 'C(=O)N(C)C',
 'C(=O)NCC',
 'Br',
 'I(C)C',
 'CC(N(=O)=O)N(=O)=O',
 'c1ccc(cc1)N',
 'CC(C)CC',
 'C(=O)C',
 'P',
 'CClC',
 'C[N+](=O)[O-]',
 'N(C)',
 'CSC',
 'CC(F)C',
 'C(=O)N(C)CC',
 'c1ccc(N(=O)=O)cc1',
 'CCP',
 'C(=O)N(C)CCCC',
 'C#N',
 'C(C)(C)C',
 'C(=O)NCCCC',
 'CCC',
 'CC(C)C#N',
 'c1ccc(O)cc1',
 'C(F)(F)C',
 'C(C=O)C',
 'C(=O)N(C)(C)C',
 'c1ccc(cc1)C',
 'N(=O)=O',
 'SC',
 'CCO',
 'CC(C)(C)O',
 'C(=O)NC',
 'OCC',
 'S',
 'F',
 'CC(C)O',
 'CCS',
 'C(Br)(Br)C',
 'C(=O)N(C)(C)CCCC',
 'C(=O)N',
 'C(C(=O)O)C(=O)O',
 'CS',
 'c1ccc(cc1)F',
 'C(=O)N(C)(C)CC',
 'c1ccc(cc1)O',
 'C(I)(I)C',
 'c1ccccc1',
 'c1ccc(C)cc1',
 'O',
 'C(=O)O',
 'CI',
 'PC',
 'C(=O)NCCC',
 'FCC',
 'I',
 'NCC',
 'Cl',
 'CCCC',
 'C(=O)',
 'FC',
 'OC',
 'C(=N)N',
 'CC(C)C',
 'CCSC',
 'C(=O)N(C)(C)CCC',
 'C(CP)C',
 'CBr',
 'NCCC',
 'S(=O)(=O)O',
 'SCC',
 'CO',
 'C(F)(F)F',
 'N(C)C',
 'C(Cl)(Cl)C',
 'CCBr',
 'CN',
 'PCC',
 'ClC',
 'C(CI)C',
 'CC',
 'BrC',
 'C(=O)N(C)CCC']
for i, j in tqdm(itertools.product(range(len(qudaiji_list)), repeat=2), total=len(qudaiji_list)**2):
# for i, j in itertools.product(range(len(qudaiji_list)), repeat=2):
    if i < j:
        smiles=f'c1({qudaiji_list[i]})cccc({qudaiji_list[j]})c1-c2c({qudaiji_list[j]})cccc2({qudaiji_list[i]})'
        try:
            mol=pcp.get_compounds(smiles, 'smiles')[0]
        except:
            try:
                mol=pcp.get_compounds(smiles, 'smiles')[0]
            except:
                continue
        if mol.cid is None:
            continue
        print(str(mol.cid)+"\n")
