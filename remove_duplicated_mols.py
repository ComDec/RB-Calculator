import os
import sys
import shutil
from openbabel import openbabel as ob

def convert_xyz_to_smiles(file_path):
    """将XYZ文件转换为SMILES编码"""
    conv = ob.OBConversion()
    conv.SetInAndOutFormats("xyz", "smi")
    mol = ob.OBMol()
    conv.ReadFile(mol, file_path)
    return conv.WriteString(mol).split()[0]

def remove_duplicates(folder_path, save_path="unique_xyz_files"):
    """检测并去除重复的分子"""
    smiles_dict = {}
    for file in os.listdir(folder_path):
        if file.endswith('.xyz'):
            file_path = os.path.join(folder_path, file)
            smiles = convert_xyz_to_smiles(file_path)
            if smiles not in smiles_dict.values():
                smiles_dict[file] = smiles
    
    # 复制保留的文件到新文件夹
    new_folder_path = os.path.join(save_path)
    os.makedirs(new_folder_path, exist_ok=True)
    for file, smiles in smiles_dict.items():
        original_file_path = os.path.join(folder_path, file)
        new_file_path = os.path.join(new_folder_path, file)
        shutil.copyfile(original_file_path, new_file_path)

# 读取文件夹中的所有xyz文件并去除重复分子
folder_path = sys.argv[1]
remove_duplicates(folder_path, sys.argv[2])
