import os
import sys
import shutil
import zipfile
from tqdm import tqdm

# 文件夹A和B的路径
folder_a = sys.argv[1]
folder_b = sys.argv[2]

# 新文件夹的路径
new_folder = sys.argv[3]

# 创建新文件夹
os.makedirs(new_folder, exist_ok=True)

# 获取文件夹A和B中的文件名
files_in_a = set(os.listdir(folder_a))
files_in_b = set(os.listdir(folder_b))

# 找出只在文件夹A中的文件
unique_files_in_a = files_in_a - files_in_b

# 将这些文件复制到新文件夹
for file in tqdm(unique_files_in_a, desc="Copying files"):
    shutil.copy(os.path.join(folder_a, file), new_folder)

# 创建一个压缩文件
os.system(f'tar -zcvf unique_files.tar.gz {new_folder}')

# 删除新文件夹
shutil.rmtree(new_folder)
