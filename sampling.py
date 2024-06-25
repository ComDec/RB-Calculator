import os
import shutil
from random import sample

def copy_random_files(source_folder, target_folder, num_files=200000):
    # 确保目标文件夹存在，如果不存在，则创建
    if not os.path.exists(target_folder):
        os.makedirs(target_folder)

    # 获取源文件夹中所有文件的列表
    all_files = [f for f in os.listdir(source_folder) if os.path.isfile(os.path.join(source_folder, f))]

    # 如果源文件夹中的文件数量少于要复制的文件数量，打印消息并退出
    if len(all_files) < num_files:
        print("Source folder does not contain enough files.")
        return

    # 随机选择文件
    selected_files = sample(all_files, num_files)

    # 复制文件到目标文件夹
    for file in selected_files:
        source_file = os.path.join(source_folder, file)
        target_file = os.path.join(target_folder, file)
        shutil.copy2(source_file, target_file)
    
    print(f"Successfully copied {num_files} files to {target_folder}")

import sys

source_folder_path = sys.argv[1]
target_folder_path = sys.argv[2]
copy_random_files(source_folder_path, target_folder_path)
