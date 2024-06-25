#!/bin/bash

# 检查是否提供了足够的参数
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 filename number_of_parts"
    exit 1
fi

# 检查文件是否存在
if [ ! -f "$1" ]; then
    echo "File not found!"
    exit 1
fi

# 检查第二个参数是否为正整数
if ! [[ "$2" =~ ^[1-9][0-9]*$ ]]; then
    echo "The number of parts must be a positive integer."
    exit 1
fi

# 计算原始文件的行数
total_lines=$(wc -l < "$1")
# 计算每个分割文件应包含的行数
((lines_per_file = (total_lines + $2 - 1) / $2))  # 使用+\$2-1来确保向上取整

# 使用split命令分割文件
split -d -l "$lines_per_file" "$1" "${1}_part_"

echo "File $1 split into $2 parts."
