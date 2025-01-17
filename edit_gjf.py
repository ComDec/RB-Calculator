import os
import sys

# Function to edit the first 6 lines of a gjf file and save to a new directory
def edit_gjf_file(file_path, output_directory):
    file_name = os.path.basename(file_path)
    id = file_name.split('.')[0]
    new_lines = [
        '%nproc=64\n',
        '%mem=128GB\n',
        f'%chk={id}.chk\n',
        '#p m062x/def2svp opt freq\n',
        '\n',
        'Title Card Required\n',
        '\n',
        '0 1\n'
    ]
    with open(file_path, 'r') as file:
        lines = file.readlines()
    output_path = os.path.join(output_directory, file_name)
    with open(output_path, 'w') as file:
        file.writelines(new_lines + lines[6:])

# Function to update nproc and mem values
def update_nproc_mem(directory, nproc, mem):
    for filename in os.listdir(directory):
        if filename.endswith('.gjf'):
            file_path = os.path.join(directory, filename)
            with open(file_path, 'r') as file:
                lines = file.readlines()
            lines[0] = f'%nproc={nproc}\n'
            lines[1] = f'%mem={mem}GB\n'
            with open(file_path, 'w') as file:
                file.writelines(lines)

# Update main function to use the new output directory
if __name__ == "__main__":
    directory = sys.argv[1]
    nproc = sys.argv[2]
    mem = sys.argv[3]
    output_directory = os.path.join(directory, 'edited_gjf_files')
    os.makedirs(output_directory, exist_ok=True)
    for filename in os.listdir(directory):
        if filename.endswith('.gjf'):
            edit_gjf_file(os.path.join(directory, filename), output_directory)
    
    # update_nproc_mem(output_directory, nproc, mem) 