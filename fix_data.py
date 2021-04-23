import re
import os
from pathlib import Path

def fix_file(file_name):
    save_file = str(Path(file_name).with_suffix('')) + '_fixed.txt'
    with open(file_name, 'r') as read_file:
        os.chdir('fixed_originals')
        with open(save_file, 'w') as write_file:
            for line in read_file:
                line = re.sub(r"\s+", " ", line)
                write_file.write(line + "\n")
        os.chdir('..')

if __name__ == '__main__':
    os.chdir('veljko_data')
    contents = os.listdir()
    print(contents)
    for file in contents:
        if os.path.isfile(file) and (file != '.DS_Store'):
            fix_file(file)

