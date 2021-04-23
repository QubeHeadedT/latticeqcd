import re

if __name__ == '__main__':
    file_name = 'ISO_128_beta0,41_Gamma3_Y_delta.txt'
    save_file = 'fixed_file.txt'

    with open(file_name, 'r') as read_file:
        with open(save_file, 'w') as write_file:
            for line in read_file:
                line = re.sub(r"\s+", " ", line)
                write_file.write(line + "\n")