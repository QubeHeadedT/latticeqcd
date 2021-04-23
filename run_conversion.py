import os
import qcd_data_converter
from pathlib import Path 

if __name__ == '__main__':
    os.chdir('veljko_data/fixed_originals')
    contents = os.listdir()
    for file in contents:
        label = str(Path(file).with_suffix('')) + '_hyperspherical.txt'
        mendi = qcd_data_converter.MendcelliConverter(file, label=label)
        mendi.save_converted_data()
    
