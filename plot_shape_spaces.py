import os
import LatticeQcdData
from pathlib import Path

if __name__ == '__main__':
    os.chdir('veljko_data/converted_hyperspherical')
    contents = os.listdir()
    for file in contents: 
        if file != '.DS_Store':
            name = str(Path(file).with_suffix(''))
            print(name)
            testData = LatticeQcdData.LatticeQcdData(data_filepath=file, name=name)
            testData.plot_shape_space()