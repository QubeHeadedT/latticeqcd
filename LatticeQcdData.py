"""Class for loading, storing, plotting and manipulating latticeQCD data sets."""

import matplotlib.pyplot as plt
import math
import sys

class LatticeQcdData:
    
    def __init__(self, name=None, beta_value=None, data_filepath=None, start_line=0, end_line=-1):
        """Initialises an object from a lattice QCD datafile."""
        lattice_data = self._load_data_file(data_filepath, start_line, end_line)
        if lattice_data:
            self._lattice_data = lattice_data
        else:
            self._lattice_data = []

        if beta_value:
            self.beta_value = beta_value
        else:
            self.beta_value = "Unknown"
        
        if name:
            self.name = f"{name} with beta value = {self.beta_value}"
        else:
            self.name = f"Untitled data with beta value = {self.beta_value}"


    def plot_shape_space(self, min_hyperradius=0, max_hyperradius=None):
        """Plots the shape space of the data points in the x-y plane."""
        if max_hyperradius:
            x = [point.x for point in self._lattice_data if min_hyperradius <= point.r <= max_hyperradius]
            y = [point.y for point in self._lattice_data if min_hyperradius <= point.r <= max_hyperradius]
            hyperradius = [point.r for point in self._lattice_data if min_hyperradius <= point.r <= max_hyperradius]
        else:
            x = [point.x for point in self._lattice_data if min_hyperradius <= point.r]
            y = [point.y for point in self._lattice_data if min_hyperradius <= point.r]
            hyperradius = [point.r for point in self._lattice_data if min_hyperradius <= point.r]

        plt.scatter(x, y, c=hyperradius, cmap='jet', marker='x')
        plt.axis([-1.0,1.0,-1.0,1.0])
        plt.xlabel('x')
        plt.ylabel('y')
        plt.title(f'{self.name}')
        plt.colorbar(mappable=None, cax=None, ax=None, label = 'Hyperradius')
        plt.show()


    def _load_data_file(self, data_filepath, start_line, end_line):
        """Loads a data file and returns data as a dictionary."""

        if not data_filepath:
            return None

        self._valid_lines(start_line, end_line)

        try:
            with open(data_filepath, 'r') as f:
                if end_line == -1:
                    lines = f.readlines()[start_line:]
                else:
                    lines = f.readlines()[start_line:end_line]
        except FileNotFoundError:
            print(f'The file {data_filepath} was not found.\n')
            return None

        lattice_data = []

        try:
            for line in lines:
                lattice_data.append(self._create_datapoint(line))
        except:
            print(f'Error reading data from {data_filepath} to python lists.')
            return None

        return lattice_data


    def _valid_lines(self, start_line, end_line):
        """Checks whether the choice of start and end read lines is valid.
        
        Params:
            start_line (int): the line to start reading data file.
            end_line (int): the line to stop reading data file.

        Returns: 
            None

        Throws: 
            Assertion error with explanation if choices are invalid.
        """
        assert isinstance(start_line, int), "start_line must be an integer."
        assert isinstance(end_line, int), "end_line must be an integer."

        if end_line == -1:
            assert start_line >= 0, "start_line must be >= 0."
        else:
            assert start_line >= 0, "start_line must be >= 0."
            assert end_line >= 0, "Chosen end_line must be >= 0."
            assert start_line < end_line, "start_line must be less than end_line."

    def _create_datapoint(self, line):
        """Constructs a latticeQCDDatapoint object from a string line of data.
        
        Data line is in the format 'r2 a b c x y z v' """
        data = line.split()
        r2 = float(data[0])
        a = float(data[1])
        b = float(data[2])
        c = float(data[3])
        x = float(data[4])
        y = float(data[5])
        z = float(data[6])
        v = float(data[7])
        return latticeQCDDataPoint(r2, a, b, c, x, y, z, v)


class latticeQCDDataPoint:
    """Represents a single data point in a lattice QCD experiment result."""

    def __init__(self, r2, a, b, c, x, y, z, v):
        
        # set input values:
        try:
            self._is_valid_data(r2, a, b, c, x, y, z, v)
            self.r2 = r2
            self.a = a
            self.b = b
            self.c = c
            self.x = x
            self.y = y
            self.z = z
            self.v = v
        except AssertionError:
            print("Error in values input as data point.")
            sys.exit(1)

        # set derived values:
        self.r = math.sqrt(self.r2)
    
    def _is_valid_data(self, r2, a, b, c, x, y, z, v):
        """Checks whether the data inputs are valid."""
        assert isinstance(r2, float), "r2 must be a float"
        assert isinstance(a, float), "a must be a float"
        assert isinstance(b, float), "b must be a float"
        assert isinstance(c, float), "c must be a float"
        assert isinstance(x, float), "x must be a float"
        assert isinstance(y, float), "y must be a float"
        assert isinstance(z, float), "z must be a float"

    def __str__(self):
        return f"{self.r2} {self.a} {self.b} {self.c} {self.x} {self.y} {self.z} {self.v}"


if __name__ == '__main__':
    """Test functions on a default example."""
    testData = LatticeQcdData(data_filepath='test_mendi2', name="Mendi data", beta_value=5.8)
    #testData.append_data('taka_data')
    testData.plot_shape_space()
    print('test finished')

