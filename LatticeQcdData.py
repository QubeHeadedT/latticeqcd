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

    def _filter_hyperradius(self, min_hyperradius=0, max_hyperradius=None):
        """Returns the poitns in the data set filtered by hyperradius."""

        if max_hyperradius:
            filtered_points = [point for point in self._lattice_data if min_hyperradius <= math.sqrt(point.r2) <= max_hyperradius]
        else:
            filtered_points = [point for point in self._lattice_data if min_hyperradius <= math.sqrt(point.r2)]
        
        return filtered_points


    def plot_shape_space(self, min_hyperradius=0, max_hyperradius=None):
        """Plots the shape space of the data points in the x-y plane."""
        
        filtered_points = self._filter_hyperradius(min_hyperradius=min_hyperradius, max_hyperradius=max_hyperradius)

        x = [point.x for point in filtered_points]
        y = [point.y for point in filtered_points]
        hyperradius = [point.r for point in filtered_points]

        plt.scatter(x, y, c=hyperradius, cmap='jet', marker='x')
        plt.axis([-1.0,1.0,-1.0,1.0])
        plt.xlabel('x')
        plt.ylabel('y')
        plt.title(f'{self.name}')
        plt.colorbar(mappable=None, cax=None, ax=None, label = 'Hyperradius')
        plt.show()

    def plot_potential(self, min_hyperradius=0, max_hyperradius=None, mode='3d'):
        """Produces a 3D plot of the potential versus hyper-radial shape space."""

        allowed_modes = ['3d', 'x5', 'y0']
        assert mode in allowed_modes, f"'mode' for plot_potential must be one of {allowed_modes}."

        filtered_points = self._filter_hyperradius(min_hyperradius=min_hyperradius, max_hyperradius=max_hyperradius)

        if mode == '3d':
            self._plot_potential_3d(filtered_points)
        elif mode == 'x5' or mode == 'y0':
            self._plot_potential_line(filtered_points, line=mode)
        
    
    def _plot_potential_3d(self, filtered_points):
        """Plots the potential values over the full shape space."""

        x = [point.x for point in filtered_points]
        y = [point.y for point in filtered_points]
        v = [point.v / math.sqrt(point.r2) for point in filtered_points]
        err = [point.err for point in filtered_points]

        fig = plt.figure()

        ax = fig.add_subplot(111, projection = '3d')
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('V / Hyper-radius')
        ax.axis([-1.0,1.0,-1.0,1.0])

        ax.scatter(x, y, v, zdir = 'z')
        ax.errorbar(x, y, v, zerr=err, fmt="o")
        plt.title(f'{self.name}')
        plt.show()

    def _plot_potential_line(self, filtered_points, line='x5', line_error=1e-15):
        """Plots the potential across a line of the shape space."""
        allowed_lines = ['x5', 'y0']
        assert line in allowed_lines, f"'line' in _plot_potential_line must be in {allowed_lines}"

        if line == 'y0':
            line_points = [point for point in filtered_points if (point.y < line_error and point.x > - 0.5)]
        elif line == 'x5':
            line_points = [point for point in filtered_points if (abs(point.x + 0.5) < line_error)]

        x = [point.x for point in line_points]
        y = [point.y for point in line_points]
        v = [point.v / math.sqrt(point.r2) for point in line_points]
        err = [point.err for point in line_points]

        if line == 'y0':
            plt.errorbar(x, v, yerr=err, fmt='x')
            plt.show()
        elif line == 'x5':
            plt.errorbar(y, v, yerr=err, fmt='x')
            plt.show()

    def potential_versus_stirnglength(self, min_hyperradius=0, max_hyperradius=None, mode='delta'):
        """Plot the potential as a function of delta or gamma string length."""
        allowed_modes = ['delta', 'gamma']
        assert mode in allowed_modes, f"'mode' for potential_versus_stringlength must be in {allowed_modes}"

        filtered_points = self._filter_hyperradius(min_hyperradius=min_hyperradius, max_hyperradius=max_hyperradius)


        string_lengths = [self._get_string_length(point, mode=mode) for point in filtered_points]
        v = [point.v for point in filtered_points]
        err = [point.err for point in filtered_points]

        plt.errorbar(string_lengths, v, yerr=err, fmt='x')
        plt.show()
    
    def _get_string_length(self, point, mode):
        """Returns the gamma or delta string length of a point."""
        allowed_modes = ['delta', 'gamma']
        assert mode in allowed_modes, f"'mode' for _get_string_length must be in {allowed_modes}"

        if mode == 'delta':
            string_length = point.a + point.b + point.c
        elif mode == 'gamma':
            pass

        return string_length
        
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

        if len(data) == 9:
            err = float(data[8])
        else:
            err = 0.0

        return latticeQCDDataPoint(r2, a, b, c, x, y, z, v, err)


class latticeQCDDataPoint:
    """Represents a single data point in a lattice QCD experiment result."""

    def __init__(self, r2, a, b, c, x, y, z, v, err=0):

        # set input values:
        try:
            self._is_valid_data(r2, a, b, c, x, y, z, v, err)
            self.r2 = r2
            self.a = a
            self.b = b
            self.c = c
            self.x = x
            self.y = y
            self.z = z
            self.v = v
            self.err = err
        except AssertionError:
            print("Error in values input as data point.")
            sys.exit(1)

        # set derived values:
        self.r = math.sqrt(self.r2)
    
    def _is_valid_data(self, r2, a, b, c, x, y, z, v, err):
        """Checks whether the data inputs are valid."""
        assert isinstance(r2, float), "r2 must be a float"
        assert isinstance(a, float), "a must be a float"
        assert isinstance(b, float), "b must be a float"
        assert isinstance(c, float), "c must be a float"
        assert isinstance(x, float), "x must be a float"
        assert isinstance(y, float), "y must be a float"
        assert isinstance(z, float), "z must be a float"
        assert isinstance(err, float), "err must be a float"

    def __str__(self):
        return f"{self.r2} {self.a} {self.b} {self.c} {self.x} {self.y} {self.z} {self.v} {self.err}"


if __name__ == '__main__':
    """Test functions on a default example."""
    testData = LatticeQcdData(data_filepath='fixed_file', name="Mendi data", beta_value=5.8)
    #testData.append_data('taka_data')
    testData.plot_shape_space()
    print('test finished')

