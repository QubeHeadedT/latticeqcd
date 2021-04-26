import csv
import math
from LatticeQcdData import latticeQCDDataPoint

class HyperSphericalConverter:
    """Handles the conversion from cartesian latticeQCD data to hyper-spherical coordinates."""


    def __init__(self, file_name=None, delimiter=",", label=None):
        self.convert(file_name=file_name, delimiter=delimiter, label=label)


    def convert(self, file_name=None, delimiter=",", label=None):
        """Sets self.converted_data to data determined from CSV."""
        self.converted_data = self.convert_data(file_name=file_name)

        if label is None:
            self.label = 'Un-named data'
        else:
            self.label = label


    def csv_to_sidelengths(self, file_name=None):
        """Interface to convert from the data file format to triangle side-lengths a, b, c."""

        with open(file_name) as file:
            csv_data=csv.reader(file, delimiter="&")
        
            sidelengths = []
            
            for line in csv_data:
                ijk = line[0].replace("(", "").replace(")", "").replace(" ", "").split(",")
                i = float(ijk[0])
                j = float(ijk[1])
                k = float(ijk[2])

                potential = line[1].split("(")[0]

                sidelength_a = math.sqrt(i*i + j*j)
                sidelength_b = math.sqrt(j*j + k*k)
                sidelength_c = math.sqrt(i*i + k*k)
                
                sidelengths.append((sidelength_a, sidelength_b, sidelength_c, potential, 0))
        
        return sidelengths


    def convert_data(self, file_name):
        """Converts and returns lattice qcd data into hyper-spherical coordinates."""

        data_points = []

        for sidelength_a, sidelength_b, sidelength_c, v, err in self.csv_to_sidelengths(file_name):

            for i in range(6):
                a, b, c = self.side_length_permutations(i, sidelength_a, sidelength_b, sidelength_c)

                x=(c*c+a*a-b*b)/a/2
                y=math.sqrt(c*c-x*x)
            
                #rho and lambda vector components
                rx=a/math.sqrt(2)
                ry=0
                lx=(a-2*x)/math.sqrt(6)
                ly=(-2*y)/math.sqrt(6)
            
                #calculate x, y, z and hyper-radius^2 values: 
                r2=rx*rx+ry*ry+lx*lx+ly*ly
                xx=2*(rx*lx+ry*ly)/r2
                yy=(rx*rx+ry*ry-lx*lx-ly*ly)/r2
                zz=2*(rx*ly-ry*lx)/r2 

                data_points.append(str(latticeQCDDataPoint(r2, a, b, c, xx, yy, zz, v, err=err)))
        
        return data_points


    def side_length_permutations(self, i,a,b,c):
        """Returns a permutation of quark positions for a three-body system."""
        perm={0:(a,b,c),1:(a,c,b),2:(b,c,a),3:(b,a,c),4:(c,a,b),5:(c,b,a)}
        return perm[i]


    def save_converted_data(self, output_file_name=None):
        """Save the converted data to file."""
        if output_file_name is None:
            output_file_name = self.label
        
        with open(output_file_name, 'w') as f:
            for point in self.converted_data:
                f.write(point)
                f.write('\n')


class MendcelliConverter(HyperSphericalConverter):
    """For reading data files from Mendicelli study."""

    def __init__(self, potential='Y', file_name=None, delimiter=",", label=None):
        assert (potential=='Y' or potential=='Delta'), 'potential parameter must be either Y or Delta'
        self._potential = potential
        super().__init__(file_name=file_name, delimiter=delimiter, label=label)


    def csv_to_sidelengths(self, file_name):

        sidelengths = []
        mode = None

        if file_name.startswith("ISO"):
            mode = "ISO"
        elif file_name.startswith("Right"):
            mode = "Right"

        with open(file_name, "r") as file:

            csv_data = csv.reader(file, delimiter=" ")

            for line in csv_data:
                #base, height, delta, Y, gamma3, Error-Jacknife
                data = [float(line[0]), float(line[1]), float(line[2]), float(line[3]), float(line[4]), float(line[5])]
                
                sidelengths.append(self._calculate_sidelengths(data, mode))
        
        return sidelengths

    def _calculate_sidelengths(self, data, mode):
        """Calculates the sidelengths from the Mendecelli data."""
        base = data[0]
        height = data[1]
        if mode == "ISO":
            hypotenuse = math.sqrt(pow((base / 2), 2) + pow(height, 2))
            length_a = base
            length_b = hypotenuse
            length_c = length_b
        elif mode == "Right":
            hypotenuse = math.sqrt(pow(base, 2) + pow(height, 2))
            length_a = base
            length_b = height
            length_c = hypotenuse
        else:
            raise ValueError("Mode for MedicelliConverter._calculate_sidelengths must be 'Right' or 'ISO'.")

        if self._potential == 'Y':
            potential = data[3]
        elif self._potential == 'Delta':
            potential = data[2]
        error = data[5]
        
        return (length_a, length_b, length_c, potential, error)

            



    



if __name__ == '__main__':
    #data = HyperSphericalConverter(file_name='beta60.csv', delimiter='&', label="Takahashi")
    #data.save_converted_data()

    mendi = MendcelliConverter('ISO_fixed_file.txt')
    mendi.save_converted_data()


