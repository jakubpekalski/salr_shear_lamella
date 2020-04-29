import ovito
from ovito.io import export_file, import_file
from ovito.modifiers import PythonScriptModifier, ComputePropertyModifier, SpatialBinningModifier, AffineTransformationModifier
import numpy as np
import os.path
from os import path

def range_step(start, end, step):
    while start <= end:
          yield start
          start += step

T = '042'

# R  
#R = ['001', '003', '005', '008', '01','02','03','05','08','1','2','4','6','8','10','12','14','16','18','20','22','24','25','26','27','28','29','30','31','32','33']
R = ['11','13','15','17','3','5','7','9','19','23']

# Setting the input and output files
input_file = "/shear_prod_npt_T"+str(T)+"_R.dump"
#input_file = "/shear_prod_npt_T045_R.dump"
input_file_dir = "/home/pekalski/Projects/SALR/lammps/LJ612/epsLJ10Y05/npt/single/T"+str(T)+"/R"
output_filepath = "/home/pekalski/Projects/OP/"

start_frame = 500
step_frame = 1
number_of_bins = 100

# Setting the modifiers
def set_the_modifiers(node):

    # Transforming the simulation box to ortogonal
    aft = AffineTransformationModifier()
    aft.relative_mode==False
    aft.target_cell = [[35.7188,0,0, -17.8594],[0,35.7188,0, -17.8594],[0,0,34.0792, -17.0396]]
    node.modifiers.append(aft)

    # Associating value 1 with all the particles
    cpm = ComputePropertyModifier()
    cpm.output_property = "indicator"
    cpm.expressions = ["1"]
    node.modifiers.append(cpm)

    # Binning particles 
    bar = SpatialBinningModifier()
    bar.property = "indicator"
    bar.bin_count_z = number_of_bins
    bar.bin_count_x = number_of_bins
    bar.bin_count_y = number_of_bins
    bar.direction = SpatialBinningModifier.Direction.YZ
    bar.Operation.SumVol
    node.modifiers.append(bar)


# Loop over the frames
#den = np.zeros(number_of_bins*number_of_bins);

for shear_rate in R:

    if path.isfile(input_file_dir+str(shear_rate)+input_file):
       node = import_file(input_file_dir+str(shear_rate)+input_file, multiple_frames = True)
       num_frames = node.source.num_frames
       set_the_modifiers(node)
       print(shear_rate)

       listofsnapshots = []
       for frame in range_step(start_frame, num_frames,step_frame):
          output = node.compute(frame)
          data = output.grids['binning[indicator]']
          data = data['indicator']
          data = np.array(data,dtype=np.bool)
          listofsnapshots = np.concatenate((listofsnapshots, data))

       np.save(output_filepath+"/snaps_T"+str(T)+"_R"+shear_rate+".npy",listofsnapshots)
    else:
       print("File: "+input_file_dir+str(shear_rate)+input_file + " not found")
