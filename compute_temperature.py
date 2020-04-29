from ovito import *
from ovito.io import export_file, import_file
from ovito.modifiers import PythonScriptModifier, ComputePropertyModifier, SpatialBinningModifier, CreateIsosurfaceModifier
import numpy as np


def range_step(start, end, step):
    while start <= end:
          yield start
          start += step

# Setting the input and output files
input_file = "shear_prod_0.16_w0.50_tg0.09_Ns2_p2.gsd"
input_file_dir = "/home/pekalski/Projects/SALR1836/MP_shear/N20300/tg009/temp/"
output_filepath = "/home/pekalski/Projects/SALR1836/analysis/temperature_YZ_50.txt"
velFile = "velzx_0.16_w0.50_tg0.09_Ns2_p2__final.npy"
vel_prof_path = "/home/pekalski/Projects/SALR1836/MP_shear/N20300/tg009/"+velFile

# import gsd file
node = import_file(input_file_dir+input_file, multiple_frames = True)
data = node.compute()

num_frames = node.source.num_frames

num_frames = 5500
start_frame = 100
step_frame = 100

# read in the velocity profile
vzx = np.load(vel_prof_path)
num_bins = np.shape(vzx)[0]
num_part = data.particles.count

#def sub_velp(frame,data):
#    data.particles_.create_property('ux', data = 

## Setting the modifiers
#
cpm = ComputePropertyModifier()
cpm.output_property = "loc_temperature"
cpm.expressions = ["(Velocity.X^2+Velocity.Y^2+Velocity.Z^2)/3"]
node.modifiers.append(cpm)
#
bar = SpatialBinningModifier()
bar.property = "loc_temperature"
bar.bin_count_x = num_bins
bar.bin_count_y = num_bins
bar.direction = SpatialBinningModifier.Direction.YZ
bar.Operation.Mean
node.modifiers.append(bar)

# Loop over the frames
temp = np.zeros(num_bins*num_bins);
for frame in range_step(start_frame, num_frames,step_frame):

    output = node.compute(frame)
    data = output.grids['binning[loc_temperature]']
    data = data['loc_temperature']
    temp = temp + data/(num_frames/step_frame)

# Save the output    
np.savetxt(output_filepath, temp)



