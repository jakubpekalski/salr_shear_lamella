import ovito
from ovito.io import export_file, import_file
from ovito.modifiers import PythonScriptModifier, WrapPeriodicImagesModifier, AffineTransformationModifier
import numpy as np
import os.path
from os import path

def range_step(start, end, step):
    while start <= end:
          yield start
          start += step

#
# Uwaga, ten kod zaklada ze interesujacy nas peak structure factor wystepuje pomiedzy 60 a 100 indeksem k_list
#
#

T = '042'

R = ['001', '003', '005', '008', '01','02','03','05','08','1','2','4','6','8','10','12','14','16','18','20','22','24','25','26','27','28','29','30','31','32','33']
#R = ['11','13','15','17','3','5','7','9','19','23']
#R = ['2','6','15','17']

#
k_start = 1.0
k_stop  = 3.0
k_len   = 500
#
start_frame = 500
step_frame = 100
#

# Setting the input and output files
input_file = "/shear_prod_npt_T"+str(T)+"_R.dump"
input_file_dir = "/home/pekalski/Projects/SALR/lammps/LJ612/epsLJ10Y05/npt/single/T"+str(T)+"/R"
output_filepath = "/home/pekalski/Projects/OP/structure_factor/"

# Setting the modifiers
def set_the_modifiers(node):

    # Transforming the simulation box to ortogonal
    aft = AffineTransformationModifier()
    aft.relative_mode=False
    aft.operate_on={'cell'}
    aft.target_cell = [[35.7188,0,0, -17.8594],[0,35.7188,0, -17.8594],[0,0,34.0792, -17.0396]]
    node.modifiers.append(aft)

    # Particle rotation so the slabs are  parallel to (z,x) plane
    aft2 = AffineTransformationModifier()
    aft2.relative_mode=True
    aft2.operate_on={'particles'}
    aft2.transformation=[[1., 0., 0., 0.],[0., 1., 0., 0.6],[0., 0., 1., 0.]]
    node.modifiers.append(aft2)

    # Wraping the particle postions to the ortogonal simulation box by applying perodic boundary conditions to the box of new shape
    wrap = WrapPeriodicImagesModifier()
    node.modifiers.append(wrap)


# output structure factor peak heights
output_file_shear = open(output_filepath+"/sf_shear_rates_T"+str(T)+".txt",'a')
output_file_peak  = open(output_filepath+"/sf_peak_T"+str(T)+".txt",'a')

# Loop over the shear rates at a given temperature T
for shear_rate in R:

    print("Computing structure factor for shear rate = "+str(shear_rate))
    # the average structure factors as a function of k are appended in the consecutive rows of this file:
    output_file_av_sf = open(output_filepath+"/sf_T"+str(T)+"_R"+shear_rate+".txt",'a')

    # Check if file is there and load it
    if path.isfile(input_file_dir+str(shear_rate)+input_file):
       node = import_file(input_file_dir+str(shear_rate)+input_file, multiple_frames = True)
       num_frames = node.source.num_frames
       set_the_modifiers(node)
       print(shear_rate)

       # Loop over the frames at a given shear rate
       av_sf = np.zeros(k_len)
       for frame in range_step(start_frame, num_frames,step_frame):
          # Apply appended modifiers to rotate particles, ortogonize the simulation box and wrap back the particles into the box
          data = node.compute(frame)
          k_list = np.linspace(k_start,k_stop,k_len)
          N = data.particles.count
          x = data.particles['Position'][:,1]
          # Compute the structure factor and find the height of the first peak
          sf = np.zeros(k_len)
          i = 0
          for k in k_list:
              sin_sum = np.sum(np.sin(k*x))
              cos_sum = np.sum(np.cos(k*x))
              sf[i] = (np.power(sin_sum,2) + np.power(cos_sum,2))/N
              i = i+1

          av_sf = av_sf + sf/k_len

       peak_h = av_sf[np.argmax(av_sf[60:120])+60]

       np.savetxt(output_file_av_sf,[av_sf])
       output_file_shear.write(shear_rate+" ")
       np.savetxt(output_file_peak,[peak_h])
       export_file(node, output_filepath+'/config_check_'+str(shear_rate)+'.dump','lammps/dump', columns = ["Position.X", "Position.Y","Position.Z"], multiple_frames = False)

    else:
       print("File: "+input_file_dir+str(shear_rate)+input_file + " not found")

output_file_shear.close()
output_file_av_sf.close()
output_file_peak.close()
