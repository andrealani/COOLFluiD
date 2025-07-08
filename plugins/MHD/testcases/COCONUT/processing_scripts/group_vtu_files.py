import pyvista as pv


# User defined
nb_proc = 180#108 # number of processors used for the simulation


# directory where the files are located
dirname = '/home/u0124639/Phd/coconut_files/gong2008_files/'


# Create filenames
filelist = []
for i in range(nb_proc):
  filelist.append('corona-flow0-P' + str(i) + '.vtu')
filenames = list(map(dirname.__add__,filelist))
# Load files
print('Reading file')
mesh = pv.read(filenames)
print('Merging file')
merged = mesh.combine()
# Create one file with all datasets
print('Saving file')
# Indicate the name of the grouped final .vtu file
merged.save(dirname+'corona-cr2207_oldbcs.vtu')
