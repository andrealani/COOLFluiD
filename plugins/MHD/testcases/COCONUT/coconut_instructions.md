09/05/25

Tinatin Baratashvili

## Install COCONUT step-by-step
-------------------------------
1. Follow the first 4 steps instructions on WIKIHOW of COOLFluiD repository:
https://github.com/andrealani/COOLFluiD/wiki/HOWTO-install-on-Genius-(for-KU-Leuven-users)
2. After step 4 (before starting the installation process) change the modules to the latest ones
- Genius

Replace the existing modules with these ones:
```
module load CMake/3.26.3-GCCcore-12.3.0
module load Boost/1.82.0-GCC-12.3.0
module load ParMETIS/4.0.3-gompi-2023a
module load PETSc/3.20.3-foss-2023a
```
Do not forget to change the PETSc also in each if statement after indicating CONF_FILE. The old PETSc module has to be replaced with
```
module load PETSc/3.20.3-foss-2023a
```

- Hortense

Replace the existing modules with these
```
module load CMake/3.26.3-GCCcore-12.3.0
module load PETSc/3.20.3-foss-2023a
```
Do not forget to change the PETSc also in each if statement after indicating CONF_FILE. The old PETSc module has to be replaced with
```
module load PETSc/3.20.3-foss-2023a
```
3. After step proceed with step 5 in the WIKIHOW. In case you get an error immediately regarding the first lines of install_COOLFluid.sh script you can locally run from the terminal
```
dos2unix install_COOLFluiD.sh
```

**Note:** The installation usually takes time and should not be interrupted. In case you can not wait for that long, or you don't trust your stable internet connection you can run it in the background with the following command (the installation options can be modified according to your preferences, the only difference is nohup and &)
```
nohup /install_COOLFluiD.sh DEBUG_NOCUDA --download=2 &
```

This will create a nohup.out file where you can monitor the installation process and you can also log out without interrupting the installation.

## Run COCONUT step-by-step

 After the installation is over, go to ~/COCONUT/Dipole folder and follow the instructions on how to run a testcase:
https://github.com/andrealani/COOLFluiD/wiki/HOWTO-run-a-testcase
If you successfully link the coolfluid solvers in your folder, you will have 3 files in your directory: _coolfluid-solver.xml, coolfluid-solver-wrapper, coolfluid-solver_.

After this the testcase can be run. The description of the necessary files are given in COCONUT_manual.pdf that can be found in this directory. You main file where you set up the boundary conditions, mesh, cfl and all the other settings is a .CFcase file. In the testcase all the paths are setup correctly. Here is an example job script to submit to the cluster to run the first testcase. These scripts are only for Genius and Hortense clusters. If you are using another cluster you can follow the example and change the paths and modules available to you.

- Genius


Change the account, jon-name, mail-user, nodes and ntasks-per-node according to the availability.  ntasks = nodes*ntasks-per-node, and time. Modify the path correctly according to your installation of COCONUT.

```
#!/bin/bash -l
#SBATCH --account="credit_name"
#SBATCH --job-name="WTD_coconut"
#SBATCH --mail-type="BEGIN,END,FAIL"
#SBATCH --mail-user="your_email"
#SBATCH --nodes="4"
#SBATCH --ntasks-per-node="36"
#SBATCH --ntasks="144"
#SBATCH --time="3:00:00"


cd /path-to-COCONUT-dir/COCONUT/Dipole

module load CMake/3.26.3-GCCcore-12.3.0
module load PETSc/3.20.3-foss-2023a
export OMP_NUM_THREADS=1

mpirun -np 144 ./path-to-COCONUT-dir/COCONUT/Dipole/map_TEST.CFcase

```

- Hortense
```
#!/bin/bash -l
#PBS -N your_credits
#PBS -l nodes=2:ppn=128
#PBS -l walltime=1:00:00
#PBS -A 2025_006
#PBS -m abe
#PBS -M your_email
cd /path-to-COCONUT-dir/COCONUT/Dipole
module load CMake/3.26.3-GCCcore-12.3.0
module load PETSc/3.20.3-foss-2023a
module load vsc-mympirun
mympirun --universe 256 ./path-to-COCONUT-dir/COCONUT/Dipole/map_TEST.CFcase
```

And then submit the job.

Good Luck!

## Pre-Processing

Before the start of the simulation, the working directory should contain the following files: CFmesh file (the computational mesh for the simulation), the input magnetic map file (.dat) format and an interactive file (.inter). The default mesh and inter files are in the COCONUT testcase folders. The default dipole.dat magnetic field configuration is also available. However, once the standard testcase is successful and you want to perform a simulation for the real solar mangetic field configuration, you should generate such an input .dat file from the mangetogram.

The magnetogram for the date of interest can be retrieved with the processing_scripts/coconut_bcfile.py. Here the user needs to define the date, the magnetogram source and the processing lmax number. This script will download the magnetogram process it and save the processed magnetic map file in the indicated directory with .dat extension. This file should be indicated in the wroking .CFcase file. For example:

```
Simulator.SubSystem.EM.Data.DistanceBased.FileNameTw = ./generated_magnetic_map_file.dat

```

## Post-processing

Once the simulation is finished you have the corona.plt, corona.CFmesh and *.vtu files. All these files contain the information about the simulations result.
- corona.plt
This file can be loaded and plotted in the Tecplot. The equations are given in ~/COCONUT/processing_scripts/equations_for_coconut.eqn, that can be loaded in Tecplot once the file is loaded to convert the solution to physical variables. The simple macros file is available in the ~/COCONUT/processing_scripts/plot_meridional_plane.mcr. If you have never used tecplot, you can load the equations, calculate them and then go to Scripting->Play Macro/Script and load this file. This will visualise the meridional plane and then you can change the variables to plot.
- corona.CFmesh
This file can be used to create the inner heliosphere file based on the output of the COCONUT model. For this you can use the script available in ~/COCONUT/processing_scripts/coupling_coconut_with_heliosphere.py. You need to indicate the paths to the files according to your setup. Modify the following lines:
```
MHDinfile = open("/path-to-your-COCONUT-output-file.CFmesh", "r")
  output_dat = '/path-to-your-temporary-boundary-file_temp.dat'
  output_dat_final = '/path-to-your-boundary-file_final.dat' # This will be the file to use later for EUHFORIA or Icarus
  mag_name = '/home/u0124639/Phd/coconut_files/test_icarus_coupling/hmi.Synoptic_Mr_small.2264.fits' # Give the path to the magnetogram file that the input file was created with
  time = '2022-11-22T01:30:00'
```
As a result, you can use the output_dat_final for starting the heliosphere model, either with Icarus or EUHFORIA.
- *.vtu
These files are partitioned. You can group them with the script ~/COCONUT/processing_scripts/group_vtu_files.py. Adjust nb_proc according to how many CPUs were used in the simulation, the directory where the *.vtu files are located (dirname) and the output file where the groupped .vtu file will be saved.
