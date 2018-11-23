#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 19 10:14:57 2018

Monitor the residuals of flow variables during a run of COOLFluiD.

Once the log of residual has fallen enough, the CFL should be modified
according to CFL*=2 in the .inter file.

Might-do list: 
    * Play an alarm sound when the CFL should be changed. For this the necessary residual
      fall-off (for instance one order of magnitude) has to be specified
    * Automatically change the CFL parameter according to CFL *= 2 by overwriting
      the appropriate line .inter on the server

@author: Peter Leitner
"""


#========================== SPECIFY RUNTIME PARAMETERS ==================================
# Specify time to wait (in minutes) until the next check for an update:
waiting_time_for_updates = 0.2   # in minutes
# Convergence file to be monitored on the server where the computation is running:
myPath = "/users/cpa/peterl/COOLFluiD/OPENMPI/optim/apps/Solver" \
         "/RESULTS_Corona_Restart/"
myFn = "convergence.plt-P0.FlowNamespace"
# SSH connection:
computation = "local" # Choose either "local" or "remote"
myHost,myUser,myPassword,myPort = "helium.esat.kuleuven.be","?","?",22
# Choose residual to monitor:
res = "T" # "rho"
# Save updated figure to eps and pdf format:
save_figs = False
#========================================================================================


import matplotlib.pylab as plt
# import StringManipulationTools as Tools # personal module of tools containing methods for string manipulation
# ... No longer needed: function find_numbers_in_string included in order to make
# the script portable

def find_numbers_in_string(a_string_containing_numbers):
    "Find either integers, floats or numbers in scientific notation."
    import re
    numeric_const_pattern = r"""
[-+]? # optional sign
(?:
(?: \d* \. \d+ ) # .1 .12 .123 etc 9.1 etc 98.1 etc
|
(?: \d+ \.? ) # 1. 12. 123. etc 1 12 123 etc
)
# followed by optional exponent part if desired
(?: [Ee] [+-]? \d+ ) ?
"""
    rx = re.compile(numeric_const_pattern, re.VERBOSE)
    listofnums = rx.findall(a_string_containing_numbers) # read the numbers out of the string
    return listofnums


# Check whether directory for figure output exists. If not, the directory is created
import os
if computation=="local":
    if not os.path.exists(myPath + "Figures-convergence/"):
        os.mkdir(myPath + "Figures-convergence")
if computation=="remote":
    import paramiko
    transport = paramiko.Transport((myHost, 22))
    transport.connect(username=myUser,password=myPassword)
    sftp = paramiko.SFTPClient.from_transport(transport)
    try:
        sftp.chdir(myPath + "Figures-convergence/")  # Test if remote_path exists
    except IOError:
        print("Path to save figures does not exist yet, I'm creating it for you.")
        sftp.mkdir(myPath + "Figures-convergence/")  # Create remote_path
        sftp.chdir(myPath + "Figures-convergence/")
    sftp.close()

    
myFile = myPath + myFn

# Initialize the convergence parameter lists
iterations = []
# Residuals
Bx_res = []; By_res = []; Bz_res = []
Ex_res = []; Ey_res = []; Ez_res = []
rho_res = []
vx_res = []; vy_res = []; vz_res = []
T_res = []
CFL = []
PhysTime = []
WallTime = []

if computation == "remote":
    client = paramiko.SSHClient()
    client.load_system_host_keys()
    client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    client.connect(myHost,myPort,myUser,myPassword)
    
    sftp = client.open_sftp()  
    fileObject = sftp.file(myFile,"rb").readlines()
    sftp.close()
elif computation == "local":
    fileObject = open(myFile,"rb").readlines()
        

lenfileObject = len(fileObject) # Number of iterations already computed

for i in range(2,len(fileObject)):   # 2 header lines
    nums = find_numbers_in_string(str(fileObject[i]))
    # There are 19 numbers in the string:
    # 0...number of iterations
    # 9...rho
    # 13...T
    # 14...CFL
    # 15...PhysTime
    # 17...WallTime
    iterations.append(int(nums[0]))
    Bx_res.append(float(nums[1])); By_res.append(float(nums[2])); Bz_res.append(float(nums[3]))
    Ex_res.append(float(nums[4])); Ey_res.append(float(nums[5])); Ez_res.append(float(nums[6]))
    rho_res.append(float(nums[9]))
    vx_res.append(float(nums[10])); vy_res.append(float(nums[11])); vz_res.append(float(nums[12]))
    T_res.append(float(nums[13]))
    CFL.append(int(nums[14]))
    PhysTime.append(float(nums[15]))
    WallTime.append(float(nums[17]))
    
WallTime = list(map(lambda x: x/60**2, WallTime)) # convert to hours


# Convert into list of floats - NO LONGER NEEDED
# iterations = [int(i) for i in iterations]
# rho_res = [float(i) for i in rho_res]
# T_res = [float(i) for i in T_res]
# CFL = [int(i) for i in CFL]
# PhysTime = [float(i) for i in PhysTime]
# WallTime = [float(i) for i in WallTime]


# Generate graph that is updated later periodically
plt.rc("text", usetex=True)
plt.rcParams["text.latex.preamble"]=[r"\usepackage{amsmath}"]
plt.close(0)
fig0, ax0 = plt.subplots(num=0)
plt.title("COOLFluiD Convergence Monitoring")
if res=="T":
    lh0, = ax0.plot(iterations,T_res,"+",color="#1f77b4")
    ax0.set_ylabel(r"$\log\, T$ Residual",color="#1f77b4")
elif res=="rho":
    lh0, = ax0.plot(iterations,rho_res,"+",color="#1f77b4")
    ax0.set_ylabel(r"$\log\, \varrho$ Residual",color="#1f77b4")    
ax0.set_xlabel(r"Number of iterations")
ax0.tick_params("y",colors="#1f77b4")
ax0.grid(True)
ax1 = ax0.twinx()
lh1, = ax1.plot(iterations,CFL,".",color="#d62728")
ax1.set_ylabel(r"CFL",color="#d62728")
ax1.tick_params("y",colors="#d62728")
plt.pause(1)   # flush the figure

plt.close(1)
fig1, ax2 = plt.subplots(num=1)
plt.title("COOLFluiD Convergence Monitoring")
lh2, = ax2.plot(iterations,Bx_res,"+",color="#ff7f0e",label=r"$\log\,B$ residual")
lh3, = ax2.plot(iterations,Ex_res,"+",color="#2ca02c",label=r"$\log\,E$ residual")
lh4, = ax2.plot(iterations,vx_res,"+",color="#7f7f7f",label=r"$\log\,U$ residual")
if res=="T":
    # Plot the residual that is not already monitored in Fig. 1
    lh5, = ax2.plot(iterations,rho_res,"+",color="#9467bd",label=r"$\log\,\varrho$ residual")
elif res=="rho":
    lh5, = ax2.plot(iterations,T_res,"+",color="#9467bd",label=r"$\log\,T$ residual")
plt.legend(loc=1)
ax2.set_xlabel(r"Number of iterations")
ax2.tick_params("y",colors="k")
ax2.grid(True)
ax3 = ax2.twinx()
lh6, = ax3.plot(iterations,CFL,".",color="#d62728")
ax3.set_ylabel(r"CFL",color="#d62728")
ax3.tick_params("y",colors="#d62728")
plt.pause(1)   # flush the figure

plt.close(2)
fig2, ax4 = plt.subplots(num=2)
plt.title("COOLFluiD Runtime")
lh7, = ax4.plot(iterations,WallTime,"+",color="#1f77b4")
ax4.set_xlabel(r"Number of iterations")
ax4.set_ylabel(r"Wall time (hrs)",color="#1f77b4")
ax4.tick_params("y",colors="#1f77b4")
ax4.grid(True)
ax5 = ax4.twinx()
lh8, = ax5.plot(iterations,CFL,".",color="#d62728")
ax5.set_ylabel(r"CFL",color="#d62728")
ax5.tick_params("y",colors="#d62728")
plt.pause(1)   # flush the figure


# Check every waiting_time_for_update minute whether the plot should be updated
import time
while True:
    time.sleep(waiting_time_for_updates*60) # check for update on the server
    print("Looking for an update...")
    
    if computation == "remote":
        sftp = client.open_sftp()
        fileObject = sftp.file(myFile).readlines()
        sftp.close()
    elif computation == "local":
        fileObject = open(myFile,"rb").readlines()
    
    if len(fileObject) != lenfileObject:
        print("Got an update, I replot.")
        nums = find_numbers_in_string(str(fileObject[len(fileObject)-1]))
        iterations.append(int(nums[0]))
        Bx_res.append(float(nums[1])); By_res.append(float(nums[2])); Bz_res.append(float(nums[3]))
        Ex_res.append(float(nums[4])); Ey_res.append(float(nums[5])); Ez_res.append(float(nums[6]))
        rho_res.append(float(nums[9]))
        vx_res.append(float(nums[10])); vy_res.append(float(nums[11])); vz_res.append(float(nums[12]))
        T_res.append(float(nums[13]))
        CFL.append(int(nums[14]))
        PhysTime.append(float(nums[15]))
        WallTime.append(float(nums[17])/60**2)

        # Replot
        lh0.set_xdata(iterations)
        lh2.set_xdata(iterations)
        if res=="T":
            lh0.set_ydata(T_res)
            lh5.set_xdata(iterations); lh5.set_ydata(rho_res)
        elif res=="rho":
            lh0.set_ydata(rho_res)
            lh5.set_xdata(iterations); lh5.set_ydata(T_res)
        lh1.set_xdata(iterations); lh1.set_ydata(CFL)
        lh2.set_xdata(iterations); lh2.set_ydata(Bx_res)
        lh3.set_xdata(iterations); lh3.set_ydata(Ex_res)
        lh4.set_xdata(iterations); lh4.set_ydata(vx_res)
        lh6.set_xdata(iterations); lh6.set_ydata(CFL)
        lh7.set_xdata(iterations); lh7.set_ydata(WallTime)
        lh8.set_xdata(iterations); lh8.set_ydata(CFL)
        ax0.relim(); ax0.autoscale_view(True,True,True)
        ax1.relim(); ax1.autoscale_view(True,True,True)
        ax2.relim(); ax2.autoscale_view(True,True,True)
        ax3.relim(); ax3.autoscale_view(True,True,True)
        ax4.relim(); ax4.autoscale_view(True,True,True)
        ax5.relim(); ax5.autoscale_view(True,True,True)
        plt.draw()
        plt.pause(1.)
        if save_figs:
            plt.savefig(myPath + "Figures-convergence/varCFL-" + res + "conv.eps",bbox_inches="tight")
            plt.savefig(myPath + "Figures-convergence/varCFL-" + res + "conv.pdf",bbox_inches="tight")

        lenfileObject = len(fileObject)
        
