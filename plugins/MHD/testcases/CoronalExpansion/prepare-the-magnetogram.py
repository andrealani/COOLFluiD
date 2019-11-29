#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  2 19:01:01 2019

Python script for preprocessing a magnetogram intended to be read from
COOLFluiD.

Python libraries needed: python-3, numpy, scipy, astropy, os, sys, urllib.request

Usage:

> python3 prepare-the-magnetogram.py magnetogram-url smoothing-sigma scaling-factor Bref

Example:

> python3 prepare-the-magnetogram.py https://gong.nso.edu/data/magmap/QR/bqj/201207/mrbqj120712/mrbqj120712t1154c2125_024.fits.gz 8 0.26 2.2e-4

@author: Peter Leitner
"""


def generate_magnetogram_BC_inputfile_for_CF(theta_vec,phi_vec,magnetogram,fn,savedir=None):
    """
    Generate a COOLFluiD input file from a magnetogram as obtained by
    calling SolarPhysics.load_magnetogram(url,date). This raw
    magnetogram is a 2-D numpy-array with Br values in Gauss.
    The output for COOLFluiD is converted to Tesla.
    
    Written by Peter Leitner on Dec 13, 2018

    Parameters
    ----------
    theta_vec: Float array
        theta vector between 0 and pi not necessarily equidistant values in rad

    phi_vec: Float array
        phi vector between 0 and 2*pi not necessarily equidistant values in rad
    
    magnetogram : Float array (Ntheta,Nphi)
        Numpy array holding the magnetogram in units of G. Use function
        SolarPhysics.load_magnetogram(url,date) to obtain the magnetogram
        from the web archive.
    
    fn : String
        Filename of the magnetogram
    
    savedir : String
        Directory where the magnetogram should be copied to.
        Typically the search path of COOLFluiD.

    Returns
    -------
    x_vec,y_vec,z_vec,Br_vec : List
        Columns of coordinate vectors and Br (T) values

    Example
    -------
    >>> x,y,z,Br = generate_magnetogram_BC_inputfile_for_CF(theta,phi,magnetogram,"m2018-04-27.dat")
    """


    import numpy as np
    import os
    
    magnetogram /= 1e4 # Convert from G to T as needed by CF
    Ntheta = len(theta_vec); Nphi = len(phi_vec)
    R = 1.0 # for the adimensional MHD module of COOLFluiD
    
    num_points = (Ntheta-2)*Nphi + 2   # due to removal of north- and south pole

    if os.path.exists(fn):
        os.remove(fn)
    CF_input_file = open(fn,"a")
    CF_input_file.write("1\n")   # number of boundaries
    CF_input_file.write("!PHOTOSPHERE " + str(num_points) + "\n")

    
    x_vec = []; y_vec = []; z_vec= []; Br_vec = []   # initialize lists
        
    for i in range(len(theta_vec)):
        for j in range(len(phi_vec)):
            
            theta = theta_vec[i]
            phi = phi_vec[j]
            
            eps = 1.0 # generous epsilon value
            if theta_vec.max() > np.pi + eps and phi_vec.max() > 2*np.pi + eps:
                # I assume theta and phi are in degrees
                print("You probably provided angles in degrees!? I converted them to rad.")
                theta = np.radians(theta_vec[i])
                phi = np.radians(phi_vec[j])
            
            x = R*np.cos(phi)*np.sin(theta)
            y = R*np.sin(phi)*np.sin(theta)
            z = R*np.cos(theta)
            
            
            # Write Cartesian coordinates and corresponding Br values to file:
            # Remove singularities at north and soutpole,
            # i.e. theta = 0, pi, corresp. to i = 0, 179
            if i == 0 and j != 0:   # South pole singularity
                break               # Break out of j-loop and continue with next i-iteration
            if i == len(theta_vec)-1 and j != 0: # South pole singularity
                break
            line = "%.16e" % x + " " + "%.16e" % y + " " + "%.16e" % z \
                   + " " + "%.16e" % magnetogram[i,j] + "\n"
            CF_input_file.write(line)
            

            x_vec.append(x); y_vec.append(y); z_vec.append(z); 
            Br_vec.append(magnetogram[i,j])
    
    CF_input_file.close()




import sys
print()
magnetogram = sys.argv[1] #"https://gong.nso.edu/data/magmap/QR/bqj/201207/mrbqj120712/mrbqj120712t1154c2125_024.fits.gz"
stddev_smoothing = float(sys.argv[2]) # 8
print("Loading magnetogram from URL " + sys.argv[1])
print("Smoothing sigma: " + str(stddev_smoothing))
scaling_factor = float(sys.argv[3]) # 0.26
print("Scaling factor: " + str(scaling_factor))
Bref = float(sys.argv[4]) # 2.2e-4
print("Magnetic-field reference value: " + str(Bref))
print()
import numpy as np
import os
####import COOLFluiD_ToolKit as CFTK
from math import radians, degrees
#import matplotlib.pylab as plt
#plt.rc("text", usetex=True)
#plt.rcParams["text.latex.preamble"]=[r"\usepackage{amsmath}"]
ans = os.system("rm -rf magnetogram.fits.gz")
ans = os.system("rm -rf magnetogram.fits")

# AL: test
import ssl
ssl._create_default_https_context = ssl._create_unverified_context

import urllib.request
fn,header = urllib.request.urlretrieve(magnetogram,"magnetogram.fits.gz")
ans = os.system("gunzip magnetogram.fits.gz")
from astropy.io import fits
image_file = "magnetogram.fits"
magnetogram = fits.getdata(image_file)
magnetogram_orig = np.copy(magnetogram) # in Gauss



#%%=== P A R T  1 ========= Interpolate from sine-lat to equidistant 1°-grid ======

# Original magnetogram 1° resolution in phi, 180 values in theta but in sine-lat
R = 1.0
th_sinelat = np.arccos(np.linspace(-1,1,180))   # in rad
th_eqdist = np.linspace(np.pi,0,180)            # in rad
th_eqdist_hres = np.linspace(np.pi,0,360)       # in rad
ph = np.linspace(0,2*np.pi,360)                 # in rad



from scipy import interpolate

f = interpolate.interp1d(th_sinelat,magnetogram[:,0])
magip0 = f(th_eqdist)



f_hres = interpolate.interp1d(th_sinelat,magnetogram[:,0])
magip0_hres = f(th_eqdist_hres)


magnetogram_ip = np.zeros((360,360),dtype=float)


for i in range(360):
    f_hres = interpolate.interp1d(th_sinelat,magnetogram[:,i])
    mag_ip_i_hres = f_hres(th_eqdist_hres)
    magnetogram_ip[:,i] = mag_ip_i_hres


##%% === P A R T  2 ======== Do the smoothing and scaling =======================

import astropy.convolution as apc
g_kernel  = apc.Gaussian2DKernel(x_stddev=stddev_smoothing/2,y_stddev=stddev_smoothing)
# Factor 1/2 because the grid resolution in theta is twice as high!

magnetogram_ip_smoothed = apc.convolve(magnetogram_ip*scaling_factor,g_kernel)

##%% ==== P A R T  3 =========== Plot the magnetograms ===============================

##%% === P A R T  4 ========= Convert to csv for plotting with Paraview =================

##%% === P A R T  5 ============= Generate the magnetogram for COOLFluiD ======

# The magnetogram is divided by Bref so that the COOLFluiD input file is already saved in
# non-dimensional units:
generate_magnetogram_BC_inputfile_for_CF(th_eqdist_hres,ph,magnetogram_ip_smoothed/Bref,
   "magnetogram_GONG.dat")

