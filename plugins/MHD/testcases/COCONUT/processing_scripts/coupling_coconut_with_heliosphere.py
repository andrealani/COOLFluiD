from scipy.io import readsav
import numpy as np
import pyvista as pv
from vtk.util.numpy_support import vtk_to_numpy
from scipy.interpolate import Rbf, griddata, RBFInterpolator, RegularGridInterpolator, LinearNDInterpolator
import matplotlib.pyplot as plt
import math
import os
import sys
import glob
import scipy
import astropy.units as u
from astropy.io import fits

from astropy.coordinates import SkyCoord
import sunpy.coordinates
from datetime import datetime
from astropy.coordinates import Angle




'''

This function creates the boundary file. It takes as input the small epsilon, number of latitudes and longitudes,
radius at which boundary is created, the path to the file for the boundary and the arrays for physical quantities.

As a result, the boundary data are stored in the output_dat file.
'''
def create_boundary_file(eps, nb_th, nb_phi,radius, output_dat, x, y, z, rho, vx, vy, vz, bx, by, bz, p, temp):

    r = np.sqrt(x*x+y*y+z*z)
    br = x/r*bx+y/r*by+z/r*bz
    vr = x/r*vx+y/r*vy+z/r*vz
    r_def = radius
    LAT = np.linspace(eps, np.pi - eps, nb_th)
    LON = np.linspace(eps, 2.0 * np.pi - eps, nb_phi)
    lat, lon = np.meshgrid(LAT, LON, indexing='ij')
    x_1 = r_def * np.sin(lat.flatten()) * np.cos(lon.flatten())
    y_1 = r_def * np.sin(lat.flatten()) * np.sin(lon.flatten())
    z_1 = r_def * np.cos(lat.flatten())

    x_1 = x_1.flatten()
    y_1 = y_1.flatten()
    z_1 = z_1.flatten()


    grid_before = np.stack((x_1, y_1, z_1), axis=1)


    grid_after = np.stack((x, y, z), axis=1)


    uni_grid_after, index_after = np.unique(grid_after, axis=0, return_index=True)


    rho_i = RBFInterpolator(uni_grid_after, rho[index_after], kernel='linear', neighbors=50)

    dens_0 = rho_i(grid_before)
    dens_int = (dens_0).reshape(len(LAT), len(LON))

    bri = RBFInterpolator(uni_grid_after, br[index_after], kernel='linear', neighbors=50)

    br_0 = bri(grid_before)
    br_int = br_0.reshape(len(LAT), len(LON))

    vri = RBFInterpolator(uni_grid_after, vr[index_after], kernel='linear', neighbors=50)

    vr_0 = vri(grid_before)
    vr_int = vr_0.reshape(len(LAT), len(LON))


    temp_i = RBFInterpolator(uni_grid_after, temp[index_after], kernel='linear', neighbors=50)

    temp_0 = temp_i(grid_before)
    temp_int = temp_0.reshape(len(LAT), len(LON))

    p_i = RBFInterpolator(uni_grid_after, p[index_after], kernel='linear', neighbors=50)

    p_0 = temp_i(grid_before)
    p_int = p_0.reshape(len(LAT), len(LON))


    if (generate_boundary_file == 'yes'):
        print('Writing boundary file...')
        with open(output_dat, 'w') as f:
            f.write('Time:\n{}\n'.format(time))
            rad_m = radius * 6.96e8
            f.write('Radius of sphere:\n{}\n'.format(rad_m))
            f.write('Number of colatitude grid points:\n{}\n'.format(nb_th))

            f.write('Colatitude grid points:\n' + '\n'.join('{:.19e}'.format(LAT[j]) for j in range(nb_th)) + '\n')

            f.write('Number of longitude grid points:\n{}\n'.format(nb_phi))
            f.write('Longitude grid points:\n' + '\n'.join('{:.19e}'.format(LON[k]) for k in range(nb_phi)) + '\n')
            f.write('vr\n' + '\n'.join(
                ['\n'.join('{:.19e}'.format(vr_int[j, k]) for j in range(nb_th)) for k in range(nb_phi)]) + '\n')
            f.write('number_density\n' + '\n'.join(
                ['\n'.join('{:.19e}'.format(dens_int[j, k]) for j in range(nb_th)) for k in range(nb_phi)]) + '\n')
            f.write('temperature\n' + '\n'.join(
                ['\n'.join('{:.19e}'.format(temp_int[j, k]) for j in range(nb_th)) for k in range(nb_phi)]) + '\n')
            f.write('Br\n' + '\n'.join(
                ['\n'.join('{:.19e}'.format(br_int[j, k]) for j in range(nb_th)) for k in range(nb_phi)]) + '\n')
        print('Done!')


# This file takes as an argument the path to the .CFmesh file and returns the np.arrays of the Quantities
# Input argument -  type (string)       : MHDinfile
# Output         -  type (numpy arrays) : x, y, z, rho, vx, vy, vz, bx, by, bz, p, temperature
def load_cfmesh(MHDinfile, radius):
    r_def = radius #*6.9551e8 #4.18*6.9551e8
    l = 1

    #Do not change lref unless you want all of the dimensions in your simulation to change:
    lref = 1.0

    #Lines below just determine the structure of the CFmesh file
    getEND = False
    doNOT = False

    string_nodes = "!LIST_NODE"
    string_data = "!LIST_STATE 1"
    string_connect = "!LIST_ELEM"


    idx1 = -1
    idx2 = -1
    idx0 = -1
    i = 1

    for MHDline in MHDlines:
        if MHDline[0:len(string_nodes)] == string_nodes:
            idx1 = i
        if MHDline[0:len(string_data)] == string_data:
            idx2 = i
        if MHDline[0:len(string_connect)] == string_connect:
            idx0 = i
        i = i + 1

    cell_centers = []
    connectivity = []
    coordinates = []

    rho = []
    center_x = []
    center_y = []
    center_z = []
    bx = []
    by = []
    bz = []
    vx = []
    vy = []
    vz = []
    p = []
    phis = []
    thetas = []
    temp = []


    spherical_coordinates = []
    spherical_data = []

    nbelements = -1


    #Here, insert a procedure to import your flux rope B field data
    #This data could have a structure of [[x1, y1, z1, Bx1, By1, Bz1], [x2, y2, z2, Bx2, By2, Bz2], ...]
    LOUIS_DATA = []
    n_LOUIS_DATA = len(LOUIS_DATA)

    for MHDline in MHDlines:


        if l == 7:
            line = MHDline.split(" ")
            # print (line)

            nbelements = int(line[1])
            # print ("number of elements: ", nbelements)
            MHDoutline = MHDline

        if getEND and MHDline[0] == "!":
            doNOT = True
            getEND = False
        if l > idx2 and MHDline[0] != "!" and MHDline[0] != "N" and len(MHDline) > 1 and doNOT == False:
            getEND = True
            vals = MHDline.split(" ")



            #We read the current cell state values
            Bx = float(vals[4])*2.2e-4
            By = float(vals[5])*2.2e-4
            Bz = float(vals[6])*2.2e-4
            # rho0 = float(vals[0])*1.67e-13/(1.67e-27*1e6)
            rho0 = float(vals[0])*1.67e-13/(1.67e-27)
            Vx0 = float(vals[1])*480248.0
            Vy0 = float(vals[2])*480248.0
            Vz0 = float(vals[3])*480248.0
            Pressure = float(vals[7])*0.03851

            cell_center = cell_centers[j]
            # print (j, cell_centers[j])


            x = cell_center[0]
            y = cell_center[1]
            z = cell_center[2]


            r = (x**2 + y**2 + z**2)**0.5
            rrho = (x**2 + y**2)**0.5

            if (x == 0 and  y ==0 ):
                x = 1e-12
                y = 1e-12

            if r < r_def * 1.05 and  r > r_def * 0.95:
                rho.append(rho0)
                center_x.append(x)
                center_y.append(y)
                center_z.append(z)
                bx.append(Bx)
                by.append(By)
                bz.append(Bz)
                vx.append(Vx0)
                vy.append(Vy0)
                vz.append(Vz0)
                p.append(Pressure)
                temp.append(Pressure/rho0/2/(1.38e-23))
            j = j + 1

        elif l > idx1 and l < idx2:
            vals = MHDline.split(" ")

            #We read the node coordinaes
            x = float(vals[0])
            y = float(vals[1])
            z = float(vals[2])

            #We manipulate the node coordinates if need be
            xMF = x * lref
            yMF = y * lref
            zMF = z * lref
            # print ("coordinates: ", l, x, y, z)


            coordinates.append([x, y, z])


        elif l == idx2:
            for i in range(0, len(connectivity)):
                n1 = connectivity[i][0]
                n2 = connectivity[i][1]
                n3 = connectivity[i][2]
                n4 = connectivity[i][3]
                n5 = connectivity[i][4]
                n6 = connectivity[i][5]

                cell_center_x = 1./6. * (coordinates[n1][0] + coordinates[n2][0] + coordinates[n3][0] + coordinates[n4][0] + coordinates[n5][0] + coordinates[n6][0])
                cell_center_y = 1./6. * (coordinates[n1][1] + coordinates[n2][1] + coordinates[n3][1] + coordinates[n4][1] + coordinates[n5][1] + coordinates[n6][1])
                cell_center_z = 1./6. * (coordinates[n1][2] + coordinates[n2][2] + coordinates[n3][2] + coordinates[n4][2] + coordinates[n5][2] + coordinates[n6][2])

                cell_centers.append(np.array([cell_center_x, cell_center_y, cell_center_z]))

            j = 0

            MHDoutline = MHDline

        elif l > idx0 and l <= idx0 + nbelements:
            vals = MHDline.split(" ")

            n1 = int(vals[0])
            n2 = int(vals[1])
            n3 = int(vals[2])
            n4 = int(vals[3])
            n5 = int(vals[4])
            n6 = int(vals[5])
            no = int(vals[6])

            connectivity.append([n1, n2, n3, n4, n5, n6])
            MHDoutline = MHDline
        else:
            MHDoutline = MHDline


        l = l + 1

    MHDinfile.close()

    cell_centers = np.array(cell_centers)
    x = np.array(center_x)
    y = np.array(center_y)
    z = np.array(center_z)
    r = np.sqrt(x*x+y*y+z*z)
    bx = np.array(bx)
    by = np.array(by)
    bz = np.array(bz)
    vx = np.array(vx)
    vy = np.array(vy)
    vz = np.array(vz)
    br = x/r*bx+y/r*by+z/r*bz
    vr = x/r*vx+y/r*vy+z/r*vz
    rho = np.array(rho)
    p = np.array(p)
    temp = np.array(temp)


    return x, y, z, rho, vx, vy, vz, bx, by, bz, p, temp


def coupling(time: str, mag_name_path: str, output_dat: str, output_dat_final: str, generate_boundary_plot: str) -> None:
    """ Create a new solar_wind.dat with the good rotation

    Args:
       time (str)             : time of the magnetogram
       mag_name (str)         : where is the magnetogram
       output_dat (str)       : where is the data
       output_dat_final (str) : where we want to save the date

    Returns:
       None

    """
    date_hmi=time
    mag_name=os.path.basename(mag_name_path)
    # Read the name of the magnetogram and find rotation angle depending on their type (GONG/ADAPT)


    # 1. Provide path that the employed magnetogram is found
    prefix = mag_name[:5]

    # GONG


    # Extract the date of interest from the name of the magnetogram
    cr = int(mag_name.split('.')[-2])
    date = datetime.strptime(date_hmi, '%Y-%m-%dT%H:%M:%S')


    # Get position of the central meridian at the corresponding date in carrington
    CM_HEEQ = SkyCoord(0 * u.deg, 0 * u.deg, frame='heliographic_stonyhurst', obstime=date,
                       observer='earth')  # position of the central meridian in stonyhurst
    CM_CAR = CM_HEEQ.transform_to('heliographic_carrington')  # position of the central meridian in carrington
    CM_CAR_value = CM_CAR.lon.value
    if CM_CAR_value < 0:
        CM_CAR_value = 360 + CM_CAR_value

    RotationAngle = CM_CAR_value + 10
    # Création des instances SkyCoord pour le méridien de Carrington et celui de Stonyhurst
    carrington = SkyCoord(0 * u.deg, 0 * u.deg, frame='heliographic_carrington', obstime='2019-07-02T00:00:00',
                          observer='earth'
                          )
    stonyhurst = SkyCoord(0 * u.deg, 0 * u.deg, frame='heliographic_stonyhurst', obstime='2019-07-02T00:00:00',
                          observer='earth'
                          )

    # Calcul de l'angle entre les deux méridiens
    angle = carrington.separation(stonyhurst)

    angle_degrees = Angle(angle).to('degree').value

    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # Give Input: the boundary.dat file from COCONUT
    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    # 2. Change paths here, assuming that the dirs from which we take the boundary data
    # and to which we will save the new EUHFORIA-formated boundary data, are different.
    print ("rotation angle ", RotationAngle)
    with open(output_dat) as f:
        lines = f.read().splitlines()

    clt_start = lines.index('Colatitude grid points:')
    lon_start = lines.index('Longitude grid points:')
    vr_start = lines.index('vr')
    n_start = lines.index('number_density')
    t_start = lines.index('temperature')
    br_start = lines.index('Br')

    clt = [float(x) for x in lines[clt_start + 1:lon_start - 2]]
    #lon = [float(x) + np.radians(angle_degrees) for x in lines[lon_start + 1:vr_start]]
    lon = [float(x) - (RotationAngle*np.pi/180.0) for x in lines[lon_start + 1:vr_start]]
    vr = [float(x) for x in lines[vr_start + 1:n_start]]  # m/s
    n = [float(x) for x in lines[n_start + 1:t_start]]  # /m3
    t = [float(x) for x in lines[t_start + 1:br_start]]  # K
    br = [float(x) for x in lines[br_start + 1:]]  # T to nT

    fp = open(output_dat_final, 'w')

    # Write time
    fp.write("Time:\n")
    date_utc = date.strftime("%Y-%m-%dT%H:%M:%S")
    fp.write(date_utc)
    fp.write("\n")

    # Write coordinates
    fp.write("Radius of sphere:\n")
    fp.write("14959787070.0\n")

    fp.write("Number of colatitude grid points:\n")
    fp.write(str(len(clt)) + "\n")
    fp.write("Colatitude grid points:\n")
    np.savetxt(fp, clt)

    fp.write("Number of longitude grid points:\n")
    fp.write(str(len(lon)) + "\n")
    fp.write("Longitude grid points:\n")
    np.savetxt(fp, lon)

    # Write variables
    # Vr
    fp.write("vr\n")
    np.savetxt(fp, vr)
    # N
    fp.write("number_density\n")
    np.savetxt(fp, n)
    # T
    fp.write("temperature\n")
    np.savetxt(fp, t)
    # Br
    fp.write("Br\n")
    np.savetxt(fp, br)
    # Close
    fp.close()




    plot_boundary_file(lon, clt, vr, n, t, br, generate_boundary_plot, 'jet')
    print('The boundary file with CR rotation is available at: '+output_dat_final)





def plot_boundary_file(lon, lat, vr, num_density, temperature, b_field, plot, cmap):
    nr_lon = len(lon)
    nr_clt=len(lat)
    vr_plot = np.empty(shape=(nr_lon,nr_clt),dtype='float')
    n_plot = np.empty(shape=(nr_lon,nr_clt),dtype='float')
    temp_plot = np.empty(shape=(nr_lon,nr_clt),dtype='float')
    br_plot = np.empty(shape=(nr_lon,nr_clt),dtype='float')

    # Assign variable values in the meshgrids
    for i in range(nr_lon):
        for j in range(nr_clt):
            vr_plot[i, j] = vr[i*nr_clt + j]
            n_plot[i, j] = num_density[i*nr_clt + j]
            temp_plot[i, j] = temperature[i*nr_clt + j]
            br_plot[i, j] = b_field[i*nr_clt + j]

    vr_plot=np.flipud(vr_plot)
    n_plot=np.flipud(n_plot)
    temp_plot=np.flipud(temp_plot)
    br_plot=np.flipud(br_plot)

    fig1, (ax5, ax6, ax7, ax8) = plt.subplots(4, 1, figsize=(8,10),sharex=True)
    # Plot overview
    if (plot == 'yes'):

        for i in range(len(lat)):
            lat[i] = 90. - lat[i]*180.0/np.pi
        for j in range(len(lon)):
            lon[j] = lon[j]*180.0/np.pi +297.6281883521901
        Lat, Long = np.meshgrid(lat, lon)



        im5 = ax5.pcolormesh(Long, Lat, np.flipud(vr_plot), cmap=cmap, shading='auto')
        fig1.colorbar(im5, ax=ax5)
        ax5.set_ylabel(r'$V_r \ (m.s^{-1})$')
        ax5.yaxis.set_label_position("right")
        ax5.set_ylim(-60,60)



        im6 = ax6.pcolormesh(Long, Lat, np.flipud(n_plot), cmap=cmap, shading='auto')
        fig1.colorbar(im6, ax=ax6)
        ax6.set_ylabel(r'$N \ (m^{-3})$')
        ax6.yaxis.set_label_position("right")
        ax6.set_ylim(-60,60)

        im7 = ax7.pcolormesh(Long, Lat, np.flipud(temp_plot), cmap=cmap, shading='auto')
        fig1.colorbar(im7, ax=ax7)
        ax7.set_ylabel(r'$T \ (K)$')
        ax7.yaxis.set_label_position("right")
        ax7.set_ylim(-60,60)

        im8 = ax8.pcolormesh(Long, Lat, np.flipud(br_plot), cmap='bwr', shading='auto')
        fig1.colorbar(im8, ax=ax8)
        ax8.set_ylabel(r'$B_r \ (T)$')
        ax8.yaxis.set_label_position("right")
        ax8.set_ylim(-60,60)

        ax8.set_xlabel("Longitudes")
        ax5.set_title("Boundary file at 0.1 AU")

    plt.show()


if __name__ == "__main__":
    generate_boundary_file = 'yes'
    generate_boundary_plot = 'yes'

    MHDinfile = open("/path-to-your-COCONUT-output-file.CFmesh", "r")
    output_dat = '/path-to-your-temporary-boundary-file_temp.dat'
    output_dat_final = '/path-to-your-boundary-file_final.dat' # This will be the file to use later for EUHFORIA or Icarus
    mag_name = '/home/u0124639/Phd/coconut_files/test_icarus_coupling/hmi.Synoptic_Mr_small.2264.fits' # Give the path to the magnetogram file that the input file was created with


    MHDlines = MHDinfile.readlines()

    time = '2022-11-22T01:30:00' # Time of the magnetogram
    nb_th = 90
    nb_phi = 180
    eps = 0.01
    radius_boundary = 21.5
    radius_tomography = 5.0
    if (generate_boundary_file == 'yes'):
        x, y, z, rho, vx, vy, vz, bx, by, bz, p, temp = load_cfmesh(MHDinfile, radius_boundary)
        create_boundary_file(eps, nb_th, nb_phi,radius_boundary, output_dat, x, y, z, rho, vx, vy, vz, bx, by, bz, p, temp)
        coupling(time, mag_name, output_dat, output_dat_final, generate_boundary_plot)
