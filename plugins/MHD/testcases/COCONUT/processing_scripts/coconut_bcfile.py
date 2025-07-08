import numpy as np
from astropy.io import fits
import scipy.special as scisp
import scipy.interpolate as interpolate
import matplotlib.pyplot as plt
import sunpy.coordinates
from datetime import datetime, timedelta
import sunpy.util.net
import requests
from bs4 import BeautifulSoup
from urllib.request import urlopen

#User Defined
date = '2020-12-07T15:00:00'
map_type ='HMI'  #'GONG' 'ADAPT', 'HMI
# NB: for the moment, only GONG, ADAPT and HMI are operational
output_dir = './'
lmax = 20 # numbers of modes
adapt_map = 6 # for ADAPT maps, between 0 and 11
amp = 1 # amplitude factor of the map
r_st = 1.0 # radius at which the magnetic field is computed
write_map = 'yes' # 'yes' or 'no'
show_map = 'yes' # 'yes' or 'no'
visu_type = 'sinlat' # 'lat' or 'sinlat'

## Name of the ouput file
year = int(date.split('-')[0])
month = int(date.split('-')[1])
day = int(date.split('-')[2].split('T')[0])
hour = int(date.split('T')[1].split(':')[0])
minute = int(date.split(':')[1])
second = int(date.split(':')[2])
date_datetime = datetime(year, month, day, hour, minute, second)
cr_number = int(sunpy.coordinates.sun.carrington_rotation_number(date_datetime))
print ("cr number ", cr_number)
# cr_number = cr_rotation
# WSO
if (map_type == 'WSO'):
  output_name = output_dir + 'map_wso_lmax' + str(lmax) + '_cr' + str(cr_number) + '.dat'
# GONG
# NB: only mrzqs for now
elif (map_type == 'GONG'):
  output_name = output_dir + 'map_gong_lmax' + str(lmax) + '_' + date[:10] + '.dat'
# ADAPT
elif (map_type == 'ADAPT'):
  output_name = output_dir + 'map_adapt_lmax' + str(lmax) + '_' + date + '.dat'
# HMI
# NB: only HMI carrington for now
elif (map_type == 'HMI'):
  output_name = output_dir + 'map_hmi_lmax' + str(lmax) + '_cr' + str(cr_number) + '.dat'
# KPVT
#output_name = '/Users/bperri/Documents/Data/Leuven/map_studies/cr1914/map_kpvt_cr1914_lmax30.dat'
# MDI
#output_name = '/Users/bperri/Documents/Data/Leuven/map_studies/cr2071/map_mdi_010808_lmax30.dat'
# SOLIS
#output_name = '/Users/bperri/Documents/Data/Leuven/map_studies/cr2193/map_solis_2193.dat'

## Name of the input magnetogram
# WSO
if (map_type == 'WSO'):
  map_name = 'wso_cr2194.txt'
#map_name = 'wso_cr2192.txt'
#map_name = 'wso_cr2219.txt'
#map_name = 'wso_cr1902.txt'
#map_name = 'wso_cr2072.txt'
# GONG
elif (map_type == 'GONG'):
  year_str = str(year)[2:]
  if (month < 10):
    month_str = '0' + str(month)
  else:
    month_str = str(month)
  if (day < 10):
    day_str = '0' + str(day)
  else:
    day_str = str(day)
  file_id = 'mrzqs'
  file_id_str = file_id[2:]
  remote_dir = 'https://gong.nso.edu/data/magmap/QR/' + file_id_str + '/' + str(year) + month_str + '/' + file_id + year_str + month_str + day_str + '/'
  # Find closest maps
  page_text = requests.get(remote_dir).text
  soup = BeautifulSoup(page_text, "html.parser")
  file_names = [node.get("href") for node in soup.find_all("a") if file_id in node.get("href")]
  time_deltas = list()
  for file_name in file_names:
    file_date = datetime.strptime(file_name.split("c")[0], file_id + "%y%m%dt%H%M")
    time_deltas.append((file_date - date_datetime).total_seconds())
  map_name = file_names[abs(np.array(time_deltas)).argmin()]
  remote_file = remote_dir + map_name
# ADAPT
elif (map_type == 'ADAPT'):
  remote_dir = 'https://gong.nso.edu/adapt/maps/gong/' + str(year) + '/'
  page_text = requests.get(remote_dir).text
  soup = BeautifulSoup(page_text, "html.parser")
  file_id = 'adapt40311'
  file_names = [node.get("href") for node in soup.find_all("a") if file_id in node.get("href")]
  time_deltas = list()
  for file_name in file_names:
    file_date = datetime.strptime(file_name.split("_")[2], "%Y%m%d%H%M")
    time_deltas.append((file_date - date_datetime).total_seconds())
  map_name = file_names[abs(np.array(time_deltas)).argmin()]
  remote_file = remote_dir + map_name
# HMI
elif (map_type == 'HMI'):
  # map_name = 'hmi.Synoptic_Mr_polfil.' + str(cr_number) + '.fits'
  map_name = 'hmi.Synoptic_Mr_small.' + str(cr_number) + '.fits'
  remote_file = 'http://jsoc.stanford.edu/data/hmi/synoptic/' + map_name
# Download file
local_file = sunpy.util.net.download_file(remote_file, directory = output_dir, overwrite = True)
# MDI
#map_name = 'synop_Mr_0.2047.fits'
#map_name = 'synop_Mr_0.2071.fits'
# SOLIS
#map_name = 'kbv7g170802t2030c2193_000_int-mas_dim-180.fits'
#map_name = 'kbv7g170802t2030c2193_000_int-mas_dim-900.fits'
# KPVT
#map_name = 'kbv7g060907t1457c2047_000_int-mas_dim-180.fits'
#map_name = 'm1914f.fits'

# Opening fits
print('Reading file')
input_file = output_dir + map_name
nb_modes_tot = int((lmax+1)*(lmax+2)/2 - 1)
# ADAPT
if (map_type == 'ADAPT'):
  input_data = fits.getdata(input_file, ext=0)
  shape = np.shape(input_data)
  nb_maps = shape[0]
  nb_th = shape[1]
  nb_phi = shape[2]
  Br_map = input_data[adapt_map,::-1,:]
  theta = np.linspace(0.,np.pi,nb_th)
  phi = np.linspace(0.,2.0*np.pi,nb_phi)
  Theta = np.tile(theta, (nb_phi,1)).T
  Phi = np.tile(phi, (nb_th,1))
  print('End of reading file')
# WSO
elif (map_type == 'wso'):
  fwso = open(input_file,'r')
  line = fwso.readline().split()
  if ('sine' in line):
    lat_type = 'sinlat'
  else:
    lat_type = 'lat'
  nb_th = int(line[1])
  nb_phi = int(360/5+1)
  nb_lines = 4*nb_phi
  nb_thplus = 4
  nb_th2 = nb_th + 2*nb_thplus
  Br_read = np.zeros((nb_th, nb_phi))
  fwso.readline()
  idx_th = nb_thplus
  idx_ph = nb_phi
  for ll in range(nb_lines):
    line = fwso.readline()
    if (line.split()[0][0] == 'C'):
      #idx_th = nb_thplus
      idx_th = 0
      idx_ph = idx_ph - 1
      for k in range(len(line.split())-1):
        Br_read[idx_th,idx_ph] = float(line.split()[k+1])
        idx_th = idx_th + 1
    else:
      for k in range(len(line.split())):
        Br_read[idx_th,idx_ph] = float(line.split()[k])
        idx_th = idx_th + 1
  fwso.close()
  if (lat_type == 'lat'):
    print('Extending Br')
    Br_ext = np.zeros((nb_th2, nb_phi))
    Br_ext[nb_thplus:nb_th2-nb_thplus,:] = Br_read
    idx_th = 0
    for k in range(nb_thplus):
      Br_ext[idx_th,:] = Br_read[0,:]
      idx_th = idx_th + 1
    idx_th = 0
    for k in range(nb_thplus):
      Br_ext[nb_th2-1 - idx_th,:] = Br_read[-1,:]
      idx_th = idx_th + 1
    Br_map = Br_ext*0.01 # from micro-tesla to gauss
    theta = (np.linspace(-90.,90.,nb_th2)+90.)*np.pi/180.
    phi = np.linspace(0.,360.,nb_phi)*np.pi/180.
    #theta = np.linspace(90.,-90.,nb_th2)
    #phi = np.linspace(0.,360.,nb_phi)
    Theta = np.tile(theta, (nb_phi,1)).T
    Phi = np.tile(phi, (nb_th2,1))
    nb_th = nb_th2
  else:
    Br_data = Br_read[::-1,:]*0.01 # from micro-tesla to gauss
    sinlat = np.linspace(-14.5/15.,14.5/15.,nb_th)
    theta_map = np.arcsin(sinlat) + np.pi/2.
    theta = np.linspace(0.,np.pi,nb_th)
    phi = np.linspace(0.,360.,nb_phi)*np.pi/180.
    Theta = np.tile(theta, (nb_phi,1)).T
    Theta_map = np.tile(theta_map, (nb_phi,1)).T
    Phi = np.tile(phi, (nb_th,1))
    #Br_data = Br_read[::-1,:]*0.01/np.cos(Theta_map) # from micro-tesla to gauss + from LOS to Br
    fbr = interpolate.RectBivariateSpline(theta_map,phi,Br_data)
    Br_map = fbr(theta,phi)
    Br_map = Br_map[::-1,:]
    #Br_map = Br_map/np.cos(Theta) # from LOS to Br
  print('End of reading file')
# GONG, HMI, MDI, KPVT, SOLIS, MWO
else:
  # input_data = fits.getdata(input_file, ext=0)
  input_data = fits.getdata(input_file)
  shape = np.shape(input_data)
  nb_th = shape[0]
  nb_phi = shape[1]
  Br_data = input_data[::-1,:]

  # print (nb_th, nb_phi, len(Br_data[0]))
  sinlat = np.linspace(-1.,1.,nb_th)
  theta = np.arcsin(sinlat) + np.pi/2.
  phi = np.linspace(0.,2.0*np.pi,nb_phi)
  Theta = np.tile(theta, (nb_phi,1)).T
  Phi = np.tile(phi, (nb_th,1))
  Br_data = np.nan_to_num(Br_data)
  Br_map = Br_data
  print('End of reading file')
Br = Br_map

print (Br)

# Decomposition of Br on spherical harmonics
print('Beginning of projection')
dtheta = np.tile(np.concatenate([np.diff(theta),[theta[-1]-theta[-2]]]),(nb_phi,1)).T
dphi = np.tile(np.concatenate([np.diff(phi),[phi[1]-phi[0]]]),(nb_th,1))
coefbr = np.zeros(nb_modes_tot, dtype=complex)
ylm = np.zeros((nb_modes_tot, nb_th, nb_phi), dtype=complex)
mod = 0
for l in range(1, lmax+1):
  print('{}'.format(l))
  for m in range(0, l+1):
    ylm[mod] = scisp.sph_harm(m, l, Phi, Theta)
    ylm_c = np.conj(ylm[mod])
    integrand_a = Br*ylm_c
    integrand_a = integrand_a*np.sin(Theta)*dphi*dtheta
    coefbr[mod] = np.sum(integrand_a)
    mod = mod+1
print('End of projection')

# Reconstruction of the field
print('Reconstructing Br')
Br_mode = np.zeros((nb_th,nb_phi))
mod = 0
# Reconstruct field up to lmax
for l in range(1,lmax+1):
  print('{}'.format(l))
  for m in range(0,l+1):
    ylm = scisp.sph_harm(m, l, Phi, Theta)
    Br_mode = Br_mode + np.real(coefbr[mod]*ylm)
    mod = mod+1
Br_mode = Br_mode/2.2 #not normalization of CF
Br_mode = Br_mode*amp # amplitude factor
print('End of reconstructing Br')

# Write boundary conditions file
if (write_map == 'yes'):
  print('Writing BC file')
  F = open(output_name,'w')
  F.write('1 \n')
  print ("nb_theta: ", nb_th, " nb_phi: ", nb_phi)
  F.write('!PHOTOSPHERE {} \n'.format((nb_th-2)*nb_phi+2))
  for j in range(nb_th):
    for k in range(nb_phi):
      xcoord = r_st*np.sin(theta[j])*np.cos(phi[k])
      ycoord = r_st*np.sin(theta[j])*np.sin(phi[k])
      zcoord = r_st*np.cos(theta[j])
      if ((j == 0) & (k != 0)):
        break
      if ((j == nb_th-1) & (k != 0)):
        break
      F.write('{:.16e} {:.16e} {:.16e} {:.16e} \n'.format(xcoord, ycoord, zcoord, Br_mode[j,k]))

  F.close()
  print('End of writing BC file')

# Show maps

if (show_map == 'yes'):
  fig, (ax1, ax2) = plt.subplots(2,1,figsize=(8,5),sharex=True)
  if (map_type == 'ADAPT'):
    visu_type = 'lat'
  # Latitudes and longitudes
  if (visu_type == 'sinlat'):
    if (map_type == 'wso'):
      if (lat_type == 'lat'):
        print('Careful, input file in latitudes, switching to lat plot!')
      else:
        longi = 180.*phi/np.pi
        Sinlat, Sinlong = np.meshgrid(sinlat, longi, indexing='ij')
    else:
      longi = 180.*phi/np.pi
      Sinlat, Sinlong = np.meshgrid(sinlat, longi, indexing='ij')
  lat = 90. - 180.*theta/np.pi
  longi = 180.*phi/np.pi
  Lat, Long = np.meshgrid(lat, longi, indexing='ij')
  # Plot original map
  if (visu_type == 'lat'):
    if (map_type == 'wso'):
      im1 = ax1.pcolormesh(Long,Lat,Br,cmap='seismic')
    else:
      im1 = ax1.pcolormesh(Long,Lat,Br,cmap='seismic')
    ax1.set_ylabel('Latitude')
  else:
    if (map_type == 'wso'):
      if (lat_type == 'lat'):
        im1 = ax1.pcolormesh(Long,Lat,Br,cmap='seismic')
      else:
        im1 = ax1.pcolormesh(Sinlong,Sinlat,Br_data,cmap='seismic')
    else:
      im1 = ax1.pcolormesh(Sinlong,Sinlat,Br_data[::-1],cmap='seismic')
      print (len(Sinlong), len(Sinlat), len(Br_data[::-1][0]))
    longi_pos = np.arange(0.,360.,60.)
    ax1.set_xticks(longi_pos)
    ax1.grid(visible=True, which='major', color='k', linestyle='-')
    ax1.set_ylabel('Sine Latitude')
  ax1.set_title('Original map')
  plt.colorbar(im1,ax=ax1)
  # Plot lmax reconstruction
  # im2 = ax2.pcolormesh(Long,Lat,Br_mode,cmap='seismic',vmin=-np.max(Br_mode),vmax=np.max(Br_mode))
  im2 = ax2.pcolormesh(Long,Lat,Br_mode,cmap='seismic')
  # print ("the resolution of the generated map: ", len(Long), len(Lat))
  ax2.set_title('BC file')
  ax2.set_ylabel('Latitude')
  ax2.set_xlabel('Longitude')
  plt.colorbar(im2,ax=ax2)
  plt.savefig(output_name+'.png')
  plt.show()
