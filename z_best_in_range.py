from astropy.io import ascii, fits
from os.path import dirname, join as pjoin
import scipy.io as sio
from scipy.io import readsav

rshift_data = readsav('./data/ceers1/ceers1_v0.2_redshift.save')
#freq_data=fits.open('./data/ceers1/CEERS_NIRCam1_v0.07.4_photom.fits')

zvals = rshift_data['z_best']
ids = rshift_data['id']
#f356 = freq_data[1].data['F356']

f=open("z_best_in_range.txt","w")
f.write("ID, Z-Index\n")
for i in range(len(zvals)) :
    if 1 < zvals[i] < 4 :
        f.write("%i, %.3f\n" %(ids[i], zvals[i]))