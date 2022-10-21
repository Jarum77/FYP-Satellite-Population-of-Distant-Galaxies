from astropy.io import ascii, fits
from os.path import dirname, join as pjoin
import scipy.io as sio
from scipy.io import readsav


sav_data=readsav('/home/AstroPhysics-Shared/DATA/JWST/CEERS/FAST/input/ceers1_v0.2.save')

print (sav_data.keys())

zvals=sav_data['z_best']

for i in zvals :
    if 1 < zvals[i] < 4 :
        print (zvals[i])