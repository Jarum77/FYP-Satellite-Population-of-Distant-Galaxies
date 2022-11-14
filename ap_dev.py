#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  8 11:20:13 2022

@author: george
"""
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.io import readsav
import numpy as np
from mpmath import mp
import NDpredict as ndp
from aperture_func import Halo_Masses

mass_data = readsav('./data/ceers1/ceers1.me_Z02_calz.fastpp.delltaugt8.5.bc03.ch.save')
rshift_data = readsav('./data/ceers1/ceers1_v0.2_redshift.save')

mass_vals = mass_data['sed_lmass']
ids = mass_data['sed_id']
zvals = rshift_data['z_best']


#prog_vals = [[],[]]
#not_prog_vals = [[],[]]
#growth_curve = [[],[],[],[]] #z_val, median, upper, lower

zi = 0.
z = 0.
M0 = 10.5   # All stellar masses are assumed to be logarithmic
N0= -2.806814232365416  #Calculated from ndp.getnum for M0 = 10.5 and z = 0

#strip any -99 fluxes and same indexs for IDs and Masses
nintynines=np.full(len(zvals),-99)

ids=ids[~np.equal(nintynines,zvals)]
mass_vals=mass_vals[~np.equal(nintynines,zvals)]
zvals=zvals[~np.equal(nintynines,zvals)]

FOURGE_MINS=[]
FOURGE_MAXS=[]
mins=[]
maxes=[]

#create curve for interpolation
x=np.linspace(0, np.ndarray.max(zvals)+1, num=71, endpoint=False)
#x=np.linspace(0, , num=71, endpoint=False)
for i in x:
    ymax=ndp.evolvingN( N0 , 0 , np.round(i, 5 ) ) + ndp.sigmaN( N0 , 0, np.round(i , 5 ) )
    ymin=ndp.evolvingN( N0 , 0 , np.round(i, 5 ) ) - ndp.sigmaN( N0 , 0, np.round(i, 5 ) )
    fourgemax=ndp.getmass(ymax, round( i , 5 ) , massfunc="zfourge")
    fourgemin=ndp.getmass(ymin, round( i , 5 ) , massfunc="zfourge")
    maxes=np.append(maxes,fourgemax)
    mins=np.append(mins,fourgemin)
    
fMax=interp1d(x,mins)
fMin=interp1d(x,maxes)

#first min determination
lower_bounds=fMin(zvals)
upper_bounds=fMax(zvals)

zeros=np.zeros(len(lower_bounds))


zvals1=zvals[np.less(zeros, (np.subtract(mass_vals,lower_bounds)))]
mass_vals1=mass_vals[np.less(zeros, (np.subtract(mass_vals,lower_bounds)))]
zeros1=zeros[np.less(zeros, (np.subtract(mass_vals,lower_bounds)))]
upper_bounds1=upper_bounds[np.less(zeros, (np.subtract(mass_vals,lower_bounds)))]
lower_bounds1=lower_bounds[np.less(zeros, (np.subtract(mass_vals,lower_bounds)))]

zvals=zvals1[np.less(zeros1, (np.subtract(upper_bounds1,mass_vals1)))]
mass_vals=mass_vals1[np.less(zeros1, (np.subtract(upper_bounds1,mass_vals1)))]


halo_masses=Halo_Masses(zvals,mass_vals)
print(halo_masses)
#for i in 




# plt.plot(x,fMin(x))
# plt.plot(x,fMax(x), label='max')
# plt.plot(zvals, mass_vals, ls='None', marker='.')
# plt.xlim(0,3)
# plt.legend()
# plt.show()



