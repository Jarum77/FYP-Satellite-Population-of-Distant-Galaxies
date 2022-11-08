#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  1 22:34:28 2022

@author: george
"""
import numpy as np
from mpmath import mp
import NDpredict as ndp
import matplotlib.pyplot as plt

zi = 0.
zf = 4.
M0 = 10.5   # All stellar masses are assumed to be logarithmic

zvals = []
Mzfourge_vals = []
Mzfourge_USD_vals = []         #Upper SD value
Mzfourge_LSD_vals = []         #Lower SD value

Millus_vals = []
Millus_USD_vals =[]
Millus_LSD_vals = []

N0= -2.806814232365416 
N0illus = ndp.getnum(M0, 0, massfunc='illustris')      #Calculated from ndp.getnum for M0 = 10.5 and z = 0
print(N0illus)
zstep = 0.25
z = 0

while (z <= zf):

    Mzfourge = ndp.newmass(M0, zi, z, massfunc='zfourge') 
    
    NDmax = ndp.evolvingN(N0,zi,z) + ndp.sigmaN(N0,zi,z)
    NDmin = ndp.evolvingN(N0,zi,z) - ndp.sigmaN(N0,zi,z)
    
    

    Mzfourge_min = ndp.getmass(NDmax, z, massfunc="zfourge")
    Mzfourge_max = ndp.getmass(NDmin, z, massfunc="zfourge")
    #print(Mzfourge_max, ", ", Mzfourge, ", ", Mzfourge_min)
    zvals.append(z)
    Mzfourge_vals.append(Mzfourge)
    Mzfourge_USD_vals.append(Mzfourge_max)
    Mzfourge_LSD_vals.append(Mzfourge_min)
 
    
 
    NDmax_illus = ndp.evolvingN(N0illus,zi,z) + ndp.sigmaN(N0illus,zi,z)
    NDmin_illus = ndp.evolvingN(N0illus,zi,z) - ndp.sigmaN(N0illus,zi,z)
 
    Millus = ndp.newmass(M0, zi, z, massfunc='illustris') 
    
    Millus_min = ndp.getmass(NDmax_illus, z, massfunc="illustris")
    Millus_max = ndp.getmass(NDmin_illus, z, massfunc="illustris")
    
    Millus_vals.append(Millus)
    Millus_USD_vals.append(Millus_max)
    Millus_LSD_vals.append(Millus_min)

    z += zstep

#Graph Plotting
fig=plt.figure()
gs=fig.add_gridspec(2,hspace=0)
ax = gs.subplots(sharex=True)
fig.set_size_inches(8,8, forward=True)
ax[0].yaxis.set_ticks(np.arange(4, 11.5, 0.5))
ax[1].yaxis.set_ticks(np.arange(4, 11.5, 0.5))

plt.xlim([0,4])
plt.ylim([7,11])

plt.xlabel("Z-Index")
plt.ylabel("Stellar Mass, $10^{10}$ $M_{☉}$")

plt.suptitle("Stellar Mass vs Redshift for Milky Way Progenitors", fontweight="bold")

ax[0].plot(zvals, Mzfourge_USD_vals, "#E71E06", label="ZFOURGE+σ", alpha=0.75)
ax[0].plot(zvals, Mzfourge_vals, "#E71E06", label="ZFOURGE")
ax[0].plot(zvals, Mzfourge_LSD_vals, "#E71E06", label="ZFOURGE-σ", alpha=0.75)

ax[1].plot(zvals, Millus_USD_vals, "#1357a6", label="ILLUSTRIS+σ", alpha=0.75)
ax[1].plot(zvals, Millus_vals, "#1357a6", label="ILLUSTRIS")
ax[1].plot(zvals, Millus_LSD_vals, "#1357a6", label="ILLUSTRIS-σ", alpha=0.75)

ax[0].legend(loc='best')
ax[1].legend(loc='best')

ax[0].vlines(x=2.0,ymin=Mzfourge_LSD_vals[8],ymax=Mzfourge_USD_vals[8],color="#251615")
ax[0].vlines(x=2.0,ymin=7,ymax=Mzfourge_LSD_vals[8],color="#251615",linestyles='dashed')
ax[0].vlines(x=2.0,ymin=Mzfourge_USD_vals[8],ymax=11,color="#251615",linestyles='dashed')
ax[0].fill_between(x=zvals, y1=Mzfourge_LSD_vals, y2=Mzfourge_USD_vals, color="#E71E06",alpha=0.5)


ax[1].vlines(x=2.0,ymin=Millus_LSD_vals[8],ymax=Millus_USD_vals[8],color="#251615")
ax[1].vlines(x=2.0,ymin=4,ymax=Millus_LSD_vals[8],color="#251615",linestyles='dashed')
ax[1].vlines(x=2.0,ymin=Millus_USD_vals[8],ymax=11,color="#251615",linestyles='dashed')
ax[1].fill_between(x=zvals, y1=Millus_LSD_vals, y2=Millus_USD_vals, color="#1357a6",alpha=0.5)

plt.show()
fig.savefig("ZFOURGE_ILLUSTRIS_WITH_SD", dpi=fig.dpi)
