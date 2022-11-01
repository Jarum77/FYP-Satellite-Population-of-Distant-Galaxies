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

N0= -2.806814232365416       #Calculated from ndp.getnum for M0 = 10.5 and z = 0

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

    z += zstep

#Graph Plotting

fig, ax = plt.subplots()
ax.yaxis.set_ticks(np.arange(4, 11.5, 0.5))

plt.xlim([0,4])
plt.ylim([4,11])

plt.xlabel("Z-Index")
plt.ylabel("Stellar Mass, $10^{10}$ $M_{☉}$")
plt.title("Stellar Mass vs Redshift for Milky Way Progenitors", fontweight="bold")

plt.plot(zvals, Mzfourge_USD_vals, "#FF7969", label="ZFOURGE+σ")
plt.plot(zvals, Mzfourge_vals, "#E71E06", label="ZFOURGE")
plt.plot(zvals, Mzfourge_LSD_vals, "#FF7969", label="ZFOURGE-σ")

plt.legend()

ax.vlines(x=2.0,ymin=Mzfourge_LSD_vals[8],ymax=Mzfourge_USD_vals[8],color="#251615")
ax.vlines(x=2.0,ymin=4,ymax=Mzfourge_LSD_vals[8],color="#251615",linestyles='dashed')
ax.vlines(x=2.0,ymin=Mzfourge_USD_vals[8],ymax=11,color="#251615",linestyles='dashed')
ax.fill_between(x=zvals, y1=Mzfourge_LSD_vals, y2=Mzfourge_USD_vals, color="#FF948C")

plt.show()

