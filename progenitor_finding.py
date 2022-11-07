import matplotlib.pyplot as plt
from scipy.io import readsav
import numpy as np
from mpmath import mp
import NDpredict as ndp

mass_data = readsav('./data/ceers1/ceers1.me_Z02_calz.fastpp.delltaugt8.5.bc03.ch.save')
rshift_data = readsav('./data/ceers1/ceers1_v0.2_redshift.save')

mass_vals = mass_data['sed_lmass']
ids = mass_data['sed_id']
zvals = rshift_data['z_best']

prog_vals = [[],[]]
not_prog_vals = [[],[]]
growth_curve = [[],[],[],[]] #z_val, median, upper, lower

zi = 0.
z = 0.
M0 = 10.5   # All stellar masses are assumed to be logarithmic
N0= -2.806814232365416  #Calculated from ndp.getnum for M0 = 10.5 and z = 0

for i in range(len(ids)):
    if ( zvals[i] != -99 ):
        NDmax = ndp.evolvingN( N0 , 0 , round( zvals[i], 5 ) ) + ndp.sigmaN( N0 , 0, round( zvals[i] , 5 ) )
        NDmin = ndp.evolvingN( N0 , 0 , round( zvals[i], 5 ) ) - ndp.sigmaN( N0 , 0, round( zvals[i] , 5 ) )

        Mzfourge_min = ndp.getmass(NDmax, round( zvals[i] , 5 ) , massfunc="zfourge")
        Mzfourge_max = ndp.getmass(NDmin, round( zvals[i] , 5 ) , massfunc="zfourge")

        if (Mzfourge_min <= mass_vals[i] <= Mzfourge_max):
            prog_vals[0].append(zvals[i])
            prog_vals[1].append(mass_vals[i])
        else:
            not_prog_vals[0].append(zvals[i])
            not_prog_vals[1].append(mass_vals[i])

while (z <= 3):

    Mzfourge = ndp.newmass(M0, zi, z, massfunc='zfourge')    
    
    NDmax = ndp.evolvingN(N0,zi,z) + ndp.sigmaN(N0,zi,z)
    NDmin = ndp.evolvingN(N0,zi,z) - ndp.sigmaN(N0,zi,z)

    Mzfourge_min = ndp.getmass(NDmax, z, massfunc="zfourge")
    Mzfourge_max = ndp.getmass(NDmin, z, massfunc="zfourge")

    growth_curve[0].append(z)
    growth_curve[1].append(Mzfourge)
    growth_curve[2].append(Mzfourge_max)
    growth_curve[3].append(Mzfourge_min)

    z += 0.25

fig, ax = plt.subplots()
ax.yaxis.set_ticks(np.arange(7, 11.5, 0.5))
plt.xlim([0,3])
plt.ylim([7,11])

plt.xlabel("Redshift")
plt.ylabel("Stellar Mass, $10^{10}$ $M_{☉}$")
plt.title("Stellar Mass vs Redshift for Milky Way Progenitors", fontweight="bold")

plt.plot(growth_curve[0], growth_curve[2], "#FF7969")
plt.plot(growth_curve[0], growth_curve[1], "#E71E06", label="ZFOURGE")
plt.plot(growth_curve[0], growth_curve[3], "#FF7969")
ax.fill_between(x=growth_curve[0], y1=growth_curve[2], y2=growth_curve[3], color="#FF948C")

plt.plot(not_prog_vals[0], not_prog_vals[1], "k", marker="x", linestyle="none", label="JWT Data")

plt.legend()

plt.show()
