import numpy as np
from mpmath import mp
import NDpredict as ndp
import matplotlib.pyplot as plt

zi = 0.
zf = 4.
M0 = np.log10(1.5*10**9)   # All stellar masses are assumed to be logarithmic

zvals = []
Mzfourge_vals = []
Millus_vals = []

zstep = 0.25
z = 0

while (z <= zf):

    Mzfourge = ndp.newmass(M0, zi, z, massfunc='zfourge')
    Millus = ndp.newmass(M0, zi, z, massfunc='illustris')

    zvals.append(z)
    Mzfourge_vals.append(Mzfourge)
    Millus_vals.append(Millus)

    z += zstep

#Graph Plotting
fig, ax = plt.subplots()
ax.yaxis.set_ticks(np.arange(8, 11.5, 0.5))
plt.plot(zvals, Mzfourge_vals, "r", label="ZFOURGE")
plt.plot(zvals, Millus_vals, "b",label="Illustris")
plt.xlim([0,4])
plt.ylim([8,11])
plt.legend()
plt.xlabel("Z-Index")
plt.ylabel("Stellar Mass, $10^{10}$ $M_{☉}$")
plt.title("Stellar Mass vs Redshift for Milky Way Progenitors", fontweight="bold")
plt.show()

