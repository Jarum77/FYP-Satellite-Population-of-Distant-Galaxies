import numpy as np
from mpmath import mp
import NDpredict as ndp
import matplotlib.pyplot as plt

zi = 0.
zf = 3.
M0 = 10.5   # All stellar masses are assumed to be logarithmic

zvals = []
mvals =[]

zstep = 0.5
z = 0

while (z <= zf):

    M0 = ndp.newmass(M0, zi, z, massfunc='zfourge')
    zvals.append(z)
    mvals.append(M0)
    z += zstep

#Graph Plotting

fig, ax = plt.subplots()
ax.yaxis.set_ticks(np.arange(0, 11.5, 0.5))
plt.plot(zvals, mvals, "r")
plt.xlim([0,3])
plt.ylim([0,11])
plt.xlabel("Z-Index")
plt.ylabel("Stellar Mass, $10^{10}$ $M_{☉}$")
plt.title("Stellar Mass vs Redshift for Milky Way Progenitors", fontweight="bold")
plt.show()

