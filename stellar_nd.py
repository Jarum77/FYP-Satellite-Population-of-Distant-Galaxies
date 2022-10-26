from mpmath import mp
import NDpredict as ndp
z0 = 2.
zf = 0.
M0 = 10.5   # All stellar masses are assumed to be logarithmic
print(ndp.newmass(M0, z0, zf, massfunc='zfourge') )