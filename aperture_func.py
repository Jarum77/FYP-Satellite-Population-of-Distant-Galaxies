import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from mpl_toolkits import mplot3d

def Halo_Masses(zvals, mvals):

    M0 = 11.339
    Mz = 0.692
    Ep0 = 0.005
    Epz = 0.689
    Beta0 = 3.344
    Betaz = -2.079
    gamma0 = 0.966

    halo_masses = []
    for i in range(len(mvals)):

        M1 = M0 + ( Mz*(zvals[i]/(zvals[i]+1)) ) 
        Beta = Beta0 + ( Betaz*(zvals[i]/(zvals[i]+1)) )
        EpN = Ep0 + ( Epz*(zvals[i]/(zvals[i]+1)) )

        Mh_vals = np.linspace( mvals[i]/2 , mvals[i]*4 ,20)     
        Ep = 2 * EpN * (( (Mh_vals/M1)**-Beta) + ( (Mh_vals/M1)**gamma0))**-1
        MG_vals = Mh_vals*10**-Ep
        f = sp.interpolate.interp1d(MG_vals, Mh_vals)
        halo_masses.append(f(mvals[i]))

    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.scatter3D(zvals, mvals, halo_masses, c=halo_masses, cmap='Reds');
    ax.set_xlabel("Redshift")
    ax.set_ylabel("Stellar Mass")
    ax.set_zlabel("Halo Mass")
    plt.show()


