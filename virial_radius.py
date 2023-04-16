import matplotlib.pyplot as plt
import numpy as np
import NDpredict as ndp
import scipy
from astropy.cosmology import Planck18 as cosmo
from astropy import units as u
from astropy.constants import M_sun, pc
from astropy.cosmology import WMAP9 as cosmo

def find_virial_radius(mass,rshift):
    
    omega_m0=cosmo.Om(0)
    omega_de0=cosmo.Ode(0)
    H0=np.divide(cosmo.H(0)*10**3, pc.value*10**6)

    G=scipy.constants.G

    hz=np.add(rshift,1)       
    hz=np.power(hz,3)    
    hz=np.multiply(hz,omega_m0)
    hz=np.add(omega_de0, hz)
    hz=np.power(hz,0.5)
    hz=np.multiply(H0, hz)
    
    mw_denom = 100*(H0*H0)*(omega_m0+omega_de0)

    mw_num = G* 10**(12.20)*M_sun

    mw_r = ((mw_num.value/mw_denom.value)**(1/3))*(1/(pc.value*1e3))


    hz_sq_100=np.multiply(100,(np.power(hz,2)))
    
    Mh_vals = np.multiply(np.power(10,mass),M_sun)

    numerator=np.multiply(Mh_vals, G)
    
    vradius_m=np.power(np.divide(numerator,hz_sq_100),(1/3))

    vradius_kpc=np.divide(vradius_m.value,pc.value*10**3)

    kpc_per_arcmin = cosmo.kpc_proper_per_arcmin(rshift)

    vradius_arcmin = np.divide(vradius_kpc, kpc_per_arcmin)
    vradius_deg = np.divide(vradius_arcmin, 60)


    # fig, ax = plt.subplots()
    # plt.rcParams["font.family"] = "serif"
    # font = {'fontname':'serif'} 
    # plt.grid(color='silver', linestyle='--', linewidth=1, zorder=1)
    # plt.scatter(mass, vradius_kpc, c=rshift,s=10,marker=".", zorder=2)

    # plt.xlabel("$Log_{10}$($M_{h}$/$M_{☉}$)", fontdict=font)
    # plt.ylabel("$R_{virial}$  $(kpc)$", **font)

    # cbar = plt.colorbar()
    # cbar.ax.set_ylabel("Redshift")
    # cbar.ax.set_yticks([np.amin(rshift),0.5,1,1.5,2,2.5,3,3.5,4])
    # plt.xlim(10.6,12.2)
    # x_ticklabels = plt.gca().get_xticklabels()
    # y_ticklabels = plt.gca().get_yticklabels()

    # for tick_label in (x_ticklabels + y_ticklabels):
    #     tick_label.set_fontname("serif")

    # plt.savefig('plots/virial_radius.png', dpi=300)
    # plt.ylim(20,200)
    # plt.show()
     
    return vradius_deg.value
