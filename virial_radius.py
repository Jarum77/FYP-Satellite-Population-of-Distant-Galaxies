import matplotlib.pyplot as plt
import numpy as np
import NDpredict as ndp
import scipy
from astropy.cosmology import Planck18 as cosmo
from astropy import units as u
from astropy.constants import M_sun, pc
from astropy.cosmology import WMAP9 as cosmo

def find_virial_radius(stell_mass, halo_mass, zvals):
    
    omega_m0=cosmo.Om(0)
    omega_de0=cosmo.Ode(0)
    H0=np.divide(cosmo.H(0)*10**3, pc.value*10**6)

    G=scipy.constants.G

    hz=np.add(zvals,1)       
    hz=np.power(hz,3)    
    hz=np.multiply(hz,omega_m0)
    hz=np.add(omega_de0, hz)
    hz=np.power(hz,0.5)
    hz=np.multiply(H0, hz)
    
    mw_denom = 100*(H0*H0)*(omega_m0+omega_de0)

    mw_num = G* 10**(12.20)*M_sun

    mw_r = ((mw_num.value/mw_denom.value)**(1/3))*(1/(pc.value*1e3))


    hz_sq_100=np.multiply(100,(np.power(hz,2)))
    
    Mh_vals = np.multiply(np.power(10,halo_mass),M_sun)

    numerator=np.multiply(Mh_vals, G)
    
    vradius_m=np.power(np.divide(numerator,hz_sq_100),(1/3))

    vradius_kpc=np.divide(vradius_m.value,pc.value*10**3)

    kpc_per_arcmin = cosmo.kpc_proper_per_arcmin(zvals)

    vradius_arcmin = np.divide(vradius_kpc, kpc_per_arcmin)
    vradius_deg = np.divide(vradius_arcmin, 60)

    '''
        Plot the Progenitor virial radius and halo halo_mass
    '''

    # fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(10, 6))
    # fig.subplots_adjust(wspace=0.3)
    # cm = plt.cm.get_cmap('viridis')
    # ## Halo Mass

    # plt.rcParams["font.family"] = "serif"
    # font = {'fontname':'serif'} 
    # ax1.set_yticks(np.arange(7.5, 11.5, 0.5))
    # ax1.set_xticks(np.arange(10.5,12.5,0.25))
    # ax1.grid(color='silver', linestyle='--', linewidth=1, zorder=1)

    # ax1.scatter(halo_mass, stell_mass, c=zvals, cmap=cm, s=15,marker=".", zorder=2,vmax=4)
    # ax1.set_xlabel("Halo Mass [$Log_{10}$($M_{DM}$/$M_{☉}$)]",**font)
    # ax1.set_ylabel("Stellar Mass [$Log_{10}$($M_{*}$/$M_{☉}$)]",**font)

    # x_ticklabels = ax1.get_xticklabels()
    # y_ticklabels = ax1.get_yticklabels()
    # for tick_label in (x_ticklabels + y_ticklabels):
    #     tick_label.set_fontname("serif")
    # ax1.set_xlim(10.75,12.5)
    # ax1.set_ylim(8,11)


    # ## Virial Radius
    # ax2.grid(color='silver', linestyle='--', linewidth=1, zorder=1)
    # sc2 = ax2.scatter(halo_mass, vradius_kpc, c=zvals,s=15, cmap=cm, marker=".", zorder=2, vmax=4)

    # ax2.set_xlabel("Halo Mass [$Log_{10}$($M_{DM}$/$M_{☉}$)]", fontdict=font)
    # ax2.set_ylabel("Virial Radius [kpc]", **font)
    # print(np.amax(vradius_kpc))
    # print(np.amin(vradius_kpc))

    # ticks_colbar = [0.5,1,1.5,2,2.5,3,3.5,4]
    # cbar = fig.colorbar(sc2, ax=[ax1, ax2], orientation="horizontal", ticks=ticks_colbar,shrink=0.7, label="Redshift", cmap=cm, pad=0.15)
    # ax2.set_xlim(10.75,12.5)
    # ax2.set_ylim(0,250)
    # x_ticklabels = ax2.get_xticklabels()
    # y_ticklabels = ax2.get_yticklabels()

    # for tick_label in (x_ticklabels + y_ticklabels):
    #     tick_label.set_fontname("serif")

    # plt.savefig('plots/virial_radius_halo_mass.png', dpi=300, bbox_inches='tight' )

    # plt.show()
     
    return vradius_deg.value
