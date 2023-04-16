import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
from astropy.coordinates import SkyCoord
from astropy.coordinates import ICRS, Galactic
from astropy.coordinates import Angle, Latitude, Longitude
from astropy import units as u
from astropy.cosmology import Planck18 as cosmo
from astropy.cosmology import WMAP9
import astropy.cosmology.units as cu
import astropy.coordinates as astropy_cords
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.ticker as mtick


def satellite_finding(prog_data, all_data, pz_vals, zgrid):

    prog_ids = prog_data[:,0]
    prog_mvals = prog_data[:,1]
    prog_ra = prog_data[:,2]
    prog_dec = prog_data[:,3]
    zvals = prog_data[:,4]
    prog_UVcol = prog_data[:,5]
    vradius = prog_data[:,6]
    prog_AU = prog_data[:,7]
    prog_AV = prog_data[:,8]

    all_ids = all_data[:,0]
    all_mvals = all_data[:,1]
    all_ra = all_data[:,2]
    all_dec = all_data[:,3]
    all_zvals = all_data[:,4]
    all_uvcol = all_data[:,5]
    
    prog_vals=SkyCoord(prog_ra,prog_dec,frame='icrs', unit='deg')

    all_vals=SkyCoord(all_ra,all_dec,frame='icrs', unit='deg')

    idx1, idx2, sep2d, dist3d = astropy_cords.search_around_sky(prog_vals, all_vals, np.amax(vradius)*u.deg)
    max = np.amax(np.amax(vradius)*u.deg)

    #Remove all values above the Virial radius and match the index for the ids
    potential_prog_ids = prog_ids[idx1[np.less(0,np.subtract(vradius[idx1], sep2d.value))]]
    potential_sats_ids = all_ids[idx2[np.less(0,np.subtract(vradius[idx1], sep2d.value))]]
    potential_prog_pz = pz_vals[idx1[np.less(0,np.subtract(vradius[idx1], sep2d.value))]]
    potential_sats_pz = pz_vals[idx2[np.less(0,np.subtract(vradius[idx1], sep2d.value))]]
    ang_sep = sep2d[np.less(0,np.subtract(vradius[idx1], sep2d.value))]
    potential_sats_indices = idx2[np.less(0,np.subtract(vradius[idx1], sep2d.value))]
    potential_prog_indices = idx1[np.less(0,np.subtract(vradius[idx1], sep2d.value))]

    #Remove all self-counted distances
    potential_prog_ids = potential_prog_ids[~np.equal(0, ang_sep)]
    potential_sats_ids = potential_sats_ids[~np.equal(0, ang_sep)]
    potential_prog_pz = potential_prog_pz[~np.equal(0, ang_sep)]
    potential_sats_pz = potential_sats_pz[~np.equal(0, ang_sep)]
    potential_sats_indices = potential_sats_indices[~np.equal(0, ang_sep)]
    potential_prog_indices = potential_prog_indices[~np.equal(0, ang_sep)]
    ang_sep = ang_sep[~np.equal(0, ang_sep)]

    # Remove pairs where a satellite has heavier mass than the main galaxy
    keep_pair = np.greater(prog_mvals[potential_prog_indices], all_mvals[potential_sats_indices])
    potential_prog_ids = potential_prog_ids[keep_pair]
    potential_sats_ids = potential_sats_ids[keep_pair]
    potential_prog_pz = potential_prog_pz[keep_pair]
    potential_sats_pz = potential_sats_pz[keep_pair]
    potential_sats_indices = potential_sats_indices[keep_pair]
    potential_prog_indices = potential_prog_indices[keep_pair]
    ang_sep = ang_sep[keep_pair] 

    ## Histogram of the separation
    # fig, ax = plt.subplots()
    # ax.grid(color='silver', linestyle='--', linewidth=1, zorder=1)
    # plt.rcParams["font.family"] = "serif"
    # font = {'fontname':'serif'} 
    # ax.set_xlabel("Angular Separation (degrees)", **font)
    # ax.set_ylabel("Count", **font)
    # ax.set_ylim(0,20000)
    # ax.yaxis.set_major_formatter(mtick.StrMethodFormatter('{x:,.0f}'))
    # # Set font family for tick labels
    # for tick in ax.get_xticklabels() + ax.get_yticklabels():
    #     tick.set_fontname("serif")
    # hist = ax.hist(ang_sep.value, bins=np.linspace(0, 0.014,35), alpha=1, color="royalblue", edgecolor="black", zorder=2)   
    # axins = inset_axes(ax, width=2, height=1.7, bbox_to_anchor=(0.66, 0.63, 0.3, 0.3), bbox_transform=ax.transAxes)
    # axins.set_xlim(0.008, 0.014)
    # axins.set_ylim(0, 100)
    # axins.grid(color='silver', linestyle='--', linewidth=1, zorder=1)
    # axins.hist(ang_sep.value, bins=np.linspace(0.008, 0.014, 15), alpha=1, color="royalblue", edgecolor="black", zorder=2)
    # plt.savefig('plots/separation_hist.png', dpi=300)
    # plt.show() 

    # Find the Progenitor - Satellite pairs with required probability
    pz_multi = np.multiply(potential_prog_pz, potential_sats_pz)
    # p_threshold = 0.6
    # Criteria is based on probability distribution integral total
    sum_pz = np.sum(pz_multi, axis=1)    

    #criteria = np.less(p_threshold, sum_pz)

    #Criteria is based on probability density peak
    #criteria = np.any(pz_multi >= p_threshold, axis = 1)

    # Combined criteria

    criteria = (np.less(0.4, sum_pz)) & (np.any(pz_multi >= 0.15, axis=1))


    '''
        Vary both thresholds and plot result
    '''

    ##Find how many per criteria

    # pcum_thresholds = np.arange(0.2, 1.0, 0.1)
    # ppeak_thresholds = np.arange(0.05, 0.45, 0.05)
    # count = []

    # for pcum in pcum_thresholds:

    #     row_count = []

    #     for ppeak in ppeak_thresholds:

    #         crit = (np.less(pcum, sum_pz)) & (np.any(pz_multi >= ppeak, axis=1))
    #         found = potential_prog_ids[crit]
    #         row_count = np.append(row_count, len(found))

    #     count.append(row_count)

    # print("Ppeak:  " + str(ppeak_thresholds) + "\n")

    # for i in range(len(pcum_thresholds)):
    #     print(f"{pcum_thresholds[i]:.2f}" + "," + str(count[i]) + "\n")


    
    found_prog_ids = potential_prog_ids[criteria]
    found_prog_indices = potential_prog_indices[criteria]
    found_prog_pz = potential_prog_pz[criteria]
    found_sats_ids = potential_sats_ids[criteria]
    found_sats_indices = potential_sats_indices[criteria]
    found_sats_pz = potential_sats_pz[criteria]
    pz_multi = pz_multi[criteria]
    ang_sep = ang_sep[criteria]


    '''
        Probability Density progenitor satellie
    '''

    # index = 10

    # multiplied = np.multiply(found_prog_pz[index],found_sats_pz[index])
    # cumulative = np.cumsum(multiplied)
    # plt.plot(zgrid, found_prog_pz[index], label=f"Progenitor, id = {found_prog_ids[index]:.0f}", color="royalblue")
    # plt.plot(zgrid, found_sats_pz[index], label=f"Satellite, id = {found_sats_ids[index]:.0f}", color="mediumorchid")
    # plt.plot(zgrid, multiplied, label="Combined Probability", color="darkslategrey")
    # plt.plot(zgrid, cumulative, label="Cumulative Distribution", color="black", ls="--")
    # plt.xlim(0,5)
    # plt.ylim(-0.1,1)
    # plt.rcParams["font.family"] = "serif"
    # font = {'fontname':'serif'} 
    # plt.grid(color='silver', linestyle='--', linewidth=1, zorder=1)
    # plt.xlabel("Redshift", **font)
    # plt.ylabel("Probability Density",**font)
    # x_ticklabels = plt.gca().get_xticklabels()
    # y_ticklabels = plt.gca().get_yticklabels()

    # for tick_label in (x_ticklabels + y_ticklabels):
    #     tick_label.set_fontname("serif")
    # plt.legend()
    # plt.savefig('plots/prob_dens_accept.png', dpi=300)

    # plt.show()


    '''
        Finding unique progenitor ids and indices through np.unique,
        then obtaining the masses, z values, and colours for our progeitors
        and their satellites
    '''
    unique_prog_indices, sat_counts = np.unique(found_prog_indices, return_counts=True)
    # unique_prog_ids, sat_counts = np.unique(found_prog_ids, return_counts=True)

    prog_ids = prog_ids[unique_prog_indices]
    prog_mvals = prog_mvals[unique_prog_indices]
    prog_zvals = zvals[unique_prog_indices]
    prog_UVcol = prog_UVcol[unique_prog_indices]
    prog_AU = prog_AU[unique_prog_indices] 
    prog_AV = prog_AV[unique_prog_indices] 


    sat_ids = all_ids[found_sats_indices]
    sat_mvals = all_mvals[found_sats_indices]
    sat_zvals = all_zvals[found_sats_indices]
    sat_UVcol = all_uvcol[found_sats_indices]


    '''
        Progenitor Mass / Greatest Satellite MASS plot.
    '''

    # prog_sat_mvals = np.split(sat_mvals, np.cumsum(sat_counts[:-1]))

    # max_func = np.vectorize(np.max)
    # max_vals = max_func(prog_sat_mvals)
    # mass_ratios = np.divide(np.power(prog_mvals,10), np.power(max_vals,10))
    # log_ratios = np.log10(mass_ratios)
    # cm = plt.cm.get_cmap('RdYlBu_r')
    # plt.scatter(prog_mvals,log_ratios, c=prog_UVcol, cmap=cm,s=3 )
    # plt.xlabel("Log Stellar Mass")
    # plt.ylabel("Log10(Mass Galaxy/ Maximum Mass Satellite)")
    # plt.colorbar(label="(U-V) rest Colour") 
    # plt.show()

    '''
        Number of satellites per kpc^3.
    '''

    # kpc_per_arcmin = cosmo.kpc_proper_per_arcmin(zvals)
    # vradius_arcmin = np.multiply(vradius, 60)
    # vradius_kpc = np.multiply(vradius_arcmin, kpc_per_arcmin)
    
    # prog_virial_vols = np.multiply((4/3)*np.pi, np.power(vradius_kpc[unique_prog_indices].value,3))

    # sats_per_vol = np.divide(sat_counts, prog_virial_vols)

    # sats_per_vol_scaled = np.log10(sats_per_vol)
    # cm = plt.cm.get_cmap('RdYlBu_r')
    # plt.scatter(prog_mvals, sats_per_vol_scaled, c=prog_UVcol, cmap=cm, s=4)
    # plt.colorbar(label="(U-V) rest Colour", ticks=np.arange(0.2, 2, 0.2))
    # plt.xlabel(" $Log_{10}$($M_{*}$/$M_{☉}$)")
    # plt.ylabel("$Log_{10}$(Satellite Count $(kPc)^{-3}$ )")
    # plt.show()

    '''
        Number of satellites for each progenitor, plotted
        as a function of stellar mass.
    '''


        
    '''
        Plotting probability density peaks for a progenitor
        and a satellite pair.
    '''
    #Probability Plotting
    # plt.plot(zgrid, found_prog_pz[index],label="Progenitor id=%i" % found_prog_ids[index], color="darkslateblue")
    # plt.plot(zgrid, found_sat_pz[index], label="Satellite id=%i" % found_sat_ids[index], color="cadetblue")
    # plt.plot(zgrid, pz_multi[index], label="Combined Probability", ls="-", color="k")
    # plt.hlines(y=p_threshold, xmin=0, xmax=4, color='black', linestyle='--', label="Threshold P = {0:3.1f}".format(p_threshold))
    # plt.xlim(0,0.5)
    # plt.ylim(-0.1,1)
    # plt.title("Probability Density vs Redshift",fontweight="bold" )
    # plt.xlabel("$Redshift$")
    # plt.ylabel("$P(z)$")
    # plt.legend()
    # plt.show()

    found_prog_data = np.c_[ prog_ids, prog_mvals , prog_zvals , prog_UVcol, prog_AU, prog_AV]
    found_sat_data = np.c_[ sat_ids, sat_mvals , sat_zvals , sat_UVcol]

    return found_prog_data, found_sat_data

    