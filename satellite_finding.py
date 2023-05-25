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
from scipy.stats import pearsonr


def satellite_finding(prog_data, all_data, pz_vals, zgrid, zfourge_min, mmf_popt):

    prog_ids = prog_data[:,0]
    prog_mvals = prog_data[:,1]
    prog_ra = prog_data[:,2]
    prog_dec = prog_data[:,3]
    zvals = prog_data[:,4]
    prog_UVcol = prog_data[:,5]
    vradius = prog_data[:,6]
    prog_sfr = prog_data[:,7]

    all_ids = all_data[:,0]
    all_mvals = all_data[:,1]
    all_ra = all_data[:,2]
    all_dec = all_data[:,3]
    all_zvals = all_data[:,4]
    all_uvcol = all_data[:,5]
    all_sfr = all_data[:,6]
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

    # prog_sfr_test = prog_sfr[found_prog_indices]
    # sat_sfr_test = all_sfr[found_sats_indices]

    # prog_zvals = zvals[found_prog_indices]
    '''
        Probability Density progenitor satellie
    '''

    # index = 4

    # multiplied = np.multiply(found_prog_pz[index],found_sats_pz[index])
    # cumulative = np.cumsum(multiplied)
    # plt.plot(zgrid, found_prog_pz[index], label=f"Progenitor, id = {found_prog_ids[index]:.0f}", color="royalblue")
    # plt.plot(zgrid, found_sats_pz[index], label=f"Satellite, id = {found_sats_ids[index]:.0f}", color="mediumorchid")
    # plt.plot(zgrid, multiplied, label="Product Probability", color="darkslategrey")
    # plt.plot(zgrid, cumulative, label="Cumulative Distribution", color="black")
    # plt.fill_between(x=zgrid, y1=0, y2=multiplied, color="darkslategrey", alpha=0.3, zorder=2)

    # plt.xlim(0.5,1.5)
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
    # plt.savefig('plots/prob_dens_accept.png', dpi=300, bbox_inches='tight' )

    # plt.show()


    '''
        Finding unique progenitor ids and indices through np.unique,
        then obtaining the masses, z values, and colours for our progeitors
        and their satellites
    '''

    sat_ids = all_ids[found_sats_indices]
    sat_mvals = all_mvals[found_sats_indices]
    sat_zvals = all_zvals[found_sats_indices]
    sat_UVcol = all_uvcol[found_sats_indices]
    sat_sfr = all_sfr[found_sats_indices]

    '''
        sSFR comparison
    '''
    zmin = .6
    zmax = 1.6
    prog_comp_sfr = prog_sfr[found_prog_indices]
    prog_comp_mvals = prog_mvals[found_prog_indices]
    prog_comp_zvals = zvals[found_prog_indices]

    unique_indices = np.unique(found_prog_indices)
    median_ssfh_vals = []
    prog_ssfh_vals = []

    mmin = 9.3
    mmax = 15
    for id in unique_indices:
        if (zvals[id] > zmin and zvals[id] < zmax):

            indices = np.where(np.array(found_prog_indices) == id)[0]

            sat_indices = [found_sats_indices[i] for i in indices]
    
            sfh_range = all_sfr[sat_indices]
            mass_range = all_mvals[sat_indices]

            sfh_range = sfh_range[np.logical_and(mass_range >= mmin, mass_range < mmax)]
            mass_range = mass_range[np.logical_and(mass_range >= mmin, mass_range < mmax)]

            ssfr_range = np.divide(np.power(10, sfh_range), np.power(10,mass_range))
            median_ssfh_vals = np.append(median_ssfh_vals, ssfr_range)

            unique_prog_sfr = prog_sfr[id]
            unique_prog_mval = prog_mvals[id]

            prog_ssfr = np.divide(np.power(10, unique_prog_sfr), np.power(10,unique_prog_mval))

            for i in range(len(ssfr_range)):
                prog_ssfh_vals = np.append(prog_ssfh_vals, prog_ssfr) 

    ##Sort data
    sorted_indices = np.argsort(prog_ssfh_vals)
    prog_ssfh_vals = prog_ssfh_vals[sorted_indices]
    median_ssfh_vals = median_ssfh_vals[sorted_indices]

    # perform correlation test
    r, p_value = pearsonr(prog_ssfh_vals, median_ssfh_vals)

    '''
        Progenitor Mass / Greatest Satellite MASS plot. Merged with Number of satellites per kpc^3.
    '''

    unique_prog_indices, sat_counts = np.unique(found_prog_indices, return_counts=True)
    # unique_prog_ids, sat_counts = np.unique(found_prog_ids, return_counts=True)

    prog_ids = prog_ids[unique_prog_indices]
    prog_mvals = prog_mvals[unique_prog_indices]
    prog_zvals = zvals[unique_prog_indices]
    prog_UVcol = prog_UVcol[unique_prog_indices]
    prog_sfr = prog_sfr[unique_prog_indices]


    # fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(10, 5))
    # fig.subplots_adjust(wspace=0.3)

    prog_sat_mvals = np.split(sat_mvals, np.cumsum(sat_counts[:-1]))

    max_func = np.vectorize(np.max)
    max_vals = max_func(prog_sat_mvals)
    mass_ratios = np.divide(np.power(prog_mvals,10), np.power(max_vals,10))
    log_ratios = np.log10(mass_ratios)

    # Find Minimum possible mass ratio
    def mmf(x,a,b,c,f):
        return a*np.log10(b*x+c)+f
    
    z_steps = np.linspace(0, 4, 1000)

    ##Single plot

    cm = plt.cm.get_cmap('seismic')
    plt.rcParams["font.family"] = "serif"
    font = {'fontname':'serif'} 
    plt.grid(color='silver', linestyle='--', linewidth=1, zorder=1)
    
    sc2 = plt.scatter(prog_mvals, log_ratios, c=prog_UVcol, cmap=cm,s=12, zorder=2)
    sc1 = plt.scatter(10.78, np.log10(np.divide(np.power(10, 10.78),np.power(10, 9.43))), marker="^", color="magenta", s=25, label="Milky Way / LMC", zorder=2)
    plt.legend(loc="upper left")
    plt.xlabel("Stellar Mass [$Log_{10}$($M_{*}$/$M_{☉}$)]", **font)
    plt.ylabel("$Log_{10}$(Progenitor Mass / Max Satellite Mass)]", **font)
    plt.xlim(8,11)
    plt.ylim(0,1.75)
    ticks_colbar = np.arange(-.2,2.2,0.2)

    cbar = plt.colorbar(sc2, ticks=ticks_colbar,label="(U-V) Colour", cmap=cm)
    plt.yticks(labels=np.arange(0,2,0.25),ticks=np.arange(0,2,0.25))
    plt.savefig("plots/mass_ratio", dpi=300, bbox_inches="tight")
    plt.show()

 
    # # perform correlation test
    # r, p_value = pearsonr(prog_mvals, log_ratios)
    # print("Pearson correlation for mass ratio: \n")
    # print(r)
    # print(np.amax(log_ratios))
    # cm = plt.cm.get_cmap('seismic')
    # plt.rcParams["font.family"] = "serif"
    # font = {'fontname':'serif'} 
    # ax1.grid(color='silver', linestyle='--', linewidth=1, zorder=1)
    
    # ax1.scatter(prog_mvals, log_ratios, c=prog_UVcol, cmap=cm,s=12 ,vmin=-0.5, vmax=2, zorder=2)
    # ax1.scatter(10.78, np.log10(np.divide(np.power(10, 10.78),np.power(10, 9.43))),color="magenta", s=25, marker="x", label="Milky Way / LMC", zorder=2)
    # ax1.legend()
    # #ax1.plot(zfourge_min(z_steps), max_ratio, label="Maximum Detectable Ratio", ls="--", color="black")
    # #ax1.fill_between(x=zfourge_min(z_steps), y1=0, y2=min_ratio, color="black", alpha=0.2, zorder=2)
    # ax1.set_xlabel("Stellar Mass [$Log_{10}$($M_{*}$/$M_{☉}$)]", **font)
    # ax1.set_ylabel("$Log_{10}$(Progenitor Mass / Max Satellite Mass)]", **font)
    # ax1.set_xlim(8,11)
    # ax1.set_ylim(0,1.75)
    # ax1.set_yticks(labels=np.arange(0,2,0.25),ticks=np.arange(0,2,0.25))


    # kpc_per_arcmin = cosmo.kpc_proper_per_arcmin(zvals)
    # vradius_arcmin = np.multiply(vradius, 60)
    # vradius_kpc = np.multiply(vradius_arcmin, kpc_per_arcmin)
    
    # prog_virial_area = np.multiply(np.pi, np.power(vradius_kpc[unique_prog_indices].value,2))

    # sats_per_area = np.divide(sat_counts, prog_virial_area)

    # sats_per_vol_scaled = np.log10(sats_per_area)
    # plt.rcParams["font.family"] = "serif"
    # font = {'fontname':'serif'} 
    # ax2.grid(color='silver', linestyle='--', linewidth=1, zorder=1)

    # ## perform correlation test
    # r, p_value = pearsonr(prog_mvals, sats_per_vol_scaled)
    # print("Pearson correlation for area: \n")
    # print(r)
    # sc2 = ax2.scatter(prog_mvals, sats_per_vol_scaled, c=prog_UVcol, cmap=cm, s=12, zorder=2, vmin=-0.5, vmax=2)  

    # ax2.set_xlabel("Stellar Mass [$Log_{10}$($M_{*}$/$M_{☉}$)]", **font)
    # ax2.set_ylabel("Satellite Count [$Log_{10}$( $(kpc)^{-2}$)]", **font)
    # ax2.set_xlim(8,11)
    # ax2.set_ylim(-5,-3)
    # x_ticklabels = ax1.get_xticklabels()
    # y_ticklabels = ax1.get_yticklabels()
    # for tick_label in (x_ticklabels + y_ticklabels):
    #     tick_label.set_fontname("serif")

    # x_ticklabels = ax2.get_xticklabels()
    # y_ticklabels = ax2.get_yticklabels()
    # for tick_label in (x_ticklabels + y_ticklabels):
    #     tick_label.set_fontname("serif")


    # ticks_colbar = np.arange(-.5,2.2,0.2)
    # cbar = fig.colorbar(sc2, ax=[ax1, ax2], fraction=0.1, orientation="horizontal", ticks=ticks_colbar,shrink=0.6, label="(U-V) Colour", cmap=cm, pad=0.15)

    # plt.savefig("plots/satellites_per_area_ratio", dpi=300, bbox_inches='tight')
    # plt.show()


    found_prog_data = np.c_[ prog_ids, prog_mvals , prog_zvals , prog_UVcol, prog_sfr]
    found_sat_data = np.c_[ sat_ids, sat_mvals , sat_zvals, sat_UVcol, sat_sfr]

    return found_prog_data, found_sat_data

    