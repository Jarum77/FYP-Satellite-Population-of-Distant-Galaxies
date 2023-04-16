import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import BoundaryNorm
from scipy.stats import ks_2samp
from statsmodels.distributions.empirical_distribution import ECDF


def colour_plot(prog_UVcol, UVcol, mass_vals, allmvals, zvals, allzvals):
    cm = plt.cm.get_cmap('RdYlBu_r')
    #colour_plot2(zvals, allzvals, UVcol, prog_UVcol,cm)
    
    #fig=plt.figure(figsize=(10,6))
#    plt.scatter(UVcol, allmvals, c='k',s=2)
#    plt.scatter(prog_UVcol, mass_vals, c='r',s=2)

    # plt.scatter(allzvals,allmvals,c='grey',s=2)
    # plt.xlabel('$z$')
    # plt.ylabel('$log M_{stellar}$')
    # plt.scatter(zvals,mass_vals,c=prog_UVcol, cmap=cm,s=3)
    # plt.colorbar()
    # plt.show()
    
def median_colour_plot(found_prog_data, found_sat_data, all_data ):

    prog_ids = found_prog_data[:,0]
    prog_mvals = found_prog_data[:,1]
    prog_zvals = found_prog_data[:,2]
    prog_UVcol = found_prog_data[:,3]

    sat_ids = found_sat_data[:,0]
    sat_mvals = found_sat_data[:,1]
    sat_zvals = found_sat_data[:,2]
    sat_UVcol = found_sat_data[:,3]

    all_ids = all_data[:,0]
    all_mvals = all_data[:,1]
    all_zvals = all_data[:,4]
    all_UVcol = all_data[:,5]

    #Remove all of the satellites from all data
    not_sats = np.in1d(all_ids, sat_ids, invert=True)
    
    all_ids = all_ids[not_sats]
    all_mvals = all_mvals[not_sats]
    all_zvals = all_zvals[not_sats]
    all_UVcol = all_UVcol[not_sats]

    prog_median_colours = []
    sat_median_colours = []
    all_median_colours = []

    prog_median_mvals = []
    sat_median_mvals = []
    all_median_mvals = []

    z_step = 0.5
    z_steps = np.arange(0.000, 4.00, z_step)

    m_step = 0.5,
    m_steps = np.arange(np.amin(all_mvals), 10.5, m_step)

    z_vals = []
    m_vals = []

    for z in z_steps:

        prog_UVcol_range = prog_UVcol[np.logical_and(prog_zvals >= z, prog_zvals <= z+z_step)]
        prog_mvals_range = prog_mvals[np.logical_and(prog_zvals >= z, prog_zvals <= z+z_step)] 
        #prog_zvals_range = prog_zvals[np.logical_and(prog_zvals >= z, prog_zvals <= z+z_step)] 

        sat_UVcol_range = sat_UVcol[np.logical_and(sat_zvals >= z, sat_zvals <= z+z_step)]
        sat_mvals_range = sat_mvals[np.logical_and(sat_zvals >= z, sat_zvals <= z+z_step)]   
        #sat_zvals_range = sat_zvals[np.logical_and(sat_zvals >= z, sat_zvals <= z+z_step)]

        all_UVcol_range = all_UVcol[np.logical_and(all_zvals >= z, all_zvals <= z+z_step)]
        all_mvals_range = all_mvals[np.logical_and(all_zvals >= z, all_zvals <= z+z_step)]   
        #all_zvals_range = all_zvals[np.logical_and(all_zvals >= z, all_zvals <= z+z_step)]

        for m in m_steps: 
            
            prog_UVcol_range2 = prog_UVcol_range[np.logical_and(prog_mvals_range >= m, prog_mvals_range <= m+m_step)]
            sat_UVcol_range2 = sat_UVcol_range[np.logical_and(sat_mvals_range >= m, sat_mvals_range <= m+m_step)]
            all_UVcol_range2 = all_UVcol_range[np.logical_and(all_mvals_range >= m, all_mvals_range <= m+m_step)]

            prog_avg_colour = np.median(prog_UVcol_range2)
            sat_avg_colour = np.median(sat_UVcol_range2)
            all_avg_colour = np.median(all_UVcol_range2)

            prog_median_colours = np.append(prog_median_colours,prog_avg_colour)
            sat_median_colours = np.append(sat_median_colours,sat_avg_colour) 
            all_median_colours = np.append(all_median_colours,all_avg_colour) 

            z_vals = np.append(z_vals, z)
            m_vals = np.append(m_vals, m)
    
    prog_z_points = z_vals[~np.isnan(prog_median_colours)]
    prog_m_points = m_vals[~np.isnan(prog_median_colours)]
    prog_median_colours = prog_median_colours[~np.isnan(prog_median_colours)]

    sat_z_points = z_vals[~np.isnan(sat_median_colours)]
    sat_m_points = m_vals[~np.isnan(sat_median_colours)]
    sat_median_colours = sat_median_colours[~np.isnan(sat_median_colours)]

    all_z_points = z_vals[~np.isnan(all_median_colours)]
    all_m_points = m_vals[~np.isnan(all_median_colours)]
    all_median_colours = all_median_colours[~np.isnan(all_median_colours)]

    fig, (ax1, ax2) = plt.subplots(ncols=2, figsize=(10, 6))
    #cm = plt.cm.get_cmap('RdYlBu_r')

    '''
        Colour values of all data and satellite
    '''
    # sc1 = ax1.scatter(all_z_points, all_m_points, c=all_median_colours, edgecolor='grey', cmap='RdYlBu_r',marker="o", s=60, linewidths=1 )
    # ax1.set_title('All Ceers Data')
    # ax1.set_ylabel("$Log_{10}$($M_{*}$/$M_{☉}$)") 
    # ax1.set_xlabel('Redshift')
    # ax1.set_yticks(np.arange(3,12,1))
    # ax1.set_yticklabels(np.arange(3,12,1))

    # ##Satellite data subplot
    # sc2 = ax2.scatter(sat_z_points, sat_m_points, c=sat_median_colours, edgecolor='grey', cmap='RdYlBu_r', marker="o", s=60, linewidths=1 )
    # ax2.set_title('Satellite Data')
    # ax2.set_ylabel("$Log_{10}$($M_{*}$/$M_{☉}$)") 
    # ax2.set_xlabel('Redshift')
    # ax2.set_yticks(np.arange(3,12,1))
    # ax2.set_yticklabels(np.arange(3,12,1))

    # ##set the range of values for the color bar based on the maximum and minimum values in both sets of data
    # vmin = min(np.min(all_median_colours), np.min(sat_median_colours))
    # vmax = max(np.max(all_median_colours), np.max(sat_median_colours))

    # ##create the shared color bar

    # cmap = plt.get_cmap('RdYlBu_r')
    # bounds = np.arange(0, vmax+0.2, 0.2)
    # norm = BoundaryNorm(bounds, cmap.N)
    # cbar = fig.colorbar(sc2, ax=[ax1, ax2], fraction=0.05, orientation="horizontal", shrink=0.6, boundaries=bounds, ticks=bounds, label="(U-V) rest Colour",extend='both')

    # plt.show()
    '''
        Plots displaying probability distribution of the satellite vs all data
    '''
    stat, p_value = ks_2samp(all_median_colours, sat_median_colours)
    print("Kolmogorov-Smirnov two-sample test statistic: ", stat)
    print("2 sample P-value: ", p_value)

    sat_counts, sat_bins = np.histogram(sat_median_colours, bins=np.arange(0, max(sat_median_colours)+0.1, 0.1))   
    sat_bin_centers = 0.5 * (sat_bins[1:] + sat_bins[:-1])
    sat_prob_density = sat_counts / np.sum(sat_counts)

    all_counts, all_bins = np.histogram(all_median_colours, bins=np.arange(0, max(all_median_colours)+0.1, 0.1))
    all_bin_centers = 0.5 * (all_bins[1:] + all_bins[:-1])
    all_prob_density = all_counts / np.sum(all_counts)


    ax1.plot(sat_bin_centers, sat_prob_density, marker="o", label="Satellite Data", color="royalblue")
    ax1.plot(all_bin_centers, all_prob_density, marker="o", label="All Data", color="indianred")
    ax1.set_xlabel("(U-V) rest Color")
    ax1.set_ylabel("Probability Density")
    ax1.set_ylim(0, 0.4)
    ax1.set_xlim(0,1.75)

    ecdf_sat = ECDF(sat_median_colours)
    ecdf_all = ECDF(all_median_colours)

    ax2.step(ecdf_sat.x, ecdf_sat.y, label='Satellite Data', color="royalblue")
    ax2.step(ecdf_all.x, ecdf_all.y, label='All Data', color="indianred")
    ax2.set_xlabel("(U-V) rest Color")
    ax2.set_ylabel("Empirical Cumulative Density")
    ax2.set_ylim(0, 1)
    ax2.set_xlim(0,1.75)

    handles, labels = ax1.get_legend_handles_labels()
    fig.legend(handles, labels,loc="upper center", ncol=2, fancybox=True)

    plt.show()
    
        