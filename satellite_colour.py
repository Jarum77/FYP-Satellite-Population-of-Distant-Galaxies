import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import BoundaryNorm
from scipy.stats import ks_2samp
from statsmodels.distributions.empirical_distribution import ECDF
from scipy.spatial.distance import cdist

    
def satellite_colour_compare(sat_data, field_data ):

    sat_mvals = sat_data[:,0]
    sat_zvals = sat_data[:,1]
    sat_UVcol = sat_data[:,2]
    sat_ids = sat_data[:,3]

    field_mvals = field_data[:,0]
    field_zvals = field_data[:,1]
    field_UVcol = field_data[:,2]  

    sat_ids, unique_indices = np.unique(sat_ids, return_index=True)

    # Use the indices to get the unique values from arrays b and c
    sat_mvals = sat_mvals[unique_indices]
    sat_zvals = sat_zvals[unique_indices]
    sat_UVcol = sat_UVcol[unique_indices]

    #Apply mass lim
    mmax = 13
    mmin = 0

    sat_UVcol = sat_UVcol[np.logical_and(sat_mvals >= mmin, sat_mvals < mmax)]
    sat_zvals = sat_zvals[np.logical_and(sat_mvals >= mmin, sat_mvals < mmax)]
    sat_mvals = sat_mvals[np.logical_and(sat_mvals >= mmin, sat_mvals < mmax)]

    # field_UVcol = field_UVcol[np.logical_and(field_mvals >= mmin, field_mvals < mmax)]
    # field_zvals = field_zvals[np.greater(field_mvals, mmin)]
    # field_mvals = field_mvals[np.greater(field_mvals, mmin)]    

    #Apply redshift lim
    zmin = 2.5
    zmax = 4
    sat_mvals = sat_mvals[np.logical_and(sat_zvals >= zmin, sat_zvals < zmax)]
    sat_UVcol = sat_UVcol[np.logical_and(sat_zvals >= zmin, sat_zvals < zmax)]
    sat_zvals = sat_zvals[np.logical_and(sat_zvals >= zmin, sat_zvals < zmax)]

    # field_mvals = field_mvals[np.less(field_zvals, zmax)]
    # field_UVcol = field_UVcol[np.less(field_zvals, zmax)]
    # field_zvals = field_zvals[np.less(field_zvals, zmax)]


    sat_d = np.column_stack((sat_zvals, sat_mvals))
    field_d = np.column_stack((field_zvals, field_mvals))

    distances = cdist(sat_d, field_d)

    closest_indices = np.argmin(distances, axis=1)

    closest_mass = field_mvals[closest_indices]
    closest_zvals = field_zvals[closest_indices]
    closest_UVcol = field_UVcol[closest_indices]

    '''
        Plots displaying probability distribution of the satellite vs all data
    '''
    stat, p_value = ks_2samp(sat_UVcol, closest_UVcol)
    print("Kolmogorov-Smirnov two-sample test statistic: ", stat)
    print("2 sample P-value: ", p_value)

    fig, ax1 = plt.subplots()

    bin_width = 0.1

    plt.rcParams["font.family"] = "serif"
    font = {'fontname':'serif'} 

    sat_counts, sat_bins = np.histogram(sat_UVcol, bins=np.arange(-1.1, 2, bin_width))  
    sat_prob_density = sat_counts / np.sum(sat_counts)

    field_counts, field_bins = np.histogram(closest_UVcol, bins=np.arange(-1.1, 2, bin_width))
    field_prob_density = field_counts / np.sum(field_counts)

    ax1.plot(field_bins[:-1], field_prob_density, drawstyle="steps-post", label="Host Dens.", color="darkturquoise") 
    ax1.plot(sat_bins[:-1], sat_prob_density, drawstyle="steps-post", label="Satellite Dens.", color="mediumorchid")   
    ax1.fill_between(x=field_bins[:-1], y1=0, y2=field_prob_density, step="post", color="darkturquoise", alpha=0.1, zorder=2)
    ax1.fill_between(x=sat_bins[:-1], y1=0, y2=sat_prob_density, step="post", color="mediumorchid", alpha=0.1, zorder=2)
    ax1.set_xlabel("(U-V) rest Color", **font)
    ax1.set_ylim(0,0.3)
    ax1.set_xlim(-1.1,2)
    x_ticklabels = ax1.get_xticklabels()
    y_ticklabels = ax1.get_yticklabels()
    for tick_label in (x_ticklabels + y_ticklabels):
        tick_label.set_fontname("serif")
    ax1.set_ylabel("Probability Density", **font)

    ax2 = ax1.twinx()

    ecdf_sat = ECDF(sat_UVcol)
    ecdf_field = ECDF(field_UVcol)

    ax2.step(ecdf_field.x, ecdf_field.y, ls="--", label='Host Dist.', color="darkturquoise")
    ax2.step(ecdf_sat.x, ecdf_sat.y, ls="--", label='Satellite Dist.', color="mediumorchid")

    ax2.set_ylim(0, 1)
    x_ticklabels = ax2.get_xticklabels()
    y_ticklabels = ax2.get_yticklabels()
    for tick_label in (x_ticklabels + y_ticklabels):
        tick_label.set_fontname("serif")
    ax2.set_ylabel("Cumulative Distribution", **font)
    handles, labels = [], []
    for ax in [ax1, ax2]:
        for h, l in zip(*ax.get_legend_handles_labels()):
            handles.append(h)
            labels.append(l)
    ax1.legend(handles, labels, loc="center right")
    #plt.savefig("plots/sat_col_compare.png", dpi=300, bbox_inches='tight')
    #plt.savefig("plots/sat_col_range_compare.png", dpi=300, bbox_inches='tight')
    plt.show()
    