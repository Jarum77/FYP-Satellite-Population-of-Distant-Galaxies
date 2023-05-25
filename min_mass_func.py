import numpy as np
from mpmath import mp
import NDpredict as ndp
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.io import readsav
from scipy.signal import medfilt
from scipy.optimize import curve_fit

def min_mass_func(zvals,mass_vals, SNR):

    zvals = zvals[np.greater(SNR, 0)]
    mass_vals = mass_vals[np.greater(SNR, 0)]
    SNR = SNR[np.greater(SNR, 0)]

    ratio = 3

    log_mmin_vals = np.add(np.log10(np.divide(ratio, SNR)), mass_vals)

    z_steps = np.arange(0.0001, 4.0001, 0.01)
    z = 0
    percentile_vals = []

    for z in z_steps:

        mass_vals_range = log_mmin_vals[np.greater(zvals, z)]
        SNR_range = SNR[np.greater(zvals, z)]  
        zvals_range = zvals[np.greater(zvals, z)]

        
        mass_vals_range = mass_vals_range[np.less(zvals_range, z+0.1)]
        SNR_range = SNR_range[np.less(zvals_range, z+0.1)]  
        zvals_range = zvals_range[np.less(zvals_range, z+0.1)]

        #mass_vals_range = medfilt(mass_vals_range)

        percentile_vals = np.append(percentile_vals, np.percentile(mass_vals_range, 10))


    # def func(x, a, b, d,c, e):
    #     return a*x**2+c+d*x**3+e*x
    
    
    def func(x,a,b,c,f):
        return a*np.log10(b*x+c)+f
    
    # def func(x,c, d, e, f):
    #     return c*x**3+d*x**2+e*x+f
    p0 = [5, 5, 1, 0]
    popt, pcov = curve_fit(func, z_steps, percentile_vals, p0=p0)
    # print(popt)
    perr = np.sqrt(np.diag(pcov))
    # print(np.log10(perr))
    # print("Minimum mass at z = 4:")
    # print(func(4,*popt))
    ##Plot
    # plt.rcParams["font.family"] = "serif"
    # font = {'fontname':'serif'} 
    # plt.grid(color='silver', linestyle='--', linewidth=1, zorder=1)

    # plt.scatter(z_steps,percentile_vals, s=10, c='k',label="90th Percentiles", zorder=2, marker="x")
    # plt.plot(z_steps,func(z_steps,*popt), label="Minimum Mass Function", zorder=2, color="orchid", linewidth=2)


    # x_ticklabels = plt.gca().get_xticklabels()
    # y_ticklabels = plt.gca().get_yticklabels()

    # for tick_label in (x_ticklabels + y_ticklabels):
    #     tick_label.set_fontname("serif")
    # plt.xlabel("Redshift", **font)
    # plt.ylabel("Minimum Stellar Mass [$Log_{10}$($M_{*}$/$M_{☉}$)]", **font)
    # plt.legend(loc="center right")
    # plt.xlim(0,4)
    # plt.ylim(4,7.5)
    # plt.savefig('plots/minimum_mass_function.png', dpi=300)
    # plt.show()

    return popt