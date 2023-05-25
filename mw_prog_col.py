import matplotlib.pyplot as plt
import numpy as np
from astropy.cosmology import FlatLambdaCDM
import scipy as sp
from scipy.misc import derivative

cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

def MW_Prog_UVcol(prog_data, sat_data, mmf_popt):

    #Assign Data
    prog_mvals = prog_data[:,1]
    prog_zvals = prog_data[:,4]
    prog_UVcol = prog_data[:,5]
    prog_sfr = prog_data[:,7]

    sat_mvals = sat_data[:,0]
    sat_zvals = sat_data[:,1]
    sat_UVcol = sat_data[:,2]

    #Bruzal Charlot SSP
    bc_ssp_data = np.loadtxt("data/SSP/ssp.txt", skiprows=1,  delimiter=',') 

    bruzal_log_age = bc_ssp_data[:,0]
    bruzal_Umag = bc_ssp_data[:,2]
    bruzal_Vmag = bc_ssp_data[:,4]   
    bruzal_Bmag = bc_ssp_data[:,3]

    bruzal_col = np.subtract(bruzal_Umag, bruzal_Vmag)
    bruzal_bvcol = np.subtract(bruzal_Bmag, bruzal_Vmag)
    bruzal_t_vals = np.divide(np.power(10, bruzal_log_age),1e9) 

    bruzal_Umag_func = sp.interpolate.interp1d(bruzal_t_vals, bruzal_Umag, bounds_error=False, fill_value=(np.min(bruzal_Umag),np.max(bruzal_Umag)))  
    bruzal_Vmag_func = sp.interpolate.interp1d(bruzal_t_vals, bruzal_Vmag, bounds_error=False, fill_value=(np.min(bruzal_Vmag),np.max(bruzal_Vmag)))   
    bruzal_Bmag_func = sp.interpolate.interp1d(bruzal_t_vals, bruzal_Bmag, bounds_error=False, fill_value=(np.min(bruzal_Bmag),np.max(bruzal_Bmag)))   


    # Milky Way Star Formation History

    MW_SFR_data = np.loadtxt("data/MW LMC SFH/MW SFR.txt", skiprows=1,  delimiter=',')

    MW_SFR_lookback_t_vals = MW_SFR_data[:,0]
    MW_SFR_vals = MW_SFR_data[:,1]
    MW_SFR_vals = np.multiply(MW_SFR_vals, 1.5610046192724059)
    MW_SFR_t_vals =np.multiply(np.subtract(MW_SFR_lookback_t_vals, np.amax(MW_SFR_lookback_t_vals)),-1)

    MW_SFR_func = sp.interpolate.interp1d(MW_SFR_t_vals, MW_SFR_vals, bounds_error=False, fill_value=(np.min(MW_SFR_vals),np.max(MW_SFR_vals))) 

    # LMC Star Formation History

    LMC_data = np.loadtxt("data/MW LMC SFH/LMC SFR.txt", skiprows=1,  delimiter=',')

    LMC_lookback_t_vals = LMC_data[:,0]
    LMC_SFR_vals = LMC_data[:,1]
    LMC_SFR_vals = np.multiply(LMC_SFR_vals, 2.660199815753284)
    LMC_t_vals =np.multiply(np.subtract(LMC_lookback_t_vals, np.amax(LMC_lookback_t_vals)),-1)

    LMC_SFR_func = sp.interpolate.interp1d(LMC_t_vals, LMC_SFR_vals, bounds_error=False, fill_value=(np.min(LMC_SFR_vals),np.max(LMC_SFR_vals))) 

    M_sun_bol = 4.74
    L_sun_bol = 3.826e26

    def luminosity(M1, M2, L2):
        power = - ( M1 - M2 ) / 2.5
        return L2 * np.power(10,power)    
    
    # Universe time
    t_universe = 13.466983947061877   

    # Buzzoni SSP Maginitude Model
    # Fe_H = .0
    # alpha1_U = 2.743
    # alpha2_U = 0.26
    # beta_U = 0.006
    # gamma_U = 2.872
    # delta_U = 0.09 * 10 ** Fe_H

    # alpha1_V = 2.163
    # alpha2_V = 0.08
    # beta_V = -0.094
    # gamma_V = 2.246

    # alpha1_B = 2.390
    # alpha2_B = 0.08
    # beta_B = 0.033
    # gamma_B = 2.909

    # offset = 2.5 * np.log10(3.23 * (1.5**2) + 6.41)
    # print(offset)

    # def mag_U(t):
    #     return (alpha1_U + alpha2_U * Fe_H ) * np.log10(t) + beta_U * Fe_H + gamma_U + delta_U + offset
    
    # def mag_V(t):
    #     return (alpha1_V + alpha2_V * Fe_H ) * np.log10(t) + beta_V * Fe_H + gamma_V + offset
      
    # def mag_B(t):
    #     return (alpha1_B + alpha2_B * Fe_H ) * np.log10(t) + beta_B * Fe_H + gamma_B 

    # Function to perform the integrations
    def mass_integral(t_max):
        h1 = 0.0001
        t_vals = np.arange(h1, t_max + h1, h1)
        
        MW_total = 0
        LMC_total = 0

        i = 0
        while ( i < len(t_vals)):
            t = t_vals[i]

            if (i == 0 or i == len(t_vals)-1):  
                MW_total += 0.5 * MW_SFR_func(t_max - t)
                LMC_total += 0.5 * LMC_SFR_func(t_max - t)
            else:
                MW_total += MW_SFR_func(t_max - t)
                LMC_total += LMC_SFR_func(t_max - t)  

            i += 1
        
        MW_total = MW_total * h1 * 1e9
        LMC_total = LMC_total * h1 * 1e9

        return MW_total, LMC_total

    # Function to perform the integrations
    def Luminosity_integral(t_max):
        # time is in Gyr
        h1 = 0.00025

        t_vals = np.arange(h1, t_max + h1, h1)

        MW_bruzal_total_u = 0
        MW_bruzal_total_v = 0

        LMC_bruzal_total_u = 0
        LMC_bruzal_total_v = 0

        i = 0

        while (i < len(t_vals)):
            t = t_vals[i]

            if (i == 0 or i == len(t_vals)-1):

                #MW
                MW_bruzal_total_u += 0.5 * ( luminosity(bruzal_Umag_func( t ), M_sun_bol, L_sun_bol)* MW_SFR_func(t_max - t))
                MW_bruzal_total_v += 0.5 * ( luminosity(bruzal_Vmag_func( t ), M_sun_bol, L_sun_bol)* MW_SFR_func(t_max - t))

                #LMC
                LMC_bruzal_total_u += 0.5 * ( luminosity(bruzal_Umag_func( t ), M_sun_bol, L_sun_bol)* LMC_SFR_func(t_max - t))
                LMC_bruzal_total_v += 0.5 * ( luminosity(bruzal_Vmag_func( t ), M_sun_bol, L_sun_bol)* LMC_SFR_func(t_max - t))


            else:
                #MW
                MW_bruzal_total_u += ( luminosity(bruzal_Umag_func( t ), M_sun_bol, L_sun_bol)* MW_SFR_func(t_max - t))
                MW_bruzal_total_v += ( luminosity(bruzal_Vmag_func( t ), M_sun_bol, L_sun_bol)* MW_SFR_func(t_max - t))

                #LMC
                LMC_bruzal_total_u += ( luminosity(bruzal_Umag_func( t ), M_sun_bol, L_sun_bol)* LMC_SFR_func(t_max - t))
                LMC_bruzal_total_v += ( luminosity(bruzal_Vmag_func( t ), M_sun_bol, L_sun_bol)* LMC_SFR_func(t_max - t))
          
            i += 1

        #MW
        MW_bruzal_total_u = MW_bruzal_total_u * h1 
        MW_bruzal_total_v = MW_bruzal_total_v * h1
        
        #LMC
        LMC_bruzal_total_u = LMC_bruzal_total_u * h1 
        LMC_bruzal_total_v = LMC_bruzal_total_v * h1

        #mw_buzz_col = -2.5 * np.log10(np.divide(MW_buzz_total_u,MW_buzz_total_v))
        mw_bruzal_col = -2.5 * np.log10(np.divide(MW_bruzal_total_u,MW_bruzal_total_v))

        #lmc_buzz_col = -2.5 * np.log10(np.divide(LMC_buzz_total_u,LMC_buzz_total_v))
        lmc_bruzal_col = -2.5 * np.log10(np.divide(LMC_bruzal_total_u,LMC_bruzal_total_v))

        return mw_bruzal_col, lmc_bruzal_col
    

    t_vals = np.linspace(0.1, 12.3, 123) #T vals are in Gyr
    t_vals = np.append(t_vals, 12.33)

    #MW
    MW_bruzal_UVcol_vals = []
    MW_mass_vals = []

    #LMC
    LMC_bruzal_UVcol_vals = []
    LMC_mass_vals = []

    # for t in t_vals:
    #     print(f"Lookback time = {t:.1f} Gyr \n")
    #     # Colour Calculation
    #     mw_bruzal_col, lmc_bruzal_col = Luminosity_integral(t)

    #     MW_bruzal_UVcol_vals = np.append(MW_bruzal_UVcol_vals, mw_bruzal_col)
    #     LMC_bruzal_UVcol_vals = np.append(LMC_bruzal_UVcol_vals, lmc_bruzal_col)

    #     #Mass Calculation
    #     MW_mass, LMC_mass = mass_integral(t)
    #     MW_mass_vals = np.append(MW_mass_vals, MW_mass)
    #     LMC_mass_vals = np.append(LMC_mass_vals, LMC_mass)
        


    ## Find the Progenitor Median Bins
    bin_width = 0.5
    z_bins = np.arange(0, 4+bin_width, bin_width)
    print(cosmo.age(0) - cosmo.age(1.35))
    print(cosmo.age(0) - cosmo.age(2.26))
    prog_col_median = []
    prog_col_upper = []
    prog_col_lower = []

    prog_mass_median = []
    prog_mass_upper = []
    prog_mass_lower = []
    
    prog_ssfr_median = []
    prog_ssfr_upper = []
    prog_ssfr_lower = []

    sat_col_median = []
    sat_col_upper = []
    sat_col_lower = []

    sat_mass_median = []
    sat_mass_upper = []
    sat_mass_lower = []

    for z in z_bins:
        # Colour Median
        prog_UVcol_range = prog_UVcol[np.logical_and(prog_zvals >= z, prog_zvals <= z+bin_width)]

        if (np.isnan(np.median(prog_UVcol_range))):
            prog_col_median = np.append(prog_col_median, prog_col_median[-1])
            prog_col_upper = np.append(prog_col_upper, prog_col_upper[-1])
            prog_col_lower = np.append(prog_col_lower, prog_col_lower[-1])
        else:
            prog_col_median = np.append(prog_col_median, np.median(prog_UVcol_range))
            prog_col_upper = np.append(prog_col_upper, np.percentile(prog_UVcol_range, 75))
            prog_col_lower = np.append(prog_col_lower, np.percentile(prog_UVcol_range, 25))

        sat_UVcol_range = sat_UVcol[np.logical_and(sat_zvals >= z, sat_zvals <= z+bin_width)]

        if (np.isnan(np.median(sat_UVcol_range))):
            sat_col_median = np.append(sat_col_median, sat_col_median[-1])
            sat_col_upper = np.append(sat_col_upper, sat_col_upper[-1])
            sat_col_lower = np.append(sat_col_lower, sat_col_lower[-1])
        else:
            sat_col_median = np.append(sat_col_median, np.median(sat_UVcol_range))
            sat_col_upper = np.append(sat_col_upper, np.percentile(sat_UVcol_range, 75))
            sat_col_lower = np.append(sat_col_lower, np.percentile(sat_UVcol_range, 25))

        # Mass median
        prog_mass_range = prog_mvals[np.logical_and(prog_zvals >= z, prog_zvals <= z+bin_width)]

        if (np.isnan(np.median(prog_mass_range))):
            prog_mass_median = np.append(prog_mass_median, prog_mass_median[-1])
            prog_mass_upper = np.append(prog_mass_upper, prog_mass_upper[-1])
            prog_mass_lower = np.append(prog_mass_lower, prog_mass_lower[-1])
        else:
            prog_mass_median = np.append(prog_mass_median, np.median(prog_mass_range))
            prog_mass_upper = np.append(prog_mass_upper, np.percentile(prog_mass_range, 75))
            prog_mass_lower = np.append(prog_mass_lower, np.percentile(prog_mass_range, 25))

        sat_mass_range = sat_mvals[np.logical_and(sat_zvals >= z, sat_zvals <= z+bin_width)]

        if (np.isnan(np.median(sat_mass_range))):
            sat_mass_median = np.append(sat_mass_median, sat_mass_median[-1])
            sat_mass_upper = np.append(sat_mass_upper, sat_mass_upper[-1])
            sat_mass_lower = np.append(sat_mass_lower, sat_mass_lower[-1])
        else:
            sat_mass_median = np.append(sat_mass_median, np.median(sat_mass_range))
            sat_mass_upper = np.append(sat_mass_upper, np.percentile(sat_mass_range, 75))
            sat_mass_lower = np.append(sat_mass_lower, np.percentile(sat_mass_range, 25))

        # sSFH median
        prog_sfh_range = prog_sfr[np.logical_and(prog_zvals >= z, prog_zvals <= z+bin_width)]
        #prog_ssfr_range = np.divide(np.power(10, prog_sfh_range), np.power(10, prog_mass_range))

        if (np.isnan(np.median(prog_sfh_range))):
            prog_ssfr_median = np.append(prog_ssfr_median, prog_ssfr_median[-1])
            prog_ssfr_upper = np.append(prog_ssfr_upper, prog_ssfr_upper[-1])
            prog_ssfr_lower = np.append(prog_ssfr_lower, prog_ssfr_lower[-1])
        else:
            prog_ssfr_median = np.append(prog_ssfr_median, np.median(prog_sfh_range))
            prog_ssfr_upper = np.append(prog_ssfr_upper, np.percentile(prog_sfh_range, 75))
            prog_ssfr_lower = np.append(prog_ssfr_lower, np.percentile(prog_sfh_range, 25))


    ##Generate the redshift interpolation
    z_range = np.arange(0,15+0.01,0.01)
    t_range = cosmo.age(z_range).value
    t_to_redshift = sp.interpolate.interp1d(t_range, z_range, bounds_error=False)
    print(t_to_redshift(13.47 - 1.5))
    # Find redshift for age values
    t_vals = np.add(np.subtract(t_universe, 12.33), t_vals)
    MW_z_vals = t_to_redshift(t_vals)

    ##Import data
    data = np.loadtxt('data/MW LMC SFH/output.txt', delimiter=",", skiprows=1)

    t_vals = data[:,0]
    MW_mass_vals = data[:,1]
    LMC_mass_vals = data[:,2]
    MW_bruzal_UVcol_vals = data[:,3]
    LMC_bruzal_UVcol_vals = data[:,4]

    '''
        MW LMC Mass vs redshift
    '''

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(11, 6), sharex=True, gridspec_kw={'height_ratios': [3, 1]})
    plt.subplots_adjust(left=0.1, right=0.75, top=0.9, bottom=0.1)

    plt.rcParams["font.family"] = "serif"
    font = {'fontname':'serif'} 
    ax1.grid(color='silver', linestyle='--', linewidth=1, zorder=1)
    ax2.grid(color='silver', linestyle='--', linewidth=1, zorder=1)

   
    ## Calculated stellar mass
    def mmf(x,a,b,c,f):
        return a*np.log10(b*x+c)+f
    ##Satellites
    ax1.plot(z_bins, sat_mass_median, label="Satellite Median", color="mediumorchid", zorder=2, alpha=0.8)
    ax1.plot(z_bins, sat_mass_upper, label="Satellite Bound", color="mediumorchid", zorder=2, ls="--", alpha=0.6)
    ax1.plot(z_bins, sat_mass_lower, color="mediumorchid", zorder=2, ls="--", alpha=0.6)
    ax1.fill_between(x=z_bins, y1=sat_mass_lower, y2=sat_mass_upper, color="mediumorchid", alpha=0.1, zorder=2)
    ##Progenitors
    ax1.plot(z_bins, prog_mass_median, label="Progenitor Median", color="royalblue", zorder=2, alpha=0.8)
    ax1.plot(z_bins, prog_mass_upper, label="Progenitor Bound", color="royalblue", zorder=2, ls="--", alpha=0.6)
    ax1.plot(z_bins, prog_mass_lower, color="royalblue", zorder=2, ls="--", alpha=0.6)
    ax1.fill_between(x=z_bins, y1=prog_mass_lower, y2=prog_mass_upper, color="royalblue", alpha=0.1, zorder=2)
    ##MMF
    ax1.plot(np.linspace(0,4,1000), mmf(np.linspace(0,4,1000),*mmf_popt), label="Minimum Mass Function", color="black", zorder=2, ls="--", alpha=0.8)
    ax1.fill_between(x=np.linspace(0,4,1000), y1=0, y2=mmf(np.linspace(0,4,1000),*mmf_popt), color="black", alpha=0.05, zorder=2)
    #LMC
    ax1.plot(MW_z_vals, np.log10(LMC_mass_vals), label="Large Magellanic Cloud", color="darkviolet", zorder=2)
    #MW
    ax1.plot(MW_z_vals, np.log10(MW_mass_vals), label="Milky Way", color="mediumblue", zorder=2)


    ax1.set_ylabel("Stellar Mass [$Log_{10}$($M_{*}$/$M_{☉}$)]", **font)
    ax1.tick_params(axis='x', length=0)
    ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5), fancybox=True, edgecolor="#474747")
    x_ticklabels = ax1.get_xticklabels()
    y_ticklabels = ax1.get_yticklabels()
    for tick_label in (x_ticklabels + y_ticklabels):
        tick_label.set_fontname("serif")

    ax1.set_xlim(0,4)
    ax1.set_ylim(7,11)


    # LMC & MW SFH
    ## MW
    MW_SFR_t_vals = np.add(MW_SFR_t_vals, np.subtract(t_universe,np.amax(MW_SFR_t_vals)))
    MW_SFH_z = t_to_redshift(MW_SFR_t_vals)
    ax2.plot(MW_SFH_z, MW_SFR_vals, color="mediumblue")
    ax2.set_ylim(0,15)
    ax2.set_ylabel("MW SFR [$M_{☉}$ / yr]", **font)
    ax2.set_xlim(0,4)
    #LMC
    LMC_SFR_vals = LMC_SFR_vals[np.less(LMC_t_vals, 12.33)]
    LMC_t_vals = LMC_t_vals[np.less(LMC_t_vals, 12.33)]
    ax2_LMC = ax2.twinx();
    LMC_t_vals = np.add(LMC_t_vals, np.subtract(t_universe,np.amax(LMC_t_vals)))
    LMC_SFH_z = t_to_redshift(LMC_t_vals)
    ax2_LMC.plot(LMC_SFH_z, LMC_SFR_vals, label="LMC", color="darkviolet")
    ax2_LMC.set_ylabel("LMC SFR [$M_{☉}$ / yr]", **font)
    ax2_LMC.set_ylim(0,0.75)
    ax2_LMC.set_yticks(ticks=[0.25,0.5,0.75],labels=[0.25,0.5,0.75])
   
    ax2.set_xlabel("Redshift",**font)
    plt.subplots_adjust(hspace=0.1)

    x_ticklabels = ax2.get_xticklabels()
    y_ticklabels = ax2.get_yticklabels()
    for tick_label in (x_ticklabels + y_ticklabels):
        tick_label.set_fontname("serif")

    #plt.savefig('plots/LMC_MW_mass.png', dpi=300, bbox_inches='tight')
    plt.show()

    # # Write Data to File    
    # t_vals = np.linspace(0.1, 12.3, 123) #T vals are in Gyr
    # t_vals = np.append(t_vals, 12.33)

    # data = np.c_[t_vals, MW_mass_vals, LMC_mass_vals, MW_bruzal_UVcol_vals, LMC_bruzal_UVcol_vals]
    # np.savetxt('data/MW LMC SFH/output.txt', data, delimiter=',', header="Age (Gyr), MW Stellar Mass, LMC Stellar Mass, MW col (Bruzal), LMC col (Bruzal)")

    '''
        Colour Evolution
    '''
    # fig, axs = plt.subplots(2, 2, figsize=(14, 7), gridspec_kw={'height_ratios': [3, 1], 'hspace':0.15, 'wspace':0.25}, sharex = "col")   
    # plt.rcParams["font.family"] = "serif"
    # font = {'fontname':'serif'} 

    # MW_col_ax = axs[0,0]
    # LMC_col_ax = axs[0,1]

    # MW_SFH_ax = axs[1,0]
    # LMC_SFH_ax = axs[1,1]

    # MW_col_ax.grid(color='silver', linestyle='--', linewidth=1)
    # LMC_col_ax.grid(color='silver', linestyle='--', linewidth=1)
    # MW_SFH_ax.grid(color='silver', linestyle='--', linewidth=1)
    # LMC_SFH_ax.grid(color='silver', linestyle='--', linewidth=1)

    # ##Milky Way Calculated colour
    # MW_col_ax.plot(z_bins, prog_col_median, label="Progenitor Median", color="royalblue", zorder=2)
    # MW_col_ax.plot(z_bins, prog_col_upper, color="royalblue", ls="--", zorder=2, alpha=0.8)
    # MW_col_ax.plot(z_bins, prog_col_lower, label="Progenitor Bound", color="royalblue", ls="--", zorder=2, alpha=0.8)
    # MW_col_ax.fill_between(x=z_bins, y1=prog_col_lower, y2=prog_col_upper, color="royalblue", alpha=0.1, zorder=2)
    # MW_col_ax.plot(MW_z_vals,MW_bruzal_UVcol_vals, label="Milky Way", color="mediumblue")
    # MW_col_ax.set_ylim(-0.5,1.75)
    # MW_col_ax.set_xlim(0,4)
    # MW_col_ax.legend()
    # MW_col_ax.set_ylabel("(U-V) Colour", **font)
    # MW_col_ax.tick_params(axis='x', length=0)

    # ## SSP Colours    
    # # mag_U_vals = mag_U(t_vals)
    # # mag_V_vals = mag_V(t_vals)
    # # mag_UV_vals = np.subtract(mag_U_vals, mag_V_vals)

    # # bc_Umag = bruzal_Umag_func(t_vals)
    # # bc_Vmag = bruzal_Vmag_func(t_vals)
    # # bc_col = np.subtract(bc_Umag, bc_Vmag) 

    # ## MW SFH
    # MW_SFR_t_vals = np.add(MW_SFR_t_vals, np.subtract(t_universe,np.amax(MW_SFR_t_vals)))
    # MW_SFH_z = t_to_redshift(MW_SFR_t_vals)

    # MW_SFH_ax.plot(MW_SFH_z, MW_SFR_vals, label="Milky Way", color="mediumblue", zorder=2)
    # MW_SFH_ax.set_xlim(0,4)
    # MW_SFH_ax.set_ylim(0,15)
    # # MW_SFH_ax.set_yticks(np.arange(0,15, 5),np.arange(0,15, 5))
    # MW_SFH_ax.set_ylabel("SFR [$M_{☉}$ / yr]", **font, labelpad=20)
    # MW_SFH_ax.set_xlabel("Redshift",**font)
    # MW_SFH_ax.legend()

    # ##MLMC Calculated colour
    # LMC_col_ax.plot(z_bins, sat_col_median, label="Satellite Median", color="mediumorchid", zorder=2)
    # LMC_col_ax.plot(z_bins, sat_col_upper, color="mediumorchid", ls="--", zorder=2, alpha=0.8)
    # LMC_col_ax.plot(z_bins, sat_col_lower, label="Satellite Bound", color="mediumorchid", ls="--", zorder=2, alpha=0.8)
    # LMC_col_ax.fill_between(x=z_bins, y1=sat_col_lower, y2=sat_col_upper, color="mediumorchid", alpha=0.1, zorder=2)
    # LMC_col_ax.plot(MW_z_vals,LMC_bruzal_UVcol_vals, label="Large Magellanic Cloud", color="darkviolet")
    # LMC_col_ax.set_ylim(-0.5,1.75)
    # LMC_col_ax.set_xlim(0,4)
    # LMC_col_ax.legend()
    # LMC_col_ax.set_ylabel("(U-V) Colour", **font)
    # LMC_col_ax.tick_params(axis='x', length=0)

    # ## LMC SFH
    # LMC_SFR_vals = LMC_SFR_vals[np.less(LMC_t_vals, 12.33)]
    # LMC_t_vals = LMC_t_vals[np.less(LMC_t_vals, 12.33)]
    # LMC_t_vals = np.add(LMC_t_vals, np.subtract(t_universe,np.amax(LMC_t_vals)))
    # LMC_SFH_z = t_to_redshift(LMC_t_vals)
    # LMC_SFH_ax.plot(LMC_SFH_z, LMC_SFR_vals, label="Large Magellanic Cloud", color="darkviolet")
    # LMC_SFH_ax.set_xlim(0,4)
    # LMC_SFH_ax.set_ylim(0,0.75)
    # LMC_SFH_ax.set_yticks(ticks=np.arange(0,1,0.25),labels=np.arange(0,1,0.25))

    # LMC_SFH_ax.set_ylabel("SFR [$M_{☉}$ / yr]", **font,labelpad=10)
    # LMC_SFH_ax.set_xlabel("Redshift",**font)
    # LMC_SFH_ax.legend()

    # ## Make tick labels serif

    # for ax in axs.flatten():
    #     xtick_labels = ax.get_xticklabels()
    #     ytick_labels = ax.get_yticklabels()
    #     for tick_label in (xtick_labels + ytick_labels):
    #         tick_label.set_fontname("serif")

    # #plt.subplots_adjust(hspace=0.1)
    # plt.savefig("plots/MC_LMC_colevolution", dpi=300, bbox_inches='tight' )
    # plt.show()
    '''
        MW sSFH vs Prog sSFH
    '''
    # plt.rcParams["font.family"] = "serif"
    # font = {'fontname':'serif'} 
    # fig, axs = plt.subplots(1, 2, figsize=(12, 5),gridspec_kw={'wspace':0.25})

    # ssfh_plot = axs[0]
    # col_mass_plot = axs[1]

    # #SFH
    # ssfh_plot.grid(color='silver', linestyle='--', linewidth=1)
    # MW_ssfh = np.divide(MW_SFR_func(t_vals), MW_mass_vals)
    # ssfh_plot.plot(z_bins, np.power(10,prog_ssfr_median), label="Progenitor Median", color="royalblue", zorder=2, alpha=0.8)
    # ssfh_plot.plot(z_bins, np.power(10,prog_ssfr_upper), label="Progenitor Bound", color="royalblue", zorder=2, ls="--", alpha=0.6)
    # ssfh_plot.plot(z_bins, np.power(10,prog_ssfr_lower), color="royalblue", zorder=2, ls="--", alpha=0.6)
    # ssfh_plot.fill_between(x=z_bins, y1=np.power(10,prog_ssfr_lower), y2=np.power(10,prog_ssfr_upper), color="royalblue", alpha=0.1, zorder=2)
    # ssfh_plot.plot(MW_z_vals, MW_SFR_func(t_vals), label="Milky Way", color="mediumblue")
    # ssfh_plot.set_xlim(0,4)
    # ssfh_plot.set_ylim(0,16)
    # ssfh_plot.set_xlabel("Redshift",**font)
    # ssfh_plot.set_ylabel("SFR [$M_{☉}$ / yr]",**font)
    # ssfh_plot.legend()

    # # Colour - Mass plot
    # prog_UVcol = np.append(prog_UVcol,0)
    # prog_mvals = np.append(prog_mvals,0)
    # prog_zvals = np.append(prog_zvals,0)
    # prog_UVcol = np.append(prog_UVcol,0)
    # prog_mvals = np.append(prog_mvals,0)
    # prog_zvals = np.append(prog_zvals,4)

    # cm = plt.cm.get_cmap('viridis')
    # col_mass_plot.grid(color='silver', linestyle='--', linewidth=1, zorder=1)
    # sc2 = col_mass_plot.scatter(prog_mvals, prog_UVcol, s=12, c=prog_zvals, cmap=cm, zorder=2)
    # col_mass_plot.scatter(10.78, 0.43, s=25, label="Milky Way", marker="^", c=[0], zorder=2)
    # col_mass_plot.set_xlim(8, 11)
    # col_mass_plot.set_ylim(-0.5, 2)
    # col_mass_plot.legend()
    # col_mass_plot.set_xlabel("Stellar Mass [$Log_{10}$($M_{*}$/$M_{☉}$)]", **font)
    # col_mass_plot.set_ylabel("(U-V$)_{rest}$ colour", **font)
    # cbar = fig.colorbar(sc2, ticks=[0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, np.amax(prog_zvals)], label="Redshift")


        # ssfh_plot.grid(color='silver', linestyle='--', linewidth=1)
    MW_ssfh = np.divide(MW_SFR_func(t_vals), MW_mass_vals)
    plt.grid(color='silver', linestyle='--', linewidth=1)
    plt.plot(z_bins, np.power(10,prog_ssfr_median), label="Progenitor Median", color="royalblue", zorder=2, alpha=0.8)
    plt.plot(z_bins, np.power(10,prog_ssfr_upper), label="Progenitor Bound", color="royalblue", zorder=2, ls="--", alpha=0.6)
    plt.plot(z_bins, np.power(10,prog_ssfr_lower), color="royalblue", zorder=2, ls="--", alpha=0.6)
    plt.fill_between(x=z_bins, y1=np.power(10,prog_ssfr_lower), y2=np.power(10,prog_ssfr_upper), color="royalblue", alpha=0.1, zorder=2)
    plt.plot(MW_z_vals, MW_SFR_func(t_vals), label="Milky Way", color="mediumblue")
    plt.xlim(0,4)
    plt.ylim(0,16)
    plt.xlabel("Redshift",**font)
    plt.ylabel("SFR [$M_{☉}$ / yr]",**font)
    plt.legend()

    plt.savefig("plots/MC_Prog_SFH", dpi=300, bbox_inches='tight' )
    plt.show()



    '''
        UV
    '''


    # t_vals = np.linspace(np.amin(bruzal_t_vals), 12.3, 100000) #T vals are in Gyr

    # # mag_U_vals = mag_U(t_vals)
    # # mag_V_vals = mag_V(t_vals)
    # # mag_B_vals = mag_B(t_vals)
    # # mag_UV_vals = np.subtract(mag_U_vals, mag_V_vals)
    # # mag_BV_vals = np.subtract(mag_B_vals, mag_V_vals)

    # bc_Umag = bruzal_Umag_func(t_vals)
    # bc_Vmag = bruzal_Vmag_func(t_vals)
    # bc_col = np.subtract(bc_Umag, bc_Vmag) 

    # # plt.semilogx(t_vals, mag_U_vals, label="Buzzoni $M_{u}$", color="blueviolet")
    # # plt.semilogx(t_vals, mag_V_vals, label="Buzzoni $M_{v}$", color="palegreen")
    # # plt.semilogx(t_vals, mag_UV_vals, label="Buzzoni (U-V)", color="dodgerblue")

    # plt.semilogx(t_vals, bc_Umag, label="SSP $M_{u}$", color="blueviolet")
    # plt.semilogx(t_vals, bc_Vmag, label="SSP $M_{v}$", color="palegreen")
    # plt.semilogx(t_vals, bc_col, label="SSP (U-V)", color="dodgerblue")
    
    # #Set the fontstyle for the plot
    # plt.rcParams["font.family"] = "serif"
    # font = {'fontname':'serif'} 
    # #Plot the grid
    # plt.grid(color='silver', linestyle='--', linewidth=1)

    # plt.ylabel("Magnitude", **font)
    # plt.xlabel("Age (Gyr)", **font)
    # # plt.title(f"(U-V) colour for the SSP with Fe/H = {Fe_H}", **font)
    # x_ticklabels = plt.gca().get_xticklabels()
    # y_ticklabels = plt.gca().get_yticklabels()
    # for tick_label in (x_ticklabels + y_ticklabels):
    #     tick_label.set_fontname("serif")

    # plt.legend()
    # plt.ylim(-2,8)
    # plt.xlim(1e-4, 12.3)
    # plt.yticks(np.arange(-2,10,2), np.arange(-2,10,2))
    # plt.savefig('plots/SSP_colours.png', dpi=300, bbox_inches='tight')
    # plt.show()



