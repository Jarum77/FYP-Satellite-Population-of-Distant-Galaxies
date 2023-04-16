import matplotlib.pyplot as plt
import numpy as np
from astropy.cosmology import FlatLambdaCDM
import scipy as sp
 
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

def MW_Prog_UVcol(prog_data, found_sat_data):

    #Assign Data
    prog_zvals = prog_data[:,0]
    prog_UVcol = prog_data[:,1]
    prog_mvals = prog_data[:,2]

    sat_zvals = found_sat_data[:,2]
    sat_UVcol = found_sat_data[:,3]
    sat_mvals = found_sat_data[:,1]

    #Bruzal Charlot SSP
    bc_ssp_data = np.loadtxt("data/SSP/ssp.txt", skiprows=1,  delimiter=',') 

    bruzal_log_age = bc_ssp_data[:,0]
    bruzal_Umag = bc_ssp_data[:,2]
    bruzal_Vmag = bc_ssp_data[:,4]   
    bruzal_Bmag = bc_ssp_data[:,3]

    bruzal_col = np.subtract(bruzal_Umag, bruzal_Vmag)
    bruzal_bvcol = np.subtract(bruzal_Bmag, bruzal_Vmag)
    bruzal_t_vals = np.divide(np.power(10, bruzal_log_age),1e9) 

    bruzal_Umag_func = sp.interpolate.interp1d(bruzal_t_vals, np.add(bruzal_Umag, 0.), bounds_error=False, fill_value=(np.min(bruzal_Umag),np.max(bruzal_Umag)))  
    bruzal_Vmag_func = sp.interpolate.interp1d(bruzal_t_vals, np.add(bruzal_Vmag, 0.), bounds_error=False, fill_value=(np.min(bruzal_Vmag),np.max(bruzal_Vmag)))   
    bruzal_Bmag_func = sp.interpolate.interp1d(bruzal_t_vals, np.add(bruzal_Bmag, -0), bounds_error=False, fill_value=(np.min(bruzal_Bmag),np.max(bruzal_Bmag)))   


    # Milky Way Star Formation History

    MW_SFR_data = np.loadtxt("data/MW LMC SFH/MW SFR.txt", skiprows=1,  delimiter=',')

    MW_SFR_lookback_t_vals = MW_SFR_data[:,0]
    MW_SFR_vals = MW_SFR_data[:,1]

    MW_SFR_t_vals =np.multiply(np.subtract(MW_SFR_lookback_t_vals, np.amax(MW_SFR_lookback_t_vals)),-1)

    MW_SFR_func = sp.interpolate.interp1d(MW_SFR_t_vals, MW_SFR_vals, bounds_error=False, fill_value=(np.min(MW_SFR_vals),np.max(MW_SFR_vals))) 

    # LMC Star Formation History

    LMC_data = np.loadtxt("data/MW LMC SFH/LMC SFR.txt", skiprows=1,  delimiter=',')

    LMC_t_vals = LMC_data[:,0]
    LMC_SFR_vals = LMC_data[:,1]

    LMC_SFR_func = sp.interpolate.interp1d(LMC_t_vals, LMC_SFR_vals, bounds_error=False, fill_value=(np.min(LMC_SFR_vals),np.max(LMC_SFR_vals))) 

    M_sun_bol = 4.74
    L_sun_bol = 3.826e26

    def luminosity(M1, M2, L2):
        power = - ( M1 - M2 ) / 2.5
        return L2 * np.power(10,power)    
    
    # Universe time
    t_universe = 13.466983947061877   

    # SSP Maginitude Model

    Fe_H = .0
    alpha1_U = 2.743
    alpha2_U = 0.26
    beta_U = 0.006
    gamma_U = 2.872
    delta_U = 0.09 * 10 ** Fe_H

    alpha1_V = 2.163
    alpha2_V = 0.08
    beta_V = -0.094
    gamma_V = 2.246

    alpha1_B = 2.390
    alpha2_B = 0.08
    beta_B = 0.033
    gamma_B = 2.909

    def mag_U(t):
        return (alpha1_U + alpha2_U * Fe_H ) * np.log10(t) + beta_U * Fe_H + gamma_U + delta_U + 0.79
    
    def mag_V(t):
        return (alpha1_V + alpha2_V * Fe_H ) * np.log10(t) + beta_V * Fe_H + gamma_V + 0.02
      
    def mag_B(t):
        return (alpha1_B + alpha2_B * Fe_H ) * np.log10(t) + beta_B * Fe_H + gamma_B  - 0.09

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

        ## Buzzoni
        buzz_total_u = 0
        buzz_total_v = 0
        buzz_total_b = 0

        ## Bruzal
        bruzal_total_u = 0
        bruzal_total_v = 0
        bruzal_total_b = 0

        i = 0
        while (i < len(t_vals)):
            t = t_vals[i]

            if (i == 0 or i == len(t_vals)-1):

                buzz_total_u += 0.5 * ( luminosity(mag_U( t ), M_sun_bol, L_sun_bol)* MW_SFR_func(t_max - t))
                buzz_total_v += 0.5 * ( luminosity(mag_V( t ), M_sun_bol, L_sun_bol)* MW_SFR_func(t_max - t))
                ##buzz_total_b += 0.5 * ( luminosity(mag_B( t ), M_sun_bol, L_sun_bol)* MW_SFR_func(t_max - t))

                bruzal_total_u += 0.5 * ( luminosity(bruzal_Umag_func( t ), M_sun_bol, L_sun_bol)* MW_SFR_func(t_max - t))
                bruzal_total_v += 0.5 * ( luminosity(bruzal_Vmag_func( t ), M_sun_bol, L_sun_bol)* MW_SFR_func(t_max - t))
                ##bruzal_total_b += 0.5 * ( luminosity(bruzal_Bmag_func( t ), M_sun_bol, L_sun_bol)* MW_SFR_func(t_max - t))

            else:
                buzz_total_u += ( luminosity(mag_U( t ), M_sun_bol, L_sun_bol)* MW_SFR_func(t_max - t))
                buzz_total_v += ( luminosity(mag_V( t ), M_sun_bol, L_sun_bol)* MW_SFR_func(t_max - t))
                ##buzz_total_b += ( luminosity(mag_B( t ), M_sun_bol, L_sun_bol)* MW_SFR_func(t_max - t))

                bruzal_total_u += ( luminosity(bruzal_Umag_func( t ), M_sun_bol, L_sun_bol)* MW_SFR_func(t_max - t))
                bruzal_total_v += ( luminosity(bruzal_Vmag_func( t ), M_sun_bol, L_sun_bol)* MW_SFR_func(t_max - t))
                ##bruzal_total_b += ( luminosity(bruzal_Bmag_func( t ), M_sun_bol, L_sun_bol)* MW_SFR_func(t_max - t))
          
            i += 1

        buzz_total_u = buzz_total_u * h1 
        buzz_total_v = buzz_total_v * h1
        buzz_total_b = buzz_total_b * h1

        bruzal_total_u = bruzal_total_u * h1 
        bruzal_total_v = bruzal_total_v * h1
        bruzal_total_b = bruzal_total_b * h1

        buzz_mw_col = -2.5 * np.log10(np.divide(buzz_total_u,buzz_total_v))
        bruzal_mw_col = -2.5 * np.log10(np.divide(bruzal_total_u,bruzal_total_v))

        return buzz_mw_col, bruzal_mw_col
    

    t_vals = np.linspace(0.1, 12.3, 123) #T vals are in Gyr
    t_vals = np.append(t_vals, 12.33)
    buzz_MW_UV_col_vals = []
    buzz_MW_BV_col_vals = []
    bruzal_MW_UV_col_vals = []
    bruzal_MW_BV_col_vals = []

    MW_mass_vals = []
    LMC_mass_vals = []

    for t in t_vals:
        print(f"Lookback time = {t:.1f} Gyr \n")
        # Colour Calculation
        buzz_mw_col, bruzal_mw_col = Luminosity_integral(t)

        buzz_MW_UV_col_vals = np.append(buzz_MW_UV_col_vals, buzz_mw_col)
        bruzal_MW_UV_col_vals = np.append(bruzal_MW_UV_col_vals, bruzal_mw_col)

        ## Mass Calculation
        # MW_mass, LMC_mass = mass_integral(t)
        # MW_mass_vals = np.append(MW_mass_vals, MW_mass)
        # LMC_mass_vals = np.append(LMC_mass_vals, LMC_mass)


    ## Find the Progenitor Median Bins
    bin_width = 0.5
    z_bins = np.arange(0, 4+bin_width, bin_width)
    print(z_bins)

    prog_col_median = []
    prog_col_upper = []
    prog_col_lower = []
    prog_mass_median = []
    prog_mass_upper = []
    prog_mass_lower = []

    sat_mass_median = []
    sat_mass_upper = []
    sat_mass_lower = []

    for z in z_bins:
        ## Colour Median
        prog_UVcol_range = prog_UVcol[np.logical_and(prog_zvals >= z, prog_zvals <= z+bin_width)]

        if (np.isnan(np.median(prog_UVcol_range))):
            prog_col_median = np.append(prog_col_median, prog_col_median[-1])
            prog_col_upper = np.append(prog_col_upper, prog_col_upper[-1])
            prog_col_lower = np.append(prog_col_lower, prog_col_lower[-1])
        else:
            prog_col_median = np.append(prog_col_median, np.median(prog_UVcol_range))
            prog_col_upper = np.append(prog_col_upper, np.percentile(prog_UVcol_range, 75))
            prog_col_lower = np.append(prog_col_lower, np.percentile(prog_UVcol_range, 25))


        ## Mass median
        # prog_mass_range = prog_mvals[np.logical_and(prog_zvals >= z, prog_zvals <= z+bin_width)]
        # prog_mass_median = np.append(prog_mass_median, np.median(prog_mass_range))
        # prog_mass_upper = np.append(prog_mass_upper, np.percentile(prog_mass_range, 75))
        # prog_mass_lower = np.append(prog_mass_lower, np.percentile(prog_mass_range, 25))

        # sat_mass_range = sat_mvals[np.logical_and(sat_zvals >= z, sat_zvals <= z+bin_width)]
        # sat_mass_median = np.append(sat_mass_median, np.median(sat_mass_range))
        # sat_mass_upper = np.append(sat_mass_upper, np.percentile(sat_mass_range, 75))
        # sat_mass_lower = np.append(sat_mass_lower np.percentile(prog_mass_range, 25))

    ##Generate the redshift interpolation
    z_range = np.arange(0,15+0.01,0.01)
    t_range = cosmo.age(z_range).value
    t_to_redshift = sp.interpolate.interp1d(t_range, z_range, bounds_error=False)

    '''
        MW Colour vs redshift
    '''

    t_vals = np.add(np.subtract(t_universe, 12.3), t_vals)
    MW_z_vals = t_to_redshift(t_vals)

    '''
        MW LMC Mass vs redshift
    '''

    plt.rcParams["font.family"] = "serif"
    font = {'fontname':'serif'} 
    plt.grid(color='silver', linestyle='--', linewidth=1, zorder=1)

    ##Milky Way
    # plt.plot(np.add(z_bins,bin_width*0.5), prog_mass_median, label="Progenitor Median", color="royalblue", zorder=2)
    # plt.fill_between(x=np.add(z_bins,bin_width*0.5), y1=prog_mass_lower, y2=prog_mass_upper, color="royalblue", alpha=0.6, label="Progenitor Median ± σ", zorder=2)
    # plt.plot(MW_z_vals, np.log10(MW_mass_vals), label="Milky Way Calculation", color="darkslategrey", zorder=2)

    ##Satellites
    # plt.plot(np.add(z_bins,bin_width*0.5), sat_mass_median, label="Satellite Median", color="firebrick", zorder=2)
    # plt.fill_between(x=np.add(z_bins,bin_width*0.5), y1=sat_mass_lower, y2=sat_mass_upper, color="firebrick", alpha=0.6, label="Satellite Median ± σ", zorder=2)
    # plt.plot(MW_z_vals, np.log10(LMC_mass_vals), label="Large Magellanic Cloud", color="darkslategrey", zorder=2)

    # plt.xlabel("Redshift",**font)
    # plt.ylabel("$Log_{10}$($M_{*}$/$M_{☉}$)", **font)
    # plt.legend(loc="lower right")
    # x_ticklabels = plt.gca().get_xticklabels()
    # y_ticklabels = plt.gca().get_yticklabels()
    # for tick_label in (x_ticklabels + y_ticklabels):
    #     tick_label.set_fontname("serif")
    # plt.xlim(0,4)
    # plt.ylim(5.5,9.5)
    # #plt.savefig('plots/MW_mass_prog.png', dpi=300)
    # plt.savefig('plots/LMC_mass_sat.png', dpi=300)
    # plt.show()


    # mag_U_vals = mag_U(t_vals)
    # mag_V_vals = mag_V(t_vals)
    # mag_UV_vals = np.subtract(mag_U_vals, mag_V_vals)

    '''
        Colour Evolution
    '''
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6, 8), sharex=True, gridspec_kw={'height_ratios': [4, 1]})
    plt.rcParams["font.family"] = "serif"
    font = {'fontname':'serif'} 
    ax1.grid(color='silver', linestyle='--', linewidth=1)
    ax2.grid(color='silver', linestyle='--', linewidth=1)

    ##Milky Way Calculated colour
    ax1.plot(MW_z_vals,buzz_MW_UV_col_vals, label="Milky Way, Buzzoni SSP", color="#284AF1")
    ax1.plot(MW_z_vals,bruzal_MW_UV_col_vals, label="Milky Way, Bruzal SSP", color="#284AF1", ls="--")
    ax1.plot(z_bins, prog_col_median, label="Progenitor Median", color="royalblue", zorder=2)
    ax1.plot(z_bins, prog_col_upper, label="Progenitor Upper", color="#28A8F1", ls="--", zorder=2)
    ax1.plot(z_bins, prog_col_lower, label="Progenitor Lower", color="#28A8F1", ls="--", zorder=2)
    ax1.fill_between(x=z_bins, y1=prog_col_lower, y2=prog_col_upper, color="#28A8F1", alpha=0.4, zorder=2)

    ax1.set_xlim(0,4)

    ax1.legend()
    ax1.set_ylabel("(U-V) Colour", **font)
    # plt.plot(MW_z_vals,mag_UV_vals, label="SSP", color="black", ls="--")
    ## MW SFH
    ax2.set_xlabel("Redshift",**font)
    MW_SFH_z = t_to_redshift(MW_SFR_t_vals)
    ax2.plot(MW_SFH_z, MW_SFR_vals, label="Milky Way SFH", color="mediumblue")
    ax2.set_xlim(0,4)
    ax2.set_ylim(0,9)
    ax2.legend()
    ax2.set_ylabel("SFR [$M_{☉}$/yr]", **font)

    # plt.plot(z_bins, median_col, label="Progenitor", color="darkslategrey")
    # plt.plot(z_bins, upper_col, label="Progenitor + σ", color="grey")
    # plt.plot(z_bins, lower_col, label="Progenitor - σ", color="grey")
    # plt.plot(z_bins, median_non_ext, label="Progenitor (No ext)", color="darkslategrey", ls="--")
    plt.subplots_adjust(hspace=0)
    plt.show()

    '''
        SSP Colour vs redshift
    '''

    # t_vals = np.linspace(np.amin(bruzal_t_vals), 12.3, 100000) #T vals are in Gyr

    # mag_U_vals = mag_U(t_vals)
    # mag_V_vals = mag_V(t_vals)
    # mag_B_vals = mag_B(t_vals)
    # mag_UV_vals = np.subtract(mag_U_vals, mag_V_vals)
    # # mag_BV_vals = np.subtract(mag_B_vals, mag_V_vals)

    # bc_Umag = bruzal_Umag_func(t_vals)
    # bc_Vmag = bruzal_Vmag_func(t_vals)
    # bc_col = np.subtract(bc_Umag, bc_Vmag) 

    '''
        SFH Milky Way & LMC
    '''
    # fig, ax = plt.subplots()
    # plt.rcParams["font.family"] = "serif"
    # font = {'fontname':'serif'} 
    # plt.grid(color='silver', linestyle='--', linewidth=1, zorder=1)
    # scaled_MW_SFR = np.divide(MW_SFR_vals, 10)
    # plt.plot(MW_SFR_t_vals, scaled_MW_SFR, label="Milky Way (scaled by 1/10)", color="mediumblue")
    # plt.plot(LMC_t_vals, LMC_SFR_vals, label="Large Magellanic Cloud", color="rebeccapurple")

    # plt.xlabel("Age [Gyr]", **font)
    # plt.ylabel("SFR [$M_{☉}$/yr]", **font)
    # x_ticklabels = plt.gca().get_xticklabels()
    # y_ticklabels = plt.gca().get_yticklabels()
    # for tick_label in (x_ticklabels + y_ticklabels):
    #     tick_label.set_fontname("serif")
    # plt.legend()
    # plt.xlim(np.amin(MW_SFR_t_vals), np.amax(MW_SFR_t_vals))
    # plt.ylim(0,1)
    # plt.savefig('plots/MW_LMC_SFH.png', dpi=300)
    # plt.show()

    '''
        Luminosity Plot
    '''
 
    # plt.plot(t_vals, MW_lu_vals, label="MW LU",color="blueviolet")
    # plt.plot(t_vals, MW_lv_vals, label="MW LV", color="palegreen")      
    # plt.rcParams["font.family"] = "serif"
    # font = {'fontname':'serif'} 
    # #Plot the grid
    # plt.grid(color='silver', linestyle='--', linewidth=1)
    # plt.legend()
    # plt.xlabel("lookback time", **font)
    # plt.ylabel("Luminosity", **font)
    # plt.xlim(0,12.3)
    # plt.show()


    '''
        BV
    '''
    # plt.semilogx(t_vals, mag_B_vals, label="Buzzoni $M_{B}$", color="blueviolet")
    # plt.semilogx(t_vals, mag_V_vals, label="Buzzoni $M_{v}$", color="palegreen")
    # plt.semilogx(t_vals, mag_BV_vals, label="Buzzoni (B-V)", color="dodgerblue")

    # plt.semilogx(bc_t_vals, bc_Bmag, ls="--", label="Bruzal $M_{B}$", color="blueviolet")
    # plt.semilogx(bc_t_vals, bc_Vmag, ls="--", label="Bruzal $M_{v}$", color="palegreen")
    # plt.semilogx(bc_t_vals, bc_bvcol, ls="--", label="Bruzal (B-V)", color="dodgerblue")

    '''
        UV
    '''

    # plt.semilogx(t_vals, mag_U_vals, label="Buzzoni $M_{u}$", color="blueviolet")
    # plt.semilogx(t_vals, mag_V_vals, label="Buzzoni $M_{v}$", color="palegreen")
    # plt.semilogx(t_vals, mag_UV_vals, label="Buzzoni (U-V)", color="dodgerblue")

    # plt.semilogx(t_vals, bc_Umag, ls="--", label="Bruzal $M_{u}$", color="blueviolet")
    # plt.semilogx(t_vals, bc_Vmag, ls="--", label="Bruzal $M_{v}$", color="palegreen")
    # plt.semilogx(t_vals, bc_col, ls="--", label="Bruzal (U-V)", color="dodgerblue")
    
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

    # plt.legend(ncol=2)
    # plt.ylim(-8,8)
    # plt.xlim(1e-4, 12.3)
    # plt.yticks(np.arange(-8,10,2), np.arange(-8,10,2))
    # plt.savefig('plots/SPP_colours.png', dpi=300)
    # plt.show()



