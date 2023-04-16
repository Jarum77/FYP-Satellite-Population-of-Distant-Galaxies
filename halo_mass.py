import numpy as np
import matplotlib.pyplot as plt
import scipy as sp

def Halo_Masses(zvals, mvals):

    M0 = 11.590
    Mz = 1.195
    Ep0 = 0.0351
    Epz = -0.0247
    Beta0 = 1.376
    Betaz = -0.826
    gamma0 = 0.608
    gammaz = 0.329

    halo_masses = []

    for i in range(len(mvals)):

        LogM1 = M0 + ( Mz*(zvals[i]/(zvals[i]+1)) )
        Beta = Beta0 + ( Betaz*(zvals[i]/(zvals[i]+1)) )
        EpN = Ep0 + ( Epz*(zvals[i]/(zvals[i]+1)) )
        gamma = gamma0 + ( gammaz*(zvals[i]/(zvals[i]+1)) )

        logMh_vals = np.linspace(10,16,10000)    
        
        LogEp = np.log10( 2 * EpN ) - np.log10(np.power(10, Beta*(LogM1-logMh_vals)) + np.power(10, gamma*(logMh_vals-LogM1)))
        
        LogMstel = logMh_vals + LogEp

        f = sp.interpolate.interp1d(LogMstel, logMh_vals, bounds_error=False, fill_value=(np.min(logMh_vals),np.max(logMh_vals)))
        LogMh_calc = f(mvals[i])
        halo_masses.append(LogMh_calc)

    # fig, ax = plt.subplots()
    # plt.rcParams["font.family"] = "serif"
    # font = {'fontname':'serif'} 
    # ax.set_yticks(np.arange(7.5, 11.5, 0.5))
    # ax.set_xticks(np.arange(10.5,12.5,0.25))
    # plt.grid(color='silver', linestyle='--', linewidth=1, zorder=1)

    # plt.scatter(halo_masses, mvals, c=zvals,s=15,marker=".", zorder=2)
    # plt.xlabel("$Log_{10}$($M_{h}$/$M_{☉}$)",**font)
    # plt.ylabel("$Log_{10}$($M_{*}$/$M_{☉}$)",**font)
    # cbar = plt.colorbar()
    # cbar.ax.set_ylabel("Redshift")
    # cbar.ax.set_yticks([np.amin(zvals),0.5,1,1.5,2,2.5,3,3.5,4])

    # x_ticklabels = plt.gca().get_xticklabels()
    # y_ticklabels = plt.gca().get_yticklabels()
    # for tick_label in (x_ticklabels + y_ticklabels):
    #     tick_label.set_fontname("serif")
    # plt.xlim(10.5,12.25)
    # plt.ylim(7.5,11)
    # plt.savefig('plots/halo_stell_prog.png', dpi=300)
    # plt.show()

    return halo_masses


def sanity_check():

    # 2013
    def thirten_calc(logMh_vals, z):
        M0 = 11.590
        Mz = 1.195
        Ep0 = 0.0351
        Epz = -0.0247
        Beta0 = 1.376
        Betaz = -0.826
        gamma0 = 0.608
        gammaz = 0.329

        LogM1 = (M0 + ( Mz*(z/(z+1))))
        Beta = Beta0 + ( Betaz*(z/(z+1)))
        EpN = Ep0 + ( Epz*(z/(z+1)))
        gamma = gamma0 + ( gammaz*(z/(z+1)))

        LogEp = np.log10( 2 * EpN ) - np.log10(np.power(10, Beta*(LogM1-logMh_vals)) + np.power(10, gamma*(logMh_vals-LogM1)))
        logMstel = LogEp + logMh_vals

        return logMstel
    
    # 2018
    def eight_calc(logMh_vals, z):
        M0 = 11.339
        Mz = 0.692
        Ep0 = 0.005
        Epz = 0.689
        Beta0 = 3.344
        Betaz = -2.079
        gamma0 = 0.966

        LogM1 = (M0 + ( Mz*(z/(z+1))))
        Beta = Beta0 + ( Betaz*(z/(z+1)))
        EpN = Ep0 + ( Epz*(z/(z+1)))
        gamma = gamma0

        LogEp = np.log10( 2 * EpN ) - np.log10(np.power(10, Beta*(LogM1-logMh_vals)) + np.power(10, gamma*(logMh_vals-LogM1)))
        logMstel = LogEp + logMh_vals

        return logMstel

    z=0
    zi = [0,1,2,3,4]
    logMh_vals = np.linspace(10,16,1000)

    Ep_vals = []
    M_max_vals = []
    thir_logMstel_vals = []
    eight_logMstel_vals = []

    for z in zi:

        thir_logMstel_vals.append(thirten_calc(logMh_vals, z)) 
        eight_logMstel_vals.append(eight_calc(logMh_vals, z))         

        # Ep_vals.append(np.power(10,LogEp))
        # logMstel_vals.append(logMstel)
    
    plt.rcParams["font.family"] = "serif"
    font = {'fontname':'serif'} 
    #Plot the grid
    plt.grid(color='silver', linestyle='--', linewidth=1)
    plt.plot(logMh_vals, thir_logMstel_vals[0], color="royalblue", label="Moster et al. 2013", marker=".", markersize="1", alpha=0.8)
    # plt.plot(logMh_vals, thir_logMstel_vals[1], color="midnightblue", label="z = 1", marker=".", markersize="1", alpha=0.8)
    # plt.plot(logMh_vals, thir_logMstel_vals[2], color="goldenrod", label="z = 2", marker=".", markersize="1", alpha=0.8)
    # plt.plot(logMh_vals, thir_logMstel_vals[3], color="lightgreen", label="z = 3", marker=".", markersize="1", alpha=0.8)
    # plt.plot(logMh_vals, thir_logMstel_vals[4], color="mediumvioletred", label="z = 4", marker=".", markersize="1", alpha=0.8)
    plt.plot(logMh_vals, eight_logMstel_vals[0], color="firebrick", label="Moster et al. 2018", marker=".", markersize="1", alpha=0.8)

    plt.xlabel(" $Log_{10}$($M_{DM}$/$M_{☉}$)",**font)
    plt.ylabel(" $Log_{10}$($M_{*}$/$M_{☉}$)",**font)
    x_ticklabels = plt.gca().get_xticklabels()
    y_ticklabels = plt.gca().get_yticklabels()

    for tick_label in (x_ticklabels + y_ticklabels):
        tick_label.set_fontname("serif")
    plt.scatter(12.2,10.81, label="Milky Way", marker="x", s=30, color="black")
    plt.xlim(10,16)
    plt.ylim(3,13)
    plt.legend()
    plt.savefig('plots/stellar_halo_relation.png', dpi=300)
    

    # plt.title("Stellar Mass vs Halo Mass At Different Redshift", fontweight="bold")

    # plt.ylim(7,13)
    # plt.plot(logMh_vals, Ep_vals[0], ls="none", marker=".", markersize="1", color="k", label="Z=0" )
    # plt.plot(logMh_vals, Ep_vals[1], ls="none", marker=".", markersize="1", color="r",label="Z=0.5")
    # plt.plot(logMh_vals, Ep_vals[2], ls="none", marker=".", markersize="1", color="g", label="Z=1")
    # plt.plot(logMh_vals, Ep_vals[3], ls="none", marker=".", markersize="1", color="y", label="Z=2")
    # plt.plot(logMh_vals, Ep_vals[4], ls="none", marker=".", markersize="1", color="b", label="Z=3")
    # plt.plot(logMh_vals, Ep_vals[5], ls="none", marker=".", markersize="1", color="c", label="Z=4")
    # plt.xlim(9.5,15.5)
    # plt.yscale("log")
    # plt.xlabel("Log(Mh/M-Sun)")
    # plt.ylabel("ε")
    # plt.ylim(2*10**-5, 1)
    # plt.xlim(10,16)
    # plt.vlines(M_max_vals[0], 2*10**-5, .3, ls="--", color="k")
    # plt.text(11.6,3*10**-5, "M1 = {:.2f}".format(M_max_vals[0]))
    # plt.legend()
    plt.show()

#sanity_check()
