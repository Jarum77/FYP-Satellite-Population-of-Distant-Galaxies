import numpy as np
from mpmath import mp
import NDpredict as ndp
import matplotlib.pyplot as plt
from read_ceers import read_in_ceers
from scipy.interpolate import interp1d
from min_mass_func import min_mass_func

zi = 0.
zf = 4.
M0 = 10.78  # All stellar masses are assumed to be logarithmic

zvals = []
Mzfourge_vals = []
Mzfourge_USD_vals = []         #Upper SD value
Mzfourge_LSD_vals = []         #Lower SD value

N0= ndp.getnum(10.75,0, massfunc="zfourge")      #Calculated from ndp.getnum for M0 = 10.5 and z = 0
print(N0)
zstep = 0.25
z = 0

ceers_data, zgrid, pz_vals = read_in_ceers('1 2 3 6',True, 4)

ceers_zvals = ceers_data[:,1]
ceers_mvals = ceers_data[:,2]
SNR = ceers_data[:,5]

all_zvals = ceers_zvals
all_mvals = ceers_mvals

popt = min_mass_func(ceers_zvals, ceers_mvals, SNR)
    
def MMF(x,a,b,c,f):
    return a*np.log10(b*x+c)+f

FOURGE_MINS=[]
FOURGE_MAXS=[]
mins=[]
maxes=[]

#create curve for interpolation
x=np.linspace(0, np.amax(ceers_zvals)+0.5, num=100, endpoint=False)

for i in x:
    ymax=ndp.evolvingN( N0 , 0 , i ) + ndp.sigmaN( N0 , 0, i )
    ymin=ndp.evolvingN( N0 , 0 , i ) - ndp.sigmaN( N0 , 0, i )
    fourgemax=ndp.getmass(ymax, i, massfunc="zfourge")
    fourgemin=ndp.getmass(ymin, i , massfunc="zfourge")
    
    maxes=np.append(maxes,fourgemax)
    mins=np.append(mins,fourgemin)
    
fMax=interp1d(x,mins)
fMin=interp1d(x,maxes)

#first min determination
lower_bounds=fMin(ceers_zvals)
upper_bounds=fMax(ceers_zvals)

zeros=np.zeros(len(lower_bounds))

ceers_zvals1=ceers_zvals[np.less(zeros, (np.subtract(ceers_mvals,lower_bounds)))]
zeros1=zeros[np.less(zeros, (np.subtract(ceers_mvals,lower_bounds)))]
upper_bounds1=upper_bounds[np.less(zeros, (np.subtract(ceers_mvals,lower_bounds)))]
lower_bounds1=lower_bounds[np.less(zeros, (np.subtract(ceers_mvals,lower_bounds)))]
ceers_mvals1=ceers_mvals[np.less(zeros, (np.subtract(ceers_mvals,lower_bounds)))]

ceers_zvals=ceers_zvals1[np.less(zeros1, (np.subtract(upper_bounds1,ceers_mvals1)))]
ceers_mvals=ceers_mvals1[np.less(zeros1, (np.subtract(upper_bounds1,ceers_mvals1)))]


#ZFourge Limits
while (z <= zf):

    Mzfourge = ndp.newmass(M0, zi, z, massfunc='zfourge')    
    
    NDmax = ndp.evolvingN(N0,zi,z) + ndp.sigmaN(N0,zi,z)
    NDmin = ndp.evolvingN(N0,zi,z) - ndp.sigmaN(N0,zi,z)

    Mzfourge_min = ndp.getmass(NDmax, z, massfunc="zfourge")
    Mzfourge_max = ndp.getmass(NDmin, z, massfunc="zfourge")
    #print(Mzfourge_max, ", ", Mzfourge, ", ", Mzfourge_min)
    zvals.append(z)
    Mzfourge_vals.append(Mzfourge)
    Mzfourge_USD_vals.append(Mzfourge_max)
    Mzfourge_LSD_vals.append(Mzfourge_min)

    z += zstep

#Graph Plotting

print(len(ceers_mvals))
plt.xlim([0,4])
plt.ylim([4,12])
plt.rcParams["font.family"] = "serif"
font = {'fontname':'serif'} 
plt.xlabel("Redshift", **font)
plt.ylabel("Stellar Mass [$Log_{10}$($M_{*}$/$M_{☉}$)]", **font)
#plt.title("Stellar Mass vs Redshift for Milky Way Progenitors", fontweight="bold")
##Plot Min Mass func
z_steps = np.arange(0.0001, 4.0001, 0.01)

#Plot the grid
plt.grid(color='silver', linestyle='--', linewidth=1, zorder=1)
# Set font family for tick labels
x_ticklabels = plt.gca().get_xticklabels()
y_ticklabels = plt.gca().get_yticklabels()

for tick_label in (x_ticklabels + y_ticklabels):
    tick_label.set_fontname("serif")
plt.scatter(all_zvals, all_mvals, marker=".", s=1.5, color="silver", label="CEERS Data")
plt.plot(zvals, Mzfourge_USD_vals, "#DE302D", linewidth=1.5, ls="--", zorder=2, alpha=0.8)
plt.plot(zvals, Mzfourge_vals, "#DE302D", label="Median",linewidth=1.5,)
plt.plot(zvals, Mzfourge_LSD_vals, "#DE302D", label="Median ± σ", ls="--",alpha=0.8,linewidth=1.5,)
plt.fill_between(x=zvals, y1=Mzfourge_LSD_vals, y2=Mzfourge_USD_vals, color="#DE302D", alpha=0.2, zorder=2)
plt.plot(z_steps,MMF(z_steps,*popt), label="MMF", zorder=2, color="orchid", linewidth=1.5, ls="--")
plt.fill_between(x=z_steps, y1=0, y2=MMF(z_steps,*popt), color="orchid", alpha=0.2, zorder=2)

plt.legend(markerscale=5,loc='upper center', bbox_to_anchor=(0.5, 1.15), ncol=4)
plt.tight_layout()

plt.savefig('plots/zfourge_mass.png', dpi=300, bbox_inches='tight' )

plt.show()

