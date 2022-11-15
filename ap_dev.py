import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.io import readsav
import numpy as np
from mpmath import mp
import NDpredict as ndp
from aperture_func import Halo_Masses
from read_ceers import read_in_ceers

ceers_data = read_in_ceers('1 2 3 6',True)
zvals = ceers_data[:,1]
mvals = ceers_data[:,2]

zi = 0.
z = 0.
M0 = 10.5   # All stellar masses are assumed to be logarithmic
N0= -2.806814232365416  #Calculated from ndp.getnum for M0 = 10.5 and z = 0

FOURGE_MINS=[]
FOURGE_MAXS=[]
mins=[]
maxes=[]

x=np.linspace(0, np.ndarray.max(zvals)+1, num=71, endpoint=False)  #Generate equally spaced Z-index from 0 to the largest in our data + 1

for i in x:
    END_max = ndp.evolvingN( N0 , 0 , i ) + ndp.sigmaN( N0 , 0, i )     #Calculate the Evolving Number Density + 1σ
    END_min=ndp.evolvingN( N0 , 0 , i ) - ndp.sigmaN( N0 , 0, i )     #Calculate the Evolving Number Density - 1σ
    fourgemax=ndp.getmass(END_min, i , massfunc="zfourge")            #Calculate the Median fourge mass + 1σ
    fourgemin=ndp.getmass(END_max, i , massfunc="zfourge")            #Calculate the Median fourge mass - 1σ

    print(f"At redshift={i}, The fourgemax={fourgemax} and the fourgemin={fourgemin}")
    maxes.append(fourgemax)
    mins.append(fourgemin)
    
fMax=interp1d(x,mins)
fMin=interp1d(x,maxes)

#first min determination
lower_bounds=fMin(zvals)
upper_bounds=fMax(zvals)

zeros=np.zeros(len(lower_bounds))


zvals1=zvals[np.less(zeros, (np.subtract(mvals,lower_bounds)))]
mass_vals1=mvals[np.less(zeros, (np.subtract(mvals,lower_bounds)))]
zeros1=zeros[np.less(zeros, (np.subtract(mvals,lower_bounds)))]
upper_bounds1=upper_bounds[np.less(zeros, (np.subtract(mvals,lower_bounds)))]
lower_bounds1=lower_bounds[np.less(zeros, (np.subtract(mvals,lower_bounds)))]

zvals=zvals1[np.less(zeros1, (np.subtract(upper_bounds1,mass_vals1)))]
mvals=mass_vals1[np.less(zeros1, (np.subtract(upper_bounds1,mass_vals1)))]


halo_masses=Halo_Masses(zvals,mvals)
print(halo_masses)
#for i in 




# plt.plot(x,fMin(x))
# plt.plot(x,fMax(x), label='max')
# plt.plot(zvals, mass_vals, ls='None', marker='.')
# plt.xlim(0,3)
# plt.legend()
# plt.show()



