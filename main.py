import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.io import readsav
import numpy as np
from mpmath import mp
import NDpredict as ndp
from halo_mass import Halo_Masses
from read_ceers import read_in_ceers
from virial_radius import find_virial_radius
from satellite_finding import satellite_finding
from min_mass_func import min_mass_func
from colourplot import median_colour_plot
from mw_prog_col import MW_Prog_UVcol

ceers_data, zgrid, pz_vals = read_in_ceers('1 2 3 6',True, 4)
ids = ceers_data[:,0]
zvals = ceers_data[:,1]
mass_vals = ceers_data[:,2]
ra_vals = ceers_data[:,3]
dec_vals = ceers_data[:,4]
SNR = ceers_data[:,5]
UVcol = ceers_data[:,6]
A_V = ceers_data[:,7]
A_U = ceers_data[:,8]


z = 0.
#M0 = np.log10(6.43e10)   # All stellar masses are assumed to be logarithmic
#N0= -3.2377700639871483   #Calculated from ndp.getnum for M0 = 10.81 and z = 0
M0 = 10.5
N0 = -2.806814232365416

FOURGE_MINS=[]
FOURGE_MAXS=[]
mins=[]
maxes=[]

#create curve for interpolation
x=np.linspace(0, np.ndarray.max(zvals)+1, num=71, endpoint=False)

for i in x:
    ymax=ndp.evolvingN( N0 , 0 , np.round(i, 5 ) ) + ndp.sigmaN( N0 , 0, np.round(i , 5 ) )
    ymin=ndp.evolvingN( N0 , 0 , np.round(i, 5 ) ) - ndp.sigmaN( N0 , 0, np.round(i, 5 ) )
    fourgemax=ndp.getmass(ymax, round( i , 5 ) , massfunc="zfourge")
    fourgemin=ndp.getmass(ymin, round( i , 5 ) , massfunc="zfourge")
    
    maxes=np.append(maxes,fourgemax)
    mins=np.append(mins,fourgemin)
    
fMax=interp1d(x,mins)
fMin=interp1d(x,maxes)

#first min determination
lower_bounds=fMin(zvals)
upper_bounds=fMax(zvals)

zeros=np.zeros(len(lower_bounds))

all_ids = ids
allra_vals = ra_vals
alldec_vals = dec_vals
allmvals = mass_vals
allzvals = zvals
alluvcol = UVcol
all_AV = A_V
all_AU = A_U

# mmf = min_mass_func(allzvals, allmvals, SNR)

ids1 = ids[np.less(zeros, (np.subtract(mass_vals, lower_bounds)))]
zvals1=zvals[np.less(zeros, (np.subtract(mass_vals,lower_bounds)))]
mass_vals1=mass_vals[np.less(zeros, (np.subtract(mass_vals,lower_bounds)))]
zeros1=zeros[np.less(zeros, (np.subtract(mass_vals,lower_bounds)))]
ra_vals1 = ra_vals[np.less(zeros, (np.subtract(mass_vals,lower_bounds)))]
dec_vals1 = dec_vals[np.less(zeros, (np.subtract(mass_vals,lower_bounds)))]
UVcol1 = UVcol[np.less(zeros, (np.subtract(mass_vals,lower_bounds)))]
AU1 = A_U[np.less(zeros, (np.subtract(mass_vals,lower_bounds)))]
AV1 = A_V[np.less(zeros, (np.subtract(mass_vals,lower_bounds)))]
upper_bounds1=upper_bounds[np.less(zeros, (np.subtract(mass_vals,lower_bounds)))]
lower_bounds1=lower_bounds[np.less(zeros, (np.subtract(mass_vals,lower_bounds)))]

ids = ids1[np.less(zeros1, np.subtract(upper_bounds1, mass_vals1))]
zvals = zvals1[np.less(zeros1, (np.subtract(upper_bounds1,mass_vals1)))]
mass_vals=mass_vals1[np.less(zeros1, (np.subtract(upper_bounds1,mass_vals1)))]
ra_vals = ra_vals1[np.less(zeros1, (np.subtract(upper_bounds1, mass_vals1)))]
dec_vals = dec_vals1[np.less(zeros1, (np.subtract(upper_bounds1, mass_vals1)))]
prog_UVcol=UVcol1[np.less(zeros1, (np.subtract(upper_bounds1, mass_vals1)))]
prog_AU = AU1[np.less(zeros1, (np.subtract(upper_bounds1, mass_vals1)))]
prog_AV = AV1[np.less(zeros1, (np.subtract(upper_bounds1, mass_vals1)))]

# Find the mean data for the progenitors, used in report

# bin_edges = [0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4]

# mean_mass = []
# mass_sdv = []
# mean_col = []
# col_sdv = []

# for z in bin_edges:
#     mass_range  = mass_vals[np.logical_and(zvals >= z, zvals < z+0.5)]
#     mean_mass = np.append(mean_mass, np.mean(mass_range))
#     mass_sdv = np.append(mass_sdv, np.std(mass_range))

#     col_range  = prog_UVcol[np.logical_and(zvals >= z, zvals < z+0.5)]
#     mean_col = np.append(mean_col, np.mean(col_range))
#     col_sdv = np.append(col_sdv, np.std(col_range))

# print(mean_col)
# print(col_sdv)
 
halo_mass_vals = Halo_Masses(zvals,mass_vals)
vradius_deg = find_virial_radius(halo_mass_vals, zvals)

#all_halo_mass = Halo_Masses(allzvals,allmvals)

prog_data = np.c_[ids,  mass_vals, ra_vals, dec_vals, zvals, prog_UVcol, vradius_deg, prog_AU, prog_AV] 
all_data = np.c_[all_ids, allmvals, allra_vals, alldec_vals, allzvals, alluvcol]

found_prog_data, found_sat_data = satellite_finding(prog_data, all_data, pz_vals, zgrid)

#median_colour_plot(found_prog_data, found_sat_data, all_data)

sat_ids = found_sat_data[:,0]
sat_mvals = found_sat_data[:,1]
sat_zvals = found_sat_data[:,2]
sat_UVcol = found_sat_data[:,3]

## Find the index of the satellite ids with all ids for Extinction correction

# Sort Sat Ids in order
sort_indices = np.argsort(sat_ids)  # get indices that would sort array1

# sort all three arrays based on the indices from the previous step
sat_ids = sat_ids[sort_indices]
sat_mvals = sat_mvals[sort_indices]
sat_zvals = sat_zvals[sort_indices]
sat_UVcol = sat_UVcol[sort_indices]

# Get unique satellites
sat_ids, unique_indices = np.unique(sat_ids, return_index=True)

# Use the indices to get the unique values from arrays b and c
sat_mvals = sat_mvals[unique_indices]
sat_zvals = sat_zvals[unique_indices]
sat_UVcol = sat_UVcol[unique_indices]

matches = np.isin(all_ids, sat_ids)  # get a boolean array where True corresponds to a match
sat_index = np.where(matches)[0]

sat_AV = all_AV[sat_index]
sat_AU = all_AU[sat_index]

# Apply extinction values
sat_UVcol_ext = np.add(np.subtract(sat_UVcol, sat_AU), sat_AV)
prog_UVcol_ext = np.add(np.subtract(prog_UVcol, prog_AU), prog_AV)

prog_data = np.c_[zvals, prog_UVcol_ext, mass_vals]
sat_data = np.c_[sat_zvals, sat_UVcol_ext, sat_mvals]

MW_Prog_UVcol(prog_data, found_sat_data)