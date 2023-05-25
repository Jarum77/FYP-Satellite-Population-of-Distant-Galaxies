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
from satellite_colour import satellite_colour_compare
from mw_prog_col import MW_Prog_UVcol
from scipy.stats import ks_2samp

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
sfr_vals = ceers_data[:,9]

# Apply Extinction
UVcol = np.add(np.subtract(UVcol, A_U), A_V)

z = 0.
#M0 = np.log10(6.43e10)   # All stellar masses are assumed to be logarithmic
#N0= -3.2377700639871483   #Calculated from ndp.getnum for M0 = 10.81 and z = 0
M0 = 10.78
N0 = ndp.getnum(10.75,0, massfunc="zfourge")  

FOURGE_MINS=[]
FOURGE_MAXS=[]
mins=[]
maxes=[]

#create curve for interpolation
x=np.linspace(0, np.ndarray.max(zvals)+1, num=100, endpoint=False)

for i in x:
    ymax=ndp.evolvingN( N0 , 0 , np.round(i, 5 ) ) + ndp.sigmaN( N0 , 0, np.round(i , 5 ) )
    ymin=ndp.evolvingN( N0 , 0 , np.round(i, 5 ) ) - ndp.sigmaN( N0 , 0, np.round(i, 5 ) )
    fourgemax=ndp.getmass(ymax, round( i , 5 ) , massfunc="zfourge")
    fourgemin=ndp.getmass(ymin, round( i , 5 ) , massfunc="zfourge")
    
    maxes=np.append(maxes,fourgemax)
    mins=np.append(mins,fourgemin)
    
zfourge_max=interp1d(x,mins)
zfourge_min=interp1d(x,maxes)

#first min determination
lower_bounds=zfourge_min(zvals)
upper_bounds=zfourge_max(zvals)

zeros=np.zeros(len(lower_bounds))

all_ids = ids
allra_vals = ra_vals
alldec_vals = dec_vals
allmvals = mass_vals
allzvals = zvals
alluvcol = UVcol
all_AV = A_V
all_AU = A_U
all_sfr = sfr_vals

mmf_popt = min_mass_func(allzvals, allmvals, SNR)

ids1 = ids[np.less(zeros, (np.subtract(mass_vals, lower_bounds)))]
zvals1=zvals[np.less(zeros, (np.subtract(mass_vals,lower_bounds)))]
mass_vals1=mass_vals[np.less(zeros, (np.subtract(mass_vals,lower_bounds)))]
zeros1=zeros[np.less(zeros, (np.subtract(mass_vals,lower_bounds)))]
ra_vals1 = ra_vals[np.less(zeros, (np.subtract(mass_vals,lower_bounds)))]
dec_vals1 = dec_vals[np.less(zeros, (np.subtract(mass_vals,lower_bounds)))]
UVcol1 = UVcol[np.less(zeros, (np.subtract(mass_vals,lower_bounds)))]
sfr_vals1 = sfr_vals[np.less(zeros, (np.subtract(mass_vals,lower_bounds)))]
upper_bounds1=upper_bounds[np.less(zeros, (np.subtract(mass_vals,lower_bounds)))]
lower_bounds1=lower_bounds[np.less(zeros, (np.subtract(mass_vals,lower_bounds)))]

ids = ids1[np.less(zeros1, np.subtract(upper_bounds1, mass_vals1))]
zvals = zvals1[np.less(zeros1, (np.subtract(upper_bounds1,mass_vals1)))]
mass_vals=mass_vals1[np.less(zeros1, (np.subtract(upper_bounds1,mass_vals1)))]
ra_vals = ra_vals1[np.less(zeros1, (np.subtract(upper_bounds1, mass_vals1)))]
dec_vals = dec_vals1[np.less(zeros1, (np.subtract(upper_bounds1, mass_vals1)))]
prog_UVcol=UVcol1[np.less(zeros1, (np.subtract(upper_bounds1, mass_vals1)))]
sfr_vals = sfr_vals1[np.less(zeros1, (np.subtract(upper_bounds1, mass_vals1)))]
 
halo_mass_vals = Halo_Masses(zvals,mass_vals)
vradius_deg = find_virial_radius(mass_vals, halo_mass_vals, zvals)

#all_halo_mass = Halo_Masses(allzvals,allmvals)


#Remove masses below MMF at z = 4
min_mass = 7.161139665164848

all_ids = all_ids[np.greater(allmvals, min_mass)]
allra_vals = allra_vals[np.greater(allmvals, min_mass)]
alldec_vals = alldec_vals[np.greater(allmvals, min_mass)]
allzvals = allzvals[np.greater(allmvals, min_mass)]
alluvcol = alluvcol[np.greater(allmvals, min_mass)]
all_sfr = all_sfr[np.greater(allmvals, min_mass)]
allmvals = allmvals[np.greater(allmvals, min_mass)]

prog_data = np.c_[ids,  mass_vals, ra_vals, dec_vals, zvals, prog_UVcol, vradius_deg, sfr_vals] 
all_data = np.c_[all_ids, allmvals, allra_vals, alldec_vals, allzvals, alluvcol, all_sfr]

found_prog_data, found_sat_data = satellite_finding(prog_data, all_data, pz_vals, zgrid, zfourge_min, mmf_popt)

sat_ids = found_sat_data[:,0]
sat_mvals = found_sat_data[:,1]
sat_zvals = found_sat_data[:,2]
sat_UVcol = found_sat_data[:,3]

hosts_ids= found_prog_data[:,0]
hosts_mvals = found_prog_data[:,1]
hosts_zvals = found_prog_data[:,2]
hosts_UVcol = found_prog_data[:,3]
# Get unique satellites
# sat_ids, unique_indices = np.unique(sat_ids, return_index=True)

# # Use the indices to get the unique values from arrays b and c
# sat_mvals = sat_mvals[unique_indices]
# sat_zvals = sat_zvals[unique_indices]
# sat_UVcol = sat_UVcol[unique_indices]

# sat_data = np.c_[sat_ids, sat_mvals, sat_zvals, sat_UVcol]


#Remove satellites from ceers galaxies
not_sats = np.in1d(all_ids, sat_ids, invert=True)

field_mvals = allmvals[not_sats]
field_UVcol = alluvcol[not_sats]
field_zvals = allzvals[not_sats]
field_ids = all_ids[not_sats]

##Find the mean data for the progenitors, used in report

bin_edges = [0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5]

prog_mean_mass = []
prog_mass_sdv = []
prog_mean_col = []
prog_col_sdv = []
prog_count = []

sat_mean_mass = []
sat_mass_sdv = []
sat_mean_col = []
sat_col_sdv = []
sat_count = []

hosts_count = []
hosts_avgsats = []
# field_mean_mass = []
# field_mass_sdv = []
# field_mean_col = []
# field_col_sdv = []
# field_count = []

for z in bin_edges:
    #Progenitors
    prog_mass_range  = mass_vals[np.logical_and(zvals >= z, zvals < z+0.5)]
    prog_mean_mass = np.append(prog_mean_mass, np.median(prog_mass_range))
    prog_mass_sdv = np.append(prog_mass_sdv, np.std(prog_mass_range))
    prog_count = np.append(prog_count, len(prog_mass_range))

    prog_col_range  = prog_UVcol[np.logical_and(zvals >= z, zvals < z+0.5)]
    prog_mean_col = np.append(prog_mean_col, np.median(prog_col_range))
    prog_col_sdv = np.append(prog_col_sdv, np.std(prog_col_range))

    #Satellites
    sat_mass_range  = sat_mvals[np.logical_and(sat_zvals >= z, sat_zvals < z+0.5)]
    sat_mean_mass = np.append(sat_mean_mass, np.median(sat_mass_range))
    sat_mass_sdv = np.append(sat_mass_sdv, np.std(sat_mass_range))
    sat_count = np.append(sat_count, len(sat_mass_range))

    sat_col_range  = sat_UVcol[np.logical_and(sat_zvals >= z, sat_zvals < z+0.5)]
    sat_mean_col = np.append(sat_mean_col, np.median(sat_col_range))
    sat_col_sdv = np.append(sat_col_sdv, np.std(sat_col_range))

    #Hosts
    host_mass_range  = hosts_mvals[np.logical_and(hosts_zvals >= z, hosts_zvals < z+0.5)]
    hosts_count = np.append(hosts_count, len(host_mass_range))
    hosts_avgsats = np.append(hosts_avgsats, np.divide(np.sum(sat_count), np.sum(hosts_count)))

    #Field
    # field_mass_range  = field_mvals[np.logical_and(field_zvals >= z, field_zvals < z+0.5)]
    # field_mean_mass = np.append(field_mean_mass, np.median(field_mass_range))
    # field_mass_sdv = np.append(field_mass_sdv, np.std(field_mass_range))
    # field_count = np.append(field_count, len(field_mass_range))

    # field_col_range  = field_UVcol[np.logical_and(field_zvals >= z, field_zvals < z+0.5)]
    # field_mean_col = np.append(field_mean_col, np.median(field_col_range))
    # field_col_sdv = np.append(field_col_sdv, np.std(field_col_range))

# print("\n Progenitor Mass and std")
# print(prog_mean_mass)
# print(prog_mass_sdv)
# print(prog_count)
# print("Progenitor col and std")
# print(prog_mean_col)
# print(prog_col_sdv)

# print("\n Satellite Mass and std")
# print(sat_mean_mass)
# print(sat_mass_sdv)
# print(sat_count)
# print("Satellite col and std")
# print(sat_mean_col)
# print(sat_col_sdv)

# print("\n Hosts data")
# print(hosts_count)
# print(hosts_avgsats)

field_data = np.c_[field_mvals, field_zvals, field_UVcol]
sat_data = np.c_[sat_mvals, sat_zvals, sat_UVcol, sat_ids]

#satellite_colour_compare(sat_data, field_data)
#MW_Prog_UVcol(prog_data, sat_data, mmf_popt)