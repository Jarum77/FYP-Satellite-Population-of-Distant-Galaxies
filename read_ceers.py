import matplotlib.pyplot as plt
from scipy.io import readsav
import numpy as np
from mpmath import mp
import NDpredict as ndp
from scipy.interpolate import interp1d

def read_in_ceers(ceers, join):
    
    mass_vals_joined=[]
    rshift_joined=[]
    ids_joined=[]
    data=[]
    data1=[]
    data2=[]
    data3=[]
    data6=[]
    first=True
    
    try:
        ceers = np.asarray([int(i) for i in ceers.split()])
    except:
        print ("Give required ceers pointing numbers as list, seperated by spaces, followed by 'True' if all data should be concatenated or 'False' if should be returned separately.")
        print ("E.g. '1 3 6, True'")
        return
    
    for i in ceers:
        mass_data = readsav('./data/ceers'+str(i)+'/ceers'+str(i)+'.me_Z02_calz.fastpp.delltaugt8.5.bc03.ch.save')
        rshift_data = readsav('./data/ceers'+str(i)+'/ceers'+str(i)+'_v0.2_redshift.save')
        
        mass_vals = mass_data['sed_lmass']
        ids = mass_data['sed_id']
        zvals = rshift_data['z_best']
                
        #strip any negative z-indices and same indexs for IDs and Masses
        ids=ids[~np.less(zvals,0)]
        mass_vals=mass_vals[~np.less(zvals,0)]
        zvals=zvals[~np.less(zvals,0)]

        c_id=int(str(i)+'00000')
        ids=ids+c_id
        
        if join == True:
            mass_vals_joined=np.append(mass_vals_joined, mass_vals)
            ids_joined=np.append(ids_joined,ids)
            rshift_joined=np.append(rshift_joined,zvals)
            data=np.c_[ids_joined,mass_vals_joined,rshift_joined]    
        
        else:
            if i == 1:
                data1=np.c_[ids,mass_vals,zvals]
            elif i == 2:
                data2=np.c_[ids,mass_vals,zvals]
            elif i == 3:
                data3=np.c_[ids,mass_vals,zvals]
            elif i == 6:
                data6=np.c_[ids,mass_vals,zvals]
        
    if len(data)!=0:
        np.savetxt("data.csv",data,delimiter=",", header="ID, Z-Index, Mass" )
        return data
    else:
        return data1, data2, data3, data6
            
 