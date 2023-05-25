import matplotlib.pyplot as plt
from scipy.io import readsav
import numpy as np
from mpmath import mp
import NDpredict as ndp
from scipy.interpolate import interp1d
from astropy.io import fits

def read_in_ceers(ceers, join, zlim):
    
    mass_vals_joined=[]
    rshift_joined=[]
    ids_joined=[]
    ra_joined=[]
    dec_joined=[]
    pzvals_joined=[]
    SNR_joined=[]
    UVcol_joined=[]
    A_V_joined=[]
    A_U_joined=[]
    sfr_joined = []
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
        hdul = fits.open('./data/ceers'+str(i)+'/CEERS_NIRCam'+str(i)+'_v0.2_photom.fits')
        rest_frame_mag = readsav('./data/ceers'+str(i)+'/CEERS'+str(i)+'.rf.save')
        prob_data = fits.open('./data/ceers'+str(i)+'/CEERS_NIRCam'+str(i)+'_v0.2_pz.fits')
  
        data=hdul[1].data
        prob_data=prob_data[1].data
        
   
        ra=data['RA']
        dec=data['DEC']
        mass_vals = mass_data['sed_lmass']
        ids = mass_data['sed_id']
        zvals = rshift_data['z_best'] 
        z_type = rshift_data['z_spec_type']
        A_V=mass_data['sed_Av']
        A_U=np.multiply(A_V,1.47)
        pzvals = prob_data['pz'][0] 
        zgrid = prob_data['zgrid'][0]
        rfU=rest_frame_mag['rfu']
        rfV=rest_frame_mag['rfv'] 
        sfr_vals = mass_data['sed_lsfr']

        #spectroscopic

        F356=data['F356']
        uF356=data['DF356']
        SNR=np.divide(F356, uF356)

        #colour
        UVcol=np.subtract(rfU,rfV)

        #Remove all negative zvals
        ids = ids[~np.less(zvals,0)]
        mass_vals = mass_vals[~np.less(zvals,0)]
        ra = ra[~np.less(zvals,0)]
        dec = dec[~np.less(zvals,0)]
        pzvals = pzvals[~np.less(zvals,0)]  
        z_type = z_type[~np.less(zvals,0)] 
        SNR=SNR[~np.less(zvals,0)]
        UVcol=UVcol[~np.less(zvals,0)]
        A_U = A_U[~np.less(zvals,0)]
        A_V = A_V[~np.less(zvals,0)]
        sfr_vals = sfr_vals[~np.less(zvals,0)] 
        zvals = zvals[~np.less(zvals,0)]    

        #Remove all negative mvals
        ids = ids[~np.less(mass_vals,0)]
        zvals = zvals[~np.less(mass_vals,0)]
        ra = ra[~np.less(mass_vals,0)]
        dec = dec[~np.less(mass_vals,0)]
        pzvals = pzvals[~np.less(mass_vals,0)]  
        z_type = z_type[~np.less(mass_vals,0)] 
        SNR=SNR[~np.less(mass_vals,0)]  
        UVcol=UVcol[~np.less(mass_vals,0)]
        A_U = A_U[~np.less(mass_vals,0)]
        A_V = A_V[~np.less(mass_vals,0)]
        sfr_vals = sfr_vals[~np.less(mass_vals,0)]
        mass_vals = mass_vals[~np.less(mass_vals,0)]     
        

        #if redshift limit exists, apply it...
        if zlim:
            mass_vals = mass_vals[~np.greater(zvals,zlim)]
            ids = ids[~np.greater(zvals,zlim)]
            ra = ra[~np.greater(zvals,zlim)]
            dec = dec[~np.greater(zvals,zlim)]
            pzvals = pzvals[~np.greater(zvals,zlim)] 
            z_type = z_type[~np.greater(zvals,zlim)]
            SNR=SNR[~np.greater(zvals,zlim)]  
            UVcol=UVcol[~np.greater(zvals,zlim)]
            A_U = A_U[~np.greater(zvals,zlim)]
            A_V = A_V[~np.greater(zvals,zlim)]
            sfr_vals = sfr_vals[~np.greater(zvals,zlim)]
            zvals = zvals[~np.greater(zvals,zlim)]

        #Normalise P(z)
        pz_norm = np.linalg.norm(pzvals, axis=1)
        pzvals = np.divide(pzvals,pz_norm[:,None])

        c_id=int(str(i)+'00000')
        ids=ids+c_id
        
        if join == True:
            mass_vals_joined = np.append(mass_vals_joined, mass_vals)
            ids_joined = np.append(ids_joined,ids)
            rshift_joined = np.append(rshift_joined,zvals)
            A_V_joined=np.append(A_V_joined,A_V)
            A_U_joined=np.append(A_U_joined,A_U)
            ra_joined = np.append(ra_joined,ra)
            dec_joined = np.append(dec_joined,dec)
            pzvals_joined = np.append(pzvals_joined, pzvals)
            SNR_joined=np.append(SNR_joined,SNR)
            UVcol_joined=np.append(UVcol_joined,UVcol)
            sfr_joined = np.append(sfr_joined, sfr_vals)
            data = np.c_[ids_joined,rshift_joined,mass_vals_joined,ra_joined,dec_joined,SNR_joined,UVcol_joined,A_V_joined,A_U_joined,sfr_joined]    
        
        else:
            if i == 1:
                data1=np.c_[ids,mass_vals,zvals,ra,dec,pzvals]
            elif i == 2:
                data2=np.c_[ids,mass_vals,zvals,ra,dec,pzvals]
            elif i == 3:
                data3=np.c_[ids,mass_vals,zvals,ra,dec,pzvals]
            elif i == 6:
                data6=np.c_[ids,mass_vals,zvals,ra,dec,pzvals]
        
    if len(data)!=0:
        pzvals_joined[pzvals_joined < 1e-6] = 0
        pzvals_joined = np.reshape(pzvals_joined, (len(ids_joined), 667))

        #np.savetxt('ceers_data.csv', data, delimiter=',', header='ids,redshift,mass,ra,dec,SNR,UV Colour,P(z)')

        return data, zgrid, pzvals_joined
    else:
        return data1, data2, data3, data6
            
 