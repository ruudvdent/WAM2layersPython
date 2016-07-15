# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 13:24:45 2016

@author: Ent00002
"""

#%% Import libraries
import numpy as np
from netCDF4 import Dataset
import scipy.io as sio
import calendar
from getconstants import getconstants
from timeit import default_timer as timer
import os

#%%BEGIN OF INPUT (FILL THIS IN)
years = np.arange(2010,2011) #fill in the years
yearpart = np.arange(0,364) # for a full (leap)year fill in np.arange(0,366)
boundary = 8 # with 8 the vertical separation is at 812.83 hPa for surface pressure = 1031.25 hPa, which corresponds to k=47 (ERA-Interim)
divt = 24 # division of the timestep, 24 means a calculation timestep of 6/24 = 0.25 hours (numerical stability purposes)
count_time = 4 # number of indices to get data from (for six hourly data this means everytime one day)

# Manage the extent of your dataset (FILL THIS IN)
# Define the latitude and longitude cell numbers to consider and corresponding lakes that should be considered part of the land
latnrs = np.arange(7,114)
lonnrs = np.arange(0,240)
isglobal = 1 # fill in 1 for global computations (i.e. Earth round), fill in 0 for a local domain with boundaries

# the lake numbers below belong to the ERA-Interim data on 1.5 degree starting at Northern latitude 79.5 and longitude 0
lake_mask_1 = np.array([9,9,9,12,12,21,21,22,22,23,24,25,23,23,25,25,53,54,61,23,24,23,24,25,27,22,23,24,25,26,27,28,22,25,26,27,28,23,23,12,18])
lake_mask_2 = np.array([120+19,120+40,120+41,120+43,120+44,120+61,120+62,120+62,120+63,120+62,120+62,120+62,120+65,120+66,120+65,120+66,142-120,142-120,143-120,152-120,152-120,153-120,153-120,153-120,153-120,154-120,154-120,154-120,154-120,154-120,154-120,154-120,155-120,155-120,155-120,155-120,155-120,159-120,160-120,144-120,120+55])
lake_mask = np.transpose(np.vstack((lake_mask_1,lake_mask_2))) #recreate the arrays of the matlab model

#END OF INPUT

#%% Datapaths (FILL THIS IN)
invariant_data = 'Interim_data/full/invariants.nc' #invariants
interdata_folder = r'C:\Users\bec\Desktop\WAM2\interdata'
input_folder = r'C:\Users\bec\Desktop\WAM2'

# other scripts use exactly this sequence, do not change it unless you change it also in the scripts
def data_path(yearnumber,a):   
    sp_data = os.path.join(input_folder, str(yearnumber) + '-sp.nc') #surface pressure
    sp_eoy_data = os.path.join(input_folder, str(yearnumber+1) + '-sp.nc') #surface pressure end of the year
    
    q_f_data = os.path.join(input_folder, str(yearnumber) + '-q_mod.nc') #specific humidity 
    q_f_eoy_data = os.path.join(input_folder, str(yearnumber+1) + '-q_mod.nc')#specific humidity end of the year

    tcw_data = os.path.join(input_folder, str(yearnumber) + '-tcw.nc') #total column water 
    tcw_eoy_data = os.path.join(input_folder, str(yearnumber+1) + '-tcw.nc') #total column water end of the year

    u_f_data = os.path.join(input_folder, str(yearnumber) + '-u_mod.nc' )
    u_f_eoy_data = os.path.join(input_folder, str(yearnumber+1) + '-u_mod.nc' )

    v_f_data = os.path.join(input_folder, str(yearnumber) + '-v_mod.nc' )
    v_f_eoy_data = os.path.join(input_folder, str(yearnumber+1) + '-v_mod.nc' )

    ewvf_data = os.path.join(input_folder, str(yearnumber) + '-ewvf.nc')
    ewvf_eoy_data = os.path.join(input_folder, str(yearnumber+1) + '-ewvf.nc')

    nwvf_data = os.path.join(input_folder, str(yearnumber) + '-nwvf.nc')
    nwvf_eoy_data = os.path.join(input_folder, str(yearnumber+1) + '-nwvf.nc')

    eclwf_data = os.path.join(input_folder, str(yearnumber) + '-eclwf.nc')
    eclwf_eoy_data = os.path.join(input_folder, str(yearnumber+1) + '-eclwf.nc')

    nclwf_data = os.path.join(input_folder, str(yearnumber) + '-nclwf.nc')
    nclwf_eoy_data = os.path.join(input_folder, str(yearnumber+1) + '-nclwf.nc')

    ecfwf_data = os.path.join(input_folder, str(yearnumber) + '-ecfwf.nc')
    ecfwf_eoy_data = os.path.join(input_folder, str(yearnumber+1) + '-ecfwf.nc')

    ncfwf_data = os.path.join(input_folder, str(yearnumber) + '-ncfwf.nc')
    ncfwf_eoy_data = os.path.join(input_folder, str(yearnumber+1) + '-ncfwf.nc')

    evaporation_precipitation_data = os.path.join(input_folder, str(yearnumber) + '-E-P.nc')
    
    save_path = os.path.join(interdata_folder, str(yearnumber) + '-' + str(a) + 'fluxes_storages.mat')
    
    return sp_data,sp_eoy_data,q_f_data,q_f_eoy_data,tcw_data,tcw_eoy_data,u_f_data,u_f_eoy_data,v_f_data,v_f_eoy_data,ewvf_data,ewvf_eoy_data,nwvf_data,nwvf_eoy_data,eclwf_data,eclwf_eoy_data,nclwf_data,nclwf_eoy_data,ecfwf_data,ecfwf_eoy_data,ncfwf_data,ncfwf_eoy_data,evaporation_precipitation_data,save_path

#%% Code (no need to look at this for running)

# The model level as downloaded from ERA-Interim are hard-defined within this code. So it has to be changed if other model levels are downloaded
def getW(latnrs,lonnrs,final_time,a,yearnumber,begin_time,count_time,
    density_water,latitude,longitude,g,A_gridcell,boundary):
    
    if a != final_time: # not the end of the year
        
        # surface pressure (state at 00.00, 06.00, 12.00, 18.00)
        sp = Dataset(datapath[0], mode = 'r').variables['sp'][begin_time:(begin_time+count_time+1),latnrs,lonnrs] #Pa
        
        # specific humidity (state at 00.00, 06.00, 12.00, 18.00)
        q = Dataset(datapath[2], mode = 'r').variables['q'][begin_time:(begin_time+count_time+1),:,latnrs,lonnrs] #kg/kg
        
        # total column water 
        tcw_ERA = Dataset(datapath[4], mode = 'r').variables['tcw'][begin_time:(begin_time+count_time+1),latnrs,lonnrs] #kg/m2
    
    else: #end of the year
        
        # surface pressure (state at 00.00, 06.00, 12.00, 18.00)
        sp_first = Dataset(datapath[0], mode = 'r').variables['sp'][begin_time:(begin_time+count_time),latnrs,lonnrs]
        sp = np.insert(sp_first,[4],(Dataset(datapath[1], mode = 'r').variables['sp'][0,latnrs,lonnrs]), axis = 0) #Pa
        
        # specific humidity (state at 00.00 06.00 12.00 18.00)
        q_first = Dataset(datapath[2], mode = 'r').variables['q'][begin_time:(begin_time+count_time),:,latnrs,lonnrs]
        q = np.insert(q_first,[4],(Dataset(datapath[3], mode = 'r').variables['q'][0,:,latnrs,lonnrs]), axis = 0) #kg/kg
        
        # total column water 
        tcw_ERA_first = Dataset(datapath[4], mode = 'r').variables['tcw'][begin_time:(begin_time+count_time),latnrs,lonnrs]
        tcw_ERA = np.insert(tcw_ERA_first,[4],Dataset(datapath[5], mode = 'r').variables['tcw'][0,latnrs,lonnrs],axis = 0)
    
    # below are the model levels k as downloaded from the ERA-Interim archive, with the pressure defined by A + B * sp
    A = np.vstack(np.array([0 , 20, 38.42534 , 63.6478 , 95.63696 , 134.4833 , 180.5844 , 234.7791 , 298.4958 , 373.9719 , 464.6181 , 
              575.651 , 713.2181 , 883.6605 , 1094.835 , 1356.475 , 1680.64 , 2082.274 , 2579.889 , 3196.422 , 3960.292 , 
              4906.708 , 6018.02 , 7306.631 , 8765.054 , 10376.13 , 12077.45 , 13775.33 , 15379.81 , 16819.47 , 18045.18 , 
              19027.70 , 19755.11 , 20222.21 , 20429.86 , 20384.48 , 20097.40 , 19584.33 , 18864.75 , 17961.36 , 16899.47 , 
              15706.45 , 14411.12 , 13043.22 , 11632.76 , 10209.50 , 8802.356 , 7438.803 , 6144.315 , 4941.778 , 3850.913 ,
              2887.697 , 2063.78 , 1385.913 , 855.3618 , 467.3336 , 210.3939 , 65.88924 , 7.367743 , 0 , 0]))

    B = np.vstack(np.array([0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 7.58*(10**-5) ,
                  0.000461 , 0.001815 , 0.005081 , 0.011143 , 0.020678 , 0.034121 , 0.05169 , 0.073534 , 0.099675 , 0.130023 , 
                  0.164384 , 0.202476 , 0.243933 , 0.288323 , 0.335155 , 0.383892 , 0.433963 , 0.484772 , 0.53571 , 0.586168 , 
                  0.635547 , 0.683269 , 0.728786 , 0.771597 , 0.811253 , 0.847375 , 0.879657 , 0.907884 , 0.93194 , 0.951822 ,
                  0.967645 , 0.979663 , 0.98827 , 0.994019 , 0.99763 , 1]))

    k = np.vstack(np.array([0 , 17 , 27 , 32 , 35 , 38 , 41 , 44 , 47 , 48 , 51 , 54 , 55 , 56 , 57 , 58 , 59 , 60]))
    
    p_swap = np.zeros((len(k), count_time+1, len(latitude), len(longitude))) #Pa
    for m in range(len(k)):
        p_swap[m] =  A[k[m]] + B[k[m]] * sp[:,:,:]

    # make cwv vector
    q_swap = np.swapaxes(q,0,1)
    cwv_swap = np.zeros((np.shape(q_swap))) #kg/m2
    for n in range(len(k)-1):
        cwv_swap[n] = (np.squeeze(q_swap[n])*np.squeeze(p_swap[n+1] - p_swap[n])) / g # column water vapor = specific humidity * pressure levels length / g [kg/m2]
    cwv = np.swapaxes(cwv_swap,0,1)
    
    # make tcwv vector
    tcwv = np.squeeze(np.sum(cwv,1)) #total column water vapor, cwv is summed over the vertical [kg/m2]
    
    # make cw vector
    cw_swap = np.zeros((np.shape(cwv_swap)))
    for n in range(len(k)-1):
        cw_swap[n] = (tcw_ERA / tcwv) * np.squeeze(cwv_swap[n])
    cw = np.swapaxes(cw_swap,0,1)
        
    # just a test, return this variable when running tests
    # tcw = np.squeeze(np.sum(cw,1)) #total column water, cw is summed over the vertical [kg/m2]    
    # testvar = tcw_ERA - tcw # should be around zero for most cells
    
    # put A_gridcell on a 3D grid
    A_gridcell2D = np.tile(A_gridcell,[1,len(longitude)])
    A_gridcell_1_2D = np.reshape(A_gridcell2D, [1,len(latitude),len(longitude)])
    A_gridcell_plus3D = np.tile(A_gridcell_1_2D,[count_time+1,1,1])
    
    # water volumes
    vapor_top = np.squeeze(np.sum(cwv[:,0:boundary,:,:],1))
    vapor_down = np.squeeze(np.sum(cwv[:,boundary:,:,:],1))
    vapor = vapor_top + vapor_down
    W_top = tcw_ERA * (vapor_top / vapor) * A_gridcell_plus3D / density_water #m3
    W_down = tcw_ERA * (vapor_down / vapor) * A_gridcell_plus3D / density_water #m3
    
    return cw, W_top, W_down

#%% Code
def getwind(latnrs,lonnrs,final_time,a,yearnumber,begin_time,count_time):
    # u stands for wind in zonal direction = west-east
    # v stands for wind in meridional direction = south-north 
    if a != final_time: # not the end of the year
        
        # read the u-wind data
        u_f = Dataset(datapath[6], mode = 'r').variables['u'][begin_time:(begin_time+count_time+1),:,latnrs,lonnrs] #m/s
        
        # read the v-wind data
        v_f = Dataset(datapath[8], mode = 'r').variables['v'][begin_time:(begin_time+count_time+1),:,latnrs,lonnrs] #m/s
        
    else: #end of year
        
        # read the u-wind data
        u_f_first = Dataset(datapath[6], mode = 'r').variables['u'][begin_time:(begin_time+count_time),:,latnrs,lonnrs]
        u_f = np.insert(u_f_first,[4],(Dataset(datapath[7], mode = 'r').variables['u'][0,:,latnrs,lonnrs]), axis = 0)
        
        # read the v-wind data
        v_f_first = Dataset(datapath[8], mode = 'r').variables['v'][begin_time:(begin_time+count_time),:,latnrs,lonnrs]
        v_f = np.insert(v_f_first,[4],(Dataset(datapath[9], mode = 'r').variables['v'][0,:,latnrs,lonnrs]), axis = 0)
    
    U = u_f
    V = v_f
    
    return U,V


#%% Code
def getFa(latnrs,lonnrs,boundary,cw,U,V,count_time,begin_time,yearnumber,a,final_time):
    
    if a != final_time: #not the end of the year
        
        #get ERA vertically integrated fluxes
        ewvf = Dataset(datapath[10], mode = 'r').variables['p71.162'][begin_time:(begin_time+count_time+1),latnrs,lonnrs]
        nwvf = Dataset(datapath[12], mode = 'r').variables['p72.162'][begin_time:(begin_time+count_time+1),latnrs,lonnrs]
        eclwf = Dataset(datapath[14], mode = 'r').variables['p88.162'][begin_time:(begin_time+count_time+1),latnrs,lonnrs]
        nclwf = Dataset(datapath[16], mode = 'r').variables['p89.162'][begin_time:(begin_time+count_time+1),latnrs,lonnrs]
        ecfwf = Dataset(datapath[18], mode = 'r').variables['p90.162'][begin_time:(begin_time+count_time+1),latnrs,lonnrs]
        ncfwf = Dataset(datapath[20], mode = 'r').variables['p91.162'][begin_time:(begin_time+count_time+1),latnrs,lonnrs]
        
    else: #end of year
        
        ewvf_first = Dataset(datapath[10], mode = 'r').variables['p71.162'][begin_time:(begin_time+count_time),latnrs,lonnrs]
        ewvf = np.insert(ewvf_first,[4],(Dataset(datapath[11], mode = 'r').variables['p71.162'][0,latnrs,lonnrs]), axis = 0)
        nwvf_first = Dataset(datapath[12], mode = 'r').variables['p72.162'][begin_time:(begin_time+count_time),latnrs,lonnrs]
        nwvf = np.insert(nwvf_first,[4],(Dataset(datapath[13], mode = 'r').variables['p72.162'][0,latnrs,lonnrs]), axis = 0)
        eclwf_first = Dataset(datapath[14], mode = 'r').variables['p88.162'][begin_time:(begin_time+count_time),latnrs,lonnrs]
        eclwf = np.insert(eclwf_first,[4],(Dataset(datapath[15], mode = 'r').variables['p88.162'][0,latnrs,lonnrs]), axis = 0)
        nclwf_first = Dataset(datapath[16], mode = 'r').variables['p89.162'][begin_time:(begin_time+count_time),latnrs,lonnrs]
        nclwf = np.insert(nclwf_first,[4],(Dataset(datapath[17], mode = 'r').variables['p89.162'][0,latnrs,lonnrs]), axis = 0)
        ecfwf_first = Dataset(datapath[18], mode = 'r').variables['p90.162'][begin_time:(begin_time+count_time),latnrs,lonnrs]
        ecfwf = np.insert(ecfwf_first,[4],(Dataset(datapath[19], mode = 'r').variables['p90.162'][0,latnrs,lonnrs]), axis = 0)
        ncfwf_first = Dataset(datapath[20], mode = 'r').variables['p91.162'][begin_time:(begin_time+count_time),latnrs,lonnrs]
        ncfwf = np.insert(ncfwf_first,[4],(Dataset(datapath[21], mode = 'r').variables['p91.162'][0,latnrs,lonnrs]), axis = 0)
        
    ewf = ewvf + eclwf + ecfwf #kg*m-1*s-1
    nwf = nwvf + nclwf + ncfwf #kg*m-1*s-1
    
    #eastward and northward fluxes
    Fa_E_p = U * cw
    Fa_N_p = V * cw
    
    # uncorrected down and top fluxes
    Fa_E_down_uncorr = np.squeeze(np.sum(Fa_E_p[:,boundary:,:,:],1)) #kg*m-1*s-1
    Fa_N_down_uncorr = np.squeeze(np.sum(Fa_N_p[:,boundary:,:,:],1)) #kg*m-1*s-1
    Fa_E_top_uncorr = np.squeeze(np.sum(Fa_E_p[:,0:boundary,:,:],1)) #kg*m-1*s-1
    Fa_N_top_uncorr = np.squeeze(np.sum(Fa_N_p[:,0:boundary,:,:],1)) #kg*m-1*s-1
    
    # correct down and top fluxes
    corr_E = np.zeros([count_time+1,len(latitude),len(longitude)])   
    corr_N = np.zeros([count_time+1,len(latitude),len(longitude)])
    for i in range(len(longitude)):
        for j in range(len(latitude)):
            for k in range(count_time+1):
                corr_E[k,j,i] = min(2,max(0,ewf[k,j,i]/(Fa_E_down_uncorr[k,j,i] + Fa_E_top_uncorr[k,j,i])))
                corr_N[k,j,i] = min(2,max(0,nwf[k,j,i]/(Fa_N_down_uncorr[k,j,i] + Fa_N_top_uncorr[k,j,i])))
    Fa_E_down = corr_E * Fa_E_down_uncorr #kg*m-1*s-1
    Fa_N_down = corr_N * Fa_N_down_uncorr #kg*m-1*s-1
    Fa_E_top = corr_E * Fa_E_top_uncorr #kg*m-1*s-1
    Fa_N_top = corr_N * Fa_N_top_uncorr #kg*m-1*s-1
    
    # make the fluxes during the timestep
    Fa_E_down = 0.5*(Fa_E_down[0:-1,:,:]+Fa_E_down[1:,:,:]);
    Fa_N_down = 0.5*(Fa_N_down[0:-1,:,:]+Fa_N_down[1:,:,:]);
    Fa_E_top = 0.5*(Fa_E_top[0:-1,:,:]+Fa_E_top[1:,:,:]);
    Fa_N_top = 0.5*(Fa_N_top[0:-1,:,:]+Fa_N_top[1:,:,:]);
    
    
    return Fa_E_top,Fa_N_top,Fa_E_down,Fa_N_down


#%% Code
def getEP(latnrs,lonnrs,yearnumber,begin_time,count_time,latitude,longitude,A_gridcell):
    
    #(accumulated after the forecast at 00.00 and 12.00 by steps of 3 hours in time
    evaporation = Dataset(datapath[22], mode = 'r').variables['e'][begin_time*2:(begin_time*2+count_time*2),latnrs,lonnrs] #m
    precipitation = Dataset(datapath[22], mode = 'r').variables['tp'][begin_time*2:(begin_time*2+count_time*2),latnrs,lonnrs] #m
    
    #values that apply for the amount fallen/evaporated during the previous three hours
    for x in np.arange(0,count_time*2,4):
        evaporation[x+3,:,:] = evaporation[x+3,:,:] - evaporation[x+2,:,:] 
        evaporation[x+2,:,:] = evaporation[x+2,:,:] - evaporation[x+1,:,:] 
        evaporation[x+1,:,:] = evaporation[x+1,:,:] - evaporation[x,:,:] 
        precipitation[x+3,:,:] = precipitation[x+3,:,:] - precipitation[x+2,:,:]
        precipitation[x+2,:,:] = precipitation[x+2,:,:] - precipitation[x+1,:,:]
        precipitation[x+1,:,:] = precipitation[x+1,:,:] - precipitation[x,:,:]
    
    #delete and transfer negative values, change sign convention to all positive
    precipitation = np.reshape(np.maximum(np.reshape(precipitation, (np.size(precipitation))) + np.maximum(np.reshape(evaporation, (np.size(evaporation))),0.0),0.0),
                        (np.int(count_time*2),len(latitude),len(longitude))) 
    evaporation = np.reshape(np.abs(np.minimum(np.reshape(evaporation, (np.size(evaporation))),0.0)),(np.int(count_time*2),len(latitude),len(longitude)))   
    
    #calculate volumes
    A_gridcell2D = np.tile(A_gridcell,[1,len(longitude)])
    A_gridcell_1_2D = np.reshape(A_gridcell2D, [1,len(latitude),len(longitude)])
    A_gridcell_max3D = np.tile(A_gridcell_1_2D,[count_time*2,1,1])

    E = evaporation * A_gridcell_max3D
    P = precipitation * A_gridcell_max3D

    return E, P

#%% Code
def getrefined(Fa_E_top,Fa_N_top,Fa_E_down,Fa_N_down,W_top,W_down,E,P,divt,count_time,latitude,longitude):
    
    #for 3 hourly information
    divt2 = divt/2.
    oddvector2 = np.zeros((1,np.int(count_time*2*divt2)))
    partvector2 = np.zeros((1,np.int(count_time*2*divt2)))
    da = np.arange(1,divt2)  
    
    for o in np.arange(0,np.int(count_time*2*divt2),12):
        for i in range(len(da)):
            oddvector2[0,o+i]    = (divt2-da[i])/divt2
            partvector2[0,o+i+1] = da[i]/divt2

    E_small = np.nan*np.zeros((np.int(count_time*2*divt2),len(latitude),len(longitude)))
    for t in range(1,np.int(count_time*2*divt2)+1):
        E_small[t-1] = (1./divt2) * E[np.int(t/divt2+oddvector2[0,t-1]-1)]
    E = E_small            

    P_small = np.nan*np.zeros((np.int(count_time*2*divt2),len(latitude),len(longitude)))
    for t in range(1,np.int(count_time*2*divt2)+1):
        P_small[t-1] = (1./divt2) * P[np.int(t/divt2+oddvector2[0,t-1]-1)]
    P = P_small
    
    # for 6 hourly info
    oddvector = np.zeros((1,np.int(count_time*divt)))
    partvector = np.zeros((1,np.int(count_time*divt)))
    da = np.arange(1,divt) 
    divt = np.float(divt)
    for o in np.arange(0,np.int(count_time*divt),np.int(divt)):
        for i in range(len(da)):
            oddvector[0,o+i]    = (divt-da[i])/divt
            partvector[0,o+i+1] = da[i]/divt    
        
    W_top_small = np.nan*np.zeros((np.int(count_time*divt+1),len(latitude),len(longitude)))
    for t in range(1,np.int(count_time*divt)+1):
        W_top_small[t-1] = W_top[np.int(t/divt+oddvector[0,t-1]-1)] + partvector[0,t-1] * (W_top[np.int(t/divt+oddvector[0,t-1])] - W_top[np.int(t/divt+oddvector[0,t-1]-1)])
        W_top_small[-1] = W_top[-1]
    W_top = W_top_small

    W_down_small = np.nan*np.zeros((np.int(count_time*divt+1),len(latitude),len(longitude)))
    for t in range(1,np.int(count_time*divt)+1):
        W_down_small[t-1] = W_down[np.int(t/divt+oddvector[0,t-1]-1)] + partvector[0,t-1] * (W_down[np.int(t/divt+oddvector[0,t-1])] - W_down[np.int(t/divt+oddvector[0,t-1]-1)])
        W_down_small[-1] = W_down[-1]
    W_down = W_down_small

    Fa_E_down_small = np.nan*np.zeros((np.int(count_time*divt),len(latitude),len(longitude)))
    Fa_N_down_small = np.nan*np.zeros((np.int(count_time*divt),len(latitude),len(longitude)))
    Fa_E_top_small = np.nan*np.zeros((np.int(count_time*divt),len(latitude),len(longitude)))
    Fa_N_top_small = np.nan*np.zeros((np.int(count_time*divt),len(latitude),len(longitude)))
    for t in range(1,np.int(count_time*divt)+1):
        Fa_E_down_small[t-1] = Fa_E_down[np.int(t/divt+oddvector[0,t-1]-1)]
        Fa_N_down_small[t-1] = Fa_N_down[np.int(t/divt+oddvector[0,t-1]-1)]
        Fa_E_top_small[t-1] = Fa_E_top[np.int(t/divt+oddvector[0,t-1]-1)]
        Fa_N_top_small[t-1] = Fa_N_top[np.int(t/divt+oddvector[0,t-1]-1)]

    Fa_E_down = Fa_E_down_small
    Fa_N_down = Fa_N_down_small
    Fa_E_top = Fa_E_top_small
    Fa_N_top = Fa_N_top_small
    
    return Fa_E_top,Fa_N_top,Fa_E_down,Fa_N_down,E,P,W_top,W_down


#%% Code
def get_stablefluxes(W_top,W_down,Fa_E_top_1,Fa_E_down_1,Fa_N_top_1,Fa_N_down_1,
                                   timestep,divt,L_EW_gridcell,density_water,L_N_gridcell,L_S_gridcell,latitude):
    
    #redefine according to units
    Fa_E_top_kgpmps = Fa_E_top_1
    Fa_E_down_kgpmps = Fa_E_down_1
    Fa_N_top_kgpmps = Fa_N_top_1
    Fa_N_down_kgpmps = Fa_N_down_1
    
    #convert to m3
    Fa_E_top = Fa_E_top_kgpmps * timestep/np.float(divt) * L_EW_gridcell / density_water # [kg*m^-1*s^-1*s*m*kg^-1*m^3]=[m3]
    Fa_E_down = Fa_E_down_kgpmps * timestep/np.float(divt) * L_EW_gridcell / density_water # [s*m*kg*m^-1*s^-1*kg^-1*m^3]=[m3]

    Fa_N_top_swap = np.zeros((len(latitude),np.int(count_time*np.float(divt)),len(longitude)))
    Fa_N_down_swap = np.zeros((len(latitude),np.int(count_time*np.float(divt)),len(longitude)))
    Fa_N_top_kgpmps_swap = np.swapaxes(Fa_N_top_kgpmps,0,1)
    Fa_N_down_kgpmps_swap = np.swapaxes(Fa_N_down_kgpmps,0,1)
    for c in range(len(latitude)):
        Fa_N_top_swap[c] = Fa_N_top_kgpmps_swap[c] * timestep/np.float(divt) * 0.5 *(L_N_gridcell[c]+L_S_gridcell[c]) / density_water # [s*m*kg*m^-1*s^-1*kg^-1*m^3]=[m3]
        Fa_N_down_swap[c] = Fa_N_down_kgpmps_swap[c] * timestep/np.float(divt) * 0.5*(L_N_gridcell[c]+L_S_gridcell[c]) / density_water # [s*m*kg*m^-1*s^-1*kg^-1*m^3]=[m3]

    Fa_N_top = np.swapaxes(Fa_N_top_swap,0,1) 
    Fa_N_down = np.swapaxes(Fa_N_down_swap,0,1) 
    
    #find out where the negative fluxes are
    Fa_E_top_posneg = np.ones(np.shape(Fa_E_top))
    Fa_E_top_posneg[Fa_E_top < 0] = -1
    Fa_N_top_posneg = np.ones(np.shape(Fa_E_top))
    Fa_N_top_posneg[Fa_N_top < 0] = -1
    Fa_E_down_posneg = np.ones(np.shape(Fa_E_top))
    Fa_E_down_posneg[Fa_E_down < 0] = -1
    Fa_N_down_posneg = np.ones(np.shape(Fa_E_top))
    Fa_N_down_posneg[Fa_N_down < 0] = -1
    
    #make everything absolute   
    Fa_E_top_abs = np.abs(Fa_E_top)
    Fa_E_down_abs = np.abs(Fa_E_down)
    Fa_N_top_abs = np.abs(Fa_N_top)
    Fa_N_down_abs = np.abs(Fa_N_down)
    
    # stabilize the outfluxes / influxes
    stab = 1./2.  # during the reduced timestep the water cannot move further than 1/x * the gridcell, 
                    #in other words at least x * the reduced timestep is needed to cross a gridcell
    Fa_E_top_stable = np.reshape(np.minimum(np.reshape(Fa_E_top_abs, (np.size(Fa_E_top_abs))), (np.reshape(Fa_E_top_abs, (np.size(Fa_E_top_abs)))  / 
                                    (np.reshape(Fa_E_top_abs, (np.size(Fa_E_top_abs)))  + np.reshape(Fa_N_top_abs, (np.size(Fa_N_top_abs))))) * stab 
                                            * np.reshape(W_top[:-1,:,:], (np.size(W_top[:-1,:,:])))),(np.int(count_time*np.float(divt)),len(latitude),len(longitude)))
    Fa_N_top_stable = np.reshape(np.minimum(np.reshape(Fa_N_top_abs, (np.size(Fa_N_top_abs))), (np.reshape(Fa_N_top_abs, (np.size(Fa_N_top_abs)))  / 
                                    (np.reshape(Fa_E_top_abs, (np.size(Fa_E_top_abs)))  + np.reshape(Fa_N_top_abs, (np.size(Fa_N_top_abs))))) * stab 
                                            * np.reshape(W_top[:-1,:,:], (np.size(W_top[:-1,:,:])))),(np.int(count_time*np.float(divt)),len(latitude),len(longitude)))
    Fa_E_down_stable = np.reshape(np.minimum(np.reshape(Fa_E_down_abs, (np.size(Fa_E_down_abs))), (np.reshape(Fa_E_down_abs, (np.size(Fa_E_down_abs)))  / 
                                    (np.reshape(Fa_E_down_abs, (np.size(Fa_E_down_abs)))  + np.reshape(Fa_N_down_abs, (np.size(Fa_N_down_abs))))) * stab 
                                             * np.reshape(W_down[:-1,:,:], (np.size(W_down[:-1,:,:])))),(np.int(count_time*np.float(divt)),len(latitude),len(longitude)))
    Fa_N_down_stable = np.reshape(np.minimum(np.reshape(Fa_N_down_abs, (np.size(Fa_N_down_abs))), (np.reshape(Fa_N_down_abs, (np.size(Fa_N_down_abs)))  / 
                                    (np.reshape(Fa_E_down_abs, (np.size(Fa_E_down_abs)))  + np.reshape(Fa_N_down_abs, (np.size(Fa_N_down_abs))))) * stab 
                                             * np.reshape(W_down[:-1,:,:], (np.size(W_down[:-1,:,:])))),(np.int(count_time*np.float(divt)),len(latitude),len(longitude)))
    
    #get rid of the nan values
    Fa_E_top_stable[np.isnan(Fa_E_top_stable)] = 0
    Fa_N_top_stable[np.isnan(Fa_N_top_stable)] = 0
    Fa_E_down_stable[np.isnan(Fa_E_down_stable)] = 0
    Fa_N_down_stable[np.isnan(Fa_N_down_stable)] = 0
    
    #redefine
    Fa_E_top = Fa_E_top_stable * Fa_E_top_posneg
    Fa_N_top = Fa_N_top_stable * Fa_N_top_posneg
    Fa_E_down = Fa_E_down_stable * Fa_E_down_posneg
    Fa_N_down = Fa_N_down_stable * Fa_N_down_posneg
    
    return Fa_E_top,Fa_E_down,Fa_N_top,Fa_N_down

#%% Code
def getFa_Vert(Fa_E_top,Fa_E_down,Fa_N_top,Fa_N_down,E,P,W_top,W_down,divt,count_time,latitude,longitude):
    
    #total moisture in the column
    W = W_top + W_down
    
    #define the horizontal fluxes over the boundaries
    # fluxes over the eastern boundary
    Fa_E_top_boundary = np.zeros(np.shape(Fa_E_top))
    Fa_E_top_boundary[:,:,:-1] = 0.5 * (Fa_E_top[:,:,:-1] + Fa_E_top[:,:,1:])
    if isglobal == 1:
        Fa_E_top_boundary[:,:,-1] = 0.5 * (Fa_E_top[:,:,-1] + Fa_E_top[:,:,0])
    Fa_E_down_boundary = np.zeros(np.shape(Fa_E_down))
    Fa_E_down_boundary[:,:,:-1] = 0.5 * (Fa_E_down[:,:,:-1] + Fa_E_down[:,:,1:])
    if isglobal == 1:
        Fa_E_down_boundary[:,:,-1] = 0.5 * (Fa_E_down[:,:,-1] + Fa_E_down[:,:,0])

    # find out where the positive and negative fluxes are
    Fa_E_top_pos = np.ones(np.shape(Fa_E_top))
    Fa_E_down_pos = np.ones(np.shape(Fa_E_down))
    Fa_E_top_pos[Fa_E_top_boundary < 0] = 0
    Fa_E_down_pos[Fa_E_down_boundary < 0] = 0
    Fa_E_top_neg = Fa_E_top_pos - 1
    Fa_E_down_neg = Fa_E_down_pos - 1

    # separate directions west-east (all positive numbers)
    Fa_E_top_WE = Fa_E_top_boundary * Fa_E_top_pos;
    Fa_E_top_EW = Fa_E_top_boundary * Fa_E_top_neg;
    Fa_E_down_WE = Fa_E_down_boundary * Fa_E_down_pos;
    Fa_E_down_EW = Fa_E_down_boundary * Fa_E_down_neg;

    # fluxes over the western boundary
    Fa_W_top_WE = np.nan*np.zeros(np.shape(P))
    Fa_W_top_WE[:,:,1:] = Fa_E_top_WE[:,:,:-1]
    Fa_W_top_WE[:,:,0] = Fa_E_top_WE[:,:,-1]
    Fa_W_top_EW = np.nan*np.zeros(np.shape(P))
    Fa_W_top_EW[:,:,1:] = Fa_E_top_EW[:,:,:-1]
    Fa_W_top_EW[:,:,0] = Fa_E_top_EW[:,:,-1]
    Fa_W_down_WE = np.nan*np.zeros(np.shape(P))
    Fa_W_down_WE[:,:,1:] = Fa_E_down_WE[:,:,:-1]
    Fa_W_down_WE[:,:,0] = Fa_E_down_WE[:,:,-1]
    Fa_W_down_EW = np.nan*np.zeros(np.shape(P))
    Fa_W_down_EW[:,:,1:] = Fa_E_down_EW[:,:,:-1]
    Fa_W_down_EW[:,:,0] = Fa_E_down_EW[:,:,-1]    

    # fluxes over the northern boundary
    Fa_N_top_boundary = np.nan*np.zeros(np.shape(Fa_N_top));
    Fa_N_top_boundary[:,1:,:] = 0.5 * ( Fa_N_top[:,:-1,:] + Fa_N_top[:,1:,:] )
    Fa_N_down_boundary = np.nan*np.zeros(np.shape(Fa_N_down));
    Fa_N_down_boundary[:,1:,:] = 0.5 * ( Fa_N_down[:,:-1,:] + Fa_N_down[:,1:,:] )

    # find out where the positive and negative fluxes are
    Fa_N_top_pos = np.ones(np.shape(Fa_N_top))
    Fa_N_down_pos = np.ones(np.shape(Fa_N_down))
    Fa_N_top_pos[Fa_N_top_boundary < 0] = 0
    Fa_N_down_pos[Fa_N_down_boundary < 0] = 0
    Fa_N_top_neg = Fa_N_top_pos - 1
    Fa_N_down_neg = Fa_N_down_pos - 1

    # separate directions south-north (all positive numbers)
    Fa_N_top_SN = Fa_N_top_boundary * Fa_N_top_pos
    Fa_N_top_NS = Fa_N_top_boundary * Fa_N_top_neg
    Fa_N_down_SN = Fa_N_down_boundary * Fa_N_down_pos
    Fa_N_down_NS = Fa_N_down_boundary * Fa_N_down_neg

    # fluxes over the southern boundary
    Fa_S_top_SN = np.nan*np.zeros(np.shape(P))
    Fa_S_top_SN[:,:-1,:] = Fa_N_top_SN[:,1:,:]
    Fa_S_top_NS = np.nan*np.zeros(np.shape(P))
    Fa_S_top_NS[:,:-1,:] = Fa_N_top_NS[:,1:,:]
    Fa_S_down_SN = np.nan*np.zeros(np.shape(P))
    Fa_S_down_SN[:,:-1,:] = Fa_N_down_SN[:,1:,:]
    Fa_S_down_NS = np.nan*np.zeros(np.shape(P))
    Fa_S_down_NS[:,:-1,:] = Fa_N_down_NS[:,1:,:]

    # check the water balance
    Sa_after_Fa_down = np.zeros([1,len(latitude),len(longitude)])
    Sa_after_Fa_top = np.zeros([1,len(latitude),len(longitude)])
    Sa_after_all_down = np.zeros([1,len(latitude),len(longitude)])
    Sa_after_all_top = np.zeros([1,len(latitude),len(longitude)])
    residual_down = np.zeros(np.shape(P)) # residual factor [m3]
    residual_top = np.zeros(np.shape(P)) # residual factor [m3]

    for t in range(np.int(count_time*divt)):
        # down: calculate with moisture fluxes:
        Sa_after_Fa_down[0,1:-1,:] = (W_down[t,1:-1,:] - Fa_E_down_WE[t,1:-1,:] + Fa_E_down_EW[t,1:-1,:] + Fa_W_down_WE[t,1:-1,:] - Fa_W_down_EW[t,1:-1,:] - Fa_N_down_SN[t,1:-1,:] 
                                     + Fa_N_down_NS[t,1:-1,:] + Fa_S_down_SN[t,1:-1,:] - Fa_S_down_NS[t,1:-1,:])

        # top: calculate with moisture fluxes:
        Sa_after_Fa_top[0,1:-1,:] = (W_top[t,1:-1,:]- Fa_E_top_WE[t,1:-1,:] + Fa_E_top_EW[t,1:-1,:] + Fa_W_top_WE[t,1:-1,:] - Fa_W_top_EW[t,1:-1,:] - Fa_N_top_SN[t,1:-1,:] 
                                     + Fa_N_top_NS[t,1:-1,:] + Fa_S_top_SN[t,1:-1,:]- Fa_S_top_NS[t,1:-1,:])
    
        # down: substract precipitation and add evaporation
        Sa_after_all_down[0,1:-1,:] = Sa_after_Fa_down[0,1:-1,:] - P[t,1:-1,:] * (W_down[t,1:-1,:] / W[t,1:-1,:]) + E[t,1:-1,:]
    
        # top: substract precipitation
        Sa_after_all_top[0,1:-1,:] = Sa_after_Fa_top[0,1:-1,:] - P[t,1:-1,:] * (W_top[t,1:-1,:] / W[t,1:-1,:])
    
        # down: calculate the residual
        residual_down[t,1:-1,:] = W_down[t+1,1:-1,:] - Sa_after_all_down[0,1:-1,:]
    
        # top: calculate the residual
        residual_top[t,1:-1,:] = W_top[t+1,1:-1,:] - Sa_after_all_top[0,1:-1,:]

    # compute the resulting vertical moisture flux
    Fa_Vert_raw = W_down[1:,:,:] / W[1:,:,:] * (residual_down + residual_top) - residual_down # the vertical velocity so that the new residual_down/W_down =  residual_top/W_top (positive downward)

    # find out where the negative vertical flux is
    Fa_Vert_posneg = np.ones(np.shape(Fa_Vert_raw))
    Fa_Vert_posneg[Fa_Vert_raw < 0] = -1

    # make the vertical flux absolute
    Fa_Vert_abs = np.abs(Fa_Vert_raw)

    # stabilize the outfluxes / influxes
    stab = 1./4. #during the reduced timestep the vertical flux can maximally empty/fill 1/x of the top or down storage
    
    Fa_Vert_stable = np.reshape(np.minimum(np.reshape(Fa_Vert_abs, (np.size(Fa_Vert_abs))), np.minimum(stab*np.reshape(W_top[1:,:,:], (np.size(W_top[1:,:,:]))), 
                                        stab*np.reshape(W_down[1:,:,:], (np.size(W_down[1:,:,:]))))),(np.int(count_time*np.float(divt)),len(latitude),len(longitude)))
                
    # redefine the vertical flux
    Fa_Vert = Fa_Vert_stable * Fa_Vert_posneg;

    return Fa_Vert_raw, Fa_Vert

# #### End of code

#%% Runtime & Results

start1 = timer()

# obtain the constants
latitude,longitude,lsm,g,density_water,timestep,A_gridcell,L_N_gridcell,L_S_gridcell,L_EW_gridcell,gridcell = getconstants(latnrs,lonnrs,lake_mask,invariant_data)

# loop through the years
for yearnumber in years:

    ly = int(calendar.isleap(yearnumber))
    final_time = 364 + ly  # number of parts-1 to divide a year in
    
    for a in yearpart:  # a > 365 (366th index) and not a leapyear
        start = timer()

        datapath = data_path(yearnumber,a) # global variable
        if a > final_time:
            pass
            # do nothing
        else:
            begin_time = a*4 # first index to get data from (netcdf is zero based) (leave at a*16)
            
            #1 integrate specific humidity to get the (total) column water (vapor)
            cw,W_top,W_down = getW(latnrs,lonnrs,final_time,a,yearnumber,begin_time,count_time,density_water,latitude,longitude,g,A_gridcell,boundary)
            
            #2 wind in between pressure levels
            U,V = getwind(latnrs,lonnrs,final_time,a,yearnumber,begin_time,count_time)
            
            #3 calculate horizontal moisture fluxes
            Fa_E_top,Fa_N_top,Fa_E_down,Fa_N_down = getFa(latnrs,lonnrs,boundary,cw,U,V,count_time,begin_time,yearnumber,a,final_time)
            
            #4 evaporation and precipitation
            E,P = getEP(latnrs,lonnrs,yearnumber,begin_time,count_time,latitude,longitude,A_gridcell)
            
            #5 put data on a smaller time step
            Fa_E_top_1,Fa_N_top_1,Fa_E_down_1,Fa_N_down_1,E,P,W_top,W_down = getrefined(Fa_E_top,Fa_N_top,Fa_E_down,Fa_N_down,W_top,W_down,E,P,divt,count_time,latitude,longitude)

            #6 stabilize horizontal fluxes and get everything in (m3 per smaller timestep)
            Fa_E_top,Fa_E_down,Fa_N_top,Fa_N_down = get_stablefluxes(W_top,W_down,Fa_E_top_1,Fa_E_down_1,Fa_N_top_1,Fa_N_down_1,
                                   timestep,divt,L_EW_gridcell,density_water,L_N_gridcell,L_S_gridcell,latitude)
            
            #7 determine the vertical moisture flux
            Fa_Vert_raw,Fa_Vert = getFa_Vert(Fa_E_top,Fa_E_down,Fa_N_top,Fa_N_down,E,P,W_top,W_down,divt,count_time,latitude,longitude)
            
            sio.savemat(datapath[23], {'Fa_E_top':Fa_E_top, 'Fa_N_top':Fa_N_top, 'Fa_E_down':Fa_E_down,'Fa_N_down':Fa_N_down, 'Fa_Vert':Fa_Vert, 'E':E, 'P':P, 
                                                                                    'W_top':W_top, 'W_down':W_down}, do_compression=True)
            
            # alternative, but slower and more spacious
            # np.savez_compressed(datapath[23],Fa_E_top,Fa_N_top,Fa_E_down,Fa_N_down,Fa_Vert,E,P,W_top,W_down)
            
        end = timer()
        print 'Runtime fluxes_and_storages for day ' + str(a+1) + ' in year ' + str(yearnumber) + ' is',(end - start),' seconds.'
end1 = timer()
print 'The total runtime is',(end1-start1),' seconds.'
