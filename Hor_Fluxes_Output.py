# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 13:24:45 2016

@author: Ent00002
"""
#%% Import libraries

import numpy as np
import scipy.io as sio
import calendar
import datetime
from getconstants import getconstants
from timeit import default_timer as timer
import os

#%% BEGIN OF INPUT (FILL THIS IN)
years = np.arange(2010,2011) #fill in the years
yearpart = np.arange(0,364) # for a full (leap)year fill in 0:366

# Manage the extent of your dataset (FILL THIS IN)
# Define the latitude and longitude cell numbers to consider and corresponding lakes that should be considered part of the land
latnrs = np.arange(7,114)
lonnrs = np.arange(0,240)

# the lake numbers below belong to the ERA-Interim data on 1.5 degree starting at Northern latitude 79.5 and longitude -180
lake_mask_1 = np.array([9,9,9,12,12,21,21,22,22,23,24,25,23,23,25,25,53,54,61,23,24,23,24,25,27,22,23,24,25,26,27,28,22,25,26,27,28,23,23,12,18])
lake_mask_2 = np.array([120+19,120+40,120+41,120+43,120+44,120+61,120+62,120+62,120+63,120+62,120+62,120+62,120+65,120+66,120+65,120+66,142-120,142-120,143-120,152-120,152-120,153-120,153-120,153-120,153-120,154-120,154-120,154-120,154-120,154-120,154-120,154-120,155-120,155-120,155-120,155-120,155-120,159-120,160-120,144-120,120+55])
lake_mask = np.transpose(np.vstack((lake_mask_1,lake_mask_2))) #recreate the arrays of the matlab model

daily = 1 # 1 for writing out daily data, 0 for only monthly data

#END OF INPUT

#%% Datapaths (FILL THIS IN)
invariant_data = r'C:\Users\bec\Desktop\WAM2/invariants_15x15.nc'#invariants
interdata_folder = r'C:\Users\bec\Desktop\WAM2\interdata'
output_folder = r'C:\Users\bec\Desktop\WAM2\output'

def data_path(y,a,years):
    load_fluxes_and_storages = os.path.join(interdata_folder, str(y) + '-' + str(a) + 'fluxes_storages.mat')

    save_path = os.path.join(output_folder, 'Hor_Fluxes_full' + str(years[0]) + '-' + str(years[-1]) + '.mat')
    
    save_path_daily = os.path.join(output_folder, 'Hor_Fluxes_daily_full' + str(y) + '.mat')
    
    return load_fluxes_and_storages,save_path,save_path_daily

#%% Runtime % Results

start1 = timer()

# obtain the constants
latitude,longitude,lsm,g,density_water,timestep,A_gridcell,L_N_gridcell,L_S_gridcell,L_EW_gridcell,gridcell = getconstants(latnrs,lonnrs,lake_mask,invariant_data)

startyear = years[0]
Fa_E_down_per_year_per_month = np.zeros((len(years),12,len(latitude),len(longitude)))
Fa_E_top_per_year_per_month = np.zeros((len(years),12,len(latitude),len(longitude)))
Fa_N_down_per_year_per_month = np.zeros((len(years),12,len(latitude),len(longitude)))
Fa_N_top_per_year_per_month = np.zeros((len(years),12,len(latitude),len(longitude)))
Fa_Vert_per_year_per_month = np.zeros((len(years),12,len(latitude),len(longitude)))

for y in years:

    ly = int(calendar.isleap(y))
    final_time = 364+ly
    
    Fa_E_down_per_day = np.zeros((365+ly,len(latitude),len(longitude)))
    Fa_E_top_per_day = np.zeros((365+ly,len(latitude),len(longitude)))
    Fa_N_down_per_day = np.zeros((365+ly,len(latitude),len(longitude)))
    Fa_N_top_per_day = np.zeros((365+ly,len(latitude),len(longitude)))
    Fa_Vert_per_day = np.zeros((365+ly,len(latitude),len(longitude)))
    
    for a in yearpart:
        start = timer()

        datapath = data_path(y,a,years)
        
        if a > final_time: # a = 365 and not a leapyear
            pass
        else:
            # load horizontal fluxes
            loading_FS = sio.loadmat(datapath[0],verify_compressed_data_integrity=False)
            Fa_E_top = loading_FS['Fa_E_top']
            Fa_N_top = loading_FS['Fa_N_top']
            Fa_E_down = loading_FS['Fa_E_down']
            Fa_N_down = loading_FS['Fa_N_down']
            Fa_Vert = loading_FS['Fa_Vert']
                        
            # save per day
            Fa_E_down_per_day[a,:,:] = np.sum(Fa_E_down, axis =0)
            Fa_E_top_per_day[a,:,:] = np.sum(Fa_E_top, axis =0)
            Fa_N_down_per_day[a,:,:] = np.sum(Fa_N_down, axis =0)
            Fa_N_top_per_day[a,:,:] = np.sum(Fa_N_top, axis =0)
            Fa_Vert_per_day[a,:,:] = np.sum(Fa_Vert, axis =0)
            
            # timer
            end = timer()
            print 'Runtime output for day ' + str(a+1) + ' in year ' + str(y) + ' is',(end - start),' seconds.'
    
    # save daily fluxes on disk
    if daily == 1:            
        sio.savemat(datapath[2],
                    {'Fa_E_down_per_day':Fa_E_down_per_day,'Fa_E_top_per_day':Fa_E_top_per_day,
                     'Fa_N_down_per_day':Fa_N_down_per_day,'Fa_N_top_per_day':Fa_N_top_per_day, 
                     'Fa_Vert_per_day':Fa_Vert_per_day},do_compression=True)    
    
    for m in range(12):
        first_day = int(datetime.date(y,m+1,1).strftime("%j"))
        last_day = int(datetime.date(y,m+1,calendar.monthrange(y,m+1)[1]).strftime("%j"))
        days = np.arange(first_day,last_day+1)-1 # -1 because Python is zero-based
       
        Fa_E_down_per_year_per_month[y-startyear,m,:,:] = (np.squeeze(np.sum(Fa_E_down_per_day[days,:,:], axis = 0)))
        Fa_E_top_per_year_per_month[y-startyear,m,:,:] = (np.squeeze(np.sum(Fa_E_top_per_day[days,:,:], axis = 0)))
        Fa_N_down_per_year_per_month[y-startyear,m,:,:] = (np.squeeze(np.sum(Fa_N_down_per_day[days,:,:], axis = 0)))
        Fa_N_top_per_year_per_month[y-startyear,m,:,:] = (np.squeeze(np.sum(Fa_N_top_per_day[days,:,:], axis = 0)))
        Fa_Vert_per_year_per_month[y-startyear,m,:,:] = (np.squeeze(np.sum(Fa_Vert_per_day[days,:,:], axis = 0)))

# save monthly fluxes on disk                         
sio.savemat(datapath[1], 
            {'Fa_E_down_per_year_per_month':Fa_E_down_per_year_per_month,'Fa_E_top_per_year_per_month':Fa_E_top_per_year_per_month, 
            'Fa_N_down_per_year_per_month':Fa_N_down_per_year_per_month,'Fa_N_top_per_year_per_month':Fa_N_top_per_year_per_month, 
            'Fa_Vert_per_year_per_month':Fa_Vert_per_year_per_month})
    
end1 = timer()
print 'The total runtime is',(end1-start1),' seconds.'

