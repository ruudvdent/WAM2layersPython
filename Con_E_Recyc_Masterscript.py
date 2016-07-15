# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 13:24:45 2016

@author: Ent00002
"""
#delayed runs
#import time
#time.sleep(500)

#%% Import libraries
import numpy as np
import scipy.io as sio
import calendar
from getconstants import getconstants
from timeit import default_timer as timer
import os

#%% BEGIN OF INPUT1 (FILL THIS IN)
years = np.arange(2009,2001,-1) # fill in the years backward
yearpart = np.arange(365,-1,-1) # for a full (leap)year fill in (365,-1,-1)
boundary = 8 # with 8 the vertical separation is at 812.83 hPa for surface pressure = 1031.25 hPa, which corresponds to k=47 (ERA-Interim)
divt = 24 # division of the timestep, 24 means a calculation timestep of 6/24 = 0.25 hours (numerical stability purposes)
count_time = 4 # number of indices to get data from (for six hourly data this means everytime one day)

# Manage the extent of your dataset (FILL THIS IN)
# Define the latitude and longitude cell numbers to consider and corresponding lakes that should be considered part of the land
latnrs = np.arange(7,114)
lonnrs = np.arange(0,240)
isglobal = 1 # fill in 1 for global computations (i.e. Earth round), fill in 0 for a local domain with boundaries

# the lake numbers below belong to the ERA-Interim data on 1.5 degree starting at Northern latitude 79.5 and longitude -180
lake_mask_1 = np.array([9,9,9,12,12,21,21,22,22,23,24,25,23,23,25,25,53,54,61,23,24,23,24,25,27,22,23,24,25,26,27,28,22,25,26,27,28,23,23,12,18])
lake_mask_2 = np.array([120+19,120+40,120+41,120+43,120+44,120+61,120+62,120+62,120+63,120+62,120+62,120+62,120+65,120+66,120+65,120+66,142-120,142-120,143-120,152-120,152-120,153-120,153-120,153-120,153-120,154-120,154-120,154-120,154-120,154-120,154-120,154-120,155-120,155-120,155-120,155-120,155-120,159-120,160-120,144-120,120+55])
lake_mask = np.transpose(np.vstack((lake_mask_1,lake_mask_2))) #recreate the arrays of the matlab model

# obtain the constants
invariant_data = 'Interim_data/full/invariants.nc'#invariants
latitude,longitude,lsm,g,density_water,timestep,A_gridcell,L_N_gridcell,L_S_gridcell,L_EW_gridcell,gridcell = getconstants(latnrs,lonnrs,lake_mask,invariant_data)

# BEGIN OF INPUT 2 (FILL THIS IN)
Region = lsm
Kvf = 3 # vertical dispersion factor (advection only is 0, dispersion the same size of the advective flux is 1, for stability don't make this more than 3)
timetracking = 1 # 0 for not tracking time and 1 for tracking time
veryfirstrun = 0 # type '1' if no run has been done before from which can be continued, otherwise type '0'

interdata_folder = r'C:\Users\bec\Desktop\WAM2\interdata' #must be an existing folder, existence is not checked

#END OF INPUT

#%% Datapaths (FILL THIS IN)

# Check if interdata folder exists:
assert os.path.isdir(interdata_folder), "Please create the interdata_folder before running the script"
# Check if sub interdata folder exists otherwise create it:
sub_interdata_folder = os.path.join(interdata_folder, 'continental_backward')
if os.path.isdir(sub_interdata_folder):
    pass
else:
    os.makedirs(sub_interdata_folder)


def data_path_ea(years,yearpart):
    save_empty_arrays_ld_track = os.path.join(sub_interdata_folder, str(years[0]+1) + '-' + str(0) + 'Sa_track.mat')
    save_empty_arrays_ld_time = os.path.join(sub_interdata_folder, str(years[0]+1) + '-' + str(0) + 'Sa_time.mat')
    
    save_empty_arrays_track = os.path.join(sub_interdata_folder, str(years[0]) + '-' + str(yearpart[0]+1) + 'Sa_track.mat')
    save_empty_arrays_time = os.path.join(sub_interdata_folder, str(years[0]) + '-' + str(yearpart[0]+1) + 'Sa_time.mat')
    
    return save_empty_arrays_ld_track,save_empty_arrays_ld_time,save_empty_arrays_track,save_empty_arrays_time

def data_path(previous_data_to_load,yearnumber,a):
    load_Sa_track = os.path.join(sub_interdata_folder, previous_data_to_load + 'Sa_track.mat')
    
    load_fluxes_and_storages = os.path.join(interdata_folder, str(yearnumber) + '-' + str(a) + 'fluxes_storages.mat')
    
    load_Sa_time = os.path.join(sub_interdata_folder, previous_data_to_load + 'Sa_time.mat')

    save_path_track = os.path.join(sub_interdata_folder, str(yearnumber) + '-' + str(a) + 'Sa_track.mat')
    save_path_time = os.path.join(sub_interdata_folder, str(yearnumber) + '-' + str(a) + 'Sa_time.mat')
    return load_Sa_track,load_fluxes_and_storages,load_Sa_time,save_path_track,save_path_time


#%% Code (no need to look at this for running)

def get_Sa_track_backward(latitude,longitude,count_time,divt,Kvf,Region,Fa_E_top,Fa_N_top,Fa_E_down,Fa_N_down,
                                       Fa_Vert,E,P,W_top,W_down,Sa_track_top_last,Sa_track_down_last):
    
    # make P_region matrix
    Region3D = np.tile(np.reshape(Region,[1,len(latitude),len(longitude)]),[len(P[:,0,0]),1,1])
    P_region = Region3D * P
    
    # Total moisture in the column
    W = W_top + W_down

    # separate the direction of the vertical flux and make it absolute
    Fa_upward = np.zeros(np.shape(Fa_Vert))
    Fa_upward[Fa_Vert <= 0 ] = Fa_Vert[Fa_Vert <= 0 ]
    Fa_downward = np.zeros(np.shape(Fa_Vert));
    Fa_downward[Fa_Vert >= 0 ] = Fa_Vert[Fa_Vert >= 0 ]
    Fa_upward = np.abs(Fa_upward)
    
    # include the vertical dispersion
    if Kvf == 0:
        pass 
        # do nothing
    else:
        Fa_upward = (1.+Kvf) * Fa_upward
        Fa_upward[Fa_Vert >= 0] = Fa_Vert[Fa_Vert >= 0] * Kvf
        Fa_downward = (1.+Kvf) * Fa_downward
        Fa_downward[Fa_Vert <= 0] = np.abs(Fa_Vert[Fa_Vert <= 0]) * Kvf
        
    # define the horizontal fluxes over the boundaries
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
        
    # defining size of output
    Sa_track_down = np.zeros(np.shape(W_down))
    Sa_track_top = np.zeros(np.shape(W_top))
    
    # assign begin values of output == last (but first index) values of the previous time slot
    Sa_track_down[-1,:,:] = Sa_track_down_last
    Sa_track_top[-1,:,:] = Sa_track_top_last
    
    # defining sizes of tracked moisture
    Sa_track_after_Fa_down = np.zeros(np.shape(Sa_track_down_last))
    Sa_track_after_Fa_P_E_down = np.zeros(np.shape(Sa_track_down_last))
    Sa_track_E_down = np.zeros(np.shape(Sa_track_down_last))
    Sa_track_W_down = np.zeros(np.shape(Sa_track_down_last))
    Sa_track_N_down = np.zeros(np.shape(Sa_track_down_last))
    Sa_track_S_down = np.zeros(np.shape(Sa_track_down_last))
    Sa_track_after_Fa_top = np.zeros(np.shape(Sa_track_top_last))
    Sa_track_after_Fa_P_E_top = np.zeros(np.shape(Sa_track_top_last))
    Sa_track_E_top = np.zeros(np.shape(Sa_track_top_last))
    Sa_track_W_top = np.zeros(np.shape(Sa_track_top_last))
    Sa_track_N_top = np.zeros(np.shape(Sa_track_top_last))
    Sa_track_S_top = np.zeros(np.shape(Sa_track_top_last))

    # define sizes of total moisture
    Sa_E_down = np.zeros(np.shape(Sa_track_down_last))
    Sa_W_down = np.zeros(np.shape(Sa_track_down_last))
    Sa_N_down = np.zeros(np.shape(Sa_track_down_last))
    Sa_S_down = np.zeros(np.shape(Sa_track_down_last))
    Sa_E_top = np.zeros(np.shape(Sa_track_top_last))
    Sa_W_top = np.zeros(np.shape(Sa_track_top_last))
    Sa_N_top = np.zeros(np.shape(Sa_track_top_last))
    Sa_S_top = np.zeros(np.shape(Sa_track_top_last))
    
    # define variables that find out what happens to the water
    north_loss = np.zeros((np.int(count_time*divt),1,len(longitude)))
    south_loss = np.zeros((np.int(count_time*divt),1,len(longitude)))
    down_to_top = np.zeros(np.shape(P))
    top_to_down = np.zeros(np.shape(P))
    water_lost = np.zeros(np.shape(P))
    water_lost_down = np.zeros(np.shape(P))
    water_lost_top = np.zeros(np.shape(P))
    
    # Sa calculation backward in time
    for t in np.arange(np.int(count_time*divt),0,-1):
        # down: define values of total moisture
        Sa_E_down[0,:,:-1] = W_down[t,:,1:] # Atmospheric storage of the cell to the east [m3]
        # to make dependent on isglobal but for now kept to avoid division by zero errors     
        Sa_E_down[0,:,-1] = W_down[t,:,0] # Atmospheric storage of the cell to the east [m3]
        Sa_W_down[0,:,1:] = W_down[t,:,:-1] # Atmospheric storage of the cell to the west [m3]
        # to make dependent on isglobal but for now kept to avoid division by zero errors      
        Sa_W_down[0,:,0] = W_down[t,:,-1] # Atmospheric storage of the cell to the west [m3]
        Sa_N_down[0,1:,:] = W_down[t,0:-1,:] # Atmospheric storage of the cell to the north [m3]
        Sa_S_down[0,:-1,:] = W_down[t,1:,:] # Atmospheric storage of the cell to the south [m3]
    
        # top: define values of total moisture
        Sa_E_top[0,:,:-1] = W_top[t,:,1:] # Atmospheric storage of the cell to the east [m3]
        # to make dependent on isglobal but for now kept to avoid division by zero errors      
        Sa_E_top[0,:,-1] = W_top[t,:,0] # Atmospheric storage of the cell to the east [m3]
        Sa_W_top[0,:,1:] = W_top[t,:,:-1] # Atmospheric storage of the cell to the west [m3]
        # to make dependent on isglobal but for now kept to avoid division by zero errors      
        Sa_W_top[0,:,0] = W_top[t,:,-1] # Atmospheric storage of the cell to the west [m3]
        Sa_N_top[0,1:,:] = W_top[t,:-1,:] # Atmospheric storage of the cell to the north [m3]
        Sa_S_top[0,:-1,:] = W_top[t,1:,:] # Atmospheric storage of the cell to the south [m3]
        
         # down: define values of tracked moisture of neighbouring grid cells
        Sa_track_E_down[0,:,:-1] = Sa_track_down[t,:,1:] # Atmospheric tracked storage of the cell to the east [m3]
        if isglobal == 1:       
            Sa_track_E_down[0,:,-1] = Sa_track_down[t,:,0] #Atmospheric tracked storage of the cell to the east [m3]
        Sa_track_W_down[0,:,1:] = Sa_track_down[t,:,:-1] # Atmospheric storage of the cell to the west [m3]
        if isglobal == 1:       
            Sa_track_W_down[0,:,0] = Sa_track_down[t,:,-1] # Atmospheric storage of the cell to the west [m3]
        Sa_track_N_down[0,1:,:] = Sa_track_down[t,:-1,:] # Atmospheric storage of the cell to the north [m3]
        Sa_track_S_down[0,:-1,:] = Sa_track_down[t,1:,:] # Atmospheric storage of the cell to the south [m3]
        
        # down: calculate with moisture fluxes
        Sa_track_after_Fa_down[0,1:-1,:] = (Sa_track_down[t,1:-1,:] 
        + Fa_E_down_WE[t-1,1:-1,:] * (Sa_track_E_down[0,1:-1,:] / Sa_E_down[0,1:-1,:]) 
        - Fa_E_down_EW[t-1,1:-1,:] * (Sa_track_down[t,1:-1,:] / W_down[t,1:-1,:]) 
        - Fa_W_down_WE[t-1,1:-1,:] * (Sa_track_down[t,1:-1,:] / W_down[t,1:-1,:]) 
        + Fa_W_down_EW[t-1,1:-1,:] * (Sa_track_W_down[0,1:-1,:] / Sa_W_down[0,1:-1,:]) 
        + Fa_N_down_SN[t-1,1:-1,:] * (Sa_track_N_down[0,1:-1,:] / Sa_N_down[0,1:-1,:]) 
        - Fa_N_down_NS[t-1,1:-1,:] * (Sa_track_down[t,1:-1,:] / W_down[t,1:-1,:]) 
        - Fa_S_down_SN[t-1,1:-1,:] * (Sa_track_down[t,1:-1,:] / W_down[t,1:-1,:]) 
        + Fa_S_down_NS[t-1,1:-1,:] * (Sa_track_S_down[0,1:-1,:] / Sa_S_down[0,1:-1,:])
        - Fa_downward[t-1,1:-1,:] * (Sa_track_down[t,1:-1,:] / W_down[t,1:-1,:]) 
        + Fa_upward[t-1,1:-1,:] * (Sa_track_top[t,1:-1,:] / W_top[t,1:-1,:]))
        
        # top: define values of tracked moisture of neighbouring grid cells
        Sa_track_E_top[0,:,:-1] = Sa_track_top[t,:,1:] # Atmospheric tracked storage of the cell to the east [m3]
        if isglobal == 1:       
            Sa_track_E_top[0,:,-1] = Sa_track_top[t,:,0] # Atmospheric tracked storage of the cell to the east [m3]
        Sa_track_W_top[0,:,1:] = Sa_track_top[t,:,:-1] # Atmospheric tracked storage of the cell to the west [m3]
        if isglobal == 1:       
            Sa_track_W_top[0,:,0] = Sa_track_top[t,:,-1] # Atmospheric tracked storage of the cell to the west [m3]
        Sa_track_N_top[0,1:,:] = Sa_track_top[t,:-1,:] # Atmospheric tracked storage of the cell to the north [m3]
        Sa_track_S_top[0,:-1,:] = Sa_track_top[t,1:,:] # Atmospheric tracked storage of the cell to the south [m3]
        
        # top: calculate with moisture fluxes 
        Sa_track_after_Fa_top[0,1:-1,:] = (Sa_track_top[t,1:-1,:]
        + Fa_E_top_WE[t-1,1:-1,:] * (Sa_track_E_top[0,1:-1,:] / Sa_E_top[0,1:-1,:]) 
        - Fa_E_top_EW[t-1,1:-1,:] * (Sa_track_top[t,1:-1,:] / W_top[t,1:-1,:]) 
        - Fa_W_top_WE[t-1,1:-1,:] * (Sa_track_top[t,1:-1,:] / W_top[t,1:-1,:]) 
        + Fa_W_top_EW[t-1,1:-1,:] * (Sa_track_W_top[0,1:-1,:] / Sa_W_top[0,1:-1,:])
        + Fa_N_top_SN[t-1,1:-1,:] * (Sa_track_N_top[0,1:-1,:] / Sa_N_top[0,1:-1,:]) 
        - Fa_N_top_NS[t-1,1:-1,:] * (Sa_track_top[t,1:-1,:]/ W_top[t,1:-1,:]) 
        - Fa_S_top_SN[t-1,1:-1,:] * (Sa_track_top[t,1:-1,:] / W_top[t,1:-1,:]) 
        + Fa_S_top_NS[t-1,1:-1,:] * (Sa_track_S_top[0,1:-1,:] / Sa_S_top[0,1:-1,:]) 
        + Fa_downward[t-1,1:-1,:] * (Sa_track_down[t,1:-1,:] / W_down[t,1:-1,:]) 
        - Fa_upward[t-1,1:-1,:] * (Sa_track_top[t,1:-1,:] / W_top[t,1:-1,:]))
        
        # losses to the north and south
        north_loss[t-1,0,:] = (Fa_N_top_NS[t-1,1,:] * (Sa_track_top[t,1,:] / W_top[t,1,:])
                                + Fa_N_down_NS[t-1,1,:] * (Sa_track_down[t,1,:] / W_down[t,1,:]))
        south_loss[t-1,0,:] = (Fa_S_top_SN[t-1,-2,:] * (Sa_track_top[t,-2,:] / W_top[t,-2,:])
                                + Fa_S_down_SN[t-1,-2,:] * (Sa_track_down[t,-2,:] / W_down[t,-2,:]))
    
        # down: add precipitation and subtract evaporation
        Sa_track_after_Fa_P_E_down[0,1:-1,:] = (Sa_track_after_Fa_down[0,1:-1,:]
        + P_region[t-1,1:-1,:] * (W_down[t,1:-1,:] / W[t,1:-1,:])   
        - E[t-1,1:-1,:] * (Sa_track_down[t,1:-1,:] / W_down[t,1:-1,:]))
    
        # top: add precipitation
        Sa_track_after_Fa_P_E_top[0,1:-1,:] = (Sa_track_after_Fa_top[0,1:-1,:] 
        + P_region[t-1,1:-1,:] * (W_top[t,1:-1,:] / W[t,1:-1,:])) 
        
        # down and top: redistribute unaccounted water that is otherwise lost from the sytem
        down_to_top[t-1,:,:] = np.reshape(np.maximum(0, np.reshape(Sa_track_after_Fa_P_E_down, (np.size(Sa_track_after_Fa_P_E_down))) - np.reshape(W_down[t-1,:,:],
                                            (np.size(W_down[t-1,:,:])))), (len(latitude),len(longitude)))
        top_to_down[t-1,:,:] = np.reshape(np.maximum(0, np.reshape(Sa_track_after_Fa_P_E_top, (np.size(Sa_track_after_Fa_P_E_top))) - np.reshape(W_top[t-1,:,:],
                                            (np.size(W_top[t-1,:,:])))), (len(latitude),len(longitude)))
        Sa_track_after_all_down = Sa_track_after_Fa_P_E_down - down_to_top[t-1,:,:] + top_to_down[t-1,:,:]
        Sa_track_after_all_top = Sa_track_after_Fa_P_E_top - top_to_down[t-1,:,:] + down_to_top[t-1,:,:]

        # down and top: water lost to the system:
        water_lost_down[t-1,:,:] = np.reshape(np.maximum(0, np.reshape(Sa_track_after_all_down, (np.size(Sa_track_after_all_down))) - np.reshape(W_down[t-1,:,:],
                                            (np.size(W_down[t-1,:,:])))), (len(latitude),len(longitude)))
        water_lost_top[t-1,:,:] = np.reshape(np.maximum(0, np.reshape(Sa_track_after_all_top, (np.size(Sa_track_after_all_top))) - np.reshape(W_top[t-1,:,:],
                                            (np.size(W_top[t-1,:,:])))), (len(latitude),len(longitude)))
        water_lost = water_lost_down + water_lost_top
    
        # down: determine Sa_region of this next timestep 100% stable
        Sa_track_down[t-1,1:-1,:] = np.reshape(np.maximum(0,np.minimum(np.reshape(W_down[t-1,1:-1,:], np.size(W_down[t-1,1:-1,:])), np.reshape(Sa_track_after_all_down[0,1:-1,:],
                                                np.size(Sa_track_after_all_down[0,1:-1,:])))), (len(latitude[1:-1]),len(longitude)))
        # top: determine Sa_region of this next timestep 100% stable
        Sa_track_top[t-1,1:-1,:] = np.reshape(np.maximum(0,np.minimum(np.reshape(W_top[t-1,1:-1,:], np.size(W_top[t-1,1:-1,:])), np.reshape(Sa_track_after_all_top[0,1:-1,:],
                                                np.size(Sa_track_after_all_top[0,1:-1,:])))), (len(latitude[1:-1]),len(longitude)))
  
    return Sa_track_top,Sa_track_down,north_loss,south_loss,down_to_top,top_to_down,water_lost

#%% Code

def get_Sa_track_backward_TIME(latitude,longitude,count_time,divt,timestep,Kvf,Region,Fa_E_top,Fa_N_top,Fa_E_down,Fa_N_down,Fa_Vert,E,P,
                                            W_top,W_down,Sa_track_top_last,Sa_track_down_last,Sa_time_top_last,Sa_time_down_last):

    # make P_region matrix
    Region3D = np.tile(np.reshape(Region,[1,len(latitude),len(longitude)]),[len(P[:,0,0]),1,1])
    P_region = Region3D * P
    
    # Total moisture in the column
    W = W_top + W_down

    # separate the direction of the vertical flux and make it absolute
    Fa_upward = np.zeros(np.shape(Fa_Vert))
    Fa_upward[Fa_Vert <= 0 ] = Fa_Vert[Fa_Vert <= 0 ]
    Fa_downward = np.zeros(np.shape(Fa_Vert));
    Fa_downward[Fa_Vert >= 0 ] = Fa_Vert[Fa_Vert >= 0 ]
    Fa_upward = np.abs(Fa_upward)
    
    # include the vertical dispersion
    if Kvf == 0:
        pass 
        # do nothing
    else:
        Fa_upward = (1.+Kvf) * Fa_upward
        Fa_upward[Fa_Vert >= 0] = Fa_Vert[Fa_Vert >= 0] * Kvf
        Fa_downward = (1.+Kvf) * Fa_downward
        Fa_downward[Fa_Vert <= 0] = np.abs(Fa_Vert[Fa_Vert <= 0]) * Kvf
        
    # define the horizontal fluxes over the boundaries
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
        
    # defining size of output
    Sa_track_down = np.zeros(np.shape(W_down))
    Sa_track_top = np.zeros(np.shape(W_top))
    Sa_time_down = np.zeros(np.shape(W_down))
    Sa_time_top = np.zeros(np.shape(W_top))
    
    # assign begin values of output == last (but first index) values of the previous time slot
    Sa_track_down[-1,:,:] = Sa_track_down_last
    Sa_track_top[-1,:,:] = Sa_track_top_last
    Sa_time_down[-1,:,:] = Sa_time_down_last
    Sa_time_top[-1,:,:] = Sa_time_top_last
    
    # defining sizes of tracked moisture
    Sa_track_after_Fa_down = np.zeros(np.shape(Sa_track_down_last))
    Sa_track_after_Fa_P_E_down = np.zeros(np.shape(Sa_track_down_last))
    Sa_track_E_down = np.zeros(np.shape(Sa_track_down_last))
    Sa_track_W_down = np.zeros(np.shape(Sa_track_down_last))
    Sa_track_N_down = np.zeros(np.shape(Sa_track_down_last))
    Sa_track_S_down = np.zeros(np.shape(Sa_track_down_last))
    Sa_track_after_Fa_top = np.zeros(np.shape(Sa_track_top_last))
    Sa_track_after_Fa_P_E_top = np.zeros(np.shape(Sa_track_top_last))
    Sa_track_E_top = np.zeros(np.shape(Sa_track_top_last))
    Sa_track_W_top = np.zeros(np.shape(Sa_track_top_last))
    Sa_track_N_top = np.zeros(np.shape(Sa_track_top_last))
    Sa_track_S_top = np.zeros(np.shape(Sa_track_top_last))

    # define sizes of total moisture
    Sa_E_down = np.zeros(np.shape(Sa_track_down_last))
    Sa_W_down = np.zeros(np.shape(Sa_track_down_last))
    Sa_N_down = np.zeros(np.shape(Sa_track_down_last))
    Sa_S_down = np.zeros(np.shape(Sa_track_down_last))
    Sa_E_top = np.zeros(np.shape(Sa_track_top_last))
    Sa_W_top = np.zeros(np.shape(Sa_track_top_last))
    Sa_N_top = np.zeros(np.shape(Sa_track_top_last))
    Sa_S_top = np.zeros(np.shape(Sa_track_top_last))
    
    # define variables that find out what happens to the water
    north_loss = np.zeros((np.int(count_time*divt),1,len(longitude)))
    south_loss = np.zeros((np.int(count_time*divt),1,len(longitude)))
    down_to_top = np.zeros(np.shape(P))
    top_to_down = np.zeros(np.shape(P))
    water_lost = np.zeros(np.shape(P))
    water_lost_down = np.zeros(np.shape(P))
    water_lost_top = np.zeros(np.shape(P))
    
    # Sa calculation backward in time
    for t in np.arange(np.int(count_time*divt),0,-1):
        # down: define values of total moisture
        Sa_E_down[0,:,:-1] = W_down[t,:,1:] # Atmospheric storage of the cell to the east [m3]
        # to make dependent on isglobal but for now kept to avoid division by zero errors      
        Sa_E_down[0,:,-1] = W_down[t,:,0] # Atmospheric storage of the cell to the east [m3]
        Sa_W_down[0,:,1:] = W_down[t,:,:-1] # Atmospheric storage of the cell to the west [m3]
        # to make dependent on isglobal but for now kept to avoid division by zero errors      
        Sa_W_down[0,:,0] = W_down[t,:,-1] # Atmospheric storage of the cell to the west [m3]
        Sa_N_down[0,1:,:] = W_down[t,0:-1,:] # Atmospheric storage of the cell to the north [m3]
        Sa_S_down[0,:-1,:] = W_down[t,1:,:] # Atmospheric storage of the cell to the south [m3]
    
        # top: define values of total moisture
        Sa_E_top[0,:,:-1] = W_top[t,:,1:] # Atmospheric storage of the cell to the east [m3]
        # to make dependent on isglobal but for now kept to avoid division by zero errors      
        Sa_E_top[0,:,-1] = W_top[t,:,0] # Atmospheric storage of the cell to the east [m3]
        Sa_W_top[0,:,1:] = W_top[t,:,:-1] # Atmospheric storage of the cell to the west [m3]
        # to make dependent on isglobal but for now kept to avoid division by zero errors      
        Sa_W_top[0,:,0] = W_top[t,:,-1] # Atmospheric storage of the cell to the west [m3]
        Sa_N_top[0,1:,:] = W_top[t,:-1,:] # Atmospheric storage of the cell to the north [m3]
        Sa_S_top[0,:-1,:] = W_top[t,1:,:] # Atmospheric storage of the cell to the south [m3]
        
         # down: define values of tracked moisture of neighbouring grid cells
        Sa_track_E_down[0,:,:-1] = Sa_track_down[t,:,1:] # Atmospheric tracked storage of the cell to the east [m3]
        if isglobal == 1:
            Sa_track_E_down[0,:,-1] = Sa_track_down[t,:,0] #Atmospheric tracked storage of the cell to the east [m3]
        Sa_track_W_down[0,:,1:] = Sa_track_down[t,:,:-1] # Atmospheric storage of the cell to the west [m3]
        if isglobal == 1:        
            Sa_track_W_down[0,:,0] = Sa_track_down[t,:,-1] # Atmospheric storage of the cell to the west [m3]
        Sa_track_N_down[0,1:,:] = Sa_track_down[t,:-1,:] # Atmospheric storage of the cell to the north [m3]
        Sa_track_S_down[0,:-1,:] = Sa_track_down[t,1:,:] # Atmospheric storage of the cell to the south [m3]
        
        # down: calculate with moisture fluxes
        Sa_track_after_Fa_down[0,1:-1,:] = (Sa_track_down[t,1:-1,:] 
        + Fa_E_down_WE[t-1,1:-1,:] * (Sa_track_E_down[0,1:-1,:] / Sa_E_down[0,1:-1,:]) 
        - Fa_E_down_EW[t-1,1:-1,:] * (Sa_track_down[t,1:-1,:] / W_down[t,1:-1,:]) 
        - Fa_W_down_WE[t-1,1:-1,:] * (Sa_track_down[t,1:-1,:] / W_down[t,1:-1,:]) 
        + Fa_W_down_EW[t-1,1:-1,:] * (Sa_track_W_down[0,1:-1,:] / Sa_W_down[0,1:-1,:]) 
        + Fa_N_down_SN[t-1,1:-1,:] * (Sa_track_N_down[0,1:-1,:] / Sa_N_down[0,1:-1,:]) 
        - Fa_N_down_NS[t-1,1:-1,:] * (Sa_track_down[t,1:-1,:] / W_down[t,1:-1,:]) 
        - Fa_S_down_SN[t-1,1:-1,:] * (Sa_track_down[t,1:-1,:] / W_down[t,1:-1,:]) 
        + Fa_S_down_NS[t-1,1:-1,:] * (Sa_track_S_down[0,1:-1,:] / Sa_S_down[0,1:-1,:])
        - Fa_downward[t-1,1:-1,:] * (Sa_track_down[t,1:-1,:] / W_down[t,1:-1,:]) 
        + Fa_upward[t-1,1:-1,:] * (Sa_track_top[t,1:-1,:] / W_top[t,1:-1,:]))
        
        # top: define values of tracked moisture of neighbouring grid cells
        Sa_track_E_top[0,:,:-1] = Sa_track_top[t,:,1:] # Atmospheric tracked storage of the cell to the east [m3]
        if isglobal == 1:
            Sa_track_E_top[0,:,-1] = Sa_track_top[t,:,0] # Atmospheric tracked storage of the cell to the east [m3]
        Sa_track_W_top[0,:,1:] = Sa_track_top[t,:,:-1] # Atmospheric tracked storage of the cell to the west [m3]
        if isglobal == 1:
            Sa_track_W_top[0,:,0] = Sa_track_top[t,:,-1] # Atmospheric tracked storage of the cell to the west [m3]
        Sa_track_N_top[0,1:,:] = Sa_track_top[t,:-1,:] # Atmospheric tracked storage of the cell to the north [m3]
        Sa_track_S_top[0,:-1,:] = Sa_track_top[t,1:,:] # Atmospheric tracked storage of the cell to the south [m3]
        
        # top: calculate with moisture fluxes 
        Sa_track_after_Fa_top[0,1:-1,:] = (Sa_track_top[t,1:-1,:]
        + Fa_E_top_WE[t-1,1:-1,:] * (Sa_track_E_top[0,1:-1,:] / Sa_E_top[0,1:-1,:]) 
        - Fa_E_top_EW[t-1,1:-1,:] * (Sa_track_top[t,1:-1,:] / W_top[t,1:-1,:]) 
        - Fa_W_top_WE[t-1,1:-1,:] * (Sa_track_top[t,1:-1,:] / W_top[t,1:-1,:]) 
        + Fa_W_top_EW[t-1,1:-1,:] * (Sa_track_W_top[0,1:-1,:] / Sa_W_top[0,1:-1,:])
        + Fa_N_top_SN[t-1,1:-1,:] * (Sa_track_N_top[0,1:-1,:] / Sa_N_top[0,1:-1,:]) 
        - Fa_N_top_NS[t-1,1:-1,:] * (Sa_track_top[t,1:-1,:]/ W_top[t,1:-1,:]) 
        - Fa_S_top_SN[t-1,1:-1,:] * (Sa_track_top[t,1:-1,:] / W_top[t,1:-1,:]) 
        + Fa_S_top_NS[t-1,1:-1,:] * (Sa_track_S_top[0,1:-1,:] / Sa_S_top[0,1:-1,:]) 
        + Fa_downward[t-1,1:-1,:] * (Sa_track_down[t,1:-1,:] / W_down[t,1:-1,:]) 
        - Fa_upward[t-1,1:-1,:] * (Sa_track_top[t,1:-1,:] / W_top[t,1:-1,:]))
        
        # losses to the north and south
        north_loss[t-1,0,:] = (Fa_N_top_NS[t-1,1,:] * (Sa_track_top[t,1,:] / W_top[t,1,:])
                                + Fa_N_down_NS[t-1,1,:] * (Sa_track_down[t,1,:] / W_down[t,1,:]))
        south_loss[t-1,0,:] = (Fa_S_top_SN[t-1,-2,:] * (Sa_track_top[t,-2,:] / W_top[t,-2,:])
                                + Fa_S_down_SN[t-1,-2,:] * (Sa_track_down[t,-2,:] / W_down[t,-2,:]))
    
        # down: add precipitation and subtract evaporation
        Sa_track_after_Fa_P_E_down[0,1:-1,:] = (Sa_track_after_Fa_down[0,1:-1,:]
                                                + P_region[t-1,1:-1,:] * (W_down[t,1:-1,:] / W[t,1:-1,:])   
                                                - E[t-1,1:-1,:] * (Sa_track_down[t,1:-1,:] / W_down[t,1:-1,:]))
    
        # top: add precipitation
        Sa_track_after_Fa_P_E_top[0,1:-1,:] = (Sa_track_after_Fa_top[0,1:-1,:] 
                                                + P_region[t-1,1:-1,:] * (W_top[t,1:-1,:] / W[t,1:-1,:])) 
        
        # down and top: redistribute unaccounted water that is otherwise lost from the sytem
        down_to_top[t-1,:,:] = np.reshape(np.maximum(0, np.reshape(Sa_track_after_Fa_P_E_down, (np.size(Sa_track_after_Fa_P_E_down))) - np.reshape(W_down[t-1,:,:],
                                            (np.size(W_down[t-1,:,:])))), (len(latitude),len(longitude)))
        top_to_down[t-1,:,:] = np.reshape(np.maximum(0, np.reshape(Sa_track_after_Fa_P_E_top, (np.size(Sa_track_after_Fa_P_E_top))) - np.reshape(W_top[t-1,:,:],
                                            (np.size(W_top[t-1,:,:])))), (len(latitude),len(longitude)))
        Sa_track_after_all_down = Sa_track_after_Fa_P_E_down - down_to_top[t-1,:,:] + top_to_down[t-1,:,:]
        Sa_track_after_all_top = Sa_track_after_Fa_P_E_top - top_to_down[t-1,:,:] + down_to_top[t-1,:,:]

        # down and top: water lost to the system:
        water_lost_down[t-1,:,:] = np.reshape(np.maximum(0, np.reshape(Sa_track_after_all_down, (np.size(Sa_track_after_all_down))) - np.reshape(W_down[t-1,:,:],
                                            (np.size(W_down[t-1,:,:])))), (len(latitude),len(longitude)))
        water_lost_top[t-1,:,:] = np.reshape(np.maximum(0, np.reshape(Sa_track_after_all_top, (np.size(Sa_track_after_all_top))) - np.reshape(W_top[t-1,:,:],
                                            (np.size(W_top[t-1,:,:])))), (len(latitude),len(longitude)))
        water_lost = water_lost_down + water_lost_top
    
        # down: determine Sa_region of this next timestep 100% stable
        Sa_track_down[t-1,1:-1,:] = np.reshape(np.maximum(0,np.minimum(np.reshape(W_down[t-1,1:-1,:], np.size(W_down[t-1,1:-1,:])), np.reshape(Sa_track_after_all_down[0,1:-1,:],
                                                np.size(Sa_track_after_all_down[0,1:-1,:])))), (len(latitude[1:-1]),len(longitude)))
        # top: determine Sa_region of this next timestep 100% stable
        Sa_track_top[t-1,1:-1,:] = np.reshape(np.maximum(0,np.minimum(np.reshape(W_top[t-1,1:-1,:], np.size(W_top[t-1,1:-1,:])), np.reshape(Sa_track_after_all_top[0,1:-1,:],
                                                np.size(Sa_track_after_all_top[0,1:-1,:])))), (len(latitude[1:-1]),len(longitude)))
        
        #############################################################
        # timetracking start
        
        # defining sizes of timed moisture
        Sa_time_after_Fa_down = np.zeros(np.shape(Sa_time_down_last))
        Sa_time_after_Fa_P_E_down = np.zeros(np.shape(Sa_time_down_last))
        Sa_time_E_down = np.zeros(np.shape(Sa_time_down_last))
        Sa_time_W_down = np.zeros(np.shape(Sa_time_down_last))
        Sa_time_N_down = np.zeros(np.shape(Sa_time_down_last))
        Sa_time_S_down = np.zeros(np.shape(Sa_time_down_last))
        Sa_time_after_Fa_top = np.zeros(np.shape(Sa_time_top_last))
        Sa_time_after_Fa_P_E_top = np.zeros(np.shape(Sa_time_top_last))
        Sa_time_E_top = np.zeros(np.shape(Sa_time_top_last))
        Sa_time_W_top = np.zeros(np.shape(Sa_time_top_last))
        Sa_time_N_top = np.zeros(np.shape(Sa_time_top_last))
        Sa_time_S_top = np.zeros(np.shape(Sa_time_top_last))

        # time increase
        ti = timestep/divt

        # down: define values of timeed moisture of neighbouring grid cells
        Sa_time_E_down[0,:,:-1] = Sa_time_down[t,:,1:] # Atmospheric timeed storage of the cell to the east [s]
        if isglobal == 1:
            Sa_time_E_down[0,:,-1] = Sa_time_down[t,:,0] # Atmospheric timeed storage of the cell to the east [s]
        Sa_time_W_down[0,:,1:] = Sa_time_down[t,:,:-1] # Atmospheric timeed storage of the cell to the west [s]
        if isglobal == 1:        
            Sa_time_W_down[0,:,0] = Sa_time_down[t,:,-1] # Atmospheric timeed storage of the cell to the west [s]
        Sa_time_N_down[0,1:,:] = Sa_time_down[t,:-1,:] # Atmospheric timeed storage of the cell to the north [s]
        Sa_time_S_down[0,:-1,:] = Sa_time_down[t,1:,:] # Atmospheric timeed storage of the cell to the south [s]

        # down: calculate with moisture fluxes
        Sa_time_after_Fa_down[0,1:-1,:] = ((Sa_track_down[t,1:-1,:] * (ti + Sa_time_down[t,1:-1,:]) 
        + Fa_E_down_WE[t-1,1:-1,:] * (ti + Sa_time_E_down[0,1:-1,:]) * (Sa_track_E_down[0,1:-1,:] / Sa_E_down[0,1:-1,:]) 
        - Fa_E_down_EW[t-1,1:-1,:] * (ti + Sa_time_down[t,1:-1,:]) * (Sa_track_down[t,1:-1,:] / W_down[t,1:-1,:]) 
        - Fa_W_down_WE[t-1,1:-1,:] * (ti + Sa_time_down[t,1:-1,:]) * (Sa_track_down[t,1:-1,:] / W_down[t,1:-1,:]) 
        + Fa_W_down_EW[t-1,1:-1,:] * (ti + Sa_time_W_down[0,1:-1,:]) * (Sa_track_W_down[0,1:-1,:] / Sa_W_down[0,1:-1,:]) 
        + Fa_N_down_SN[t-1,1:-1,:] * (ti + Sa_time_N_down[0,1:-1,:]) * (Sa_track_N_down[0,1:-1,:] / Sa_N_down[0,1:-1,:]) 
        - Fa_N_down_NS[t-1,1:-1,:] * (ti + Sa_time_down[t,1:-1,:]) * (Sa_track_down[t,1:-1,:] / W_down[t,1:-1,:]) 
        - Fa_S_down_SN[t-1,1:-1,:] * (ti + Sa_time_down[t,1:-1,:]) * (Sa_track_down[t,1:-1,:] / W_down[t,1:-1,:]) 
        + Fa_S_down_NS[t-1,1:-1,:] * (ti + Sa_time_S_down[0,1:-1,:]) * (Sa_track_S_down[0,1:-1,:] / Sa_S_down[0,1:-1,:]) 
        - Fa_downward[t-1,1:-1,:] * (ti + Sa_time_down[t,1:-1,:]) * (Sa_track_down[t,1:-1,:] / W_down[t,1:-1,:]) 
        + Fa_upward[t-1,1:-1,:] * (ti + Sa_time_top[t,1:-1,:]) * (Sa_track_top[t,1:-1,:] / W_top[t,1:-1,:])
        ) / Sa_track_after_Fa_down[0,1:-1,:])

        where_are_NaNs = np.isnan(Sa_time_after_Fa_down)
        Sa_time_after_Fa_down[where_are_NaNs] = 0 

        # top: define values of timeed moisture of neighbouring grid cells
        Sa_time_E_top[0,:,:-1] = Sa_time_top[t,:,1:] # Atmospheric storage of the cell to the east [m3]
        if isglobal == 1:
            Sa_time_E_top[0,:,-1] = Sa_time_top[t,:,0] # Atmospheric storage of the cell to the east [m3]
        Sa_time_W_top[0,:,1:] = Sa_time_top[t,:,:-1] # Atmospheric storage of the cell to the west [m3]
        if isglobal == 1:
            Sa_time_W_top[0,:,0] = Sa_time_top[t,:,-1] # Atmospheric storage of the cell to the west [m3]
        Sa_time_N_top[0,1:,:] = Sa_time_top[t,:-1,:] # Atmospheric storage of the cell to the north [m3]
        Sa_time_S_top[0,:-1,:] = Sa_time_top[t,1:,:] # Atmospheric storage of the cell to the south [m3]

        # top: calculate with moisture fluxes 
        Sa_time_after_Fa_top[0,1:-1,:] = ((Sa_track_top[t,1:-1,:] * (ti + Sa_time_top[t,1:-1,:]) 
        + Fa_E_top_WE[t-1,1:-1,:] * (ti + Sa_time_E_top[0,1:-1,:]) * (Sa_track_E_top[0,1:-1,:] / Sa_E_top[0,1:-1,:]) 
        - Fa_E_top_EW[t-1,1:-1,:] * (ti + Sa_time_top[t,1:-1,:]) * (Sa_track_top[t,1:-1,:] / W_top[t,1:-1,:])  
        - Fa_W_top_WE[t-1,1:-1,:] * (ti + Sa_time_top[t,1:-1,:]) * (Sa_track_top[t,1:-1,:] / W_top[t,1:-1,:]) 
        + Fa_W_top_EW[t-1,1:-1,:] * (ti + Sa_time_W_top[0,1:-1,:]) * (Sa_track_W_top[0,1:-1,:] / Sa_W_top[0,1:-1,:])
        + Fa_N_top_SN[t-1,1:-1,:] * (ti + Sa_time_N_top[0,1:-1,:]) * (Sa_track_N_top[0,1:-1,:] / Sa_N_top[0,1:-1,:]) 
        - Fa_N_top_NS[t-1,1:-1,:] * (ti + Sa_time_top[t,1:-1,:]) * (Sa_track_top[t,1:-1,:] / W_top[t,1:-1,:])
        - Fa_S_top_SN[t-1,1:-1,:] * (ti + Sa_time_top[t,1:-1,:]) * (Sa_track_top[t,1:-1,:] / W_top[t,1:-1,:]) 
        + Fa_S_top_NS[t-1,1:-1,:] * (ti + Sa_time_S_top[0,1:-1,:]) * (Sa_track_S_top[0,1:-1,:] / Sa_S_top[0,1:-1,:]) 
        + Fa_downward[t-1,1:-1,:] * (ti + Sa_time_down[t,1:-1,:]) * (Sa_track_down[t,1:-1,:] / W_down[t,1:-1,:]) 
        - Fa_upward[t-1,1:-1,:] * (ti + Sa_time_top[t,1:-1,:]) * (Sa_track_top[t,1:-1,:] / W_top[t,1:-1,:]) 
        ) / Sa_track_after_Fa_top[0,1:-1,:])

        where_are_NaNs = np.isnan(Sa_time_after_Fa_top)
        Sa_time_after_Fa_top[where_are_NaNs] = 0

        # down: add precipitation and substract evaporation
        Sa_time_after_Fa_P_E_down[0,1:-1,:] = ((Sa_track_after_Fa_down[0,1:-1,:] * Sa_time_after_Fa_down[0,1:-1,:] 
        + P_region[t-1,1:-1,:] * ti/2. * (W_down[t,1:-1,:] / W[t,1:-1,:]) 
        - E[t-1,1:-1,:] * (ti + Sa_time_down[t,1:-1,:]) * (Sa_track_down[t,1:-1,:] / W_down[t,1:-1,:]) 
        ) / Sa_track_after_Fa_P_E_down[0,1:-1,:])

        where_are_NaNs = np.isnan(Sa_time_after_Fa_P_E_down)
        Sa_time_after_Fa_P_E_down[where_are_NaNs] = 0

        # top: add precipitation
        Sa_time_after_Fa_P_E_top[0,1:-1,:] = ((Sa_track_after_Fa_top[0,1:-1,:] * Sa_time_after_Fa_top[0,1:-1,:]
       + P_region[t-1,1:-1,:] * ti/2 * (W_top[t,1:-1,:] / W[t,1:-1,:])
        ) / Sa_track_after_Fa_P_E_top[0,1:-1,:])

        where_are_NaNs = np.isnan(Sa_time_after_Fa_P_E_top)
        Sa_time_after_Fa_P_E_top[where_are_NaNs] = 0

        # down: redistribute water
        Sa_time_after_all_down = ((Sa_track_after_Fa_P_E_down * Sa_time_after_Fa_P_E_down 
        - down_to_top[t-1,:,:] * Sa_time_after_Fa_P_E_down
        + top_to_down[t-1,:,:] * Sa_time_after_Fa_P_E_top
        ) / Sa_track_after_all_down)

        where_are_NaNs = np.isnan(Sa_time_after_all_down)
        Sa_time_after_all_down[where_are_NaNs] = 0

        # top: redistribute water
        Sa_time_after_all_top = ((Sa_track_after_Fa_P_E_top * Sa_time_after_Fa_P_E_top
        - top_to_down[t-1,:,:] * Sa_time_after_Fa_P_E_top 
        + down_to_top[t-1,:,:] * Sa_time_after_Fa_P_E_down
        ) / Sa_track_after_all_top)

        where_are_NaNs = np.isnan(Sa_time_after_all_top)
        Sa_time_after_all_top[where_are_NaNs] = 0

        # down: determine Sa_region of this next timestep 100% stable
        Sa_time_down[t-1,1:-1,:] = Sa_time_after_all_down[0,1:-1,:]

        # top: determine Sa_region of this next timestep 100% stable
        Sa_time_top[t-1,1:-1,:] = Sa_time_after_all_top[0,1:-1,:]
        #############################################################
                                                        
    return Sa_time_top,Sa_time_down,Sa_track_top,Sa_track_down,north_loss,south_loss,down_to_top,top_to_down,water_lost

#%% create empty array for track and time

def create_empty_array(count_time,divt,latitude,longitude,yearpart,years):
    Sa_time_top = np.zeros((np.int(count_time*divt)+1,len(latitude),len(longitude)))
    Sa_time_down = np.zeros((np.int(count_time*divt)+1,len(latitude),len(longitude)))
    Sa_track_top = np.zeros((np.int(count_time*divt)+1,len(latitude),len(longitude)))
    Sa_track_down = np.zeros((np.int(count_time*divt)+1,len(latitude),len(longitude)))
    if yearpart[0] == 365:
        sio.savemat(datapathea[0],{'Sa_track_top':Sa_track_top,'Sa_track_down':Sa_track_down,},do_compression=True)
        sio.savemat(datapathea[1],{'Sa_time_top':Sa_time_top,'Sa_time_down':Sa_time_down},do_compression=True)
    else:
        sio.savemat(datapathea[2], {'Sa_track_top':Sa_track_top,'Sa_track_down':Sa_track_down},do_compression=True)
        sio.savemat(datapathea[3], {'Sa_time_top':Sa_time_top,'Sa_time_down':Sa_time_down},do_compression=True)
    
    return

#%% Runtime & Results

start1 = timer()

# The two lines below create empty arrays for first runs/initial values are zero. 
datapathea = data_path_ea(years,yearpart) #define paths for empty arrays
if veryfirstrun == 1:
    create_empty_array(count_time,divt,latitude,longitude,yearpart,years) #creates empty arrays for first day run

# loop through the years
for i in range(len(years)):
    yearnumber = years[i]
    ly = int(calendar.isleap(yearnumber))
    
    if (yearpart[0] == 365) & (calendar.isleap(yearnumber) == 0):
        thisyearpart = yearpart[1:]
    else: # a leapyear
        thisyearpart = yearpart
        
    for j in range(len(thisyearpart)):
        start = timer()
        a = thisyearpart[j]
        if a == (364 + calendar.isleap(yearnumber)): # a == 31 December
            previous_data_to_load = (str(yearnumber+1) + '-0')
        else: # a != 31 December
            previous_data_to_load = (str(yearnumber) + '-' + str(a+1))
            
        datapath = data_path(previous_data_to_load,yearnumber,a)
        
        loading_ST = sio.loadmat(datapath[0],verify_compressed_data_integrity=False)
        Sa_track_top = loading_ST['Sa_track_top']
        Sa_track_down = loading_ST['Sa_track_down']
        Sa_track_top_last_scheef = Sa_track_top[0,:,:]
        Sa_track_down_last_scheef = Sa_track_down[0,:,:]
        Sa_track_top_last =  np.reshape(Sa_track_top_last_scheef, (1,len(latitude),len(longitude)))
        Sa_track_down_last =  np.reshape(Sa_track_down_last_scheef, (1,len(latitude),len(longitude)))
        
        loading_FS = sio.loadmat(datapath[1],verify_compressed_data_integrity=False)
        Fa_E_top = loading_FS['Fa_E_top']
        Fa_N_top = loading_FS['Fa_N_top']
        Fa_E_down = loading_FS['Fa_E_down']
        Fa_N_down = loading_FS['Fa_N_down']
        Fa_Vert = loading_FS['Fa_Vert']
        E = loading_FS['E']
        P = loading_FS['P']
        W_top = loading_FS['W_top']
        W_down = loading_FS['W_down']
        
        # call the backward tracking function
        if timetracking == 0:
            Sa_track_top,Sa_track_down,north_loss,south_loss,down_to_top,top_to_down,water_lost = get_Sa_track_backward(latitude,longitude,count_time,divt,Kvf,Region,Fa_E_top,Fa_N_top,Fa_E_down,Fa_N_down,
                                       Fa_Vert,E,P,W_top,W_down,Sa_track_top_last,Sa_track_down_last)
        elif timetracking == 1:
            loading_STT = sio.loadmat(datapath[2],verify_compressed_data_integrity=False)
            Sa_time_top = loading_STT['Sa_time_top'] # [seconds]
            Sa_time_down = loading_STT['Sa_time_down']            
            Sa_time_top_last_scheef = Sa_time_top[0,:,:]
            Sa_time_down_last_scheef = Sa_time_down[0,:,:]
            Sa_time_top_last =  np.reshape(Sa_time_top_last_scheef, (1,len(latitude),len(longitude)))
            Sa_time_down_last =  np.reshape(Sa_time_down_last_scheef, (1,len(latitude),len(longitude)))
            
            Sa_time_top,Sa_time_down,Sa_track_top,Sa_track_down,north_loss,south_loss,down_to_top,top_to_down,water_lost = get_Sa_track_backward_TIME(latitude,longitude,count_time,divt,
                                            timestep,Kvf,Region,Fa_E_top,Fa_N_top,Fa_E_down,Fa_N_down,Fa_Vert,E,P,W_top,W_down,Sa_track_top_last,Sa_track_down_last,Sa_time_top_last,Sa_time_down_last)
        # save this data 
        sio.savemat(datapath[3], {'Sa_track_top':Sa_track_top,'Sa_track_down':Sa_track_down,'north_loss':north_loss, 'south_loss':south_loss,'down_to_top':down_to_top,'top_to_down':top_to_down,'water_lost':water_lost},do_compression=True)
        if timetracking == 1:
            sio.savemat(datapath[4], {'Sa_time_top':Sa_time_top,'Sa_time_down':Sa_time_down},do_compression=True)
        
        end = timer()       
        print 'Runtime Sa_track for day ' + str(a+1) + ' in year ' + str(yearnumber) + ' is',(end - start),' seconds.'

end1 = timer()
print 'The total runtime of Con_E_Recyc_Masterscript is',(end1-start1),' seconds.'


