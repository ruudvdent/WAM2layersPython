#!/home/users/lguo/anaconda2/bin/python
#BSUB -o %J.o
#BSUB -e %J.e
#BSUB -q short-serial
#BSUB -W 24:00
#BSUB -R "rusage[mem=16000]"
#BSUB -M 16000

import numpy as np
import scipy.io as sio
import calendar
from timeit import default_timer as timer
import os
from netCDF4 import Dataset
import datetime
import cf
import pdb

#%% BEGIN OF INPUT1 (FILL THIS IN)
years = np.arange(1991,1992) #fill in the years
yearpart = np.arange(0,366) # for a full (leap)year fill in np.arange(0,366)
boundary = 8 # with 8 the vertical separation is at 812.83 hPa for surface pressure = 1031.25 hPa, which corresponds to k=47 (ERA-Interim)
divt = 24 # division of the timestep, 24 means a calculation timestep of 6/24 = 0.25 hours (numerical stability purposes)
count_time = 4 # number of indices to get data from (for six hourly data this means everytime one day)

# Manage the extent of your dataset (FILL THIS IN)
# Define the latitude and longitude cell numbers to consider and corresponding lakes that should be considered part of the land
latnrs = np.arange(7,114)
lonnrs = np.arange(0,240)
isglobal = 1 # fill in 1 for global computations (i.e. Earth round), fill in 0 for a local domain with boundaries

# obtain the constants
invariant_data = '/home/users/lguo/perchance/WAM_input/erai/invariants.nc' #invariants
#latitude,longitude,lsm,g,density_water,timestep,A_gridcell,L_N_gridcell,L_S_gridcell,L_EW_gridcell,gridcell = getconstants(latnrs,lonnrs,lake_mask,invariant_data)
# 
latitude = Dataset(invariant_data,mode='r').variables['latitude'][latnrs]
longitude = Dataset(invariant_data,mode='r').variables['longitude'][lonnrs]
lsm = np.squeeze(Dataset(invariant_data,mode='r').variables['lsm'][0,latnrs,lonnrs]) # 0 = sea, 1 = land
lsm[0,:] = 0 # the northern boundary is always oceanic = 0
lsm[-1,:] = 0 # the southern boundary is always oceanic = 0
# Constants 
g = 9.80665 # [m/s2] from ERA-interim archive
density_water = 1000 # [kg/m3]
dg = 111089.56 # [m] length of 1 degree latitude
timestep = 6*3600 # [s] timestep in the ERA-interim archive (watch out! P & E have 3 hour timestep)
Erad = 6.371e6 # [m] Earth radius
# Semiconstants
gridcell = np.abs(longitude[1]-longitude[0]) # [degrees] grid cell size, it is 1.5 degree
lat_n_bound = np.minimum(90.0,latitude+0.5*gridcell)	# Finding north and south boundaries of each gridcell
lat_s_bound = np.maximum(-90.0,latitude-0.5*gridcell)
A_gridcell = np.zeros([len(latitude),1])
A_gridcell[:,0] = (np.pi/180.0)*Erad**2*abs(np.sin(lat_s_bound*np.pi/180.0)-np.sin(lat_n_bound*np.pi/180.0))*gridcell
L_N_gridcell = gridcell*np.cos((latitude+gridcell/2.0)*np.pi/180.0)*dg # [m] length northern boundary of a cell, [107]
L_S_gridcell = gridcell*np.cos((latitude-gridcell/2.0)*np.pi/180.0)*dg # [m] length southern boundary of a cell, [107]
L_EW_gridcell = gridcell*dg # [m] length eastern/western boundary of a cell, [1]
#
source_region = Dataset('/home/users/lguo/perchance/WAM_input/erai_filter/source_cn1_7915_montly.nc',mode='r').variables['sr']
Kvf = 3 # vertical dispersion factor (advection only is 0, dispersion the same size of the advective flux is 1, for stability don't make this more than 3)
timetracking = 1 # 0 for not tracking time and 1 for tracking time
veryfirstrun = 1 # type '1' if no run has been done before from which can be continued, otherwise type '0'
#
interdata_folder = '/home/users/lguo/perchance/WAM_inter/interdata_erai' # must be an existing folder, existence is not checked
#END OF INPUT

# Check if interdata folder exists:
assert os.path.isdir(interdata_folder), "Please create the interdata_folder before running the script"
# Check if sub interdata folder exists otherwise create it:
sub_interdata_folder = os.path.join(interdata_folder,'cn1_source_forward')
if os.path.isdir(sub_interdata_folder):
    pass
else:
    os.makedirs(sub_interdata_folder)

def data_path_ea(years,yearpart):
    save_empty_arrays_ly_track = os.path.join(sub_interdata_folder,str(years[0]-1)+'-'+str(364+calendar.isleap(years[0]-1)).zfill(3)+'Sa_track.nc')
    save_empty_arrays_track = os.path.join(sub_interdata_folder,str(years[0])+'-'+str(yearpart[0]-1).zfill(3)+'Sa_track.nc')
    return save_empty_arrays_ly_track,save_empty_arrays_track

def data_path(previous_data_to_load,yearnumber,a):
    load_Sa_track = os.path.join(sub_interdata_folder,previous_data_to_load+'Sa_track.nc')
    load_fluxes_and_storages = os.path.join(interdata_folder,str(yearnumber)+'-'+str(a).zfill(3)+'fluxes_storages.nc')
    save_path_track = os.path.join(sub_interdata_folder,str(yearnumber)+'-'+str(a).zfill(3)+'Sa_track.nc')
    return load_Sa_track,load_fluxes_and_storages,save_path_track

#%% Code (no need to look at this for running)
def get_Sa_track_forward(latitude,longitude,count_time,divt,Kvf,Region,Fa_E_top,Fa_N_top,Fa_E_down,Fa_N_down,Fa_Vert,E,P,W_top,W_down,Sa_track_top_last,Sa_track_down_last):
    # make E_region matrix
    Region3D = np.tile(np.matlib.reshape(Region,[1,len(latitude),len(longitude)]),[len(P[:,0,0]),1,1])
    E_region = Region3D * E # Ocean is 0, Continent is 1.

    # Total moisture in the column
    W = W_top + W_down

    # separate the direction of the vertical flux and make it absolute (downward movement is positive)
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
    
    # assign begin values of output == last values of the previous time slot
    Sa_track_down[0,:,:] = Sa_track_down_last
    Sa_track_top[0,:,:] = Sa_track_top_last

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
    
    # Sa calculation forward in time
    for t in range(np.int(count_time*divt)):
        # down: define values of total moisture
        Sa_E_down[0,:,:-1] = W_down[t,:,1:] # Atmospheric storage of the cell to the east [m3]
        # to make dependent on isglobal but for now kept to avoid division by zero errors      
        Sa_E_down[0,:,-1] = W_down[t,:,0] # Atmospheric storage of the cell to the east [m3]
        Sa_W_down[0,:,1:] = W_down[t,:,:-1] # Atmospheric storage of the cell to the west [m3]
        # to make dependent on isglobal but for now kept to avoid division by zero errors        
        Sa_W_down[0,:,0] = W_down[t,:,-1] # Atmospheric storage of the cell to the west [m3]
        Sa_N_down[0,1:,:] = W_down[t,:-1,:] # Atmospheric storage of the cell to the north [m3]
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
            Sa_track_E_down[0,:,-1] = Sa_track_down[t,:,0] # Atmospheric tracked storage of the cell to the east [m3]
        Sa_track_W_down[0,:,1:] = Sa_track_down[t,:,:-1] # Atmospheric tracked storage of the cell to the west [m3]
        if isglobal == 1:        
            Sa_track_W_down[0,:,0] = Sa_track_down[t,:,-1] # Atmospheric tracked storage of the cell to the west [m3]
        Sa_track_N_down[0,1:,:] = Sa_track_down[t,:-1,:] # Atmospheric tracked storage of the cell to the north [m3]
        Sa_track_S_down[0,:-1,:] = Sa_track_down[t,1:,:] # Atmospheric tracked storage of the cell to the south [m3]

        # down: calculate with moisture fluxes
        Sa_track_after_Fa_down[0,1:-1,:] = (Sa_track_down[t,1:-1,:] 
        - Fa_E_down_WE[t,1:-1,:] * (Sa_track_down[t,1:-1,:] / W_down[t,1:-1,:])
        + Fa_E_down_EW[t,1:-1,:] * (Sa_track_E_down[0,1:-1,:] / Sa_E_down[0,1:-1,:]) 
        + Fa_W_down_WE[t,1:-1,:] * (Sa_track_W_down[0,1:-1,:] / Sa_W_down[0,1:-1,:]) 
        - Fa_W_down_EW[t,1:-1,:] * (Sa_track_down[t,1:-1,:] / W_down[t,1:-1,:]) 
        - Fa_N_down_SN[t,1:-1,:] * (Sa_track_down[t,1:-1,:] / W_down[t,1:-1,:]) 
        + Fa_N_down_NS[t,1:-1,:] * (Sa_track_N_down[0,1:-1,:] / Sa_N_down[0,1:-1,:])
        + Fa_S_down_SN[t,1:-1,:] * (Sa_track_S_down[0,1:-1,:] / Sa_S_down[0,1:-1,:]) 
        - Fa_S_down_NS[t,1:-1,:] * (Sa_track_down[t,1:-1,:] / W_down[t,1:-1,:]) 
        + Fa_downward[t,1:-1,:] * (Sa_track_top[t,1:-1,:] / W_top[t,1:-1,:])
        - Fa_upward[t,1:-1,:] * (Sa_track_down[t,1:-1,:] / W_down[t,1:-1,:]))

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
        - Fa_E_top_WE[t,1:-1,:] * (Sa_track_top[t,1:-1,:] / W_top[t,1:-1,:])  
        + Fa_E_top_EW[t,1:-1,:] * (Sa_track_E_top[0,1:-1,:] / Sa_E_top[0,1:-1,:]) 
        + Fa_W_top_WE[t,1:-1,:] * (Sa_track_W_top[0,1:-1,:] / Sa_W_top[0,1:-1,:]) 
        - Fa_W_top_EW[t,1:-1,:] * (Sa_track_top[t,1:-1,:] / W_top[t,1:-1,:]) 
        - Fa_N_top_SN[t,1:-1,:] * (Sa_track_top[t,1:-1,:] / W_top[t,1:-1,:]) 
        + Fa_N_top_NS[t,1:-1,:] * (Sa_track_N_top[0,1:-1,:] / Sa_N_top[0,1:-1,:]) 
        + Fa_S_top_SN[t,1:-1,:] * (Sa_track_S_top[0,1:-1,:] / Sa_S_top[0,1:-1,:]) 
        - Fa_S_top_NS[t,1:-1,:] * (Sa_track_top[t,1:-1,:] / W_top[t,1:-1,:]) 
        - Fa_downward[t,1:-1,:] * (Sa_track_top[t,1:-1,:] / W_top[t,1:-1,:]) 
        + Fa_upward[t,1:-1,:] * (Sa_track_down[t,1:-1,:] / W_down[t,1:-1,:]))

        # losses to the north and south
        north_loss[t,0,:] = (Fa_N_top_SN[t,1,:] * (Sa_track_top[t,1,:] / W_top[t,1,:])
                            + Fa_N_down_SN[t,1,:] * (Sa_track_down[t,1,:] / W_down[t,1,:]))
        south_loss[t,0,:] = (Fa_S_top_NS[t,-2,:] * (Sa_track_top[t,-2,:] / W_top[t,-2,:])
                            + Fa_S_down_NS[t,-2,:] * (Sa_track_down[t,-2,:] / W_down[t,-2,:]))

        # down: substract precipitation and add evaporation
        Sa_track_after_Fa_P_E_down[0,1:-1,:] = (Sa_track_after_Fa_down[0,1:-1,:]
                                                     - P[t,1:-1,:] * (Sa_track_down[t,1:-1,:] / W[t,1:-1,:]) 
                                                    + E_region[t,1:-1,:])

        # top: substract precipitation
        Sa_track_after_Fa_P_E_top[0,1:-1,:] = (Sa_track_after_Fa_top[0,1:-1,:] 
                                                - P[t,1:-1,:] * (Sa_track_top[t,1:-1,:] / W[t,1:-1,:])) 
        
        # down and top: redistribute unaccounted water that is otherwise lost from the sytem
        down_to_top[t,:,:] = np.reshape(np.maximum(0, np.reshape(Sa_track_after_Fa_P_E_down, (np.size(Sa_track_after_Fa_P_E_down))) - np.reshape(W_down[t+1,:,:],
                                            (np.size(W_down[t+1,:,:])))), (len(latitude),len(longitude)))
        top_to_down[t,:,:] = np.reshape(np.maximum(0, np.reshape(Sa_track_after_Fa_P_E_top, (np.size(Sa_track_after_Fa_P_E_top))) - np.reshape(W_top[t+1,:,:],
                                            (np.size(W_top[t+1,:,:])))), (len(latitude),len(longitude)))
        Sa_track_after_all_down = Sa_track_after_Fa_P_E_down - down_to_top[t,:,:] + top_to_down[t,:,:]
        Sa_track_after_all_top = Sa_track_after_Fa_P_E_top - top_to_down[t,:,:] + down_to_top[t,:,:]

        # down and top: water lost to the system: 
        water_lost_down[t,:,:] = np.reshape(np.maximum(0, np.reshape(Sa_track_after_all_down, (np.size(Sa_track_after_all_down))) - np.reshape(W_down[t+1,:,:],
                                            (np.size(W_down[t+1,:,:])))), (len(latitude),len(longitude)))
        water_lost_top[t,:,:] = np.reshape(np.maximum(0, np.reshape(Sa_track_after_all_top, (np.size(Sa_track_after_all_top))) - np.reshape(W_top[t+1,:,:],
                                            (np.size(W_top[t+1,:,:])))), (len(latitude),len(longitude)))
        water_lost[t,:,:] = water_lost_down[t,:,:] + water_lost_top[t,:,:]

        # down: determine Sa_region of this next timestep 100% stable
	#LG: eliminate negative values?
        Sa_track_down[t+1,1:-1,:] = np.reshape(np.maximum(0,np.minimum(np.reshape(W_down[t+1,1:-1,:], np.size(W_down[t+1,1:-1,:])), np.reshape(Sa_track_after_all_down[0,1:-1,:],
                                                np.size(Sa_track_after_all_down[0,1:-1,:])))), (len(latitude[1:-1]),len(longitude)))
        # top: determine Sa_region of this next timestep 100% stable
        Sa_track_top[t+1,1:-1,:] = np.reshape(np.maximum(0,np.minimum(np.reshape(W_top[t+1,1:-1,:], np.size(W_top[t+1,1:-1,:])), np.reshape(Sa_track_after_all_top[0,1:-1,:],
                                                np.size(Sa_track_after_all_top[0,1:-1,:])))), (len(latitude[1:-1]),len(longitude)))
    
    return Sa_track_top,Sa_track_down,north_loss,south_loss,down_to_top,top_to_down,water_lost

def get_Sa_track_forward_TIME(latitude,longitude,count_time,divt,timestep,Kvf,Region,Fa_E_top,Fa_N_top,Fa_E_down,Fa_N_down,Fa_Vert,E,P,W_top,W_down,Sa_track_top_last,Sa_track_down_last,Sa_time_top_last,Sa_time_down_last,Sa_dist_top_last,Sa_dist_down_last,L_N_gridcell,L_S_gridcell,L_EW_gridcell):
    
    # make E_region matrix
    Region3D = np.tile(np.reshape(Region,[1,len(latitude),len(longitude)]),[len(P[:,0,0]),1,1])
    E_region = Region3D * E

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
    Sa_dist_down = np.zeros(np.shape(W_down))
    Sa_dist_top = np.zeros(np.shape(W_top))
    
    # assign begin values of output == last values of the previous time slot
    Sa_track_down[0,:,:] = Sa_track_down_last
    Sa_track_top[0,:,:] = Sa_track_top_last
    Sa_time_down[0,:,:] = Sa_time_down_last
    Sa_time_top[0,:,:] = Sa_time_top_last
    Sa_dist_down[0,:,:] = Sa_dist_down_last
    Sa_dist_top[0,:,:] = Sa_dist_top_last

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
    
    # Sa calculation forward in time
    for t in range(np.int(count_time*divt)):
        # down: define values of total moisture
        Sa_E_down[0,:,:-1] = W_down[t,:,1:] # Atmospheric storage of the cell to the east [m3]
        # to make dependent on isglobal but for now kept to avoid division by zero errors      
        Sa_E_down[0,:,-1] = W_down[t,:,0] # Atmospheric storage of the cell to the east [m3]
        Sa_W_down[0,:,1:] = W_down[t,:,:-1] # Atmospheric storage of the cell to the west [m3]
        # to make dependent on isglobal but for now kept to avoid division by zero errors      
        Sa_W_down[0,:,0] = W_down[t,:,-1] # Atmospheric storage of the cell to the west [m3]
        Sa_N_down[0,1:,:] = W_down[t,:-1,:] # Atmospheric storage of the cell to the north [m3]
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
            Sa_track_E_down[0,:,-1] = Sa_track_down[t,:,0] # Atmospheric tracked storage of the cell to the east [m3]
        Sa_track_W_down[0,:,1:] = Sa_track_down[t,:,:-1] # Atmospheric tracked storage of the cell to the west [m3]
        if isglobal == 1:
            Sa_track_W_down[0,:,0] = Sa_track_down[t,:,-1] # Atmospheric tracked storage of the cell to the west [m3]
        Sa_track_N_down[0,1:,:] = Sa_track_down[t,:-1,:] # Atmospheric tracked storage of the cell to the north [m3]
        Sa_track_S_down[0,:-1,:] = Sa_track_down[t,1:,:] # Atmospheric tracked storage of the cell to the south [m3]

        # down: calculate with moisture fluxes
        Sa_track_after_Fa_down[0,1:-1,:] = (Sa_track_down[t,1:-1,:] 
        - Fa_E_down_WE[t,1:-1,:] * (Sa_track_down[t,1:-1,:] / W_down[t,1:-1,:])
        + Fa_E_down_EW[t,1:-1,:] * (Sa_track_E_down[0,1:-1,:] / Sa_E_down[0,1:-1,:]) 
        + Fa_W_down_WE[t,1:-1,:] * (Sa_track_W_down[0,1:-1,:] / Sa_W_down[0,1:-1,:]) 
        - Fa_W_down_EW[t,1:-1,:] * (Sa_track_down[t,1:-1,:] / W_down[t,1:-1,:]) 
        - Fa_N_down_SN[t,1:-1,:] * (Sa_track_down[t,1:-1,:] / W_down[t,1:-1,:]) 
        + Fa_N_down_NS[t,1:-1,:] * (Sa_track_N_down[0,1:-1,:] / Sa_N_down[0,1:-1,:])
        + Fa_S_down_SN[t,1:-1,:] * (Sa_track_S_down[0,1:-1,:] / Sa_S_down[0,1:-1,:]) 
        - Fa_S_down_NS[t,1:-1,:] * (Sa_track_down[t,1:-1,:] / W_down[t,1:-1,:]) 
        + Fa_downward[t,1:-1,:] * (Sa_track_top[t,1:-1,:] / W_top[t,1:-1,:])
        - Fa_upward[t,1:-1,:] * (Sa_track_down[t,1:-1,:] / W_down[t,1:-1,:]))

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
        - Fa_E_top_WE[t,1:-1,:] * (Sa_track_top[t,1:-1,:] / W_top[t,1:-1,:])  
        + Fa_E_top_EW[t,1:-1,:] * (Sa_track_E_top[0,1:-1,:] / Sa_E_top[0,1:-1,:]) 
        + Fa_W_top_WE[t,1:-1,:] * (Sa_track_W_top[0,1:-1,:] / Sa_W_top[0,1:-1,:]) 
        - Fa_W_top_EW[t,1:-1,:] * (Sa_track_top[t,1:-1,:] / W_top[t,1:-1,:]) 
        - Fa_N_top_SN[t,1:-1,:] * (Sa_track_top[t,1:-1,:] / W_top[t,1:-1,:]) 
        + Fa_N_top_NS[t,1:-1,:] * (Sa_track_N_top[0,1:-1,:] / Sa_N_top[0,1:-1,:]) 
        + Fa_S_top_SN[t,1:-1,:] * (Sa_track_S_top[0,1:-1,:] / Sa_S_top[0,1:-1,:]) 
        - Fa_S_top_NS[t,1:-1,:] * (Sa_track_top[t,1:-1,:] / W_top[t,1:-1,:]) 
        - Fa_downward[t,1:-1,:] * (Sa_track_top[t,1:-1,:] / W_top[t,1:-1,:]) 
        + Fa_upward[t,1:-1,:] * (Sa_track_down[t,1:-1,:] / W_down[t,1:-1,:]))

        # losses to the north and south
        north_loss[t,0,:] = (Fa_N_top_SN[t,1,:] * (Sa_track_top[t,1,:] / W_top[t,1,:])
                            + Fa_N_down_SN[t,1,:] * (Sa_track_down[t,1,:] / W_down[t,1,:]))
        south_loss[t,0,:] = (Fa_S_top_NS[t,-2,:] * (Sa_track_top[t,-2,:] / W_top[t,-2,:])
                            + Fa_S_down_NS[t,-2,:] * (Sa_track_down[t,-2,:] / W_down[t,-2,:]))

        # down: substract precipitation and add evaporation
        Sa_track_after_Fa_P_E_down[0,1:-1,:] = (Sa_track_after_Fa_down[0,1:-1,:]
                                                - P[t,1:-1,:] * (Sa_track_down[t,1:-1,:] / W[t,1:-1,:]) 
                                                + E_region[t,1:-1,:])

        # top: substract precipitation
        Sa_track_after_Fa_P_E_top[0,1:-1,:] = (Sa_track_after_Fa_top[0,1:-1,:] 
                                                - P[t,1:-1,:] * (Sa_track_top[t,1:-1,:] / W[t,1:-1,:])) 
        
        # down and top: redistribute unaccounted water that is otherwise lost from the sytem
        down_to_top[t,:,:] = np.reshape(np.maximum(0, np.reshape(Sa_track_after_Fa_P_E_down, (np.size(Sa_track_after_Fa_P_E_down))) - np.reshape(W_down[t+1,:,:],
                                            (np.size(W_down[t+1,:,:])))), (len(latitude),len(longitude)))
        top_to_down[t,:,:] = np.reshape(np.maximum(0, np.reshape(Sa_track_after_Fa_P_E_top, (np.size(Sa_track_after_Fa_P_E_top))) - np.reshape(W_top[t+1,:,:],
                                            (np.size(W_top[t+1,:,:])))), (len(latitude),len(longitude)))
        Sa_track_after_all_down = Sa_track_after_Fa_P_E_down - down_to_top[t,:,:] + top_to_down[t,:,:]
        Sa_track_after_all_top = Sa_track_after_Fa_P_E_top - top_to_down[t,:,:] + down_to_top[t,:,:]

        # down and top: water lost to the system: 
        water_lost_down[t,:,:] = np.reshape(np.maximum(0, np.reshape(Sa_track_after_all_down, (np.size(Sa_track_after_all_down))) - np.reshape(W_down[t+1,:,:],
                                            (np.size(W_down[t+1,:,:])))), (len(latitude),len(longitude)))
        water_lost_top[t,:,:] = np.reshape(np.maximum(0, np.reshape(Sa_track_after_all_top, (np.size(Sa_track_after_all_top))) - np.reshape(W_top[t+1,:,:],
                                            (np.size(W_top[t+1,:,:])))), (len(latitude),len(longitude)))
        water_lost[t,:,:] = water_lost_down[t,:,:] + water_lost_top[t,:,:]

        # down: determine Sa_region of this next timestep 100% stable
        Sa_track_down[t+1,1:-1,:] = np.reshape(np.maximum(0,np.minimum(np.reshape(W_down[t+1,1:-1,:], np.size(W_down[t+1,1:-1,:])), np.reshape(Sa_track_after_all_down[0,1:-1,:],
                                                np.size(Sa_track_after_all_down[0,1:-1,:])))), (len(latitude[1:-1]),len(longitude)))
        # top: determine Sa_region of this next timestep 100% stable
        Sa_track_top[t+1,1:-1,:] = np.reshape(np.maximum(0,np.minimum(np.reshape(W_top[t+1,1:-1,:], np.size(W_top[t+1,1:-1,:])), np.reshape(Sa_track_after_all_top[0,1:-1,:],
                                                np.size(Sa_track_after_all_top[0,1:-1,:])))), (len(latitude[1:-1]),len(longitude)))
        
        #########################################
        #time tracking start

        
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
        - Fa_E_down_WE[t,1:-1,:] * (ti + Sa_time_down[t,1:-1,:]) * (Sa_track_down[t,1:-1,:] / W_down[t,1:-1,:]) 
        + Fa_E_down_EW[t,1:-1,:] * (ti + Sa_time_E_down[0,1:-1,:]) * (Sa_track_E_down[0,1:-1,:] / Sa_E_down[0,1:-1,:]) 
        + Fa_W_down_WE[t,1:-1,:] * (ti + Sa_time_W_down[0,1:-1,:]) * (Sa_track_W_down[0,1:-1,:] / Sa_W_down[0,1:-1,:]) 
        - Fa_W_down_EW[t,1:-1,:] * (ti + Sa_time_down[t,1:-1,:]) * (Sa_track_down[t,1:-1,:] / W_down[t,1:-1,:]) 
        - Fa_N_down_SN[t,1:-1,:] * (ti + Sa_time_down[t,1:-1,:]) * (Sa_track_down[t,1:-1,:] / W_down[t,1:-1,:]) 
        + Fa_N_down_NS[t,1:-1,:] * (ti + Sa_time_N_down[0,1:-1,:]) * (Sa_track_N_down[0,1:-1,:] / Sa_N_down[0,1:-1,:]) 
        + Fa_S_down_SN[t,1:-1,:] * (ti + Sa_time_S_down[0,1:-1,:]) * (Sa_track_S_down[0,1:-1,:] / Sa_S_down[0,1:-1,:]) 
        - Fa_S_down_NS[t,1:-1,:] * (ti + Sa_time_down[t,1:-1,:]) * (Sa_track_down[t,1:-1,:] / W_down[t,1:-1,:]) 
        + Fa_downward[t,1:-1,:] * (ti + Sa_time_top[t,1:-1,:]) * (Sa_track_top[t,1:-1,:] / W_top[t,1:-1,:]) 
        - Fa_upward[t,1:-1,:] * (ti + Sa_time_down[t,1:-1,:]) * (Sa_track_down[t,1:-1,:] / W_down[t,1:-1,:]) 
        ) / Sa_track_after_Fa_down[0,1:-1,:])

        where_are_NaNs = np.isnan(Sa_time_after_Fa_down)
        Sa_time_after_Fa_down[where_are_NaNs] = 0 

        # top: define values of total moisture
        Sa_time_E_top[0,:,:-1] = Sa_time_top[t,:,1:] # Atmospheric storage of the cell to the east [s]
        if isglobal == 1:        
            Sa_time_E_top[0,:,-1] = Sa_time_top[t,:,0] # Atmospheric storage of the cell to the east [s]
        Sa_time_W_top[0,:,1:] = Sa_time_top[t,:,:-1] # Atmospheric storage of the cell to the west [s]
        if isglobal == 1:
            Sa_time_W_top[0,:,0] = Sa_time_top[t,:,-1] # Atmospheric storage of the cell to the west [s]
        Sa_time_N_top[0,1:,:] = Sa_time_top[t,:-1,:] # Atmospheric storage of the cell to the north [s]
        Sa_time_S_top[0,:-1,:] = Sa_time_top[t,1:,:] # Atmospheric storage of the cell to the south [s]

        # top: calculate with moisture fluxes 
        Sa_time_after_Fa_top[0,1:-1,:] = ((Sa_track_top[t,1:-1,:] * (ti + Sa_time_top[t,1:-1,:]) 
        - Fa_E_top_WE[t,1:-1,:] * (ti + Sa_time_top[t,1:-1,:]) * (Sa_track_top[t,1:-1,:] / W_top[t,1:-1,:])  
        + Fa_E_top_EW[t,1:-1,:] * (ti + Sa_time_E_top[0,1:-1,:]) * (Sa_track_E_top[0,1:-1,:] / Sa_E_top[0,1:-1,:])
        + Fa_W_top_WE[t,1:-1,:] * (ti + Sa_time_W_top[0,1:-1,:]) * (Sa_track_W_top[0,1:-1,:] / Sa_W_top[0,1:-1,:]) 
        - Fa_W_top_EW[t,1:-1,:] * (ti + Sa_time_top[t,1:-1,:]) * (Sa_track_top[t,1:-1,:] / W_top[t,1:-1,:]) 
        - Fa_N_top_SN[t,1:-1,:] * (ti + Sa_time_top[t,1:-1,:]) * (Sa_track_top[t,1:-1,:] / W_top[t,1:-1,:]) 
        + Fa_N_top_NS[t,1:-1,:] * (ti + Sa_time_N_top[0,1:-1,:]) * (Sa_track_N_top[0,1:-1,:] / Sa_N_top[0,1:-1,:]) 
        + Fa_S_top_SN[t,1:-1,:] * (ti + Sa_time_S_top[0,1:-1,:]) * (Sa_track_S_top[0,1:-1,:] / Sa_S_top[0,1:-1,:]) 
        - Fa_S_top_NS[t,1:-1,:] * (ti + Sa_time_top[t,1:-1,:]) * (Sa_track_top[t,1:-1,:] / W_top[t,1:-1,:]) 
        - Fa_downward[t,1:-1,:] * (ti + Sa_time_top[t,1:-1,:]) * (Sa_track_top[t,1:-1,:] / W_top[t,1:-1,:]) 
        + Fa_upward[t,1:-1,:] * (ti + Sa_time_down[t,1:-1,:]) * (Sa_track_down[t,1:-1,:] / W_down[t,1:-1,:]) 
        ) / Sa_track_after_Fa_top[0,1:-1,:])

        where_are_NaNs = np.isnan(Sa_time_after_Fa_top)
        Sa_time_after_Fa_top[where_are_NaNs] = 0

        # down: substract precipitation and add evaporation
	#
        Sa_time_after_Fa_P_E_down[0,1:-1,:] = ((Sa_track_after_Fa_down[0,1:-1,:] * Sa_time_after_Fa_down[0,1:-1,:] 
        - P[t,1:-1,:] * (ti + Sa_time_down[t,1:-1,:]) * (Sa_track_down[t,1:-1,:] / W[t,1:-1,:])
        + E_region[t,1:-1,:] * ti/2 
        ) / Sa_track_after_Fa_P_E_down[0,1:-1,:])

        where_are_NaNs = np.isnan(Sa_time_after_Fa_P_E_down)
        Sa_time_after_Fa_P_E_down[where_are_NaNs] = 0

        # top: substract precipitation (does not change time)
        Sa_time_after_Fa_P_E_top = Sa_time_after_Fa_top

        # down: redistribute water
        Sa_time_after_all_down = ((Sa_track_after_Fa_P_E_down * Sa_time_after_Fa_P_E_down 
        - down_to_top[t,:,:] * Sa_time_after_Fa_P_E_down
        + top_to_down[t,:,:] * Sa_time_after_Fa_P_E_top
        ) / Sa_track_after_all_down)

        where_are_NaNs = np.isnan(Sa_time_after_all_down)
        Sa_time_after_all_down[where_are_NaNs] = 0

        # top: redistribute water
        Sa_time_after_all_top = ((Sa_track_after_Fa_P_E_top * Sa_time_after_Fa_P_E_top 
        - top_to_down[t,:,:] * Sa_time_after_Fa_P_E_top 
        + down_to_top[t,:,:] * Sa_time_after_Fa_P_E_down
        ) / Sa_track_after_all_top)

        where_are_NaNs = np.isnan(Sa_time_after_all_top)
        Sa_time_after_all_top[where_are_NaNs] = 0

        # down: determine Sa_region of this next timestep 100% stable
        Sa_time_down[t+1,1:-1,:] = Sa_time_after_all_down[0,1:-1,:]

        # top: determine Sa_region of this next timestep 100% stable
        Sa_time_top[t+1,1:-1,:] = Sa_time_after_all_top[0,1:-1,:]
        #############################################################
        # distance tracking start
        
        # defining sizes of timed moisture (1,107,240)
        Sa_dist_after_Fa_down = np.zeros(np.shape(Sa_dist_down_last))
        Sa_dist_after_Fa_P_E_down = np.zeros(np.shape(Sa_dist_down_last))
        Sa_dist_E_down = np.zeros(np.shape(Sa_dist_down_last))
        Sa_dist_W_down = np.zeros(np.shape(Sa_dist_down_last))
        Sa_dist_N_down = np.zeros(np.shape(Sa_dist_down_last))
        Sa_dist_S_down = np.zeros(np.shape(Sa_dist_down_last))
        Sa_dist_after_Fa_top = np.zeros(np.shape(Sa_dist_top_last))
        Sa_dist_after_Fa_P_E_top = np.zeros(np.shape(Sa_dist_top_last))
        Sa_dist_E_top = np.zeros(np.shape(Sa_dist_top_last))
        Sa_dist_W_top = np.zeros(np.shape(Sa_dist_top_last))
        Sa_dist_N_top = np.zeros(np.shape(Sa_dist_top_last))
        Sa_dist_S_top = np.zeros(np.shape(Sa_dist_top_last))

        # distance increase [m]
	# delta_x is the mean of north and south boundaries [1,107,240]
	L_SN_gridcell = 0.5*(L_N_gridcell+L_S_gridcell)
	xi2d = np.tile(L_SN_gridcell,[1,len(longitude)])
        xi = np.reshape(xi2d,[1,len(latitude),len(longitude)])
	# delta_y [1]
	yi = L_EW_gridcell
	# delta_d is weighted mean of delta_x and delta_y [1,107,240]
	# Fa_E and Fa_N can be either positive or negative
	di = (xi*np.abs(Fa_E_down[t-1,:,:])+yi*np.abs(Fa_N_down[t-1,:,:]))/(np.abs(Fa_E_down[t-1,:,:])+np.abs(Fa_N_down[t-1,:,:]))
	#di = 0.5*(xi+yi)

        # down: define values of distance tracking moisture of neighbouring grid cells (1,107,240)
	Sa_dist_E_down[0,:,:-1] = Sa_dist_down[t,:,1:] # Atmospheric storage of the cell to the east [units: m]
        if isglobal == 1:
            Sa_dist_E_down[0,:,-1] = Sa_dist_down[t,:,0] # Atmospheric storage of the cell to the east [m]
        Sa_dist_W_down[0,:,1:] = Sa_dist_down[t,:,:-1] # Atmospheric storage of the cell to the west [m]
        if isglobal == 1:        
            Sa_dist_W_down[0,:,0] = Sa_dist_down[t,:,-1] # Atmospheric storage of the cell to the west [m]
        Sa_dist_N_down[0,1:,:] = Sa_dist_down[t,:-1,:] # Atmospheric storage of the cell to the north [m]
        Sa_dist_S_down[0,:-1,:] = Sa_dist_down[t,1:,:] # Atmospheric storage of the cell to the south [m]

	# down: calculate with moisture fluxes. LG: calculate the first, second and third terms of the numerator of Eq. 2.10
	Sa_dist_after_Fa_down[0,1:-1,:] = ((Sa_track_down[t,1:-1,:]  * (               Sa_dist_down[t,1:-1,:])
	- Fa_E_down_WE[t,1:-1,:] * (xi[0,1:-1,:] + Sa_dist_down[t,1:-1,:]  ) * (Sa_track_down[t,1:-1,:] / W_down[t,1:-1,:]) 
	+ Fa_E_down_EW[t,1:-1,:] * (               Sa_dist_E_down[0,1:-1,:]) * (Sa_track_E_down[0,1:-1,:] / Sa_E_down[0,1:-1,:])
	+ Fa_W_down_WE[t,1:-1,:] * (               Sa_dist_W_down[0,1:-1,:]) * (Sa_track_W_down[0,1:-1,:] / Sa_W_down[0,1:-1,:])
	- Fa_W_down_EW[t,1:-1,:] * (xi[0,1:-1,:] + Sa_dist_down[t,1:-1,:]  ) * (Sa_track_down[t,1:-1,:] / W_down[t,1:-1,:])
	- Fa_N_down_SN[t,1:-1,:] * (yi           + Sa_dist_down[t,1:-1,:]  ) * (Sa_track_down[t,1:-1,:] / W_down[t,1:-1,:])
	+ Fa_N_down_NS[t,1:-1,:] * (               Sa_dist_N_down[0,1:-1,:]) * (Sa_track_N_down[0,1:-1,:] / Sa_N_down[0,1:-1,:])
	+ Fa_S_down_SN[t,1:-1,:] * (               Sa_dist_S_down[0,1:-1,:]) * (Sa_track_S_down[0,1:-1,:] / Sa_S_down[0,1:-1,:])
	- Fa_S_down_NS[t,1:-1,:] * (yi           + Sa_dist_down[t,1:-1,:]  ) * (Sa_track_down[t,1:-1,:] / W_down[t,1:-1,:])
	+ Fa_downward[t,1:-1,:]  * (               Sa_dist_top[t,1:-1,:]   ) * (Sa_track_top[t,1:-1,:] / W_top[t,1:-1,:])
	- Fa_upward[t,1:-1,:]    * (               Sa_dist_down[t,1:-1,:]  ) * (Sa_track_down[t,1:-1,:] / W_down[t,1:-1,:])
        ) / Sa_track_after_Fa_down[0,1:-1,:])

        where_are_NaNs = np.isnan(Sa_dist_after_Fa_down)
        Sa_dist_after_Fa_down[where_are_NaNs] = 0 

        # top: define values of distance tracking moisture of neighbouring grid cells
	Sa_dist_E_top[0,:,:-1] = Sa_dist_top[t,:,1:] # Atmospheric storage of the cell to the east [units: s]
        if isglobal == 1:
            Sa_dist_E_top[0,:,-1] = Sa_dist_top[t,:,0] # Atmospheric storage of the cell to the east [m]
        Sa_dist_W_top[0,:,1:] = Sa_dist_top[t,:,:-1] # Atmospheric storage of the cell to the west [m]
        if isglobal == 1:
            Sa_dist_W_top[0,:,0] = Sa_dist_top[t,:,-1] # Atmospheric storage of the cell to the west [m]
        Sa_dist_N_top[0,1:,:] = Sa_dist_top[t,:-1,:] # Atmospheric storage of the cell to the north [m]
        Sa_dist_S_top[0,:-1,:] = Sa_dist_top[t,1:,:] # Atmospheric storage of the cell to the south [m]

        # top: calculate with moisture fluxes 
	Sa_dist_after_Fa_top[0,1:-1,:] = ((Sa_track_top[t,1:-1,:] * (               Sa_dist_top[t,1:-1,:]) 
	- Fa_E_top_WE[t-1,1:-1,:] * (xi[0,1:-1,:] + Sa_dist_top[t,1:-1,:]  ) * (Sa_track_top[t,1:-1,:] / W_top[t,1:-1,:]) 
	+ Fa_E_top_EW[t-1,1:-1,:] * (               Sa_dist_E_top[0,1:-1,:]) * (Sa_track_E_top[0,1:-1,:] / Sa_E_top[0,1:-1,:])
	+ Fa_W_top_WE[t-1,1:-1,:] * (               Sa_dist_W_top[0,1:-1,:]) * (Sa_track_W_top[0,1:-1,:] / Sa_W_top[0,1:-1,:])
	- Fa_W_top_EW[t-1,1:-1,:] * (xi[0,1:-1,:] + Sa_dist_top[t,1:-1,:]  ) * (Sa_track_top[t,1:-1,:] / W_top[t,1:-1,:])
	- Fa_N_top_SN[t-1,1:-1,:] * (yi           + Sa_dist_top[t,1:-1,:]  ) * (Sa_track_top[t,1:-1,:] / W_top[t,1:-1,:])
	+ Fa_N_top_NS[t-1,1:-1,:] * (               Sa_dist_N_top[0,1:-1,:]) * (Sa_track_N_top[0,1:-1,:] / Sa_N_top[0,1:-1,:])
	+ Fa_S_top_SN[t-1,1:-1,:] * (               Sa_dist_S_top[0,1:-1,:]) * (Sa_track_S_top[0,1:-1,:] / Sa_S_top[0,1:-1,:])
	- Fa_S_top_NS[t-1,1:-1,:] * (yi           + Sa_dist_top[t,1:-1,:]  ) * (Sa_track_top[t,1:-1,:] / W_top[t,1:-1,:])
	- Fa_downward[t-1,1:-1,:] * (               Sa_dist_top[t,1:-1,:]  ) * (Sa_track_top[t,1:-1,:] / W_top[t,1:-1,:])
	+ Fa_upward[t-1,1:-1,:]   * (               Sa_dist_down[t,1:-1,:] ) * (Sa_track_down[t,1:-1,:] / W_down[t,1:-1,:])
        ) / Sa_track_after_Fa_top[0,1:-1,:])

        where_are_NaNs = np.isnan(Sa_dist_after_Fa_top)
        Sa_dist_after_Fa_top[where_are_NaNs] = 0

        # down: substract precipitation but no adding evaporation (as dist=0).
	# note that calcualtions for the forth and fifth terms are swapped as this is backward calculation.
        Sa_dist_after_Fa_P_E_down[0,1:-1,:] = ((Sa_track_after_Fa_down[0,1:-1,:] * Sa_dist_after_Fa_down[0,1:-1,:] 
	- P[t,1:-1,:] * Sa_dist_down[t,1:-1,:] * (Sa_track_down[t,1:-1,:] / W_down[t,1:-1,:])
        ) / Sa_track_after_Fa_P_E_down[0,1:-1,:])

        where_are_NaNs = np.isnan(Sa_dist_after_Fa_P_E_down)
        Sa_dist_after_Fa_P_E_down[where_are_NaNs] = 0

        # top: add precipitation
	# LG: in Con_P_Recyc_Masterscript.py, Sa_dist_after_Fa_P_E_top is unchanged. But, I think it should.
        Sa_dist_after_Fa_P_E_top[0,1:-1,:] = ((Sa_track_after_Fa_top[0,1:-1,:] * Sa_dist_after_Fa_top[0,1:-1,:]
	- P[t,1:-1,:] * Sa_dist_top[t,1:-1,:] * (Sa_track_top[t,1:-1,:] / W_top[t,1:-1,:])
        ) / Sa_track_after_Fa_P_E_top[0,1:-1,:])

        where_are_NaNs = np.isnan(Sa_dist_after_Fa_P_E_top)
        Sa_dist_after_Fa_P_E_top[where_are_NaNs] = 0

        # down: redistribute water
        Sa_dist_after_all_down = ((Sa_track_after_Fa_P_E_down * Sa_dist_after_Fa_P_E_down 
        - down_to_top[t,:,:] * Sa_dist_after_Fa_P_E_down
        + top_to_down[t,:,:] * Sa_dist_after_Fa_P_E_top
        ) / Sa_track_after_all_down)

        where_are_NaNs = np.isnan(Sa_dist_after_all_down)
        Sa_dist_after_all_down[where_are_NaNs] = 0

        # top: redistribute water
        Sa_dist_after_all_top = ((Sa_track_after_Fa_P_E_top * Sa_dist_after_Fa_P_E_top
        - top_to_down[t,:,:] * Sa_dist_after_Fa_P_E_top 
        + down_to_top[t,:,:] * Sa_dist_after_Fa_P_E_down
        ) / Sa_track_after_all_top)

        where_are_NaNs = np.isnan(Sa_dist_after_all_top)
        Sa_dist_after_all_top[where_are_NaNs] = 0

        # down: determine Sa_region of this next timestep 100% stable (97,107,240)
        Sa_dist_down[t+1,1:-1,:] = Sa_dist_after_all_down[0,1:-1,:]

        # top: determine Sa_region of this next timestep 100% stable
        Sa_dist_top[t+1,1:-1,:] = Sa_dist_after_all_top[0,1:-1,:]
        #############################################################
                                                        
    return Sa_dist_top,Sa_dist_down,Sa_time_top,Sa_time_down,Sa_track_top,Sa_track_down,north_loss,south_loss,down_to_top,top_to_down,water_lost

def wrap_netcdf(year_o,yearpart_o,var_o,standard_name_o,units_o):
    # Define Coordinaties
    start_date = (datetime.datetime(year_o,1,1)+datetime.timedelta(yearpart_o)).strftime('%Y-%m-%d')
    dim0 = cf.DimensionCoordinate(properties={'standard_name':'time'},data=cf.Data(0.,cf.Units('days since '+start_date,calendar='standard')))
    dim1 = cf.DimensionCoordinate(data=cf.Data(latitude,'degrees_north'),properties={'standard_name':'latitude'})
    dim2 = cf.DimensionCoordinate(data=cf.Data(longitude,'degrees_east'),properties={'standard_name':'longitude'})
    # Define cf.Field then insert variable and coordinates
    f = cf.Field(properties={'standard_name':standard_name_o})
    f.insert_dim(dim0)
    f.insert_dim(dim1)
    f.insert_dim(dim2)
    data = cf.Data(var_o,units_o)
    f.insert_data(data)
    return f

def create_empty_array(datapathea,count_time,divt,latitude,longitude,yearpart,years):
    Sa_track_down_last = np.zeros((1,len(latitude),len(longitude)))
    Sa_track_down_day = np.zeros((1,len(latitude),len(longitude)))
    Sa_track_top_last = np.zeros((1,len(latitude),len(longitude)))
    Sa_track_top_day = np.zeros((1,len(latitude),len(longitude)))
    #
    Sa_time_down_last = np.zeros((1,len(latitude),len(longitude)))
    Sa_time_down_day = np.zeros((1,len(latitude),len(longitude)))
    Sa_time_top_last = np.zeros((1,len(latitude),len(longitude)))
    Sa_time_top_day = np.zeros((1,len(latitude),len(longitude)))
    #
    Sa_dist_down_last = np.zeros((1,len(latitude),len(longitude)))
    Sa_dist_down_day = np.zeros((1,len(latitude),len(longitude)))
    Sa_dist_top_last = np.zeros((1,len(latitude),len(longitude)))
    Sa_dist_top_day = np.zeros((1,len(latitude),len(longitude)))

    # Build cf.field here.
    if yearpart[0] == 0:
        year_o = years[0]-1
	yearpart_o = 364 + calendar.isleap(year_o)
    else:
        year_o = years[0]
        yearpart_o = yearpart[0]-1

    f0l = wrap_netcdf(year_o,yearpart_o,Sa_track_down_last,'Sa_track_down_last','m3')
    f0d = wrap_netcdf(year_o,yearpart_o,Sa_track_down_day,'Sa_track_down','m3')
    f1l = wrap_netcdf(year_o,yearpart_o,Sa_track_top_last,'Sa_track_top_last','m3')
    f1d = wrap_netcdf(year_o,yearpart_o,Sa_track_top_day,'Sa_track_top','m3')
    #
    f2l = wrap_netcdf(year_o,yearpart_o,Sa_time_down_last,'Sa_time_down_last','s')
    f2d = wrap_netcdf(year_o,yearpart_o,Sa_time_down_day,'Sa_time_down','s')
    f3l = wrap_netcdf(year_o,yearpart_o,Sa_time_top_last,'Sa_time_top_last','s')
    f3d = wrap_netcdf(year_o,yearpart_o,Sa_time_top_day,'Sa_time_top','s')
    #
    f4l = wrap_netcdf(year_o,yearpart_o,Sa_dist_down_last,'Sa_dist_down_last','m')
    f4d = wrap_netcdf(year_o,yearpart_o,Sa_dist_down_day,'Sa_dist_down','m')
    f5l = wrap_netcdf(year_o,yearpart_o,Sa_dist_top_last,'Sa_dist_top_last','m')
    f5d = wrap_netcdf(year_o,yearpart_o,Sa_dist_top_day,'Sa_dist_top','m')

    # Write out netcdf
    if yearpart[0] == 0:
        datapathnc = datapathea[0]
    else:
        datapathnc = datapathea[1]
    f = cf.FieldList([f0l,f0d,f1l,f1d,f2l,f2d,f3l,f3d,f4l,f4d,f5l,f5d])
    cf.write(f,datapathnc,single=True,unlimited='time')
    return

#%% Runtime & Results
start1 = timer()

datapathea = data_path_ea(years,yearpart)
if veryfirstrun == 1:
    create_empty_array(datapathea,count_time,divt,latitude,longitude,yearpart,years)

# loop through the years
for yearnumber in years:
    if (yearpart[-1] == 365) & (calendar.isleap(yearnumber) == 0):
        thisyearpart = yearpart[:-1]
    else: # a leapyear
        thisyearpart = yearpart
        
    for a in thisyearpart:
        start = timer()

        if a == 0:
            previous_data_to_load = (str(yearnumber-1)+'-'+str(364+calendar.isleap(yearnumber-1)).zfill(3))
        else:
            previous_data_to_load = (str(yearnumber)+'-'+str(a-1).zfill(3))
        
        datapath = data_path(previous_data_to_load,yearnumber,a)

        # choose monthly source region mask
	# find corresponding month using year and days of the year information
        sr_mn = int((datetime.datetime(yearnumber,1,1)+datetime.timedelta(a)).strftime('%m')) - 1
	print yearnumber,a,sr_mn
	Region = source_region[sr_mn,:,:]
        
        ST = cf.read(datapath[0])
        Sa_track_top_last = ST.select('Sa_track_top_last')[0].array
        Sa_track_down_last = ST.select('Sa_track_down_last')[0].array
        if timetracking == 1:
            Sa_time_top_last = ST.select('Sa_time_top_last')[0].array
            Sa_time_down_last = ST.select('Sa_time_down_last')[0].array
            Sa_dist_top_last = ST.select('Sa_dist_top_last')[0].array
            Sa_dist_down_last = ST.select('Sa_dist_down_last')[0].array
        
#       loading_FS = sio.loadmat(datapath[1],verify_compressed_data_integrity=False)
#       Fa_E_top = loading_FS['Fa_E_top']
#       Fa_N_top = loading_FS['Fa_N_top']
#       Fa_E_down = loading_FS['Fa_E_down']
#       Fa_N_down = loading_FS['Fa_N_down']
#       Fa_Vert = loading_FS['Fa_Vert']
#       E = loading_FS['E']
#       P = loading_FS['P']
#       W_top = loading_FS['W_top']
#       W_down = loading_FS['W_down']

        FS = cf.read(datapath[1])
	Fa_E_top = FS.select('Fa_E_top')[0].array
        Fa_N_top = FS.select('Fa_N_top')[0].array
        Fa_E_down = FS.select('Fa_E_down')[0].array
        Fa_N_down = FS.select('Fa_N_down')[0].array
        Fa_Vert = FS.select('Fa_Vert')[0].array
        E = FS.select('E')[0].array
        P = FS.select('P')[0].array
        W_top = FS.select('W_top')[0].array
        W_down = FS.select('W_down')[0].array
        
        # call the forward tracking function
        if timetracking == 0:
            Sa_track_top,Sa_track_down,north_loss,south_loss,down_to_top,top_to_down,water_lost = get_Sa_track_forward(latitude,longitude,count_time,divt,Kvf,Region,Fa_E_top,Fa_N_top,Fa_E_down,Fa_N_down,
                                       Fa_Vert,E,P,W_top,W_down,Sa_track_top_last,Sa_track_down_last)   
        elif timetracking == 1:
            Sa_dist_top,Sa_dist_down,Sa_time_top,Sa_time_down,Sa_track_top,Sa_track_down,north_loss,south_loss,down_to_top,top_to_down,water_lost = get_Sa_track_forward_TIME(latitude,longitude,count_time,divt,timestep,Kvf,Region,Fa_E_top,Fa_N_top,Fa_E_down,Fa_N_down,Fa_Vert,E,P,W_top,W_down,Sa_track_top_last,Sa_track_down_last,Sa_time_top_last,Sa_time_down_last,Sa_dist_top_last,Sa_dist_down_last,L_N_gridcell,L_S_gridcell,L_EW_gridcell)
            
        # write out as netCDF format
	f0l = wrap_netcdf(yearnumber,a,Sa_track_down[[-1],:,:],'Sa_track_down_last','m3') # use [] to keep the degrading dimension
	f1l = wrap_netcdf(yearnumber,a,Sa_track_top[[-1],:,:],'Sa_track_top_last','m3')
	f2l = wrap_netcdf(yearnumber,a,Sa_time_down[[-1],:,:],'Sa_time_down_last','s')
	f3l = wrap_netcdf(yearnumber,a,Sa_time_top[[-1],:,:],'Sa_time_top_last','s')
 	f4l = wrap_netcdf(yearnumber,a,Sa_dist_down[[-1],:,:],'Sa_dist_down_last','m')
 	f5l = wrap_netcdf(yearnumber,a,Sa_dist_top[[-1],:,:],'Sa_dist_top_last','m')

        
        Sa_track_down_day = np.mean(Sa_track_down[:-1,:,:],axis=0,keepdims=True)
	f0d = wrap_netcdf(yearnumber,a,Sa_track_down_day,'Sa_track_down','m3')
        Sa_track_top_day = np.mean(Sa_track_top[:-1,:,:],axis=0,keepdims=True)
	f1d = wrap_netcdf(yearnumber,a,Sa_track_top_day,'Sa_track_top','m3')
        Sa_time_down_day = np.mean(Sa_time_down[:-1,:,:]*Sa_track_down[:-1,:,:],axis=0,keepdims=True)/Sa_track_down_day
	f2d = wrap_netcdf(yearnumber,a,Sa_time_down_day,'Sa_time_down','s')
        Sa_time_top_day = np.mean(Sa_time_top[:-1,:,:]*Sa_track_top[:-1,:,:],axis=0,keepdims=True)/Sa_track_top_day
	f3d = wrap_netcdf(yearnumber,a,Sa_time_top_day,'Sa_time_top','s')
        Sa_dist_down_day = np.mean(Sa_dist_down[:-1,:,:]*Sa_track_down[:-1,:,:],axis=0,keepdims=True)/Sa_track_down_day
 	f4d = wrap_netcdf(yearnumber,a,Sa_dist_down_day,'Sa_dist_down','m')
        Sa_dist_top_day = np.mean(Sa_dist_top[:-1,:,:]*Sa_track_top[:-1,:,:],axis=0,keepdims=True)/Sa_track_top_day
 	f5d = wrap_netcdf(yearnumber,a,Sa_dist_top_day,'Sa_dist_top','m')

        W = W_top + W_down
        Sa_track = Sa_track_top + Sa_track_down
        P_track = P*(Sa_track[:-1,:,:]/W[:-1,:,:])
	P_track_day = np.nansum(P_track,axis =0,keepdims=True)	# using np.nansum rather than np.sum, in which NaNs are treated as zero.
	f6d = wrap_netcdf(yearnumber,a,P_track_day,'P_track','m3')

        P_track_down = P*(Sa_track_down[:-1,:,:]/W[:-1,:,:])
        P_track_top = P*(Sa_track_top[:-1,:,:]/W[:-1,:,:])
        P_time_down = 0.5*(Sa_time_down[:-1,:,:]+Sa_time_down[1:,:,:])
        P_time_top = 0.5*(Sa_time_top[:-1,:,:]+Sa_time_top[1:,:,:])
        P_time_day = np.nansum((P_time_down*P_track_down+P_time_top*P_track_top),axis=0,keepdims=True)/P_track_day
        where_are_NaNs = np.isnan(P_time_day)
        P_time_day[where_are_NaNs] = 0
        f7d = wrap_netcdf(yearnumber,a,P_time_day,'P_time','s')

        P_dist_down = 0.5*(Sa_dist_down[:-1,:,:]+Sa_dist_down[1:,:,:])
        P_dist_top = 0.5*(Sa_dist_top[:-1,:,:]+Sa_dist_top[1:,:,:])
        P_dist_day = np.nansum((P_dist_down*P_track_down+P_dist_top*P_track_top),axis=0,keepdims=True)/P_track_day
        where_are_NaNs = np.isnan(P_dist_day)
        P_dist_day[where_are_NaNs] = 0
 	f8d = wrap_netcdf(yearnumber,a,P_dist_day,'P_dist','m')

	f9d = wrap_netcdf(yearnumber,a,np.nansum(E,axis=0,keepdims=True),'E','m3')
	f10d = wrap_netcdf(yearnumber,a,np.nansum(P,axis=0,keepdims=True),'P','m3')

	f11d = wrap_netcdf(yearnumber,a,np.mean(W_down[1:,:,:],axis=0,keepdims=True),'W_down','m3')
	f12d = wrap_netcdf(yearnumber,a,np.mean(W_top[1:,:,:],axis=0,keepdims=True),'W_top','m3')

	f13d = wrap_netcdf(yearnumber,a,np.nansum(Fa_E_down,axis=0,keepdims=True),'Fa_E_down','m3')
	f14d = wrap_netcdf(yearnumber,a,np.nansum(Fa_E_top,axis=0,keepdims=True),'Fa_E_top','m3')
	f15d = wrap_netcdf(yearnumber,a,np.nansum(Fa_N_down,axis=0,keepdims=True),'Fa_N_down','m3')
	f16d = wrap_netcdf(yearnumber,a,np.nansum(Fa_N_top,axis=0,keepdims=True),'Fa_N_top','m3')
	f17d = wrap_netcdf(yearnumber,a,np.nansum(Fa_Vert,axis=0,keepdims=True),'Fa_Vert','m3')

 	f = cf.FieldList([f0l,f0d,f1l,f1d,f2l,f2d,f3l,f3d,f4l,f4d,f5l,f5d,f6d,f7d,f8d,f9d,f10d,f11d,f12d,f13d,f14d,f15d,f16d,f17d])
	cf.write(f,datapath[2],single=True,unlimited='time')
            
        end = timer()
        print 'Runtime Sa_track for day '+str(a+1)+' in year '+str(yearnumber)+' is',(end - start),' seconds.'
end1 = timer()
print 'The total runtime of Con_P_Recyc_Masterscript is',(end1-start1),' seconds.'
