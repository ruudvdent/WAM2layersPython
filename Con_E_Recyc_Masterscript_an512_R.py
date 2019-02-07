#!/home/users/lguo/anaconda2/bin/python
#BSUB -o %J.o
#BSUB -e %J.e
#BSUB -q long-serial
#BSUB -W 72:00
#BSUB -R "rusage[mem=32000]"
#BSUB -M 32000

#%% Import libraries
import sys
sys.path.append('/home/users/lguo/perchance/WAM')
import numpy as np
from netCDF4 import Dataset
import os
from timeit import default_timer as timer
import scipy.io as sio
import WAM_functions as WAM
import cf
import datetime
#import pdb

#%% BEGIN OF INPUT1 (FILL THIS IN)
years = np.arange(0003,0001,-1)
yearpart0 = np.arange(29,-1,-1)
yearpart = np.arange(359,-1,-1)
divt = 240
count_time = 1

# Manage the extent of your dataset (FILL THIS IN)
latnrs = np.arange(0,469) # 79.1667N ~ 30.2778S
lonnrs = np.arange(0,540) # -9.58334 ~ 184.583
isglobal = 0 # fill in 1 for global computations (i.e. Earth round), fill in 0 for a local domain with boundaries

# obtain the constants
invariant_data = '/home/users/lguo/perchance/WAM_input/cn512/invariants.nc'

latitude = Dataset(invariant_data,mode='r').variables['latitude'][latnrs]
longitude = Dataset(invariant_data,mode='r').variables['longitude'][lonnrs]

# Constants 
g = 9.80665 # [m/s2] from ERA-interim archive
density_water = 1000 # [kg/m3]
dg = 111089.56 # [m] length of 1 degree latitude
timestep = 24*3600 # [s] timestep for um daily output (watch out! P & E have 3 hour timestep)
Erad = 6.371e6 # [m] Earth radius
    
# Semiconstants
# note that, in models, delta_lat is not equal to delta_lon
gridcell = np.abs(longitude[1]-longitude[0]) # [degrees E] grid cell size, 1.875 degree
gridcella = np.abs(latitude[1]-latitude[0]) # [degrees N] grid cell size, 1.25  degree

# new area size calculation:
lat_n_bound = np.minimum(90.0,latitude+0.5*gridcella) # Finding north and south boundaries of each gridcell
lat_s_bound = np.maximum(-90.0,latitude-0.5*gridcella)
A_gridcell = np.zeros([len(latitude),1]) # area of each gridcell in units: m2 with size (y). see
A_gridcell[:,0] = Erad**2*(np.pi/180.0)*gridcell*abs(np.sin(lat_s_bound*np.pi/180.0)-np.sin(lat_n_bound*np.pi/180.0))

L_N_gridcell = gridcell*np.cos((latitude+gridcella*0.5)*np.pi/180.0)*dg # [m] length northern boundary of a cell, [107]
L_S_gridcell = gridcell*np.cos((latitude-gridcella*0.5)*np.pi/180.0)*dg # [m] length southern boundary of a cell, [107]
L_EW_gridcell = np.float32(gridcella*dg) # [m] length eastern/western boundary of a cell, [1]

# read in mask files of China regions
cn_mask = Dataset('/home/users/lguo/perchance/WAM_input/masks/mask_china_simplified2_n512_R.nc',mode='r').variables['msk_china'][latnrs,lonnrs]

# BEGIN OF INPUT 2 (FILL THIS IN)
Region = np.zeros((np.shape(cn_mask)),dtype=np.float32)
Region[cn_mask==1] = 1
Kvf = 3 # vertical dispersion factor (advection only is 0, dispersion the same size of the advective flux is 1, for stability don't make this more than 3)
timetracking = 1 # 0 for not tracking time and 1 for tracking time
veryfirstrun = 1 # type '1' if no run has been done before from which can be continued, otherwise type '0'

interdata_folder = '/home/users/lguo/ncas_climate_v1/WAM_inter/interdata_an512'

# Number of CPU cores to use (OMP threads)
num_threads = 1
#END OF INPUT

#%% Datapaths (FILL THIS IN)
# Check if interdata folder exists:
assert os.path.isdir(interdata_folder),"Please create the interdata_folder before running the script"

# Check if sub interdata folder exists otherwise create it:
sub_interdata_folder = os.path.join(interdata_folder,'cn1_backward')
if os.path.isdir(sub_interdata_folder):
    pass
else:
    os.makedirs(sub_interdata_folder)

# create sub interdata folder for unwanted data:
sub_interdata_tmp_folder = os.path.join(interdata_folder,'cn1_tmp')
if os.path.isdir(sub_interdata_tmp_folder):
    pass
else:
    os.makedirs(sub_interdata_tmp_folder)

#
def wrap_netcdf(year_o,yearpart_o,var_o,standard_name_o,units_o):
    var_shape = var_o.shape

    # Define Coordinaties
    start_date = (datetime.datetime(year_o,1,1)+datetime.timedelta(yearpart_o)).strftime('%Y-%m-%d')
    if var_shape[0] == 1:
        dim0 = cf.DimensionCoordinate(properties={'standard_name':'time'},data=cf.Data(0.,cf.Units('days since '+start_date,calendar='standard')))
    elif var_shape[0] == 240:
        nc_time = (86400.0/count_time/divt)*np.arange(count_time*divt)
	dim0 = cf.DimensionCoordinate(properties={'standard_name':'time'},data=cf.Data(nc_time,cf.Units('seconds since '+start_date+' 0:3:0',calendar='standard')))
    elif var_shape[0] == 241:
        nc_time = (86400.0/count_time/divt)*np.arange(count_time*divt+1)
	dim0 = cf.DimensionCoordinate(properties={'standard_name':'time'},data=cf.Data(nc_time,cf.Units('seconds since '+start_date+' 0:0:0',calendar='standard')))
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
    if yearpart[0] == 365:
        year_o = years[0]+1
	yearpart_o = 0
    else:
        year_o = years[0]
        yearpart_o = yearpart[0]+1

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
    if yearpart[0] == 365:
        datapathnc = datapathea[0]
    else:
        datapathnc = datapathea[1]
    f_list = cf.FieldList([f0l,f0d,f1l,f1d,f2l,f2d,f3l,f3d,f4l,f4d,f5l,f5d])
    cf.write(f_list,datapathnc,single=True,unlimited='time',compress=3)

    return

def data_path_ea(years,yearpart,sub_interdata_tmp_folder):
    save_empty_arrays_ld_track = os.path.join(sub_interdata_tmp_folder,str(years[0]+1)+'-'+str(0).zfill(3)+'Sa_track.nc')
    save_empty_arrays_track = os.path.join(sub_interdata_tmp_folder,str(years[0])+'-'+str(yearpart[0]+1).zfill(3)+'Sa_track.nc')
    return save_empty_arrays_ld_track,save_empty_arrays_track

def data_path(previous_data_to_load,yearnumber,a,i,veryfistrun,interdata_folder,sub_interdata_folder,sub_interdata_tmp_folder):
    if i == 0 and veryfirstrun == 1:
        load_Sa_track = os.path.join(sub_interdata_tmp_folder,previous_data_to_load+'Sa_track.nc')
        load_fluxes_and_storages = os.path.join(interdata_folder,str(yearnumber)+'-'+str(a).zfill(3)+'fluxes_storages.nc')
        save_path_track = os.path.join(sub_interdata_tmp_folder,str(yearnumber)+'-'+str(a).zfill(3)+'Sa_track.nc')
    elif i == 1 and a == 359:
        load_Sa_track = os.path.join(sub_interdata_tmp_folder,previous_data_to_load+'Sa_track.nc')
        load_fluxes_and_storages = os.path.join(interdata_folder,str(yearnumber)+'-'+str(a).zfill(3)+'fluxes_storages.nc')
        save_path_track = os.path.join(sub_interdata_folder,str(yearnumber)+'-'+str(a).zfill(3)+'Sa_track.nc')
    else:
        load_Sa_track = os.path.join(sub_interdata_folder,previous_data_to_load+'Sa_track.nc')
        load_fluxes_and_storages = os.path.join(interdata_folder,str(yearnumber)+'-'+str(a).zfill(3)+'fluxes_storages.nc')
        save_path_track = os.path.join(sub_interdata_folder,str(yearnumber)+'-'+str(a).zfill(3)+'Sa_track.nc')
    return load_Sa_track,load_fluxes_and_storages,save_path_track

#%% Runtime & Results
start1 = timer()

# The two lines below create empty arrays for first runs/initial values are zero.
datapathea = data_path_ea(years,yearpart0,sub_interdata_tmp_folder)
if veryfirstrun == 1:
    create_empty_array(datapathea,count_time,divt,latitude,longitude,yearpart0,years)

# loop through the years
for i in range(len(years)):
    yearnumber = years[i]

    # the first year (30 days) data is unwanted
    if i == 0 and veryfirstrun == 1:
        thisyearpart = yearpart0
    else:
        thisyearpart = yearpart

    # Loop through one year
    for j in range(len(thisyearpart)):

        start = timer()

        a = thisyearpart[j]
        if a == 359:
            previous_data_to_load = (str(yearnumber+1)+'-000')
        else:
            previous_data_to_load = (str(yearnumber)+'-'+str(a+1).zfill(3))

        datapath = data_path(previous_data_to_load,yearnumber,a,i,veryfirstrun,interdata_folder,sub_interdata_folder,sub_interdata_tmp_folder)

        ST = cf.read(datapath[0])
        Sa_track_top_last = ST.select('Sa_track_top_last')[0].array
        Sa_track_down_last = ST.select('Sa_track_down_last')[0].array
        if timetracking == 1:
            Sa_time_top_last = ST.select('Sa_time_top_last')[0].array
            Sa_time_down_last = ST.select('Sa_time_down_last')[0].array
            Sa_dist_top_last = ST.select('Sa_dist_top_last')[0].array
            Sa_dist_down_last = ST.select('Sa_dist_down_last')[0].array

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

        # call the backward tracking function
        if timetracking == 0:
            Sa_track_top,Sa_track_down,north_loss,south_loss,down_to_top,top_to_down,water_lost = WAM.get_Sa_track_backward(latitude,longitude,count_time,divt,Kvf,Region,Fa_E_top,Fa_N_top,Fa_E_down,Fa_N_down,Fa_Vert,E,P,W_top,W_down,Sa_track_top_last,Sa_track_down_last,isglobal)
        elif timetracking == 1:
            Sa_dist_top,Sa_dist_down,Sa_time_top,Sa_time_down,Sa_track_top,Sa_track_down,north_loss,south_loss,down_to_top,top_to_down,water_lost = WAM.get_Sa_track_backward_TIME(latitude,longitude,count_time,divt,timestep,Kvf,Region,Fa_E_top,Fa_N_top,Fa_E_down,Fa_N_down,Fa_Vert,E,P,W_top,W_down,Sa_track_top_last,Sa_track_down_last,Sa_time_top_last,Sa_time_down_last,Sa_dist_top_last,Sa_dist_down_last,L_N_gridcell,L_S_gridcell,L_EW_gridcell,isglobal,num_threads=num_threads)

        # save this data
        f0l = wrap_netcdf(yearnumber,a,Sa_track_down[[0],:,:],'Sa_track_down_last','m3')
	f1l = wrap_netcdf(yearnumber,a,Sa_track_top[[0],:,:],'Sa_track_top_last','m3')
	f2l = wrap_netcdf(yearnumber,a,Sa_time_down[[0],:,:],'Sa_time_down_last','s')
	f3l = wrap_netcdf(yearnumber,a,Sa_time_top[[0],:,:],'Sa_time_top_last','s')
	f4l = wrap_netcdf(yearnumber,a,Sa_dist_down[[0],:,:],'Sa_dist_down_last','m')
	f5l = wrap_netcdf(yearnumber,a,Sa_dist_top[[0],:,:],'Sa_dist_top_last','m')

	f0d = wrap_netcdf(yearnumber,a,np.mean(Sa_track_down[1:,:,:],axis =0,keepdims=True),'Sa_track_down','m3')
	f1d = wrap_netcdf(yearnumber,a,np.mean(Sa_track_top[1:,:,:],axis =0,keepdims=True),'Sa_track_top','m3')
	f2d = wrap_netcdf(yearnumber,a,np.mean(Sa_time_down[:-1,:,:],axis=0,keepdims=True),'Sa_time_down','s')
	f3d = wrap_netcdf(yearnumber,a,np.mean(Sa_time_top[:-1,:,:],axis=0,keepdims=True),'Sa_time_top','s')
	f4d = wrap_netcdf(yearnumber,a,np.mean(Sa_dist_down[:-1,:,:],axis=0,keepdims=True),'Sa_dist_down','m')
	f5d = wrap_netcdf(yearnumber,a,np.mean(Sa_dist_top[:-1,:,:],axis=0,keepdims=True),'Sa_dist_top','m')

        # This is a difference between P and E, E_track is calculated from Sa_track_down, while P_track is caculated from Sa_track
        # Sa_track = Sa_track_top+Sa_track_down
        E_track = E*(Sa_track_down[1:,:,:]/W_down[1:,:,:])
	E_track_day = np.nansum(E_track,axis =0,keepdims=True)
        f6d = wrap_netcdf(yearnumber,a,E_track_day,'E_track','m3')

        E_time = 0.5*(Sa_time_down[:-1,:,:]+Sa_time_down[1:,:,:])
        E_time_day = np.nansum((E_time*E_track),axis=0,keepdims=True)/E_track_day
        where_are_NaNs = np.isnan(E_time_day)
        E_time_day[where_are_NaNs] = 0
        f7d = wrap_netcdf(yearnumber,a,E_time_day,'E_time','s')

        E_dist = 0.5*(Sa_dist_down[:-1,:,:]+Sa_dist_down[1:,:,:])
        E_dist_day = np.nansum((E_dist*E_track),axis=0,keepdims=True)/E_track_day
        where_are_NaNs = np.isnan(E_dist_day)
        E_dist_day[where_are_NaNs] = 0
	f8d = wrap_netcdf(yearnumber,a,E_dist_day,'E_dist','m')

        f9d = wrap_netcdf(yearnumber,a,np.nansum(E,axis=0,keepdims=True),'E','m3')
        f9s = wrap_netcdf(yearnumber,a,E,'E','m3')
        f10d = wrap_netcdf(yearnumber,a,np.nansum(P,axis=0,keepdims=True),'P','m3')
	f10s = wrap_netcdf(yearnumber,a,P,'P','m3')

        f11d = wrap_netcdf(yearnumber,a,np.mean(W_down[1:,:,:],axis=0,keepdims=True),'W_down','m3')
        f12d = wrap_netcdf(yearnumber,a,np.mean(W_top[1:,:,:],axis=0,keepdims=True),'W_top','m3')

        f13d = wrap_netcdf(yearnumber,a,np.nansum(Fa_E_down,axis=0,keepdims=True),'Fa_E_down','m3')
        f14d = wrap_netcdf(yearnumber,a,np.nansum(Fa_E_top,axis=0,keepdims=True),'Fa_E_top','m3')
        f15d = wrap_netcdf(yearnumber,a,np.nansum(Fa_N_down,axis=0,keepdims=True),'Fa_N_down','m3')
        f16d = wrap_netcdf(yearnumber,a,np.nansum(Fa_N_top,axis=0,keepdims=True),'Fa_N_top','m3')
        f17d = wrap_netcdf(yearnumber,a,np.nansum(Fa_Vert,axis=0,keepdims=True),'Fa_Vert','m3')

        f = cf.FieldList([f0l,f0d,f1l,f1d,f2l,f2d,f3l,f3d,f4l,f4d,f5l,f5d,f6d,f7d,f8d,f9d,f10d,f11d,f12d,f13d,f14d,f15d,f16d,f17d])
        cf.write(f,datapath[2],single=True,unlimited='time',compress=3)

        end = timer()
        print 'Runtime Sa_track for day ' + str(a+1) + ' in year ' + str(yearnumber) + ' is',(end - start),' seconds.'

end1 = timer()
print 'The total runtime of Con_E_Recyc_Masterscript is',(end1-start1),' seconds.'
