#!/home/users/lguo/anaconda2/bin/python
#BSUB -o %J.o
#BSUB -e %J.e
#BSUB -q short-serial
#BSUB -W 24:00
##BSUB -R "rusage[mem=16000]"
##BSUB -M 16000
##!/usr/bin/python2.7

import numpy as np
import os
import cf
import pdb

region = 'cn1'
year = '1991'
band = 'b001'

interdata_folder = '/home/users/lguo/perchance/WAM_inter/interdata_erai_'+band
sub_interdata_folder = os.path.join(interdata_folder,region+'_source_forward')
output_folder = '/home/users/lguo/klingaman/output/daily_erai_band'

f = cf.read(sub_interdata_folder+'/'+year+'-???Sa_track.nc')
f.pop(21)
f.pop(19)
f.pop(17)
f.pop(15)
f.pop(13)
f.pop(11)

f_out = output_folder+'/'+region+'_source_forward_'+band+'_'+year+'_daily.nc'
print f_out
cf.write(f,f_out,single=True)

if region == 'cn1':
	E = f.select_field('E',exact=True).collapse('sum',axes='time',group=cf.M())
	P = f.select_field('P',exact=True).collapse('sum',axes='time',group=cf.M())
	W_down = f.select_field('W_down',exact=True).collapse('mean',axes='time',group=cf.M())
	W_top = f.select_field('W_top',exact=True).collapse('mean',axes='time',group=cf.M())

	Fa_E_down = f.select_field('Fa_E_down',exact=True).collapse('sum',axes='time',group=cf.M())
	Fa_E_top = f.select_field('Fa_E_top',exact=True).collapse('sum',axes='time',group=cf.M())
	Fa_N_down = f.select_field('Fa_N_down',exact=True).collapse('sum',axes='time',group=cf.M())
	Fa_N_top = f.select_field('Fa_N_top',exact=True).collapse('sum',axes='time',group=cf.M())
	Fa_Vert = f.select_field('Fa_Vert',exact=True).collapse('sum',axes='time',group=cf.M())

Sa_track_down = f.select_field('Sa_track_down',exact=True).collapse('mean',axes='time',group=cf.M())
Sa_track_top = f.select_field('Sa_track_top',exact=True).collapse('mean',axes='time',group=cf.M())
P_track = f.select_field('P_track',exact=True).collapse('sum',axes='time',group=cf.M())

Sa_time_down = (f.select_field('Sa_time_down',exact=True) * f.select_field('Sa_track_down',exact=True)).collapse('mean',axes='time',group=cf.M()) / Sa_track_down
Sa_time_top = (f.select_field('Sa_time_top',exact=True) * f.select_field('Sa_track_top',exact=True)).collapse('mean',axes='time',group=cf.M()) / Sa_track_top
P_time = (f.select_field('P_time',exact=True) * f.select_field('P_track',exact=True)).collapse('sum',axes='time',group=cf.M()) / P_track
P_time.where(P_time!=P_time,0,i=True)
P_time.setprop('standard_name','P_time')

Sa_dist_down = (f.select_field('Sa_dist_down',exact=True) * f.select_field('Sa_track_down',exact=True)).collapse('mean',axes='time',group=cf.M()) / Sa_track_down
Sa_dist_top = (f.select_field('Sa_dist_top',exact=True) * f.select_field('Sa_track_top',exact=True)).collapse('mean',axes='time',group=cf.M()) / Sa_track_top
P_dist = (f.select_field('P_dist',exact=True) * f.select_field('P_track',exact=True)).collapse('sum',axes='time',group=cf.M()) / P_track
P_dist.where(P_dist!=P_dist,0,i=True)
P_dist.setprop('standard_name','P_dist')

f_out = output_folder+'/'+region+'_source_forward_'+band+'_'+year+'_monthly.nc'
print f_out
if region == 'cn1':
	cf.write([E,P,W_down,W_top,Sa_track_down,Sa_track_top,P_track,Sa_time_down,Sa_time_top,P_time,Sa_dist_down,Sa_dist_top,P_dist,Fa_E_down,Fa_E_top,Fa_N_down,Fa_N_top,Fa_Vert],f_out,single=True)
else:
	cf.write([Sa_track_down,Sa_track_top,P_track,Sa_time_down,Sa_time_top,P_time,Sa_dist_down,Sa_dist_top,P_dist],f_out,single=True)
