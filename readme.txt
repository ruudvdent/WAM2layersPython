RELEASE NOTES

Main developer of the code Ruud van der Ent (r.j.vanderent@uu.nl). Please let me know when you intend to use this code for research, 
I might be able to offer some kind of support if it is in my interest as well, but I cannot guarantee to do any kind of troubleshooting for free.

Released under the GNU General Public License

WAM-2layers v2.4.05 | 6-7-2016
- solve another bug in backward time tracking
- more accurate area size calculation

WAM-2layers v2.4.05 | 5-7-2016
- solve bug in backward time tracking
- fix the area size of the gridcells at 90deg latitude in getconstants

WAM-2layers v2.4.04
- added the backtracking scripts. 
- The way time-output is analysed may still contain some bugs compared to the Matlab version. Will be updated in future release

WAM-2layers v2.4.03 | 28-6-2016
- added another 'verify_compressed_data_integrity=False' statement. Without one could encounter errors

WAM-2layers v2.4.02 | 28-6-2016
- fixed a small bug to "lake_mask" so that lakes are identified correctly

WAM-2layers v2.4.01 | 24-6-2016
- Python version of the model with credits to Tolga Cömert and Niek van de Koppel
- forward tracking only, backward tracking will follow soon
- structure of the code is clearer, datapaths and input need to be provided at the top of the scripts and output will follow automatically

WAM-2layers v2.3.03 - non-released version
- changed leapyear function to isleap

WAM-2layers v2.3.02
- changed non problematic mistype in Con_P_Recyc_Output_2layers_sigma.m

WAM-2layers v2.3.01
- Included timetracking (beta)
- Improved alignment of 'end' in a 'for-loop' in the 'Output' files.
- Small improvements in the comment lines

WAM-2layers v2.2.0.1
- Major bugfix in getwind_2layers_sigma.m
