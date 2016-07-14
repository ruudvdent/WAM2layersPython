RELEASE NOTES

Main developer of the code Ruud van der Ent (r.j.vanderent@uu.nl). Please let me know when you intend to use this code for research, 
I might be able to offer some kind of support if it is in my interest as well, but I cannot guarantee to do any kind of troubleshooting for free
as I have to work on my own projects as well.

The file "Readme WAM2 - not fully up-to-date.pdf" was written by Tolga Cömert and Niek van de Koppel and provides a general guideline, 
but as the name says is not fully up-to-date with the latest changes from v2.4.01 onward. One thing to clarify is that I consider the definitions:
"get_Sa_track_forward'_TIME'", "get_Sa_track_backward'_TIME'" within the "Masterscripts" to be the CORE module of "WAM-2layers". Everything else is 
either pre-processing or post-processing and merely provided here as an example for the user of the CORE module of "WAM-2layers".

Released under the GNU General Public License. Please cite this Github and 
"Van der Ent, R. J. (2014), A new view on the hydrological cycle over continents, Ph.D. thesis, 96 pp, Delft University of Technology, Delft. 
http://dx.doi.org/10.4233/uuid:0ab824ee-6956-4cc3-b530-3245ab4f32be."

WAM-2layers v2.4.07 | 7-7-2016
- added the Hor_fluxes_output file
- small changes to the example plotting files

WAM-2layers v2.4.06 | 6-7-2016
- solve another bug in backward time tracking
- more accurate area size calculation

WAM-2layers v2.4.05 | 5-7-2016
- solve bug in backward time tracking
- fix the area size of the gridcells at 90deg latitude in getconstants

WAM-2layers v2.4.04
- added the backtracking scripts. 
- The way time-output is analysed may still contain some bugs compared to the Matlab version. Will be updated in future release

WAM-2layers v2.4.03 | 28-6-2016
- added another 'verify_compressed_data_integrity=False' statement. Without, one could encounter errors

WAM-2layers v2.4.02 | 28-6-2016
- fixed a small bug to "lake_mask" so that lakes are identified correctly

WAM-2layers v2.4.01 | 24-6-2016
- Python version of the model with credits to Tolga Cömert and Niek van de Koppel
- forward tracking only, backward tracking will follow soon
- structure of the code is clearer, datapaths and input need to be provided at the top of the scripts and output will follow automatically

WAM-2layers v2.4.00 - non-released version
- Python/Jupyter version of the model by Tolga Cömert and Niek van de Koppel

WAM-2layers v2.3.03 - non-released version
- changed leapyear function to isleap

WAM-2layers v2.3.02 | 2014: latest distributed Matlab version
- changed non problematic mistype in Con_P_Recyc_Output_2layers_sigma.m

WAM-2layers v2.3.01
- Included timetracking (beta)
- Improved alignment of 'end' in a 'for-loop' in the 'Output' files.
- Small improvements in the comment lines

WAM-2layers v2.2.0.1
- Major bugfix in getwind_2layers_sigma.m

-----------------------------------------------------------------------------------------------------------------------------------------------------------
List of papers that use WAM-2layers or earlier versions of the WAM:

Please see:
"Van der Ent, R. J. (2014), A new view on the hydrological cycle over continents, Ph.D. thesis, 96 pp, Delft University of Technology, Delft. 
http://dx.doi.org/10.4233/uuid:0ab824ee-6956-4cc3-b530-3245ab4f32be."

and the upcoming PhD thesis of Patrick W. Keys at Stockholm University


