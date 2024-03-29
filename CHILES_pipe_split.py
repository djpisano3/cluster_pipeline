# CHILES_pipe_split.py
# This module is run after the user is satisfied with the results of the CHILES 
# pipeline.
# 9/21/16 DJP
# 1/10/17 DJP: Include html, logs, and plots in FINAL directory
# 10/16/19 DJP:  Included code to convert html to PDF and copy these files to FINAL directory.  WARNING: THIS ONLY WORKS on Spruce KNOB!!!
# 08/06/21 DJP:  Added code to make plot only for spws with data (not always 15 spws).

logprint ("Starting CHILES_pipe_split.py", logfileout='logs/split.log')
time_list=runtiming('split', 'start')

# If running this module right after starting casa, must run CHILES_pipe_restore first.
# Step 1, load needed functions.
import copy
import numpy as np
import pylab as pylab
import re as re
import sys

# Set up directory to save results
if os.path.exists("FINAL")==False:
    os.system("mkdir FINAL")

logprint ('Split calibrated deepfield uv data', logfileout='logs/split.log')

# Do final flag summary and compile statistics.
# Summary of flagging, after final flagging (for testing purposes only)
logprint ("Summary of flags after flagging target", logfileout='logs/split.log')
default('flagdata')
vis=ms_active
mode='summary'
spw='0~14'
field='deepfield'           # Only flagged deepfield, no need to check others here.
correlation='RR,LL'
spwchan=True
spwcorr=True
basecnt = True
action='calculate'
s_t=flagdata()

target_flag=s_t['flagged']/s_t['total']

logprint ("Final percentage of all data flagged: "+str(target_flag*100)+'%', logfileout='logs/split.log')

# Save final version of flags
logprint("Saving flags before split",logfileout='logs/target.log')

default('flagmanager')
vis=ms_active
mode='save'
versionname='finalflags'
comment='Final flags saved before split'
merge='replace'
flagmanager()
logprint ("Flag column saved to "+versionname, logfileout='logs/split.log')


# Step 3: Diagnostic Plots
logprint("Making diagnostic plots",logfileout='logs/split.log')

# Make plot of flagging statistics
# s_t is the output of flagdata run (above)
# Get information for flagging percentage vs. uvdistance
#gantdata = get_antenna_data(ms_active)
#create adictionary with flagging info
#base_dict = create_baseline_dict(ms_active, gantdata)
#gantdata and base_dict are already in the initial module so no need to retrieve that information again.
#match flagging data to dictionary entry
datamatch = flag_match_baseline(s_t['baseline'], base_dict)
#bin the statistics
binned_stats = bin_statistics(datamatch, 'B', 25)  # 25 is the number of uvdist bins such that there is minimal error in uvdist.

#Plot flagging % vs. uvdist
### Plot the Data
barwidth = binned_stats[0][1]
totflagged = 'Target Flagging: '+ str(target_flag*100) + '% Data Flagged'
pylab.close()
pylab.bar(binned_stats[0],binned_stats[1], width=barwidth, color='grey', align='edge')
pylab.title(totflagged)
pylab.grid()
pylab.ylabel('flagged data [%]')
pylab.xlabel('average UV distance [m]')
pylab.savefig('finaltarget_flag_uvdist.png')
pylab.close()
os.system("mv finaltarget_flag_uvdist.png plots/.") 

# Make plot of percentage of data flagged as a function of channel (for both correlations combined):
flag_frac=[]
ct=-1
chan=[]
freq=[]
flagged=[]
totals=[]
# Extract frequency of first channel of spw=0 from listobs output
nu0=reference_frequencies[0]/1.e6 #get reference frequency in MHz
dnu=0.015625 # channel width in MHz
freq=[]

###------------------------------ NL Addition -------------------------------###
# Make list of all spws with valid data
numalllinespw = len(s_t['spw']) # Number of line spws, either 14 or 15
seq_list = []
# Loop through the possible spectral windows.
for ii in range(0,numalllinespw):
    perc = s_t['spw'][str(ii)]['flagged']/s_t['spw'][str(ii)]['total']
    if perc < 1.: # Only retrieve the unflagged ones
        seq_list.append(int(ii))
# Make seq the only good values.
seq = np.array(seq_list)
###------------------------------ NL Addition -------------------------------###


for s in seq:
    for c in range(2048):
        ct+=1
        chan.append(ct)
        freq.append(nu0+dnu*ct)
        flagged.append(s_t['spw:channel'][str(s)+':'+str(c)]['flagged'])
        totals.append(s_t['spw:channel'][str(s)+':'+str(c)]['total'])
        flag_frac.append(flagged[ct]/totals[ct])

# Make updated plot of fraction of data flagged.
fig=pylab.figure()
pylab.plot(freq,flag_frac,'k-')
pylab.xlim(940.,1445.)
pylab.ylim(0.,1.)
pylab.xlabel('Frequency [MHz]')
pylab.ylabel('Fraction of data flagged')
pylab.savefig("finaltarget_flag.png")
pylab.close(fig)


# Write flagging statistics to file
flagstatname=msname.rstrip('ms') + 'flag_stats.txt'

np.savetxt(flagstatname,np.column_stack((freq,flagged,totals,flag_frac)),fmt='%s',header='freq [MHz] #flags #total %flagged')


# Split with no averaging.
outputms=ms_active[:-3]+'_calibrated_deepfield.ms'
targetfile=outputms

# Delete final MS if already present
if os.path.exists(targetfile):
    os.system("rm -rf "+targetfile)
    os.system("rm -rf FINAL/"+targetfile)

# Split full spectral, temporal resolution of Target:        
default('oldsplit')
vis=ms_active
datacolumn='corrected'
outputvis=targetfile
field='deepfield'
spw ='0~14'
width=1
timebin=''
oldsplit()

os.system("mv "+targetfile+" FINAL/")

# Split with temporal and spectral smoothing.  
outputms=ms_active[:-3]+'_calibrated_deepfield_SMOOTH.ms'
targetfile=outputms

# Delete final MS if already present
if os.path.exists(targetfile):
    os.system("rm -rf "+targetfile)
    os.system("rm -rf FINAL/"+targetfile)

# Smooth to 16s time resolution and 62.4 kHz channels       
default('oldsplit')
vis=ms_active
datacolumn='corrected'
outputvis=targetfile
field='deepfield'
spw ='0~14'
width=4
timebin='16s'
oldsplit()

os.system("mv "+targetfile+" FINAL/")


# Save calibration tables
if os.path.exists('antposcal.p')==True:
    shutil.move('antposcal.p','FINAL')

shutil.move('gain_curves.g','FINAL')
shutil.move('finaldelay.k','FINAL')
shutil.move('finalBPcal.b','FINAL')
shutil.move('finalphase_scan.gcal','FINAL')
shutil.move('finalphase_int.gcal','FINAL')
shutil.move('finalamp.gcal','FINAL')
shutil.move('finalflux.gcal','FINAL')

# Save flagging commands & statistics

os.system("cp manual_flags_*.txt FINAL/.")
os.system("cp *flag_stats.txt FINAL/.")

# Save all Flagversions tables

shutil.move(ms_active+'.flagversions','FINAL')

#Move plots, images to sub-directory, and save those.  

os.system("mv *.png plots")

#Create webpage with results

if os.path.exists("split.html"):
    os.system("rm split.html")
wlog = open("split.html","w")
wlog.write('<html>\n')
wlog.write('<head>\n')
wlog.write('<title>CHILES Pipeline Web Log</title>\n')
wlog.write('</head>\n')
wlog.write('<body>\n')
wlog.write('<br>\n')
wlog.write('<hr>\n')
wlog.write('<center>CHILES_pipe_split results:</center>')
wlog.write('<li> Session: '+SDM_name+'</li>\n')
wlog.write('<li><a href="logs/split.log">Split Log</a></li>\n')
wlog.write('<li> Percentage of data flagged as a function of frequency: \n')
wlog.write('<br><img src="plots/finaltarget_flag.png">\n')
wlog.write('<li> Percentage of data flagged as a function of uvdistance: \n')
wlog.write('<br><img src="plots/finaltarget_flag_uvdist.png">\n')
wlog.write('<br>')
wlog.write('<li><a href="logs/timing.log">Timing Log</a></li>\n')
wlog.write('<hr>\n')
wlog.write('</body>\n')
wlog.write('</html>\n')
wlog.close()


logprint("Convert html files to PDF files",logfileout='logs/split.log')
# Convert html files into PDF files using wkhtmltopdf (on Spruce Knob only!)
for i in ('initial','bandpass','phasecal','target','QA','split'):
    cmd_string='singularity exec /shared/software/containers/wkhtmltox.simg wkhtmltopdf --zoom 0.6 --page-size letter '+i+'.html '+i+'.pdf'
#    cmd_string='wkhtmltopdf --zoom 0.6 --page-size letter '+i+'.html '+i+'.pdf'  # If not on Spruce Knob, assuming wkhtmltopdf is in path
    os.system(cmd_string)
    shutil.move(i+'.pdf','FINAL')

# Copy html files and plots to FINAL directory
#os.system("cp -r *.pdf FINAL/.")
shutil.move('plots','FINAL/.')

logprint ("Finished CHILES_pipe_split.py", logfileout='logs/split.log')
time_list=runtiming('split', 'end')

pipeline_save()

# Save variable values
shutil.move('pipeline_shelf.restore','FINAL')

# Copy logs to FINAL directory (needs to be after final "save" to preserve all information
os.system("cp *.log logs")
shutil.move('logs','FINAL')
