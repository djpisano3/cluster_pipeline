# CHILES_pipe_QA1.py
# This module is run after the target module in the CHILES pipeline and makes some
# small test cubes to verify that the pipeline calibration and flagging is of
# sufficient quality.  
# 3/3/16 DJP
# 6/28/16, DJP:  Reducing number of iterations to 100 from 1000.  
#                Testing splitting data first. 
# 9/21/16 DJP:  Reduced number of iterations to 1.  As with original code, but running on each spw separately.
# 12/8/16 DJP: Added try, except to spw loop in case one is entirely flagged.
# 12/8/16 DJP: Calculated, print noise for each spectrum.  Changed output layout.
# 3/5/18 DJP: Changed clean to tclean, imaging in topocentric frequency to help with RFI identification.
# 4/22/18 DJP: Changing split to oldsplit
# 11/6/18 DJP: Changed name from testcubes to QA.  Does all imaging/plotting for QA purposes. 
# 12/19/18 DJP:  Including plots showing what's been flagged by masks. 
# 01/11/19 DJP:  Tried running tclean with parallel=True using mpicasa - didn't work for cubes
# 01/13/19 DJP:  Removed parallel=True.  May split QA module into QA1 and QA2 depending on run times
# 01/18/19 DJP:  Had found bug from mpicasa run.  Trying "parallel=True" again.
# 01/21/19 DJP:  Make iterations=2048 for cubes.  Designed to run with mpicasa (should work okay in casa too).
# 01/25/19 DJP:  Removed all phasecal cube imaging/measurements
# 01/27/19 DJP:  Had trouble with plotms and imview running under mpicasa.  Trying with casa and splitting into QA1, QA2
# 08/04/19 DJP:  Fixed html code (for converting to PDF).  Left out use of /dev/shm 

  

logprint ("Starting CHILES_pipe_QA1.py", logfileout='logs/QA1.log')
time_list=runtiming('QA1', 'start')

# If running this module right after starting casa, must run CHILES_pipe_restore first.
# Step 1, load needed functions.
import copy
import numpy as np
import pylab as pylab
import matplotlib.pyplot as plt
import re as re
import sys
import shutil

# Make list of all spws from 0-14
seq=range(15)

# First step involves doing all of the needed runs of "split" for plotting
# purposes.

logprint ('Split calibrated uv data as continuum data', logfileout='logs/QA1.log')

# Flux Calibrator, Frequency-Averaging
ms_name=ms_active[:-3]
output_ms_flux=ms_name+'_flux_averaged.ms'

# Remove old averaged MS if it exists.
if os.path.exists(output_ms_flux):
    os.system("rm -rf "+output_ms_flux)
    
#18: average the data
default('oldsplit')
vis=ms_active
outputvis=output_ms_flux
datacolumn='corrected'
field='1331+305=3C286'
spw='0~14:128~1920'    # Extend the region used for imaging/plots, only excluding the edges.
width='1793'
antenna=''
timebin=''
timerange=''
scan=''
intent=''
array=''
uvrange=''
correlation=''
observation=''
keepflags=False
keepmms=False
oldsplit()

# Phase Calibrator, Frequency-Averaging
ms_name=ms_active[:-3]
output_ms_phase=ms_name+'_phasecal_averaged.ms'

# Remove old averaged MS if it exists.
if os.path.exists(output_ms_phase):
    os.system("rm -rf "+output_ms_phase)
    
#Split the phase cal.
default('oldsplit')
vis=ms_active
outputvis=output_ms_phase
datacolumn='corrected'
field='J0943-0819'
spw='0~14:128~1920'  # Average over all channels, except the very edges
width='1793'
antenna=''
timebin=''
timerange=''
scan=''
intent=''
array=''
uvrange=''
correlation=''
observation=''
keepflags=False
keepmms=False
oldsplit()

#Split the target.
ms_name=ms_active[:-3]
output_ms_target=ms_name+'_target_averaged.ms'

# Remove old averaged MS if it exists.
if os.path.exists(output_ms_target):
    os.system("rm -rf "+output_ms_target)

default('oldsplit')
vis=ms_active
outputvis=output_ms_target
datacolumn='corrected'
field='deepfield'
spw='0~14:128~1920'  # Extend averaging to include all but 128 edge channels
width='1793'
antenna=''
timebin=''
timerange=''
scan=''
intent=''
array=''
uvrange=''
correlation=''
observation=''
keepflags=False
keepmms=False
oldsplit()

fluxms=ms_active[:-3]+'_calibrated_fluxcal'
phasems=ms_active[:-3]+'_calibrated_phasecal'
targetms=ms_active[:-3]+'_calibrated_deepfield'

# Split spectral data for fluxcal, phasecal, and target
for ii in seq:

    logprint ('Split calibrated uv data, Spw='+str(ii), logfileout='logs/QA1.log')

    fluxfile=fluxms+'_spw'+str(ii)+'.ms'
    phasefile=phasems+'_spw'+str(ii)+'.ms'
    targetfile=targetms+'_spw'+str(ii)+'.ms'
    if os.path.exists(fluxfile):
        os.system("rm -rf "+fluxfile)
        os.system("rm -rf fluxcube_spw"+str(ii)+".*")
        os.system("rm -rf images/fluxcube_spw"+str(ii)+".*")
        os.system("rm -rf images/fluxcalibrator_spw"+str(ii)+".*")
    if os.path.exists(phasefile):
        os.system("rm -rf "+phasefile)
        os.system("rm -rf phasecube_spw"+str(ii)+".*")
        os.system("rm -rf images/phasecube_spw"+str(ii)+".*")
        os.system("rm -rf images/phasecalibrator_spw"+str(ii)+".*")
    if os.path.exists(targetfile):
        os.system("rm -rf "+targetfile)
        os.system("rm -rf targetcube_10mJy_spw"+str(ii)+".*")
        os.system("rm -rf targetcube_HIdet_spw"+str(ii)+".*")
        os.system("rm -rf images/targetcube_10mJy_spw"+str(ii)+".*")
        os.system("rm -rf images/targetcube_HIdet_spw"+str(ii)+".*")
        os.system("rm -rf images/target_spw"+str(ii)+".*")

    try:

# Split off flux calibrator
        default('oldsplit')
        vis=ms_active
        datacolumn='corrected'
        outputvis=fluxfile
        field='1331+305=3C286'
        spw =str(ii)
        width=1
        timebin=''
        oldsplit()
  
# Split off phase calibrator
        default('oldsplit')
        vis=ms_active
        datacolumn='corrected'
        outputvis=phasefile
        field='J0943-0819'
        spw =str(ii)
        width=1
        timebin=''
        oldsplit()
                
# Split off deepfield                
        default('oldsplit')
        vis=ms_active
        datacolumn='corrected'
        outputvis=targetfile
        field='deepfield'
        spw =str(ii)
        width=1
        timebin=''
        oldsplit()
    except:
        logprint ('Unable to split uv data, Spw='+str(ii), logfileout='logs/QA1.log')

################################################################################


# Create subdirectory in /dev/shm to put MS files for quick, efficient access

#sharedir='/dev/shm/chiles/'

#if os.path.exists(sharedir):
#    shutil.rmtree(sharedir)

#os.mkdir(sharedir)

# Make plot of short baseline of phase calibrator
# Make plots showing effect of masking RFI on phase calibrator
# The following methold plots the flagged and unflagged data for the shortest 
# baselines used in the calibration.
for ii in seq:
    print('Plot Amplitude vs. Frequency for Spw '+str(ii))
    phasefile=phasems+'_spw'+str(ii)+'.ms'
    p1ax = mask_check_plot_v2(phasefile)
    plt.scatter(p1ax[0], p1ax[1][1][0],color='red',s=0.5)
    plt.scatter(p1ax[0], p1ax[1][1][1],color='red',s=0.5)
    plt.scatter(p1ax[0], p1ax[1][0][0],color='blue',s=0.5)
    plt.scatter(p1ax[0], p1ax[1][0][1],color='blue',s=0.5)
    plt.xlabel('Frequency [GHz]')
    plt.ylabel('Amplitude [Jy]')
    plt.ylim(-0.01,20.0)
    plt.title('Flagged Phase Calibrator Data for Spw'+str(ii)+': 1500-1600m')
    plt.grid()
    plt.savefig('flag1500m_Spw'+str(ii)+'.png')
    plt.close()

os.system('mv flag*.png plots')

# Make plots for Bandpass QA:
#--------------------------------------------------------------------
#Part V: Data inspection
logprint ("Create data inspection plots for Bandpass module", logfileout='logs/QA1.log')
logprint ("Plotting flagging percentage as a function of uvdistance for phase calibrator", logfileout='logs/QA1.log')

# Make plot of flagging statistics
# s_b is the output of flagdata run (above)
# Get information for flagging percentage vs. uvdistance
#gantdata = get_antenna_data(ms_active)
#create adictionary with flagging info
#base_dict = create_baseline_dict(ms_active, gantdata)
#gantdata and base_dict are already in the initial module so no need to retrieve that information again.
#match flagging data to dictionary entry
datamatch = flag_match_baseline(s_b['baseline'], base_dict)
#bin the statistics
binned_stats = bin_statistics(datamatch, 'B', 25)  # 25 is the number of uvdist bins such that there is minimal error in uvdist.

#Plot flagging % vs. uvdist
### Plot the Data
barwidth = binned_stats[0][1]
totflagged = 'BP Flagging: '+ str(flux_flag*100) + '% Data Flagged'
pylab.close()
pylab.bar(binned_stats[0],binned_stats[1], width=barwidth, color='grey', align='edge')
pylab.title(totflagged)
pylab.grid()
pylab.ylabel('Flagged data [%]')
pylab.xlabel('Average UV distance [m]')
pylab.savefig('bp_flag_uvdist.png')
os.system("mv bp_flag_uvdist.png plots/.") 

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
for s in seq:
    for c in range(2048):
        ct+=1
        chan.append(ct)
        freq.append(nu0+dnu*ct)
        flagged.append(s_b['spw:channel'][str(s)+':'+str(c)]['flagged'])
        totals.append(s_b['spw:channel'][str(s)+':'+str(c)]['total'])
        flag_frac.append(flagged[ct]/totals[ct])

fig=pylab.figure()
pylab.plot(freq,flag_frac,'k-')
pylab.xlim(940.,1445.)
pylab.xlabel('Frequency [MHz]')
pylab.ylabel('Fraction of data flagged')
pylab.savefig("bp_flag.png")
pylab.close(fig)

logprint ("Plot delays", logfileout='logs/QA1.log')

#16: plot delays
nplots=int(numAntenna/3)

if ((numAntenna%3)>0):
    nplots = nplots + 1

for ii in range(nplots):
    filename='finaldelay'+str(ii)+'.png'
    syscommand='rm -rf '+filename
    os.system(syscommand)
    
    antPlot=str(ii*3)+'~'+str(ii*3+2)
    
    default('plotcal')
    caltable='finaldelay.k'
    xaxis='freq'
    yaxis='delay'
    poln=''
    field=''
    antenna=antPlot
    spw=''
    timerange=''
    subplot=311
    overplot=False
    clearpanel='Auto'
    iteration='antenna'
    plotrange=[]
    showflags=False
    plotsymbol='o'
    plotcolor='blue'
    markersize=5.0
    fontsize=10.0
    showgui=False
    figfile=filename
    plotcal()

#17: Bandpass solutions
logprint ("Plotting bandpass solutions", logfileout='logs/QA1.log')

tb.open('finalBPcal.b')

# Gets all channel frequencies from BP table
freq_table=tb.taql('select CHAN_FREQ from finalBPcal.b/SPECTRAL_WINDOW')

for ii in seq:
    # Extract frequency for given spw
    freq=freq_table.getcol('CHAN_FREQ',ii,1)[:,0]
    freq/=1e6  # Convert frequencies to MHz
    # Extract BP data for given spectral window
    soltable=tb.taql('select CPARAM from finalBPcal.b where SPECTRAL_WINDOW_ID='+str(ii))
    # Extract complex gain values separately for each corr; returns array of channels, antennas
    g0=soltable.getcol('CPARAM')[0]
    g1=soltable.getcol('CPARAM')[1]
    # Convert complex gain to amplitude & phase
    a0=np.absolute(g0)
    a1=np.absolute(g1)
    p0=np.angle(g0,deg=True)
    p1=np.angle(g1,deg=True)
    # Convert a=1 data and p=0 data to NaNs
    a0[a0==1.]=np.nan
    a1[a1==1.]=np.nan
    p0[p0==0.]=np.nan
    p1[p1==0.]=np.nan
    
    # Plot results for Amplitudes/Gains
    plt.close()
    for ant in range(26):
        plt.plot(freq,a0[:,ant])
        plt.plot(freq,a1[:,ant])
    plt.xlabel('Frequency [MHz]')
    plt.ylabel('Gain')
    plt.title('BP Gain Solutions')
    plt.savefig('bpamp_Spw'+str(ii)+'.png')
    plt.close()
    # Phases
    for ant in range(26):
        plt.plot(freq,p0[:,ant])
        plt.plot(freq,p1[:,ant])
    plt.xlabel('Frequency [MHz]')
    plt.ylabel('Phase')
    plt.title('BP Phase Solutions')
    plt.savefig('bpphase_Spw'+str(ii)+'.png')
    plt.close()
    
tb.close()

#Plot UV spectrum (averaged over all baselines & time) of flux calibrator
logprint ("Plotting UVSPEC of flux calibrator", logfileout='logs/QA1.log')

temp_spec = []
for ii in seq:
    fluxfile=fluxms+'_spw'+str(ii)+'.ms'
    temp_spec.append(uvspec_v2(fluxfile))

# Amplitude Spectrum
fig=plt.figure()
for ii in seq:
    if ii % 2 == 0:
        plt.plot(temp_spec[ii][0],temp_spec[ii][1][0],'.',markersize=0.5)
        plt.plot(temp_spec[ii][0],temp_spec[ii][1][1],'.',markersize=0.5)
    else:
        plt.plot(temp_spec[ii][0],temp_spec[ii][1][0],'.',markersize=0.5)
        plt.plot(temp_spec[ii][0],temp_spec[ii][1][1],'.',markersize=0.5)
plt.xlabel('Frequency [GHz]')
plt.ylabel('Amplitude [Jy]')
plt.xlim(0.95,1.43)
plt.draw()
plt.savefig('fluxcal_amp_spectrum_full.png')
plt.ylim(14.5,18.5)
plt.draw()
plt.savefig('fluxcal_amp_spectrum_zoom.png')
plt.close()

# Phase Spectrum
for ii in seq:
    if ii % 2 == 0:
        plt.plot(temp_spec[ii][0],temp_spec[ii][2][0],'.',markersize=0.5)
        plt.plot(temp_spec[ii][0],temp_spec[ii][2][1],'.',markersize=0.5)
    else:
        plt.plot(temp_spec[ii][0],temp_spec[ii][2][0],'.',markersize=0.5)
        plt.plot(temp_spec[ii][0],temp_spec[ii][2][1],'.',markersize=0.5)

plt.xlabel('Frequency [GHz]')
plt.ylabel('Phase [deg]')
plt.xlim(0.95,1.43)
plt.draw()
plt.savefig('fluxcal_phase_spectrum_full.png')
plt.ylim(-0.5,0.5)
plt.draw()
plt.savefig('fluxcal_phase_spectrum_zoom.png')
plt.close()

# The following method plots the amplitude and phase information as a function 
# of each other, and of time.
logprint ("Plotting Amplitude vs. Phase vs. Time for flux calibrator", logfileout='logs/QA1.log')

for ii in seq:
    # Retrieve Data
    f2ax = amp_phase_time_v2(output_ms_flux,str(ii))
    # Amp vs. Phase
    plt.scatter(f2ax[2][0][0],f2ax[1][0][0],color='blue',s=0.5)
    plt.scatter(f2ax[2][1][0],f2ax[1][1][0],color='red',s=0.5)
    plt.xlabel('Phase [deg]')
    plt.ylabel('Amplitude [Jy]')
    plt.title('Amp. vs. Phase, Spw='+str(ii))
    plt.grid()
    plt.savefig('fluxcal_ampphase_Spw'+str(ii)+'.png')
    plt.close()
    # Amplitude vs. Time
    thr=(f2ax[0]-min(f2ax[0]))*60.
    plt.scatter(f2ax[0],f2ax[1][0][0],color='blue',s=0.5)
    plt.scatter(f2ax[0],f2ax[1][1][0],color='red',s=0.5)
    plt.xlabel('Time [hr]')
    plt.ylabel('Amplitude [Jy]')
    plt.title('Amp. vs. Time, Spw='+str(ii))
    plt.grid()
    plt.savefig('fluxcal_amptime_Spw'+str(ii)+'.png')
    plt.close()
    # Phase vs. Time
    plt.scatter(f2ax[0],f2ax[2][0][0],color='blue',s=0.5)
    plt.scatter(f2ax[0],f2ax[2][1][0],color='red',s=0.5)
    plt.xlabel('Time [hr]')
    plt.ylabel('Phase [deg]')
    plt.title('Phase vs. Time, Spw='+str(ii))
    plt.grid()
    plt.savefig('fluxcal_phasetime_Spw'+str(ii)+'.png')
    plt.close()



#19: 
logprint ("Imaging Flux Calibrator", logfileout='logs/QA1.log')

for ii in seq:
    #19: Imaging flux calibrator in continuum
    print 'STARTS IMAGING FLUX CALIBRATOR OF SPW='+str(ii)
    default('tclean')
    image_name='fluxcalibrator_spw'+str(ii)
    fieldid='1331*'
    grid_mode=''
    number_w=1
    image_size=[512,512]
    iteration=1000
    mask_name=['']
    
    vis=output_ms_flux
    imagename=image_name
    selectdata=True
    datacolumn='data'
    field=fieldid
    spw=str(ii)
    specmode='mfs'
    nterms=1
    niter=iteration
    gain=0.1
    gridmode=grid_mode
    wprojplanes=number_w
    threshold='0.0mJy'
    deconvolver='clark'
    imagermode='csclean'
    cyclefactor=1.5
    cyclespeedup=-1
    multiscale=[]
    interactive=False
    mask=mask_name
    imsize=image_size
    cell=['1.5arcsec','1.5arcsec']
    stokes='I'
    weighting='briggs'
    robust=0.8
    uvtaper=[]
    modelimage=''
    restoringbeam=['']
    pblimit=-0.2
    pbcor=False
    usescratch=False
    allowchunk=False
    async=False
#    parallel=True
    parallel=False   
    tclean()

    
#20:  Calculate bmaj, bmin, peak, and rms for images of flux calibrator in each spw.
#     Output results as text file and plots of these values.
logprint ("Plotting image properties for flux calibrator", logfileout='logs/QA1.log')

box_flux='300,30,460,200'

bmaj_flux=[]
bmin_flux=[]
max_flux=[]
rms_flux=[]

for ii in seq:
    
    image_fluxcal='fluxcalibrator_spw'+str(ii)+'.image'
    bmaj1=imhead(image_fluxcal,mode='get',hdkey='beammajor')
    if bmaj1==None:
        bmaj_flux.append(0.0)
        bmin_flux.append(0.0)
        max_flux.append(0.0)
        rms_flux.append(0.0)
    else:
        bmaj11=bmaj1['value']
        bmaj_flux.append(bmaj11)
        bmin1=imhead(image_fluxcal,mode='get',hdkey='beamminor')
        bmin11=bmin1['value']
        bmin_flux.append(bmin11)
        imst1=imstat(image_fluxcal)
        max1=imst1['max']
        max_flux.append(max1[0])
        imst1=imstat(image_fluxcal,box=box_flux)
        rms1=imst1['rms']
        rms_flux.append(rms1[0]*1e3)

# Output Image statistics in text file
f=open('statisticsFlux.txt','w')

print >> f, "Flux calibrator"
print >> f, "major axis [\"] \t", "".join(["%12.3f \t" %x for x in bmaj_flux])
print >> f, "minor axis [\"] \t", "".join(["%12.3f \t" %x for x in bmin_flux])
print >> f, "peak value [Jy] \t", "".join(["%12.3f \t" %x for x in max_flux])
print >> f, "RMS noise [mJy] \t", "".join(["%12.3f \t" %x for x in rms_flux])
f.close()

#Make plots of bmaj, bmin, peak, and rms
fig=pylab.figure()
pylab.plot(seq,bmaj_flux,'r--x',label='Bmaj')
pylab.plot(seq,bmin_flux,'b--x',label='Bmin')
pylab.xlabel('Spectral Window')
pylab.ylabel('Beam Size ["]')
pylab.legend()
pylab.savefig('fluxcal_beamsize.png')
pylab.close(fig)

fig=pylab.figure()
pylab.plot(seq,max_flux,'k-x')
pylab.xlabel('Spectral Window')
pylab.ylabel('Peak Flux [Jy]')
#pylab.yscale('log')
pylab.savefig('fluxcal_peak.png')
pylab.close(fig)

fig=pylab.figure()
pylab.plot(seq,rms_flux,'k-x')
pylab.xlabel('Spectral Window')
pylab.ylabel('RMS [mJy]')
pylab.yscale('log')
pylab.savefig('fluxcal_rms.png')
pylab.close(fig)

#Move plots, images to sub-directory

os.system("mv *.png plots")
if os.path.exists('images')==False:
    os.system("mkdir images")
os.system("mv fluxcalibrator_spw*.* images")

#Create webpage with BP results

if os.path.exists('bandpass.html'):
    os.system("rm bandpass.html")
wlog = open("bandpass.html","w")
wlog.write('<html>\n')
wlog.write('<head>\n')
wlog.write('<title>CHILES Pipeline Web Log</title>\n')
wlog.write('<style>table tr {page-break-inside: avoid}</style>\n')
wlog.write('</head>\n')
wlog.write('<body>\n')
wlog.write('<br>\n')
wlog.write('<hr>\n')
wlog.write('<center>CHILES_pipe_bandpass results:</center>')
wlog.write('<li> Session: '+SDM_name+'</li>\n')
wlog.write('<li><a href="logs/bandpass.log">Bandpass Log</a></li>\n')
wlog.write('<li>Delay Solutions: \n')
wlog.write('<br><img src="plots/finaldelay0.png">\n')
wlog.write('<br><img src="plots/finaldelay1.png">\n')
wlog.write('<br><img src="plots/finaldelay2.png">\n')
wlog.write('<br><img src="plots/finaldelay3.png">\n')
wlog.write('<br><img src="plots/finaldelay4.png">\n')
wlog.write('<br><img src="plots/finaldelay5.png">\n')
wlog.write('<br><img src="plots/finaldelay6.png">\n')
wlog.write('<br><img src="plots/finaldelay7.png">\n')
wlog.write('<br><img src="plots/finaldelay8.png"></li>\n')
wlog.write('<li> Bandpass solutions (amplitude and phase) for reference antenna: \n')
wlog.write('<li> Color coded by antenna, both polarizations shown \n')
wlog.write('<table> \n')
for ii in seq:
    wlog.write('<tr><td><img src="plots/bpamp_Spw'+str(ii)+'.png"></td>\n')
    wlog.write('<td><img src="plots/bpphase_Spw'+str(ii)+'.png"></td></tr>\n')
wlog.write('</table> \n')
wlog.write('<br>')
wlog.write('<li> Amp vs. Phase (averaged over all channels in a spw): \n')
for ii in seq:
    wlog.write('<br><img src="plots/fluxcal_ampphase_Spw'+str(ii)+'.png">\n')
wlog.write('<br>')
wlog.write('<li> Spectrum of Flux calibrator (both LL & RR, averaged over all time & baselines): \n')
wlog.write('<table> \n')
wlog.write('<tr><td><br><img src="plots/fluxcal_amp_spectrum_full.png"></td>\n')
wlog.write('<td><img src="plots/fluxcal_amp_spectrum_zoom.png"></td></tr>\n')
wlog.write('<tr><td><br><img src="plots/fluxcal_phase_spectrum_full.png"></td>\n')
wlog.write('<td><img src="plots/fluxcal_phase_spectrum_zoom.png"></td></tr>\n')
wlog.write('</table> \n')
wlog.write('<br>')
wlog.write('<li> Amp. & Phase vs. time for Flux Calibrator (averaged over frequency): \n')
wlog.write('<table> \n')
for ii in seq:
    wlog.write('<tr><td><img src="plots/fluxcal_amptime_Spw'+str(ii)+'.png"></td>\n')
    wlog.write('<td><img src="plots/fluxcal_phasetime_Spw'+str(ii)+'.png"></td></tr>\n')
wlog.write('</table> \n')
wlog.write('</li>')
wlog.write('<li> Measured properties of flux calibrator: \n')
wlog.write('<br><img src="plots/fluxcal_beamsize.png">\n')
wlog.write('<br><img src="plots/fluxcal_peak.png">\n')
wlog.write('<br><img src="plots/fluxcal_rms.png"></li>\n')
wlog.write('<br>\n')
wlog.write('<br> Flagging percentage vs. frequency:\n')
wlog.write('<li><br><img src="./plots/bp_flag.png"></li>\n')
wlog.write('<br>\n')
wlog.write('<br> Flagging percentage vs. uvdist\n')
wlog.write('<li><br><img src="./plots/bp_flag_uvdist.png"></li>\n')
wlog.write('<br>\n')
wlog.write('<br> Percentage of 3C286 data flagged: '+str(flux_flag*100)+'\n')
wlog.write('<br>')
wlog.write('<hr>\n')
wlog.write('</body>\n')
wlog.write('</html>\n')
wlog.close()

# Make QA plots for phasecal module
# Step 10: Makes diagnostic plots for assessment
# Look at calibration tables for phase cal (amp, phase vs. time, frequency)
# Make images of phase cal and look at flux,beam vs. spw

logprint ("Making diagnostic plots for phase calibrator", logfileout='logs/QA1.log')

# Make plot of flagging statistics
# s_p is the output of flagdata run (above)
# Get information for flagging percentage vs. uvdistance
#gantdata = get_antenna_data(ms_active)
#create adictionary with flagging info
#base_dict = create_baseline_dict(ms_active, gantdata)
#gantdata and base_dict are already in the initial module so no need to retrieve that information again.
#match flagging data to dictionary entry
datamatch = flag_match_baseline(s_p['baseline'], base_dict)
#bin the statistics
binned_stats = bin_statistics(datamatch, 'B', 25)  # 25 is the number of uvdist bins such that there is minimal error in uvdist.

logprint ("Plot flagging percentage vs. uvdistance for phasecal", logfileout='logs/QA1.log')
#Plot flagging % vs. uvdist
### Plot the Data
barwidth = binned_stats[0][1]
totflagged = 'Phase Cal Flagging: '+ str(phase_flag*100) + '% Data Flagged'
pylab.close()
pylab.bar(binned_stats[0],binned_stats[1], width=barwidth, color='grey', align='edge')
pylab.title(totflagged)
pylab.grid()
pylab.ylabel('flagged data [%]')
pylab.xlabel('average UV distance [m]')
pylab.savefig('phase_flag_uvdist.png')
pylab.close()

os.system("mv phase_flag_uvdist.png plots/.") 

# Make plot of percentage of data flagged as a function of channel (for both correlations combined):
logprint ("Plot flagging percentage vs. channel for phasecal", logfileout='logs/QA1.log')

flag_frac=[]
ct=-1
chan=[]
freq=[]
flagged=[]
totals=[]
# Extract frequency of first channel of spw=0 from listobs output
nu0=reference_frequencies[0]/1.e6 #get reference frequency in MHz
dnu=0.015625 # channel width in MHz
for s in seq:
    for c in range(2048):
        ct+=1
        chan.append(ct)
        freq.append(nu0+dnu*ct)
        flagged.append(s_p['spw:channel'][str(s)+':'+str(c)]['flagged'])
        totals.append(s_p['spw:channel'][str(s)+':'+str(c)]['total'])
        flag_frac.append(flagged[ct]/totals[ct])

fig=pylab.figure()
pylab.plot(freq,flag_frac,'k-')
pylab.xlim(940.,1445.)
pylab.xlabel('Frequency [MHz]')
pylab.ylabel('Fraction of data flagged')
pylab.savefig("phase_flag.png")
pylab.close(fig)

#Plot UV spectrum (averaged over all baselines & time) of phase calibrator
logprint ("Plot UVSPEC of phasecalibrator", logfileout='logs/QA1.log')

temp_spec = []
for ii in seq:
    phasefile=phasems+'_spw'+str(ii)+'.ms'
    temp_spec.append(uvspec_v2(phasefile))

# Amplitude Spectrum
fig=plt.figure()
for ii in seq:
    if ii % 2 == 0:
        plt.scatter(temp_spec[ii][0],temp_spec[ii][1][0],s=0.5)
        plt.scatter(temp_spec[ii][0],temp_spec[ii][1][1],s=0.5)
    else:
        plt.scatter(temp_spec[ii][0],temp_spec[ii][1][0],s=0.5)
        plt.scatter(temp_spec[ii][0],temp_spec[ii][1][1],s=0.5)
plt.xlabel('Frequency [GHz]')
plt.ylabel('Amplitude [Jy]')
plt.xlim(0.95,1.43)
plt.draw()
plt.savefig('phasecal_amp_spectrum_full.png')
plt.ylim(2.7,3.5)
plt.draw()
plt.savefig('phasecal_amp_spectrum_zoom.png')
plt.close()

# Phase Spectrum
for ii in seq:
    if ii % 2 == 0:
        plt.scatter(temp_spec[ii][0],temp_spec[ii][2][0],s=0.5)
        plt.scatter(temp_spec[ii][0],temp_spec[ii][2][1],s=0.5)
    else:
        plt.scatter(temp_spec[ii][0],temp_spec[ii][2][0],s=0.5)
        plt.scatter(temp_spec[ii][0],temp_spec[ii][2][1],s=0.5)

plt.xlabel('Frequency [GHz]')
plt.ylabel('Phase [deg]')
plt.xlim(0.95,1.43)
plt.draw()
plt.savefig('phasecal_phase_spectrum_full.png')
plt.ylim(-0.5,0.5)
plt.draw()
plt.savefig('phasecal_phase_spectrum_zoom.png')
plt.close()


# The following method plots the amplitude and phase information as a function 
# of each other, and of time.
logprint ("Plot Amplitude vs. Phase vs. Time for phasecal", logfileout='logs/QA1.log')

for ii in seq:
    # Retrieve Data
    p2ax = amp_phase_time_v2(output_ms_phase,str(ii))
    # Amp vs. Phase
    plt.scatter(p2ax[2][0][0],p2ax[1][0][0],color='blue',s=0.5)
    plt.scatter(p2ax[2][1][0],p2ax[1][1][0],color='red',s=0.5)
    plt.xlabel('Phase [deg]')
    plt.ylabel('Amplitude [Jy]')
    plt.title('Amp. vs. Phase, Spw='+str(ii))
    plt.grid()
    plt.savefig('phasecal_ampphase_Spw'+str(ii)+'.png')
    plt.close()
    # Amplitude vs. Time
    plt.scatter(p2ax[0],p2ax[1][0][0],color='blue',s=0.5)
    plt.scatter(p2ax[0],p2ax[1][1][0],color='red',s=0.5)
    plt.xlabel('Time [hr]')
    plt.ylabel('Amplitude [Jy]')
    plt.title('Amp. vs. Time, Spw='+str(ii))
    plt.grid()
    plt.savefig('phasecal_amptime_Spw'+str(ii)+'.png')
    plt.close()
    # Phase vs. Time
    plt.scatter(p2ax[0],p2ax[2][0][0],color='blue',s=0.5)
    plt.scatter(p2ax[0],p2ax[2][1][0],color='red',s=0.5)
    plt.xlabel('Time [hr]')
    plt.ylabel('Phase [deg]')
    plt.title('Phase vs. Time, Spw='+str(ii))
    plt.grid()
    plt.savefig('phasecal_phasetime_Spw'+str(ii)+'.png')
    plt.close()


#Image phase cal: 
for ii in seq:
    print 'STARTS IMAGING PHASE CALIBRATOR OF SPW='+str(ii)
    default('tclean')
    image_name='phasecalibrator_spw'+str(ii)
    fieldid='J0943*'
    grid_mode=''
    number_w=1
    image_size=[512,512]
    iteration=1000
    mask_name=['']
    
    vis=output_ms_phase
    imagename=image_name
    selectdata=True
    datacolumn='data'
    field=fieldid
    spw=str(ii)
    specmode='mfs'
    nterms=1
    niter=iteration
    gain=0.1
    gridmode=grid_mode
    wprojplanes=number_w
    threshold='0.0mJy'
    deconvolver='clark'
    imagermode='csclean'
    cyclefactor=1.5
    cyclespeedup=-1
    multiscale=[]
    interactive=False
    mask=mask_name
    imsize=image_size
    cell=['1.5arcsec','1.5arcsec']
    stokes='I'
    weighting='briggs'
    robust=0.8
    uvtaper=[]
    modelimage=''
    restoringbeam=['']
    pblimit=-0.2
    pbcor=False
    usescratch=False
    allowchunk=False
    async=False
#    parallel=True
    parallel=False
    tclean()
    

logprint ("Measure properties of phasecal continuum images", logfileout='logs/QA1.log')

#     Calculate bmaj, bmin, peak, and rms for images of flux calibrator in each spw.
#     Output results as text file and plots of these values.
box_phase='300,30,460,200'

bmaj_phase=[]
bmin_phase=[]
max_phase=[]
rms_phase=[]

for ii in seq:
    
    image_phasecal='phasecalibrator_spw'+str(ii)+'.image'
    bmaj1=imhead(image_phasecal,mode='get',hdkey='beammajor')
    if bmaj1==None:
        bmaj_phase.append(0.0)
        bmin_phase.append(0.0)
        max_phase.append(0.0)
        rms_phase.append(0.0)
    else:
        bmaj11=bmaj1['value']
        bmaj_phase.append(bmaj11)
        bmin1=imhead(image_phasecal,mode='get',hdkey='beamminor')
        bmin11=bmin1['value']
        bmin_phase.append(bmin11)
        imst1=imstat(image_phasecal)
        max1=imst1['max']
        max_phase.append(max1[0])
        imst1=imstat(image_phasecal,box=box_phase)
        rms1=imst1['rms']
        rms_phase.append(rms1[0]*1e3)

# Output Image statistics in text file
f=open('statisticsPhase.txt','w')

print >> f, "Phase calibrator"
print >> f, "major axis [\"] \t", "".join(["%12.3f \t" %x for x in bmaj_phase])
print >> f, "minor axis [\"] \t", "".join(["%12.3f \t" %x for x in bmin_phase])
print >> f, "peak value [Jy] \t", "".join(["%12.3f \t" %x for x in max_phase])
print >> f, "RMS noise [mJy] \t", "".join(["%12.3f \t" %x for x in rms_phase])
f.close()

#Make plots of bmaj, bmin, peak, and rms
fig=pylab.figure()
pylab.plot(seq,bmaj_phase,'r--x',label='Bmaj')
pylab.plot(seq,bmin_phase,'b--x',label='Bmin')
pylab.xlabel('Spectral Window')
pylab.ylabel('Beam Size ["]')
pylab.legend()
pylab.savefig('phasecal_beamsize.png')
pylab.close(fig)

fig=pylab.figure()
pylab.plot(seq,max_phase,'k-x')
pylab.xlabel('Spectral Window')
pylab.ylabel('Peak Flux [Jy]')
#pylab.yscale('log')
pylab.savefig('phasecal_peak.png')
pylab.close(fig)

fig=pylab.figure()
pylab.plot(seq,rms_phase,'k-x')
pylab.xlabel('Spectral Window')
pylab.ylabel('RMS [mJy]')
pylab.yscale('log')
pylab.savefig('phasecal_rms.png')
pylab.close(fig)

logprint ("Plot calibration tables", logfileout='logs/QA1.log')
#Plot calibration tables
default('plotcal')
caltable='finalamp.gcal'
xaxis='time'
yaxis='amp'
showgui=False
figfile='caltable_finalamp_amp.png'
plotcal()

yaxis='phase'
figfile='caltable_finalamp_phase.png'
plotcal()

#Move plots, images to sub-directory

os.system("mv *.png plots")
if os.path.exists('images')==False:
    os.system("mkdir images")
os.system("mv phasecalibrator_spw*.* images")

#Create webpage with phasecal results

if os.path.exists('phasecal.html'):
    os.system("rm phasecal.html")
wlog = open("phasecal.html","w")
wlog.write('<html>\n')
wlog.write('<head>\n')
wlog.write('<title>CHILES Pipeline Web Log</title>\n')
wlog.write('<style>table tr {page-break-inside: avoid}</style>\n')
wlog.write('</head>\n')
wlog.write('<body>\n')
wlog.write('<br>\n')
wlog.write('<hr>\n')
wlog.write('<center>CHILES_pipe_phasecal results:</center>')
wlog.write('<li> Session: '+SDM_name+'</li>\n')
wlog.write('<li><a href="logs/phasecal.log">Phasecal Log</a></li>\n')
wlog.write('<li>finalamp.gcal: Amp vs. Time: \n')
wlog.write('<br><img src="plots/caltable_finalamp_amp.png"></li>\n')
wlog.write('<li>finalamp.gcal: Phase vs. Time: \n')
wlog.write('<br><img src="plots/caltable_finalamp_phase.png"></li>\n')
wlog.write('<li> Amp vs. Phase: \n')
for ii in seq:
    wlog.write('<br><img src="plots/phasecal_ampphase_Spw'+str(ii)+'.png">\n')
wlog.write('<br>')
wlog.write('<li> Spectrum of Phase calibrator (both LL & RR, averaged over all time & baselines): \n')
wlog.write('<table> \n')
wlog.write('<tr><td><br><img src="plots/phasecal_amp_spectrum_full.png"></td>\n')
wlog.write('<td><img src="plots/phasecal_amp_spectrum_zoom.png"></td></tr>\n')
wlog.write('<tr><td><br><img src="plots/phasecal_phase_spectrum_full.png"></td>\n')
wlog.write('<td><img src="plots/phasecal_phase_spectrum_zoom.png"></td></tr>\n')
wlog.write('</table> \n')
wlog.write('<br>')
wlog.write('<li> Amp. & Phase vs. time for Phase Calibrator (averaged over frequency): \n')
wlog.write('<table> \n')
for ii in seq:
    wlog.write('<tr>\n')
    wlog.write('<td><img src="plots/phasecal_amptime_Spw'+str(ii)+'.png"></td>\n')
    wlog.write('<td><img src="plots/phasecal_phasetime_Spw'+str(ii)+'.png"></td>\n')
    wlog.write('</tr>\n')
wlog.write('</table> \n')
wlog.write('</li>')
wlog.write('<li> Measured properties of phase calibrator: \n')
wlog.write('<br><img src="plots/phasecal_beamsize.png">\n')
wlog.write('<br><img src="plots/phasecal_peak.png">\n')
wlog.write('<br><img src="plots/phasecal_rms.png"></li>\n')
wlog.write('<br>\n')
wlog.write('<br> Flagging percentage vs. frequency:\n')
wlog.write('<li><br><img src="./plots/phase_flag.png"></li>\n')
wlog.write('<br>\n')
wlog.write('<br> Flagging percentage vs. uvdist\n')
wlog.write('<li><br><img src="./plots/phase_flag_uvdist.png"></li>\n')
wlog.write('<br>\n')
wlog.write('<br> Percentage of J0943 data flagged: '+str(phase_flag*100)+'\n')
wlog.write('<br>')
wlog.write('<hr>\n')
wlog.write('</body>\n')
wlog.write('</html>\n')
wlog.close()

# Make QA plots for target module

# Step 3: Diagnostic Plots
logprint("Making diagnostic plots for deepfield",logfileout='logs/QA1.log')

logprint ("Plot flagging percentage vs. frequency for deepfield", logfileout='logs/QA1.log')

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

for s in seq:
    for c in range(2048):
        ct+=1
        chan.append(ct)
        freq.append(nu0+dnu*ct)
        flagged.append(s_t['spw:channel'][str(s)+':'+str(c)]['flagged'])
        totals.append(s_t['spw:channel'][str(s)+':'+str(c)]['total'])
        flag_frac.append(flagged[ct]/totals[ct])

fig=pylab.figure()
pylab.plot(freq,flag_frac,'k-')
pylab.xlim(940.,1445.)
pylab.xlabel('Frequency [MHz]')
pylab.ylabel('Fraction of data flagged')
pylab.savefig("target_flag.png")
pylab.close(fig)

# Use plotms to make plot of amplitude vs. frequency, divided by spw

for ii in seq:
    logprint ("Plot Amplitude vs. Frequency/time for deepfield, Spw "+str(ii), logfileout='logs/QA1.log')
    
    targetfile=targetms+'_spw'+str(ii)+'.ms'

#    shutil.copytree(targetfile,sharedir+targetfile)

#    tt = target_plots_v2(sharedir+targetfile,ii)
    tt = target_plots_v2(targetfile,ii)

    fig, (ax1, ax2, ax3) = plt.subplots(3, sharex = True)
    ax1.set_title('Deepfield : U-V Data Mean/Max and standard deviation')
    ax1.plot(tt[0][:,0], tt[0][:,3],color='red',linewidth=1.0)
    ax1.set_ylabel('Max Amplitude [Jy]')
    ax1.grid()
    ax2.plot(tt[0][:,0], tt[0][:,1],color='black',linewidth=1.0)
    ax2.set_ylabel('Mean Amplitude [Jy]')
    ax2.set_xlabel('Scan number')
    ax2.grid()
    ax3.plot(tt[0][:,0], tt[0][:,2],color='green',linewidth=1.0)
    ax3.set_ylabel('Std Dev [Jy]')
    ax3.set_xlabel('Scan number')
    ax3.grid()
    fig.subplots_adjust(hspace=0)
    plt.savefig('target_amptime_Spw'+str(ii)+'.png')
    plt.close()

    fig, (ax1, ax2, ax3) = plt.subplots(3, sharex = True)
    ax1.set_title('Deepfield: U-V Data Mean/Max and standard deviation')
    ax1.plot(tt[1][:,0], tt[1][:,3],color='red',linewidth=1.0)
    ax1.set_ylabel('Max Amplitude [Jy]')
    ax1.ticklabel_format(useOffset=False)
    ax1.grid()
    ax2.plot(tt[1][:,0], tt[1][:,1],color='black',linewidth=1.0)
    ax2.set_ylabel('Mean [Jy]')
    ax2.set_xlabel('Freq [MHz]')
    ax2.grid()
    ax3.plot(tt[1][:,0], tt[1][:,2],color='green',linewidth=1.0)
    ax3.set_ylabel('Std Dev [Jy]')
    ax3.set_xlabel('Freq [MHz]')
    ax3.grid()
    fig.subplots_adjust(hspace=0)
    ax3.ticklabel_format(useOffset=False)
    plt.savefig('target_ampfreq_Spw'+str(ii)+'.png')
    plt.close()

#    shutil.rmtree(sharedir+targetfile)
    

# Remove shared subdirectory 
#shutil.rmtree(sharedir)


#Image target: 
for ii in seq:
    print 'STARTS IMAGING Deepfield OF SPW='+str(ii)
    default('tclean')
    image_name='target_spw'+str(ii)
    fieldid='deepfield'
    grid_mode=''
    number_w=1
    image_size=[2048,2048]
    iteration=1000
    mask_name=['']
    vis=output_ms_target
    imagename=image_name
    selectdata=True
    datacolumn='data'
    field=fieldid
    spw=str(ii)
    specmode='mfs'
    nterms=1
    niter=iteration
    gain=0.1
    gridmode=grid_mode
    wprojplanes=number_w
    threshold='0.0mJy'
    deconvolver='clark'
    imagermode='csclean'
    cyclefactor=1.5
    cyclespeedup=-1
    multiscale=[]
    interactive=False
    mask=mask_name
    imsize=image_size
    cell=['1.5arcsec','1.5arcsec']
    stokes='I'
    weighting='briggs'
    robust=0.8
    uvtaper=[]
    modelimage=''
    restoringbeam=['']
    pblimit=-0.2
    pbcor=False
    usescratch=False
    allowchunk=False
    async=False
#    parallel=True
    parallel=False
    tclean()

# Measure statistics of deepfield image:
box_target='1300,1100,1900,1600'
bmaj_target=[]
bmin_target=[]
max_target=[]
rms_target=[]


logprint("START READING STATISTICS FROM IMHEAD() AND IMSTAT()",logfileout='logs/QA1.log')
for ii in seq:
	image_target='target_spw'+str(ii)+'.image'
	bmaj1=imhead(image_target,mode='get',hdkey='beammajor')
	if bmaj1==None :
		bmaj_target.append(0.0)
		bmin_target.append(0.0)
		max_target.append(0.0)
		rms_target.append(0.0)
	else:
		bmaj11=bmaj1['value']
		bmaj_target.append(bmaj11)
		bmin1=imhead(image_target,mode='get',hdkey='beamminor')
		bmin11=bmin1['value']
		bmin_target.append(bmin11)
		imst1=imstat(image_target)
		max1=imst1['max']
		max_target.append(max1[0]*1e3)
		imst1=imstat(image_target,box=box_target)
		rms1=imst1['rms']
		rms_target.append(rms1[0]*1e6)

logprint("FINISH READING STATISTICS FROM IMHEAD() AND IMSTAT()",logfileout='logs/QA1.log')

f=open('statisticsTarget.txt','w')
print >> f, "Deep field"
print >> f, "major axis [\"] \t", "".join(["%12.3f \t" %x for x in bmaj_target])
print >> f, "minor axis [\"] \t", "".join(["%12.3f \t" %x for x in bmin_target])
print >> f, "peak value [mJy] \t", "".join(["%12.3f \t" %x for x in max_target])
print >> f, "RMS noise [uJy] \t", "".join(["%12.3f \t" %x for x in rms_target])

f.close()

logprint ("Plot continuum image properties for deepfield", logfileout='logs/QA1.log')

#Make plots of bmaj, bmin, peak, and rms
fig=pylab.figure()
pylab.plot(seq,bmaj_target,'r--x',label='Bmaj')
pylab.plot(seq,bmin_target,'b--x',label='Bmin')
pylab.xlabel('Spectral Window')
pylab.ylabel('Beam Size ["]')
pylab.legend()
pylab.savefig('target_beamsize.png')
pylab.close(fig)

fig=pylab.figure()
pylab.plot(seq,max_target,'k-x')
pylab.xlabel('Spectral Window')
pylab.ylabel('Peak Flux [mJy]')
#pylab.yscale('log')
pylab.savefig('target_peak.png')
pylab.close(fig)

fig=pylab.figure()
pylab.plot(seq,rms_target,'k-x')
pylab.xlabel('Spectral Window')
pylab.ylabel('RMS [uJy]')
pylab.yscale('log')
pylab.savefig('target_rms.png')
pylab.close(fig)

#Want to plot image of deepfield in each spw.  Use "imview"

for ii in seq:
    image_target='target_spw'+str(ii)+'.image'
    kntr_levels=[-2*rms_target[ii]/1000.,2*rms_target[ii]/1000.,0.1*max_target[ii],0.3*max_target[ii],0.5*max_target[ii],0.7*max_target[ii],0.9*max_target[ii]]
    imview(raster={'file':image_target, 'colorwedge':True},contour={'file':image_target,'levels':kntr_levels},out='target_spw'+str(ii)+'.png')


#Move plots, images to sub-directory

os.system("mv *.png plots")
os.system("mv target_spw*.* images")

#Create webpage with results

if os.path.exists("target.html"):
    os.system("rm target.html")
wlog = open("target.html","w")
wlog.write('<html>\n')
wlog.write('<head>\n')
wlog.write('<title>CHILES Pipeline Web Log</title>\n')
wlog.write('<style>table tr {page-break-inside: avoid}</style>\n')
wlog.write('</head>\n')
wlog.write('<body>\n')
wlog.write('<br>\n')
wlog.write('<hr>\n')
wlog.write('<center>CHILES_pipe_target results:</center>')
wlog.write('<li> Session: '+SDM_name+'</li>\n')
wlog.write('<li><a href="logs/target.log">Target Log</a></li>\n')
wlog.write('<li> Amp vs. Frequency (time-averaged, all channels shown) and Amp vs. Time (all channels averaged): \n')
wlog.write('<table> \n')
for ii in seq:
    wlog.write('<tr>\n')
    wlog.write('<td><img src="plots/target_ampfreq_Spw'+str(ii)+'.png"></td>\n')
    wlog.write('<td><img src="plots/target_amptime_Spw'+str(ii)+'.png"></td>\n')
    wlog.write('</tr>\n')
wlog.write('</table> \n')
# wlog.write('<li> Spectrum of Deepfield (both LL & RR, averaged over all time & baselines): \n')
# wlog.write('<br><img src="plots/target_spectrum_full.png">\n')
# wlog.write('<br><img src="plots/target_spectrum_zoom.png">\n')
wlog.write('<li> Images of Deepfield: \n')
for ii in seq:
    wlog.write('<br><img src="plots/target_spw'+str(ii)+'.png">\n')
wlog.write('</li>')
wlog.write('<li> Measured properties of deepfield: \n')
wlog.write('<br><img src="plots/target_beamsize.png">\n')
wlog.write('<br><img src="plots/target_peak.png">\n')
wlog.write('<br><img src="plots/target_rms.png"></li>\n')
wlog.write('<br>\n')
wlog.write('<br> Percentage of deepfield data flagged: '+str(target_flag*100)+'\n')
wlog.write('<br> Flagging percentage vs. frequency (before removing channels that are more than 90% flagged)\n')
wlog.write('<br><img src="plots/target_flag.png">\n')
wlog.write('<br>')
wlog.write('<hr>\n')
wlog.write('</body>\n')
wlog.write('</html>\n')
wlog.close()


os.system('cp *.html FINAL/.')

logprint ("Finished CHILES_pipe_QA1.py", logfileout='logs/QA1.log')
time_list=runtiming('QA1', 'end')

pipeline_save()
