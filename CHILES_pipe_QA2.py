# CHILES_pipe_QA2.py
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
# 08/02/19 DJP:  Including use of /dev/shm to avoid large amounts of I/O.  Hopefully it runs faster too.  
# 08/03/21 DJP:  Added code to only loop through spws that have valid data.

  

logprint ("Starting CHILES_pipe_QA2.py", logfileout='logs/QA2.log')
time_list=runtiming('QA2', 'start')

# If running this module right after starting casa, must run CHILES_pipe_restore first.
# Step 1, load needed functions.
import copy
import numpy as np
import pylab as pylab
import matplotlib.pyplot as plt
import re as re
import sys
import shutil

# Define some variables
ms_name=ms_active[:-3]
targetms=ms_active[:-3]+'_calibrated_deepfield'

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


# Make cubes for QA purposes
# Split data with no averaging (will slow down imaging, but avoids other issues)
# Skipping imaging cube of flux calibrator 

logprint('Making QA cubes of target field', logfileout='logs/QA2.log')

# Create subdirectory in /dev/shm to put MS files for quick, efficient access

# sharedir='/dev/shm/chiles/'
# 
# if os.path.exists(sharedir):
#     shutil.rmtree(sharedir)
# 
# os.mkdir(sharedir)

# Now make cubes
for ii in seq:
    logprint ('Image calibrated uv data, Spw='+str(ii), logfileout='logs/QA2.log')

    try:



# Make test cube of 10 mJy source, all channels
        logprint ('Make cube of 10 mJy source, all channels, Spw='+str(ii), logfileout='logs/QA2.log')

        default('tclean')
    
        targetfile=targetms+'_spw'+str(ii)+'.ms'
        image_name='targetcube_10mJy_spw'+str(ii)
        fieldid='deepfield'
        phasecenter='J2000 10h01m31.4 +02d26m40'
        iteration=2048   # Total number of iterations for the entire cube (per spw), ~1 per channel 
        grid_mode=''
        number_w=1
        mask_name=['']    
        image_size=[64,64]

        #shutil.copytree(targetfile,sharedir+targetfile)

        #vis=sharedir+targetfile
        vis=targetfile
        imagename=image_name
        selectdata=False
        field=fieldid
        spw=''
        specmode='cubedata'
        nterms=1
        niter=iteration
        gain=0.1
        gridmode=grid_mode
        wprojplanes=number_w
        threshold='10.0mJy'  # Adding threshold
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
        parallel=True
#        parallel=False        
        tclean()

# Make small cube of lowest z detection

        default('tclean')
    
        image_name='targetcube_HIdet_spw'+str(ii)
        fieldid='deepfield'
        phasecenter='J2000 10h01m15.2 +02d18m24'
        iteration=2048   # Total number of iterations for the entire spw cube, ~1 per channel, only used for HI source spw 
        
        #vis=sharedir+targetfile
        vis=targetfile
        imagename=image_name
        selectdata=False
        field=fieldid
        spw=''
        specmode='cubedata'
        nterms=1
        niter=iteration
        gain=0.1
        gridmode=grid_mode
        wprojplanes=number_w
        threshold='10.0mJy'
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
        parallel=True
 #       parallel=False
# Only run if spw=13
        if ii==13:
            logprint ('Make cube of strongest pilot detection, Spw='+str(ii), logfileout='logs/QA2.log')
            tclean()
          
    except:
        logprint("Unable to image "+str(field)+", spw "+str(ii), logfileout='logs/QA2.log')

    #shutil.rmtree(sharedir+targetfile)

logprint("Finished making all cubes", logfileout='logs/QA2.log')

# Remove /dev/shm/ subdirectory since we don't need it anymore.  
#shutil.rmtree(sharedir)



# Prepare to extract information about cubes.
# Initialize variables

sigma=np.zeros(15,float)

# bmaj
#bmaj_f=[]
#bmaj_p=[]
bmaj_t10=[]
bmaj_tHI=[]

# bmin
#bmin_f=[]
#bmin_p=[]
bmin_t10=[]
bmin_tHI=[]

# rms
#fsigma=[]
#psigma=[]
t10sigma=[]
tHIsigma=[]

# flux
#bp_flux=[]
#p_flux=[]
t10_flux=[]
tHI_flux=[]

# Frequencies
#freq_f=[]
#freq_p=[]
freq_t10=[]
freq_tHI=[]


for ii in seq:       
    logprint ('Extract Data from Cubes, Spw='+str(ii), logfileout='logs/QA2.log')
    
#Extract spectra for 10 mJy cube
    try:
        logprint ('Extract Spectra from 10 mJy cube, Spw='+str(ii), logfileout='logs/QA2.log')
        image_name='targetcube_10mJy_spw'+str(ii)+'.image'
        header=imhead(image_name)
        default('imval')
        stokes='I'
        imagename=image_name
        results=imval()
        for i in range(2048):
            bmaj_t10.append(header['perplanebeams']['beams']['*'+str(i)]['*0']['major']['value'])
            bmin_t10.append(header['perplanebeams']['beams']['*'+str(i)]['*0']['minor']['value'])
            freq_t10.append(results['coords'][i,3]/1e6)
            t10_flux.append(results['data'][i])
            
        default('imval')
        stokes='I'
        box='0,0,16,16'
        imagename=image_name
        results=imval()
        for i in range(2048):
            t10sigma.append(np.std(results['data'][:,:,i]))

        sigma[ii]=np.std(results['data'][:,:,:])  # Take noise averaged over all channels of 10 mJy cube
              
    except:
        logprint("Unable to extract data from 10 mJy cube, spw "+str(ii), logfileout='logs/QA2.log')

#Extract spectrum for deepfield cubes
    if ii==13:
        try:    
            logprint ('Extract Spectra from HI source cube, Spw='+str(ii), logfileout='logs/QA2.log')    
            image_name='targetcube_HIdet_spw'+str(ii)+'.image'
            header=imhead(image_name)
            default('imval')
            stokes='I'
            box='16,16,48,48'
            imagename=image_name
            results=imval()
            for i in range(2048):
                bmaj_tHI.append(header['perplanebeams']['beams']['*'+str(i)]['*0']['major']['value'])
                bmin_tHI.append(header['perplanebeams']['beams']['*'+str(i)]['*0']['minor']['value'])
                freq_tHI.append(results['coords'][0,0,i,3]/1e6)
                tHI_flux.append(np.mean(results['data'][:,:,i]))
                
            default('imval')
            stokes='I'
            box='0,0,16,16'
            imagename=image_name
            results=imval()
            for i in range(2048):
                tHIsigma.append(np.std(results['data'][:,:,i]))
                

    
        except:
            logprint("Unable to extract data from HI source cube, spw "+str(ii), logfileout='logs/QA2.log')

# Convert lists to arrays and replace zeros with nan (for plotting)
bmaj_t10=np.array(bmaj_t10,float)
bmaj_tHI=np.array(bmaj_tHI,float)
bmin_t10=np.array(bmin_t10,float)
bmin_tHI=np.array(bmin_tHI,float)
t10sigma=np.array(t10sigma,float)
tHIsigma=np.array(tHIsigma,float)
t10_flux=np.array(t10_flux,float)
tHI_flux=np.array(tHI_flux,float)

# Replace zeros with NaNs
bmaj_t10[bmaj_t10==0.0]=np.nan
bmaj_tHI[bmaj_tHI==0.0]=np.nan
bmin_t10[bmin_t10==0.0]=np.nan
bmin_tHI[bmin_t10==0.0]=np.nan
t10sigma[t10sigma==0.0]=np.nan
tHIsigma[tHIsigma==0.0]=np.nan
t10_flux[t10_flux==0.0]=np.nan
tHI_flux[tHI_flux==0.0]=np.nan

# Plot spectra for 10 mJy source
fig=pylab.figure()
pylab.plot(freq_t10,bmaj_t10,'b-')
pylab.plot(freq_t10,bmin_t10,'r-')
pylab.xlim(min(freq_t10),max(freq_t10))
pylab.ylim(4.,12.)
pylab.xlabel('Frequency [MHz]')
pylab.ylabel('Beam size ["]')
pylab.title('Beamsize for 10 mJy Source')
pylab.savefig('image_10mJy_beamsize.png')
pylab.close(fig)

fig=pylab.figure()
pylab.plot(freq_t10,t10_flux,'k-')
pylab.xlim(min(freq_t10),max(freq_t10))
pylab.xlabel('Frequency [MHz]')
pylab.ylabel('Flux [Jy/bm]')
pylab.title('Flux of 10 mJy source')
pylab.savefig('image_10mJy_spectrum.png')
pylab.close(fig)

fig=pylab.figure()
pylab.plot(freq_t10,t10sigma,'k-')
pylab.xlim(min(freq_t10),max(freq_t10))
pylab.xlabel('Frequency [MHz]')
pylab.ylabel('RMS [Jy/bm]')
pylab.title('Noise of 10 mJy source')
pylab.savefig('image_10mJy_rms.png')
pylab.close(fig)

# Plot spectra for HI field
fig=pylab.figure()
pylab.plot(freq_tHI,bmaj_tHI,'b-')
pylab.plot(freq_tHI,bmin_tHI,'r-')
pylab.xlim(min(freq_tHI),max(freq_tHI))
pylab.ylim(4.,12.)
pylab.xlabel('Frequency [MHz]')
pylab.ylabel('Beam size ["]')
pylab.title('Beamsize for HI field')
pylab.savefig('image_HIdet_beamsize.png')
pylab.close(fig)

fig=pylab.figure()
pylab.plot(freq_tHI,tHI_flux,'k-')
pylab.xlim(min(freq_tHI),max(freq_tHI))
pylab.xlabel('Frequency [MHz]')
pylab.ylabel('Flux [Jy/bm]')
pylab.title('Flux of HI field')
pylab.savefig('image_HIdet_spectrum.png')
pylab.close(fig)

fig=pylab.figure()
pylab.plot(freq_tHI,tHIsigma,'k-')
pylab.xlim(min(freq_tHI),max(freq_tHI))
pylab.xlabel('Frequency [MHz]')
pylab.ylabel('RMS [Jy/bm]')
pylab.title('Noise of HI field')
pylab.savefig('image_HIdet_rms.png')
pylab.close(fig)

#Move plots, images to sub-directory
os.system("mv *.png plots")
os.system("mv *cube*.* images")

#Make output webpage
if os.path.exists("QA.html"):
    os.system("rm QA.html")
wlog = open("QA.html","w")
wlog.write('<html>\n')
wlog.write('<head>\n')
wlog.write('<title>CHILES Pipeline Web Log</title>\n')
wlog.write('<style>table tr {page-break-inside: avoid}</style>\n')
wlog.write('</head>\n')
wlog.write('<body>\n')
wlog.write('<br>\n')
wlog.write('<hr>\n')
wlog.write('<center>CHILES_pipe_QA results:</center>')
wlog.write('<li> Session: '+SDM_name+'</li>\n')
wlog.write('<li><a href="logs/QA2.log">QA Log</a></li>\n')
wlog.write('<hr>\n')
wlog.write('<br>\n')
wlog.write('<li> Amp. vs. Frequency for phase calibrator showing uncalibrated masking: \n')
wlog.write('<br>\n')
wlog.write('<table> \n')
for ii in seq:
    wlog.write('<tr><img src="plots/flag1500m_Spw'+str(ii)+'.png"></tr>\n')
wlog.write('</table> \n')
wlog.write('</li>')
wlog.write('<br>')
wlog.write('<br>\n')
wlog.write('<li> Beam Size vs. Frequency for 10 mJy source</li>\n')
wlog.write('<li><img src="plots/image_10mJy_beamsize.png"></li><br>\n')
wlog.write('<li> Flux vs. Frequency for 10 mJy source</li>\n')
wlog.write('<li><img src="plots/image_10mJy_spectrum.png"></li><br>\n')
wlog.write('<li> RMS Noise vs. Frequency for 10 mJy source</li>\n')
wlog.write('<li><img src="plots/image_10mJy_rms.png"></li><br>\n')
wlog.write('<li> Beam Size vs. Frequency for HI source</li>\n')
wlog.write('<li><img src="plots/image_HIdet_beamsize.png"></li><br>\n')
wlog.write('<li> Flux vs. Frequency for HI source</li>\n')
wlog.write('<li><img src="plots/image_HIdet_spectrum.png"></li><br>\n')
wlog.write('<li> RMS Noise vs. Frequency for HI source</li>\n')
wlog.write('<li><img src="plots/image_HIdet_rms.png"></li><br>\n')
wlog.write('<br>\n')
wlog.write('<li> RMS Noise per Spw for 10 mJy source cube</li>\n')
wlog.write('<table>\n')
for ii in seq:
    wlog.write('<tr>\n')
    wlog.write('<td> Noise in spw '+str(ii)+'= '+str(sigma[ii]*1e3)+' mJy</td>\n')
    wlog.write('</tr>\n')
wlog.write('</table>\n')
wlog.write('<br>\n')
wlog.write('<hr>\n')
wlog.write('</body>\n')
wlog.write('</html>\n')
wlog.close()

#os.system('cp *.html FINAL/.')
#os.system("cp -r plots FINAL/.")



logprint ("Finished CHILES_pipe_QA2.py", logfileout='logs/QA2.log')
time_list=runtiming('QA2', 'end')

pipeline_save()

# Save variable values
#os.system("cp -r pipeline_shelf.restore FINAL/.")

# Copy logs to FINAL directory (needs to be after final "save" to preserve all information
#os.system("cp *.log logs")
#os.system("cp -r logs FINAL/.")
