"""
Apply manual flags to the deepfield.
@Nick Luber
"""

i = int(var[0]) # A session variable.

import numpy as np

# Find the measurement set and the flag file.
cparams = np.loadtxt('/scratch/nml0011/cluster_pipeline/chiles_parameters.txt', delimiter=';', dtype=str)
visib1 = '/scratch/nml0011/'+cparams[i][1]+'/FINAL/'+cparams[i][2]+'_calibrated_deepfield_SMOOTH.ms'
visib2 = '/scratch/nml0011/'+cparams[i][1]+'/FINAL/'+cparams[i][2]+'_calibrated_deepfield.ms'
flgfile = '/scratch/nml0011/extra_flags/'+cparams[i][0]+'_extraflag.txt'

# Flag the data.
default('flagdata')
vis=visib1
mode='list'
inpfile=flgfile
action='apply'
flagbackup=False
overwrite=False
savepars=True
flagdata()

default('flagdata')
vis=visib2
mode='list'
inpfile=flgfile
action='apply'
flagbackup=False
overwrite=False
savepars=True
flagdata()
