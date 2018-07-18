"""Runs tweakreg and creates plot

This script sets tweakreg prameter and outputs png 
residual plots for each filter.


Author:

--------- Myles McKay, October 20, 2016

Use:

python mytweakreg_quad_main.py --path='/grp/hst/wfc3v/martlin/nov_4_2016_Myles' --refim='filename.foo' --filter='F410M'


Output:

Residual plots in png format

Dependencies:

Notes:
"""


import glob
from stsci.tools import teal
from drizzlepac import tweakreg
import os
from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from stwcs import updatewcs
import argparse

def mytweakreg_quad_main(path, refim):


	os.chdir(path)
	os.system('pwd')
	updatewcs.updatewcs('*flc.fits')
	

	teal.unlearn('tweakreg')
	findpars = {'computesig':True,'threshold':100.,'conv_width':3.5,'peakmin':100.,'peakmax':90000.,'nsigma':1.5,
					'ratio':1.0, 'theta':0.0,'dqbits':None, 'use_sharp_round':False, 'ylimit':0.1 }
	
	tweakreg.TweakReg(glob.glob('*flc.fits'), expand_refcat=False, enforce_user_order=False,
	                    updatehdr=False, shiftfile=True, writecat=True, clean=False, interactive= True, verbose=False,
	                    headerlet=False, minobj=15, searchrad=250.0, searchunits='pixels', use2dhist=True, see2dplot=True,
	                    separation=0.5, fitgeometry= 'rscale', residplot='both', labelsize=8,nclip=5, sigma=3.0,
	                    refimage=refim,
	                    imagefindcfg=findpars, 
	                    refimagefindcfg=findpars)
	  
		
# -------------------------------------------------------------------
# For command line execution
# -------------------------------------------------------------------

def parse_args():
    """Parses command line arguments.

    Parameters:
        nothing

    Returns:
        args : argparse.Namespace object
            An argparse object containing all of the added arguments.

    Outputs:
        nothing
    """

    #Create help string:
    path_help = 'Path to the folder with files to run tweakreg.'
    refim_help = 'Reference Image for tweakreg'
    
    # Add arguments:
    parser = argparse.ArgumentParser()
    parser.add_argument('--path', '-path', dest = 'path', action = 'store',
                        type = str, required = True, help = path_help)
    parser.add_argument('--refim', '-refim', dest = 'refim', action = 'store',
                        type = str, required = True, help = refim_help)
    # Parse args:
    args = parser.parse_args()

    return args
# -------------------------------------------------------------------
if __name__ == '__main__':
	args = parse_args()

	mytweakreg_quad_main(args.path, args.refim)

		
		
		
		
		
		