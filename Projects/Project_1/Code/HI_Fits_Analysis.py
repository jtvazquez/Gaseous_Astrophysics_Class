import numpy as np 
import matplotlib.pyplot as plt 
from astropy.io import fits
import os

# Location of the GASS HI datacube. Defaults to a Data/ directory alongside this
# script; override with the GASS_DATA_DIR environment variable.
directory = os.environ.get('GASS_DATA_DIR',
                           os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'Data'))

fname = os.path.join(directory, 'gass_data.fits')

fits_data = fits.open(fname)[0]



data = fits_data.data

hdr = fits_data.header


