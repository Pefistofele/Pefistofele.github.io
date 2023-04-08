# Plot channel map of outflow
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from astropy.wcs import WCS as wcs
import os

def axes_setting(ax):
	lon = ax.coords[0]
	lat = ax.coords[1]
		
	# lon.set_ticks_visible(False)
	lon.set_ticklabel_visible(False)
	# lat.set_ticks_visible(False)
	lat.set_ticklabel_visible(False)
	lon.set_axislabel('')
	lat.set_axislabel('')



# Change working directory
cwd =os.getcwd()
os.chdir(cwd)

# Read CO outflow datatcube
fits_file = "/run/media/pefistofele/TOSHIBA EXT/CENSUS/S29/lines/CO_2_1/cygxs29_CO_2_1.image.fits"
co_data, co_header = fits.getdata(fits_file, header = True)
naxis, cdelt, crval, crpix = co_header["NAXIS3"], co_header["CDELT3"], co_header["CRVAL3"], co_header["CRPIX3"]

vlsr_list = [(crval+cdelt*(i-crpix+1))/1000 for i in range(naxis)] # velocity in unit km/s

vsys_index = 43 # vsys ~ 5 km/s

# Set figure size
row, col =3, 4
fig = plt.figure(figsize=(col*3, row*3))
co_wcs = wcs(co_header)


for index in range(row*col):
	ax = plt.subplot(row, col, index+1, projection = co_wcs[0,:,:])

	axes_setting(ax)
	xshape, yshape = co_data[0,:,:].shape
	ax.set_xlim(0.2*xshape, 0.8*xshape)
	ax.set_ylim(0.2*yshape, 0.8*yshape)

	chan_data_blue = co_data[vsys_index-8-index,:,:]
	rms = np.nanstd(chan_data_blue, ddof=0)
	ax.contour(chan_data_blue, levels = [i*rms for i in [3,6,9,12,15]], colors = 'blue')

	chan_data_red = co_data[vsys_index+8+index,:,:]
	rms = np.nanstd(chan_data_red, ddof =1)
	ax.contour(chan_data_red, levels = [i*rms for i in [3,6,9,12,15]], colors = 'red')

	# Plot source name
	textstr = "%.1f km/s"%(vlsr_list[vsys_index-8-index])
			
	# these are matplotlib.patch.Patch properties
	props = dict(boxstyle='round', facecolor='white', alpha=0.5)
	ax.text(0.05, 0.95, textstr, fontsize=8, transform=ax.transAxes, verticalalignment='top', bbox=props, color = 'blue')

	textstr = " %.1f km/s"%(vlsr_list[vsys_index+8+index])
			
	# these are matplotlib.patch.Patch properties
	props = dict(boxstyle='round', facecolor='white', alpha=0.5)
	ax.text(0.05, 0.85, textstr, fontsize=8, transform=ax.transAxes, verticalalignment='top', bbox=props, color = 'red')

plt.savefig("S29OutflowChanMap.pdf", format='pdf', inches_bbox = 'tight')
