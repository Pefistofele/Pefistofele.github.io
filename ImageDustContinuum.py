from cgitb import text
from  astropy.io import fits
from matplotlib import transforms
import matplotlib.pyplot as plt
import os
from numpy.lib.utils import source
import pandas as pd
import numpy as np
from astropy.wcs import WCS
import astropy.coordinates as coord
from astropy import units as u
import pyregion
from pvextractor import Path
from pvextractor import extract_pv_slice
from matplotlib import colors
from matplotlib.colors import LogNorm
import astropy.units as u

from astropy.visualization import (MinMaxInterval, SqrtStretch,
								   ImageNormalize)

from matplotlib.patches import Ellipse
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredAuxTransformBox
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar

def header_editor(header):
	header.remove('CTYPE3')
	header.remove('CRVAL3')
	header.remove('CDELT3')
	header.remove('CRPIX3')
	header.remove('CTYPE4')
	header.remove('CRVAL4')
	header.remove('CDELT4')
	header.remove('CRPIX4')
	header.remove('NAXIS3')
	header.update(NAXIS=2)
	header.remove('NAXIS4')
	return 0

def axes_setting(ax, show = False, datashape = (512, 512), pixsize = 0.2):
	"""
	Set axes for combined maps
	pixsize: unit arcsec
	"""
	lon = ax.coords[0]
	lat = ax.coords[1]
	# lon.set_ticklabel_visible(False)
	# lat.set_ticklabel_visible(False)

	lon.set_ticks(spacing=10. * u.arcsec)
	lat.set_ticks(spacing=10. * u.arcsec)

	if not(show):
		lon.set_axislabel(' ', minpad = 0.1)
		lat.set_axislabel(' ', minpad = 0.1)
		# lon.set_ticklabel_visible(False)
		# lat.set_ticklabel_visible(False)
	else:
		# print(1)
		lon.set_axislabel('RA', minpad = 0.1)
		lat.set_axislabel('DEC', minpad = 0.1)
	
	lon.set_ticks(number=4)	
	lat.set_ticks(number=4)
	lon.display_minor_ticks(True)
	lat.display_minor_ticks(True)
		# lon.set_ticklabel_visible(False)
		# lat.set_ticklabel_visible(False)
		# lon.set_ticklabel_visible(False)

def Add_beam_size(ax, contin_header):
	"""
	Set beam size and scale bar
	"""
	#Add beam sizes
	box = AnchoredAuxTransformBox(ax.transData, loc='lower left')
	# Read beam info
	try:
		beam_major_deg, beam_minor_deg, beam_angle_deg = contin_header["BMAJ"], contin_header["BMIN"], contin_header["BPA"] # in unit degree\
	except KeyError:
		beam_major_deg, beam_minor_deg, beam_angle_deg = 4.85e-4, 3.69e-4, -85.25
	
	beam_major, beam_minor = beam_major_deg*3600, beam_minor_deg*3600 # in unit arcsec
	# print(beam_major)

	pix_scale = contin_header['CDELT1']*3600 # Convert to arcsec

	el = Ellipse((0, 0), width = beam_minor/pix_scale, height = beam_major/pix_scale, angle = beam_angle_deg, facecolor = 'none', edgecolor = 'black')
	box.drawing_area.add_artist(el)
	
	ax.add_artist(box)

	pixel_length = (pix_scale*distance*u.au).to(u.pc) # in unit AU
	set_scale = 0.1*u.pc # in unit AU
	scalebar = AnchoredSizeBar(ax.transData,
							(set_scale/pixel_length).value, '%.1f pc'%set_scale.value, 'lower right', 
							pad=0.1,
							color='black',
							frameon=False,
							size_vertical=0.8)
	ax.add_artist(scalebar)

	return 0


def condenpos(conden_info, contin_wcs):
	"""
	Get the x-y position of target condensation
	"""
	spatial_ra = conden_info.iloc[0]["Ra"]
	spatial_dec = conden_info.iloc[0]["Dec"]
		
	radius = round(float(conden_info['Size_FWHM'])/distance/np.pi*180*3600/2, 4)

	regionx, regiony = round(float(spatial_ra), 4), round(float(spatial_dec), 4)

	xcen_raw, ycen_raw = contin_wcs.all_world2pix(regionx, regiony, 1)
	xcen, ycen = xcen_raw-1, ycen_raw-1
	
	return xcen, ycen


def plot_cont_map_part(source_list, fig, row, col):
	"""
	Plot the continuum emission partially
	"""
	for source_name, subplot_index in zip(source_list, range(len(source_list))):
	# Read continuum data
		source_info= sample_info[sample_info.source_name == source_name]
		contin_file_name = source_info.iloc[0]["continuum_file"]
		contin_data, contin_header = fits.getdata(contin_path+contin_file_name, header = True)

		header_editor(contin_header)
		contin_wcs = WCS(contin_header)
		#####################################################
		# Set axes for outflow maps
		ax = fig.add_subplot(row, col, subplot_index+1, projection = contin_wcs)
	
		if subplot_index != col*(row-1):
			axes_setting(ax)
		# Set relative coordinate ticks
		elif subplot_index == col*(row-1):
			axes_setting(ax, show = True, datashape = contin_data[0,0].shape, pixsize = 0.2)
		
		# Zoom in
		xmax, ymax = contin_data[0, 0].shape
		if not("DR" in source_name):
			ax.set_xlim(0.2*xmax, 0.8*xmax)
			ax.set_ylim(0.2*ymax, 0.8*ymax)

		# Plot SMA 1.3 mm continuum map
		continuum_max, continuum_min = np.nanmax(contin_data), np.nanmin(contin_data)
		# Create interval object
		interval = MinMaxInterval()
		vmin, vmax = interval.get_limits(contin_data[0, 0])
		norm = ImageNormalize(vmin=vmin, vmax=vmax, stretch=SqrtStretch())

		# Add beam size
		Add_beam_size(ax, contin_header=contin_header)

		# Create an ImageNormalize object using a SqrtStretch object
		ax.imshow(contin_data[0,0], transform = ax.get_transform(contin_wcs), cmap = plt.cm.gray_r, norm=norm)

		# Continuum contour
		contin_sigma = np.nanstd(contin_data[0,0])/2
		contin_levels = [contin_sigma*i for i in [-4, 3, 6, 9, 12, 18, 24, 36, 48, 60]]
		ax.contour(contin_data[0, 0], levels = contin_levels, colors = 'black', linewidths = 0.5)

		# Annonate condensation name
		target_name = list(filter(lambda x: source_info.iloc[0] ["other_name"] in x, conden_cata["Name_core"]))

		conden_info_list = [conden_cata[conden_cata.Name_core == cond_name] for cond_name in target_name]

		conden_pos_list = [condenpos(conden_info, contin_wcs) for conden_info in conden_info_list]

		# Plot condensation name
		for conden_pos, conden_index in zip(conden_pos_list, range(len(conden_pos_list))):
			# ax.text(conden_pos[0], conden_pos[1], "MM%d"%(conden_index+1), size= 6)
			ax.scatter(conden_pos[0], conden_pos[1], marker = "+", color = "red", s = 16, linewidth = 0.5)
			# if conden_pos[0]>xlim_min and conden_pos[0]<xlim_max and conden_pos[1]>ylim_min and conden_pos[1]<ylim_max:
			ax.text(conden_pos[0], conden_pos[1], "MM%d"%(conden_index+1), size= 6, color = 'white')
		
		# Plot source name
		textstr = '\n'.join((
						source_name,
					))
		
		# these are matplotlib.patch.Patch properties
		props = dict(boxstyle='round', facecolor='white', alpha=0.5)
		ax.text(0.05, 0.95, textstr, fontsize=12, transform=ax.transAxes, verticalalignment='top', bbox=props)

# Change working directory
cur_dir = "/run/media/pefistofele/TOSHIBA EXT/CENSUS/Condensation_sum/Outflow_Sum/"
os.chdir(cur_dir)

# Path for associated file
contin_path =  "/run/media/pefistofele/TOSHIBA EXT/CENSUS/Condensation_sum/Conden_Sum/"
outflow_path = "/run/media/pefistofele/TOSHIBA EXT/CENSUS/Condensation_sum/Outflow_Sum"

conden_cata = pd.read_table('/run/media/pefistofele/TOSHIBA EXT/CENSUS/Condensation_sum/Conden_Sum/condensations_edition.txt',delim_whitespace=True)


# Observation settings
pix_scale = 0.2
distance = 1.4e3 # unit pc

sample_info = pd.read_csv('/run/media/pefistofele/TOSHIBA EXT/CENSUS/Condensation_sum/sample_contin_info.csv')

source_list_all = sample_info.iloc[:]["source_name"]
source_list = source_list_all[:-2]
# source_list = ["N12"]
# Set figure size
row0, col0 = 5, 3
row1, col1 = 5, 4

fig0 = plt.figure(num=0, figsize = (4*col0, 3.5*row0))
fig1 = plt.figure(num=1, figsize = (4*col1, 3.5*row1))
# fig2 = plt.figure(num=2, figsize = (4*col, 3.5*row))
# fig0.tight_layout()

plot_cont_map_part(source_list[:row0*col0], fig0, row0, col0)
plot_cont_map_part(source_list[row0*col0:], fig1, row1, col1)
# plot_cont_map_part(source_list_all[2*row*col:], fig2)
# fig0.subplots_adjust(wspace=0.0, hspace=0.0)
fig0.savefig("SMA1.3mmCont001.pdf", dpi = 300, format ="pdf", bbox_inches = 'tight')
fig1.savefig("SMA1.3mmCont002.pdf", dpi = 300, format ="pdf", bbox_inches = 'tight')
# fig2.savefig("SMA1.3mmCont003.pdf", dpi = 300, format ="pdf", bbox_inches = 'tight')

plt.show()