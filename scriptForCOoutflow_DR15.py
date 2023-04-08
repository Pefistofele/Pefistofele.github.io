# Plot Outflow map
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

from astropy.visualization import (MinMaxInterval, SqrtStretch,
								   ImageNormalize)
from matplotlib.patches import Ellipse
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredAuxTransformBox
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar

def axes_setting(ax, show = False, datashape = (512, 512), pixsize = 0.2):
	"""
	Set axes for combined maps
	pixsize: unit arcsec
	"""
	lon = ax.coords[0]
	lat = ax.coords[1]
	# lon.set_ticklabel_visible(False)
	# lat.set_ticklabel_visible(False)

	# lon.set_ticks(spacing=10. * u.arcsec)
	# lat.set_ticks(spacing=10. * u.arcsec)

	if not(show):
		lon.set_axislabel(' ', minpad = 0.1)
		lat.set_axislabel(' ', minpad = 0.1)
		lon.set_ticks(number=4)
		lat.set_ticks(number=4)
		# lon.set_ticklabel_visible(False)
		# lat.set_ticklabel_visible(False)
	else:
		# print(1)
		lon.set_axislabel('RA', minpad = 0.1)
		lat.set_axislabel('DEC', minpad = 0.1)
		lon.set_ticks(number=4)
		lat.set_ticks(number=4)
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
							(set_scale/pixel_length).value, '%.2f pc'%set_scale.value, 'lower right', 
							pad=0.1,
							color='black',
							frameon=False,
							size_vertical=0.5)
	ax.add_artist(scalebar)

	return 0

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


def arrow_plot(ax, pa, color, index, conden_pos_list):
	"""
	Plot arrow to represent outflow orientation 
	"""

	xcen, ycen = conden_pos_list[index]
	if index == 4:
		l = 100
	else:
		l = 50
	
	dx, dy = int(l*np.sin(pa/180*np.pi)), int(l*np.cos(pa/180*np.pi))
	# offset from the center
	r = 6
	rx, ry = r*np.sin(pa/180*np.pi), r*np.cos(pa/180*np.pi)

	ax.arrow(float(xcen-rx), float(ycen+ry), -dx, dy, length_includes_head = True, head_width = 1, color = color)

# Change working directory
cur_dir = "/Office_PC/Outflow_Sum/"
os.chdir(cur_dir)

# Path for associated file
contin_path =  "/Office_PC/Conden_Sum/"
outflow_path = "/Office_PC/Outflow_Sum"

conden_cata = pd.read_table('/Office_PC/Conden_Sum/condensations_edition.txt',delim_whitespace=True)


# Observation settings
pix_scale = 0.2
distance = 1.4e3 # unit pc

pd_value = [np.nan]

sample_info = pd.DataFrame(
    {
        "source_name": ["DR15"],
        "other_name": ["DR15"],
        "continuum_file": ["DR15.cont.comb.lms.fits"],
        "outflow_index0": [0],
		"outflow_angle0_red":[-20],
		"outflow_angle0_blue":[160],
		"disk_angle0":[148],
		"disk_pvsigma0": pd_value,
		"kepler_central_mass": pd_value,
		"outflow_index1": [4],
		"outflow_angle1_red":[-75],
		"outflow_angle1_blue":[105]
    }
)

source_list = ["DR15"]
# source_list = ["N12"]
# Set figure size
row, col = 3, 3
fig0 = plt.figure(num=0, figsize = (2*col, 2*row))
fig0.tight_layout()

for source_name, subplot_index in zip(source_list, range(len(source_list))):
	# Read continuum data
	source_info= sample_info
	contin_file_name = source_info.iloc[0]["continuum_file"]
	contin_data, contin_header = fits.getdata(contin_path+contin_file_name, header = True)

	header_editor(contin_header)
	contin_wcs = WCS(contin_header)
	#####################################################
	# Set axes for outflow maps
	ax = fig0.add_subplot(111, projection = contin_wcs)

	axes_setting(ax, show = True)
	# Set relative coordinate ticks
	# Plot SMA 1.3 mm continuum map
	continuum_max, continuum_min = np.nanmax(contin_data), np.nanmin(contin_data)
	# Create interval object
	interval = MinMaxInterval()
	vmin, vmax = interval.get_limits(contin_data[0, 0])
	norm = ImageNormalize(vmin=vmin, vmax=vmax, stretch=SqrtStretch())
	
	# Create an ImageNormalize object using a SqrtStretch object
	ax.imshow(contin_data[0,0], transform = ax.get_transform(contin_wcs), cmap = plt.cm.gray_r, norm=norm)

	# Continuum contour
	contin_sigma = np.nanstd(contin_data[0,0])
	contin_levels = [contin_sigma*i for i in [3, 6, 9, 18, 36, 48, 60]]
	ax.contour(contin_data[0, 0], levels = contin_levels, colors = 'black', linewidths = 0.5)

	# Annonate condensation name
	target_name = list(filter(lambda x: source_info.iloc[0] ["other_name"] in x, conden_cata["Name_core"]))

	conden_info_list = [conden_cata[conden_cata.Name_core == cond_name] for cond_name in target_name]

	conden_pos_list = [condenpos(conden_info, contin_wcs) for conden_info in conden_info_list]

	# Plot condensation name
	for conden_pos, conden_index in zip(conden_pos_list, range(len(conden_pos_list))):
		# ax.text(conden_pos[0], conden_pos[1], "MM%d"%(conden_index+1), size= 6)
		ax.scatter(conden_pos[0], conden_pos[1], marker = "+", color = "red", s = 16, linewidth = 0.8)
		# if conden_pos[0]>xlim_min and conden_pos[0]<xlim_max and conden_pos[1]>ylim_min and conden_pos[1]<ylim_max:
		ax.text(conden_pos[0], conden_pos[1], "MM%d"%(conden_index+1), size= 6, color = 'white')
	
	# Plot source name
	textstr = '\n'.join((
					source_name,
				))
	
	# these are matplotlib.patch.Patch properties
	props = dict(boxstyle='round', facecolor='white', alpha=0.5)
	ax.text(0.05, 0.95, textstr, fontsize=12, transform=ax.transAxes, verticalalignment='top', bbox=props)

	Add_beam_size(ax, contin_header=contin_header)
	fig0.subplots_adjust(wspace=0.0, hspace=0.0)
	fig0.savefig("SMA1.3mmCont_DR15.pdf", dpi = 360, format ="pdf", bbox_inches = 'tight')

	# Read outflow data(including mosaic mode)
	
	outflow_list = [
		outflow_path + "/Red_lobe/cygx"+ source_name.lower()+"_CO_2_1.red.fits",
		outflow_path + "/Blue_lobe/cygx"+ source_name.lower()+"_CO_2_1.blue.fits",
		outflow_path + "/Red_lobe/cygx"+ source_name.lower()+"_SiO_5_4.red.fits",
		outflow_path + "/Blue_lobe/cygx"+ source_name.lower()+"_SiO_5_4.blue.fits",
	]

	# for filename in filename_list:
	color_list = []; filename_list = []
	linestyle_list = []; alpha_list = []
	if os.path.exists(outflow_list[0]):
		color_list.append('red')
		linestyle_list.append("solid")
		alpha_list.append(0.8)
		filename_list.append(outflow_list[0])

	if os.path.exists(outflow_list[1]):
		color_list.append('blue')
		linestyle_list.append("solid")
		alpha_list.append(0.8)
		filename_list.append(outflow_list[1])
	
	if os.path.exists(outflow_list[2]):
		color_list.append('yellow')
		linestyle_list.append("solid")
		alpha_list.append(0.8)
		filename_list.append(outflow_list[2])

	if os.path.exists(outflow_list[3]):
		color_list.append('purple')
		linestyle_list.append("solid")
		alpha_list.append(0.8)
		filename_list.append(outflow_list[3])

	# Plot outflow contour  
		
	for filename, color, linestyle, alpha in zip(filename_list, color_list, linestyle_list, alpha_list):
		imgdata, header = fits.getdata(filename, header = True)
		outflow_wcs = WCS(header)
		vmax = np.nanmax(imgdata)
	
		upperpercent = 0.95; lowerpercent = 0.20
		num_bin = 5
		rms_level = [i*vmax for i in np.linspace(lowerpercent, upperpercent, num_bin)]

		cntr = ax.contour(imgdata, transform= ax.get_transform(outflow_wcs), levels = rms_level, colors = color, linewidths = 1, linestyles = linestyle, alpha = alpha)

	# # Zoom in the 1.3 mm continuum map
	# xmax, ymax = contin_data[0, 0].shape
	# xlim_min, xlim_max, ylim_min, ylim_max = xmax/5, xmax/5*4, ymax/5, ymax/5*4 
	# ax.set_xlim(xlim_min, xlim_max)
	# ax.set_ylim(ylim_min, ylim_max)

	# Add arrow to show outflow orietation
	arrow_color_list = ["red", "blue"]
	for num in range(2):
		outflow_index = source_info.iloc[0]["outflow_index%d"%num]
		if not(np.isnan(outflow_index)):
			outflow_index = int(outflow_index)
			angle_list = [source_info.iloc[0]["outflow_angle%d_%s"%(num, color)] for color in arrow_color_list]
			for color, angle in zip(arrow_color_list, angle_list):
				if not(np.isnan(angle)):
					arrow_plot(ax, angle, color, outflow_index,conden_pos_list)

	# plt.savefig("CO/%s_CO_outflow.pdf"%source_name, format="pdf", dpi = 360)
	

fig0.subplots_adjust(wspace=0.0, hspace=0.0)
fig0.savefig("outflowVScontinuum_DR15.pdf", dpi = 400, format ="pdf", bbox_inches = 'tight')
# fig0.savefig("N12outflowVScontinuum.pdf", dpi = 400, format ="pdf", bbox_inches = 'tight')