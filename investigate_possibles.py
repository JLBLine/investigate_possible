#!/usr/bin/python
import numpy as np
import statsmodels.api as sm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Ellipse
from itertools import combinations
import atpy
import optparse
import pyfits as pf
from astLib import astWCS, astPlots
import os
from matplotlib import patches
from sys import path
#path.append('/home/jline/Documents/cataloguing/PUMA/scripts')
#path.append('/home/jline/Documents/cataloguing/FHD/carroll/paper/complete')
import plot_outcomes_lib as pol
import make_table_lib as mkl
from pylab import MaxNLocator

dr= np.pi/180.0

parser = optparse.OptionParser()

parser.add_option('-f', '--cat_freqs', 
	help='Enter the frequencies of each catalogue')

parser.add_option('-p', '--pref_cats',
	help='Enter names of matched cataloges, in order of preferable position')

parser.add_option('-m', '--matched_cats', 
	help='Enter the names of each catalogue, base catalogue first')

parser.add_option('-i', '--input_bayes', 
	help='Enter name of eyeball bayes file')

parser.add_option('-o', '--output_name', 
	help='Enter name for output catalogue')

options, args = parser.parse_args()

closeness = 1.5/60.0

cat_fs = options.cat_freqs.split(',')
cat_freqs= []
for fs in cat_fs:
	split = fs.split('~')
	split = map(float,split)
	cat_freqs.append(split)

matched_cats = options.matched_cats.split(',')

pref_cats = ['weighted']
for pref in options.pref_cats.split(','): pref_cats.append(pref)

num_freqs = []
for freq in cat_fs: num_freqs.append(len(freq.split('~')))

mkl.closeness = 1.5/60.0
mkl.high_prob = 0.95
mkl.low_prob = 0.8
mkl.chi_thresh = 10.0
mkl.jstat_thresh = 0.1
mkl.num_freqs = num_freqs
mkl.split = '00:02:00'

pol.closeness = closeness
pol.high_prob = 0.95
pol.low_prob = 0.8
pol.chi_thresh = 10.0
pol.jstat_thresh = 0.1
pol.num_freqs = num_freqs
pol.split = '00:02:00'
pol.matched_cats = matched_cats

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def deg_to_degmins(x,style):	    #converts angle degrees form in to dd:mm:ss.ss
	x=float(x)
	deg=abs(x)
	degr=deg-int(deg)
	mins=(degr)*60.00
	secs=(mins-int(mins))*60.0
	if mins!=0:
		if -1e-5<=(secs-60)<1e-5:
			mins=mins+1
			secs=0.0
	if style == 'info':
		if x>0:
			return '+%02d %02d %04.1f' %(int(deg),int(mins),secs)
		if x<0:
			return '-%02d %02d %04.1f' %(int(deg),int(mins),secs)
	elif style == 'name':
		if x>0:
			return '+%02d%02d%04.1f' %(int(deg),int(mins),secs)
		if x<0:
			return '-%02d%02d%04.1f' %(int(deg),int(mins),secs)

def deg_to_hour(x,style):    #converts angle in degrees in to hh:mm:ss.ss, must input as a string
	x=float(x)
	deg=abs(x)
	hr=deg/15.0
	mins=(hr-int(hr))*60.0
	secs=(mins-int(mins))*60.0
	if mins!=0:
		if -1e-5<=(secs-60)<1e-5:
			mins=mins+1
			secs=0.0
	if style == 'info':
		if x>0:
			return '%02d %02d %04.1f' %(int(hr),int(mins),secs)
		if x<0:
			return '-%02d %02d %04.1f' %(int(hr),int(mins),secs)
	elif style == 'name':
		if x>0:
			return '%02d%02d%04.1f' %(int(hr),int(mins),secs)
		if x<0:
			return '-%02d%02d%04.1f' %(int(hr),int(mins),secs)
	
	
def xtick_RA(x):    #converts angle in degrees in to hh:mm:ss.ss, must input as a string
	x=float(x)
	deg=abs(x)
	hr=deg/15.0
	mins=(hr-int(hr))*60.0
	secs=(mins-int(mins))*60.0
	if mins!=0:
		if -1e-5<=(secs-60)<1e-5:
			mins=mins+1
			secs=0.0
	if x>0:
		return '%02d:%02d:%02d' %(int(hr),int(mins),secs)
	if x<0:
		return '-%02d:%02d:%02d' %(int(hr),int(mins),secs)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

		
def plot_all_image(cat,name,ra,rerr,dec,derr,major,minor,PA,ax,image_wcs):
	try:
		ra_img,dec_img = image_wcs.wcs2pix(ra,dec)
		img_maj,img_min = image_wcs.wcs2pix(ra+float(major),dec+float(minor))
		ra_img_err,dec_img_err = image_wcs.wcs2pix(ra+rerr,dec+derr)
	except AttributeError:
		ra_img,dec_img = image_wcs.wcs_world2pix(ra,dec,0)
		img_maj,img_min = image_wcs.wcs_world2pix(ra+float(major),dec+float(minor),0)
		ra_img_err,dec_img_err = image_wcs.wcs_world2pix(ra+rerr,dec+derr,0)
	
	img_maj = img_maj - ra_img
	img_min = img_min - dec_img
	
	ra_img_err = ra_img_err - ra_img
	dec_img_err = dec_img_err - dec_img
	
	'''Plots a position and an ellipse of specific colours/labels for known catalogues'''
	
	##Set up a load of markers, colours and alpha values
	markers = ['o','*','s','^','D','8']
	marker_sizes = np.array([8,11,8,8,7,10])

	marker_colours = other_colours = ['m', ##fhd
	"b",  ##vlssr
	"#FFEF00",  ##mrc  #AA6600
	"r",  ##sumss
	"g",  ##nvss
	"#c0531f"]  ##mwacs

	ell_colours1 = ["#C370C8", ##fhd
	"#698ACD",  ##vlssr
	"#AA6600",  ##mrc
	"#D84C77",  ##sumss
	"#5FC34F",  ##nvss
	"#D46733"]  ##mwacs

	ell_colours2 = ['m', ##fhd
	"b",  ##vlssr
	"#8B5A00",  ##mrc  #AA6600
	"r",  ##sumss
	"g",  ##nvss
	"#c0531f"]  ##mwacs

	alphas = [0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4]
	
	##Plot the colour by index of catalogue in matched_cats 
	ind = matched_cats.index(cat)
	pol.plot_pos(markers[ind],marker_colours[ind],ra_img,dec_img,ra_img_err,dec_img_err,name,8,ax,1.0)

	if float(minor)!=-100000.0:
		if float(major)!=-100000.0:
			pol.plt_ell(ra_img,dec_img,float(img_maj),float(img_min),-float(PA),ax,ell_colours1[ind],ell_colours2[ind],alphas[ind],1.0)
	
	##HAVE PUT NEGATIVE PA BECAUSE IMAGES ARE AUTOMATICALLY BACKWARDS???
		
def find_extra(cat_name,print_cat,cat_data,image_cats,over_plots,wcss,source_names,ax_spectral,colour1,colour2,freq,handles,plotted_names,present_already,all_fluxs):
	'''Finds any sources in the specified catalogue that wasn't in the original search, and plots
	them to the relevant over_plot subplot'''
	plot_ind = image_cats.index(cat_name)
	##Select the correct plot and wsc
	extra_plot = over_plots[plot_ind]
	image_wcs = wcss[plot_ind]
	##Work out the search boundaries based on the size of the 
	limits = extra_plot.axis()
	ra_max,dec_max = image_wcs.pix2wcs(limits[0],limits[3])
	ra_min,dec_min = image_wcs.pix2wcs(limits[1],limits[2])
	##Search the catalogue and plot the found sources
	for line in cat_data[:-1]:
		info = line.split('|')[1:-1]
		#print info
		name,ra,rerr,dec,derr,flux,ferr,major,minor,PA,flag,ID = info
		name = name.split()[0]
		if (name not in source_names) and (ra_min<float(ra)<ra_max) and (dec_min<float(dec)<dec_max):
			#if present_already=='no': print "extra ['%s', %.5f, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f]" %(name.split()[0],float(ra), float(rerr), float(dec), float(derr), float(freq), float(flux), float(ferr))
			#print "extra ['%s', %.5f, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f]" %(name.split()[0],float(ra), float(rerr), float(dec), float(derr), float(freq), float(flux), float(ferr))
			ra = float(ra)
			dec = float(dec)
			rerr = float(rerr)
			derr = float(derr)
			if rerr==-1.0: rerr=0.0
			if derr==-1.0: derr=0.0
			ra_img,dec_img = image_wcs.wcs2pix(ra,dec)
			img_maj,img_min = image_wcs.wcs2pix(ra+float(major),dec+float(minor))
			img_maj = img_maj - ra_img
			img_min = img_min - dec_img
			ra_img_err,dec_img_err = image_wcs.wcs2pix(ra+rerr,dec+derr)
			ra_img_err = ra_img_err - ra_img
			dec_img_err = dec_img_err - dec_img
			#plot_pos('s',colour1,ra_img,dec_img,ra_img_err,dec_img_err,print_cat+' '+name,6,extra_plot,'')
			
			p = extra_plot.errorbar(ra_img,dec_img,ra_img_err,dec_img_err,marker='s',ms=6,mfc=colour1,mec=colour1,ecolor=colour2,markeredgewidth=1,label=name,linestyle='None')
			
			if print_cat=='mrc':
				if PA!='-100000.0': pol.plt_ell(ra_img,dec_img,float(img_maj),float(img_min),-float(PA),extra_plot,colour2,colour2,0.15,1.0)
			elif print_cat=='comp_v5':
				pass
			else:
				if PA!='-100000.0': pol.plt_ell(ra_img,dec_img,float(img_maj),float(img_min),-float(PA),extra_plot,colour2,colour2,0.3,1.0)
				
			pol.plot_errors('s',colour1,freq,float(flux),float(ferr),print_cat+' '+name,6,ax_spectral)
			all_fluxs.append(float(flux))
			
			if name not in plotted_names:
				handles.append(p,)
				plotted_names.append(name)
				
def find_extra_fhd(print_cat,ra_min,ra_max,dec_min,dec_max,extra_plot,cat_data,source_names,colour1,colour2,handles,plotted_names,present_already,KGS_wcs):

	##Search the catalogue and plot the found sources
	for line in cat_data[:-1]:
		info = line.split('|')[1:-1]
		#print info
		name,ra,rerr,dec,derr,flux,ferr,major,minor,PA,flag,ID = info
		name = name.split()[0]
		if (name not in source_names) and (ra_min<float(ra)<ra_max) and (dec_min<float(dec)<dec_max):
			#if present_already=='no': print "extra ['%s', %.5f, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f]" %(name.split()[0],float(ra), float(rerr), float(dec), float(derr), float(freq), float(flux), float(ferr))
			#print "extra ['%s', %.5f, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f]" %(name.split()[0],float(ra), float(rerr), float(dec), float(derr), float(freq), float(flux), float(ferr))
			ra = float(ra)
			dec = float(dec)
			rerr = float(rerr)
			derr = float(derr)
			if rerr==-1.0: rerr=0.0
			if derr==-1.0: derr=0.0
			ra_img,dec_img = KGS_wcs.wcs_world2pix(ra,dec,0)
			img_maj,img_min = KGS_wcs.wcs_world2pix(ra+float(major),dec+float(minor),0)
			img_maj = img_maj - ra_img
			img_min = img_min - dec_img
			ra_img_err,dec_img_err = KGS_wcs.wcs_world2pix(ra+rerr,dec+derr,0)
			ra_img_err = ra_img_err - ra_img
			dec_img_err = dec_img_err - dec_img
			#plot_pos('s',colour1,ra_img,dec_img,ra_img_err,dec_img_err,print_cat+' '+name,6,extra_plot,'')
			
			p = extra_plot.errorbar(ra_img,dec_img,ra_img_err,dec_img_err,marker='s',ms=6,mfc=colour1,mec=colour1,ecolor=colour2,markeredgewidth=1,label=name,linestyle='None')
			
			if print_cat=='mrc':
				if PA!='-100000.0': pol.plt_ell(ra_img,dec_img,float(img_maj),float(img_min),float(PA),extra_plot,colour2,colour2,0.15,1.0)
			elif print_cat=='comp_v5':
				pass
			else:
				if PA!='-100000.0': pol.plt_ell(ra_img,dec_img,float(img_maj),float(img_min),float(PA),extra_plot,colour2,colour2,0.3,1.0)
				
			if name not in plotted_names:
				handles.append(p,)
				plotted_names.append(name)

def plot_ind(match,ax_spectral,ra_bottom,ra_top,dec_bottom,dec_top,dom_crit,comb_crit):
	'''Takes a string of information of a particular combination and uses it
	to create a plot of the single combination, and fit and plot a line to
	the spectral information. Returns the positional probability, fitted paramaters
	and the residuals of the fit'''
	
	##Get the information from the particular match given
	info = match.split()
	indexes = [(14+((i-1)*3)) for i in num_freqs]
	starts = [0]
	for i in xrange(len(indexes)-1): starts.append(sum(indexes[:i+1]))
	fluxs = []
	freqs = []
	ferrs = []
	for j in xrange(len(starts)): 
		ind = starts[j]
		cat = info[ind]
		if cat!='-100000.0':
			freq = num_freqs[j]
			name = info[ind+1]
			ra = float(info[ind+2])
			rerr = float(info[ind+3])
			dec = float(info[ind+4])
			derr = float(info[ind+5])
			nu = float(info[ind+6])
			flux = float(info[ind+7])
			ferr = float(info[ind+8])/flux
			major = info[ind+9+((freq-1)*3)]
			minor = info[ind+10+((freq-1)*3)]
			PA = info[ind+11+((freq-1)*3)]
			fluxs.append(flux)
			freqs.append(nu)
			ferrs.append(ferr)
			##Plot each source on the individual combo plot
	##Sort the frequencies, fluxes and log them
	log_fluxs = np.log([flux for (freq,flux) in sorted(zip(freqs,fluxs),key=lambda pair: pair[0])])
	sorted_ferrs = np.array([ferr for (freq,ferr) in sorted(zip(freqs,ferrs),key=lambda pair: pair[0])])
	log_freqs = np.log(sorted(freqs))
	ferrs = np.array(ferrs)
	prob = info[-1]
	
	##Fit a line using weighted least squares and plot it
	lin_fit,jstat,bse,chi_red = mkl.fit_line(log_freqs,log_fluxs,sorted_ferrs)
	
	##Plot the fitted spectral line, and return the plot object so we can create a legend from it
	if dom_crit=='No dom. source':
		if 'split' in comb_crit:
			spec_plot = 'na'
		else:
			spec_plot, = ax_spectral.plot(np.exp(log_freqs),np.exp(lin_fit.fittedvalues),linestyle='-',linewidth=1,alpha=0.3)
	else:
		spec_plot, = ax_spectral.plot(np.exp(log_freqs),np.exp(lin_fit.fittedvalues),linestyle='-',linewidth=1,alpha=1)
	
	return spec_plot,lin_fit.params

def do_plot_image(comp,all_info,image_cats,image_files,present_cats,source_names,src_g,matches):
	##Work out how many plots across we need
	num_plots = len(image_files)+2
	
	fig = plt.figure(figsize = (20,10))

	##Find out the ra,dec of the base catalogue source
	info_base = all_info[0].split()
	ra_main = float(info_base[2])
	dec_main = float(info_base[4])
	
	main_dims = [0.08, 0.5, (1.0/num_plots)-0.03, 0.35]
	spec_dims = [0.08,0.1,(1.0/num_plots)-0.03,0.35]
	
	ax_main,ax_spectral,tr_fk5,wcs = pol.make_left_plots(fig,main_dims,spec_dims,ra_main,dec_main)
	
	##Find the limits out to search area - have to do each edge individual,
	##because of the arcdistance projection malarky
	delta_RA = np.arccos((np.cos((2*closeness)*dr)-np.sin(dec_main*dr)**2)/np.cos(dec_main*dr)**2)/dr
	
	plot_lim = (2*closeness) + (0.1/60.0)
	ra_up_lim = ra_main + delta_RA + (0.1/60.0)
	ra_down_lim = ra_main - delta_RA - (0.1/60.0)
	dec_up_lim = dec_main + plot_lim
	dec_down_lim = dec_main - plot_lim
	
	##Need the axes and world coordinate systems for each catalogue image to
	##plot on later
	over_plots = []
	wcss = []
	for cat in image_cats:
		im_ind = image_cats.index(cat)
		
		##Create two images, one to build a colorbar on, one to
		##plot the sources on
		ax_overplot = fig.add_subplot(2,num_plots,im_ind+3)
		ax_image = fig.add_subplot(2,num_plots,num_plots+im_ind+3)
		
		over_plots.append(ax_overplot)
		#ax_image.set_title('%s image' %cat, fontsize=18)
		#ax_overplot.set_title('%s sources' %cat, fontsize=18)
		
		ax_image.text( 0.02, 0.96, '%s image' %cat, transform=ax_image.transAxes,verticalalignment='center',fontsize=18,color='w')
		ax_overplot.text( 0.02, 0.96, '%s sources' %cat, transform=ax_overplot.transAxes,verticalalignment='center',fontsize=18,color='w')
		
		##Open the fits file - fits from nvss/vlssr website are
		##formatted differently
		image_file = pf.open(image_files[im_ind]) 
		image = image_file[0].data
		#if cat in ['vlssr','nvss']: image = image[0][0]
		
		##Get wcs and append to be used later
		image_wcs = astWCS.WCS(image_files[im_ind])
		wcss.append(image_wcs)
		image_file.close()
		
		##Plot images - plot the over_plot image with vmax of 1/4 of max flux - should show any faint/
		##sidelobe structure #vmin=0,vmax=np.max(image)/10.0
		overplot = ax_overplot.imshow(image,cmap='gist_gray',vmin=0,vmax=np.max(image)/4,interpolation='none',origin='lower')
		#overplot = ax_overplot.imshow(image,cmap='gray',vmin=0,vmax=0.05,interpolation='none',origin='lower')
		source_image = ax_image.imshow(image,cmap='gnuplot',interpolation='none',origin='lower')
		##source_image = ax_image.imshow(image,cmap='spectral',vmin=0,vmax=0.05,interpolation='none',origin='lower')
		
		##Make room for a colorbar and plot it with labels on top
		divider = make_axes_locatable(ax_image)
		cax_1 = divider.append_axes("top", size="5%", pad=0.1)
		plt.colorbar(source_image, cax=cax_1, orientation='horizontal')
		cax_1.xaxis.set_tick_params(labeltop='on')
		cax_1.xaxis.set_tick_params(labelbottom='off')
		cax_1.xaxis.set_major_locator(MaxNLocator(5))
		#cax_1 = divider.append_axes("right", size="5%", pad=0.1)
		#plt.colorbar(source_image, cax=cax_1)
		
		divider = make_axes_locatable(ax_overplot)
		cax_2 = divider.append_axes("top", size="5%", pad=0.1)
		plt.colorbar(overplot, cax=cax_2, orientation='horizontal')
		cax_2.xaxis.set_tick_params(labeltop='on')
		cax_2.xaxis.set_tick_params(labelbottom='off')
		cax_2.xaxis.set_major_locator(MaxNLocator(5))
		#cax_2 = divider.append_axes("right", size="5%", pad=0.1)
		#plt.colorbar(overplot, cax=cax_2)
		
		#label_colorbar = [[cax_1,num_plots+im_ind+2],[cax_2,im_ind+2]]
		
		#for cax,plot_num in label_colorbar:
			#if plot_num % num_plots == 0: cax.set_ylabel("Flux [Jy/Beam]")
		
		for ax in [ax_image,ax_overplot]:
			##Find limits of the axis for the plots, create a range of a smaller number,
			##and then use wcs to translate the image coords into RA and Dec labels
			limits = ax.axis()
			
			big_ra_high = image_wcs.pix2wcs(limits[0],(limits[3]-limits[2])/2)[0]
			big_ra_low = image_wcs.pix2wcs(limits[1],(limits[3]-limits[2])/2)[0]
			
			big_dec_low = image_wcs.pix2wcs((limits[1]-limits[0])/2,limits[2])[1]
			big_dec_high = image_wcs.pix2wcs((limits[1]-limits[0])/2,limits[3])[1]
			
			ax.set_xlim(limits[0],limits[1])
			ax.set_ylim(limits[2],limits[3])
			
			#x_range = np.arange(limits[0],limits[1],abs(limits[1]-limits[0])/4)
			#x_labels = [xtick_RA(image_wcs.pix2wcs(x,0)[0]) for x in x_range[1:]]
			#ax.set_xticks(x_range[1:])
			##ax.set_xticklabels(x_labels)
			
			#y_range = np.arange(limits[2],limits[3],abs(limits[3]-limits[2])/7)
			#y_labels = ['%.2f' %image_wcs.pix2wcs(0,y)[1] for y in y_range[1:]]
			#ax.set_yticks(y_range[1:])
			##ax.set_yticklabels(y_labels)
			
			ax.set_xticklabels([])
			ax.set_yticklabels([])
		#if im_ind==0: ax_overplot.set_ylabel('Dec [deg]')
		
		#ax_image.set_xlabel("RA [hh:mm:ss]")
		
		
	##Plot Pattis catalogue here
	KGS_file = pf.open("eor0_Restored_I_med.fits") 
	KGS_img = KGS_file[1].data
	
	from astropy.wcs import WCS
	
	header = {'BITPIX'  : -32 ,                         
	'NAXIS'   :                    2,
	'NAXIS1'  :                 3072,
	'NAXIS2'  :                 3072,
	'CTYPE1'  : 'RA---SIN',                                                 
	'CTYPE2'  : 'DEC--SIN',
	'EQUINOX' :               2000.0,
	'CD1_1'   :            0.0260417,
	'CD2_1'   :                  0.0,
	'CD1_2'   :                  0.0,
	'CD2_2'   :            0.0260417,
	'CRPIX1'  :               1537.0,
	'CRPIX2'  :               1537.0,
	'CRVAL1'  :         359.84942627,
	'CRVAL2'  :       -26.7836437225,
	'LONPOLE' :                180.0,                                                 
	'LATPOLE' :                  0.0,                                                
	'PV2_1'   :    1.70107273334E-06,
	'PV2_2'   :   -9.51798142523E-08,
	'PV1_0'   :                  0.0,
	'PV1_1'   :                  0.0,
	'PV1_2'   :                 90.0,
	'PV1_3'   :                180.0,
	'PV1_4'   :                  0.0,
	'RADESYS' : 'FK5     ',
	'BUNIT'   : 'Jy/beam ',
	}
	
	KGS_wcs = WCS(header=header)
	
	KGS_RA_up = KGS_wcs.wcs_world2pix(big_ra_high,dec_main,0)[0]
	KGS_RA_down = KGS_wcs.wcs_world2pix(big_ra_low,dec_main,0)[0]
	KGS_DEC_up = KGS_wcs.wcs_world2pix(ra_main,big_dec_high,0)[1]
	KGS_DEC_down = KGS_wcs.wcs_world2pix(ra_main,big_dec_low,0)[1]
	
	ax_KGS_overplot = fig.add_subplot(2,num_plots,2)
	ax_KGS_image = fig.add_subplot(2,num_plots,num_plots+2)

	KGS_scale = KGS_img[KGS_DEC_down:KGS_DEC_up+1,KGS_RA_down:KGS_RA_up+1]
	
	KGS_overplot = ax_KGS_overplot.imshow(KGS_img,cmap='gist_gray',vmin=0,vmax=np.max(KGS_scale)/4,interpolation='none',origin='lower')
	KGS_image = ax_KGS_image.imshow(KGS_img,cmap='gnuplot',vmin=0,vmax=np.max(KGS_scale),interpolation='none',origin='lower') #vmin=0,vmax=np.max(KGS_img)/4
	
	ra_cent,dec_cent = KGS_wcs.wcs_world2pix(ra_main,dec_main,0)
	
	for ax,img in zip([ax_KGS_image,ax_KGS_overplot],[KGS_image, KGS_overplot]):
		
		ax.set_xlim(KGS_RA_down,KGS_RA_up)
		ax.set_ylim(KGS_DEC_down,KGS_DEC_up)
		
		divider = make_axes_locatable(ax)
		cax = divider.append_axes("top", size="5%", pad=0.1)
		plt.colorbar(img, cax=cax,orientation='horizontal')
		cax.xaxis.set_tick_params(labeltop='on')
		cax.xaxis.set_tick_params(labelbottom='off')
		
		cax.xaxis.set_major_locator(MaxNLocator(5))
		
		ax.set_xticklabels([])
		ax.set_yticklabels([])
		
		ax.invert_xaxis()
		
		ax.text( 0.02, 0.94, 'KGS image', transform=ax.transAxes,verticalalignment='center',fontsize=18,color='w')

	##Get the information and plot the positions on the over_plots
	ras = []
	decs = []
	all_freqs = []
	all_fluxs = []
	all_ferrs = []
	
	for i in xrange(len(all_info)):
		info=all_info[i].split()
		cat = info[0]
		name = info[1]
		ra = float(info[2])
		ras.append(ra)
		rerr = float(info[3])
		dec = float(info[4])
		decs.append(dec)
		derr = float(info[5])
		major = info[-5]
		minor = info[-4]
		PA = info[-3]
		ID = info[-1]
		
		if len(info)==14:
			all_fluxs.append(float(info[7]))
			all_freqs.append(float(info[6]))
		##If the catalogue has more than one frequency:
		else:
			extra = (len(info)-14) / 3
			for i in xrange(extra+1):
				all_freqs.append(float(info[6+(3*i)]))
				all_fluxs.append(float(info[7+(3*i)]))
		
		##Plot positions and elliptical fits on bottom left plot
		#plot_all(cat,name,ra,rerr,dec,derr,major,minor,PA,ax_main)
		##If catalogue in one of the images, plot the informaion on the applicaable over_plot
		#if cat=='mrc':
			#for i in xrange(len(over_plots)): plot_all_image(cat,name,ra,rerr,dec,derr,major,minor,PA,over_plots[i],wcss[i])
		if cat=='comp_v5':
			for i in xrange(len(over_plots)): 
				plot_all_image(cat,name,ra,rerr,dec,derr,major,minor,PA,over_plots[i],wcss[i])
				plot_all_image(cat,name,ra,rerr,dec,derr,major,minor,PA,ax_KGS_overplot,KGS_wcs)
				
		if cat in image_cats:
			cat_ind = image_cats.index(cat)
			plot_all_image(cat,name,ra,rerr,dec,derr,major,minor,PA,over_plots[cat_ind],wcss[cat_ind])
			
			
	##Use these to put in the legend later
	handles = []
	plotted_names = []
	
	###If already exist as a HACK, don't print out loads of extra info
	present_already = 'no'
	#print info_base[1]
	#if info_base[1] in HACK_sources: present_already = 'yes'

	
	find_extra_fhd('mwacs',big_ra_low,big_ra_high,big_dec_low,big_dec_high,ax_KGS_overplot,mwacs_data,source_names,'c','c',handles,plotted_names,present_already,KGS_wcs)
	find_extra_fhd('comp_v5',big_ra_low,big_ra_high,big_dec_low,big_dec_high,ax_KGS_overplot,fhd_data,source_names,'#8B008B','#8B008B',handles,plotted_names,present_already,KGS_wcs)
	
	##Plot cat entries that weren't in the original match
	if 'sumss' in image_cats:
		find_extra('sumss','mwacs',mwacs_data,image_cats,over_plots,wcss,source_names,ax_spectral,'c','c',180.0,handles,plotted_names,present_already,all_fluxs)
		find_extra('sumss','comp_v5',fhd_data,image_cats,over_plots,wcss,source_names,ax_spectral,'#8B008B','#8B008B',182.435,handles,plotted_names,present_already,all_fluxs)
		find_extra('sumss','sumss',sumss_data,image_cats,over_plots,wcss,source_names,ax_spectral,'#FF0033','#990000',843.0,handles,plotted_names,present_already,all_fluxs)
		#find_extra('sumss','mrc',mrc_data,image_cats,over_plots,wcss,source_names,ax_spectral,'y','y',408.0,handles,plotted_names,present_already,all_fluxs)
		#find_extra('sumss','A154',A154_data,image_cats,over_plots,wcss,source_names,ax_spectral,'c','c',154.0,handles,plotted_names,present_already,all_fluxs)
		#find_extra('sumss','A182',A182_data,image_cats,over_plots,wcss,source_names,ax_spectral,'#CCFF00','#CCFF00',182.0,handles,plotted_names,present_already,all_fluxs)
		
		
	if 'vlssr' in image_cats:
		find_extra('vlssr','mwacs',mwacs_data,image_cats,over_plots,wcss,source_names,ax_spectral,'c','c',180.0,handles,plotted_names,present_already,all_fluxs)
		find_extra('vlssr','comp_v5',fhd_data,image_cats,over_plots,wcss,source_names,ax_spectral,'#8B008B','#8B008B',182.435,handles,plotted_names,present_already,all_fluxs)
		find_extra('vlssr','vlssr',vlssr_data,image_cats,over_plots,wcss,source_names,ax_spectral,'#0099FF','#000099',74.0,handles,plotted_names,present_already,all_fluxs)
		#find_extra('vlssr','mrc',mrc_data,image_cats,over_plots,wcss,source_names,ax_spectral,'y','y',408.0,handles,plotted_names,present_already,all_fluxs)
		#find_extra('vlssr','A154',A154_data,image_cats,over_plots,wcss,source_names,ax_spectral,'c','c',154.0,handles,plotted_names,present_already,all_fluxs)
		#find_extra('vlssr','A182',A182_data,image_cats,over_plots,wcss,source_names,ax_spectral,'#CCFF00','#CCFF00',182.0,handles,plotted_names,present_already,all_fluxs)
		
	if 'nvss' in image_cats:
		find_extra('nvss','mwacs',mwacs_data,image_cats,over_plots,wcss,source_names,ax_spectral,'c','c',180.0,handles,plotted_names,present_already,all_fluxs)
		find_extra('nvss','comp_v5',fhd_data,image_cats,over_plots,wcss,source_names,ax_spectral,'#8B008B','#8B008B',182.435,handles,plotted_names,present_already,all_fluxs)
		find_extra('nvss','nvss',nvss_data,image_cats,over_plots,wcss,source_names,ax_spectral,'#339900','#00CC66',1400.0,handles,plotted_names,present_already,all_fluxs)
		#find_extra('nvss','mrc',mrc_data,image_cats,over_plots,wcss,source_names,ax_spectral,'y','y',408.0,handles,plotted_names,present_already,all_fluxs)
		#find_extra('nvss','A154',A154_data,image_cats,over_plots,wcss,source_names,ax_spectral,'c','c',154.0,handles,plotted_names,present_already,all_fluxs)
		#find_extra('nvss','A182',A182_data,image_cats,over_plots,wcss,source_names,ax_spectral,'#CCFF00','#CCFF00',182.0,handles,plotted_names,present_already,all_fluxs)
		
	
	
	meh,num_matches,accept_matches,accepted_inds,accept_type,stage = stats.split()
	accepted_inds = map(int,accepted_inds.split(','))
	
	fig2 = pol.create_plot_fig(comp,accepted_inds,stage,accept_type)
	
	bbox_props = dict(boxstyle="round", fc="w", ec="k", lw=1, alpha=0.7)
	ax_main.text( 0.025, 0.045, 'PUMA: %s %s' %(accept_type,stage), transform=ax_main.transAxes,verticalalignment='center',fontsize=18,color='k', bbox=bbox_props)
	
	spec_labels = []
	SIs = []
	
	comb_sources = 'blah'
	
	##If no repeated catalogues to combine, skip
	if num_matches==0 or num_matches==1:
		pass
	##Otherwise, plot the combined fluxes
	else:
		##Calculate and plot the combined fluxes of the two sources, even if source one or two has been accepted
		##just as a guide
		src_all = mkl.get_allinfo(all_info)
		
		if accepted_inds=='Nope':
			pass
		else:
			comb_crit, ra_ws, rerr_ws, dec_ws, derr_ws, temp_freqs, comb_freqs, comb_fluxs, comb_ferrs, comb_fit, comb_jstat, comb_chi_red, combined_names, set_freqs, set_fluxs, set_fits, comb_sources = mkl.combine_flux(src_all,src_g,accepted_inds,'plot=yes',len(matches))
		
		##If the criteria sent the double to be combined, actually plot the fitted line
		if 'combine' in stage or 'split' in stage:
			
			split_colors = ['#AE70ED','#FFB60B','#62A9FF','#59DF00']
			for fit in set_fits:
				ind = set_fits.index(fit)
				ax_spectral.plot(set_freqs[ind],set_fluxs[ind],linestyle='--',linewidth=1.0,color=split_colors[ind],alpha=0.7)
				split_p, = ax_spectral.plot(temp_freqs,np.exp(fit.params[1] + np.log(temp_freqs)*fit.params[0]),linestyle='-',linewidth=1.5,color=split_colors[ind])
				spec_labels.append(split_p)
				SIs.append([fit.params[0],'split %d' %(ind+1)])
			
			bright_colours = ['#FF6600','#33FF33','#FF47A3']
			
			for freq in xrange(len(comb_freqs)):
				pol.plot_errors('*',bright_colours[freq],comb_freqs[freq],comb_fluxs[freq],comb_ferrs[freq],'combo',9,ax_spectral)
			comb_p, = ax_spectral.plot(temp_freqs,np.exp(comb_fit.fittedvalues),linestyle='--',linewidth=1.5,color='k')
			spec_labels.append(comb_p)
			SIs.append([comb_fit.params[0],'comb'])
			
			##Send the combined fluxes to the all_fluxs so that ax_spectral is scaled appropriately
			for flux in comb_fluxs:
				all_fluxs.append(flux)
			
			for pos in xrange(len(ra_ws)):
				patch = pol.plot_pos('*',bright_colours[pos],ra_ws[pos],dec_ws[pos],rerr_ws[pos],derr_ws[pos],combined_names[pos],10,ax_main,ax_main.get_transform("fk5"))
		
		if len(matches)<11:
			for match in matches:
				ind = matches.index(match)
				
				spec_plot,params = plot_ind(match,ax_spectral,ra_down_lim,ra_up_lim,dec_down_lim,dec_up_lim,'No dom. source',comb_crit)
				if spec_plot=='na':
					pass
				else:
					SIs.append([params[0],str(ind+1)])
					spec_labels.append(spec_plot)
		else:
			pass
			
		
	
	##==============================================================

	##Fill the left hand plots with information goodness
	
	pol.fill_left_plots(all_info,ra_main,dec_main,ax_main,ax_spectral,tr_fk5,wcs,all_fluxs,ra_down_lim,ra_up_lim,dec_down_lim,dec_up_lim,delta_RA)
	
	##Make room at the top of the plot for a legend for ax_main, make the legend
	fig.subplots_adjust(top=0.85)
	
	leg_labels = [r'$\alpha_{%s}$ = %.2f' %(SI[1],SI[0]) for SI in SIs]
	main_handles,main_labels = ax_main.get_legend_handles_labels()
	
	for hand,lab in zip(handles,plotted_names):
		main_handles.append(hand)
		main_labels.append(lab)
	
	main_leg = fig.add_axes([0.05,0.87,0.9,0.05])
	main_leg.axis('off')
	main_leg.legend(main_handles,main_labels,loc='lower center',prop={'size':11},ncol=9) #,bbox_to_anchor=(0,1.02),
	
	if len(spec_labels)<4:
		ax_spectral.legend(spec_labels,leg_labels,loc='best',prop={'size':14},fancybox=True)
	else:
		ax_spectral.legend(spec_labels[0:1],leg_labels[0:1],loc='best',prop={'size':14},fancybox=True)
		
	#fig.tight_layout()
	
	return fig,fig2,comb_sources
		
##======================================================================================
##Check whether the source is already a HACK

##Catalogue data
mrc_data = open('mrc_simple.txt' ).read().split('\n')
sumss_data = open('sumss_simple.txt' ).read().split('\n')
vlssr_data = open('vlssr_simple.txt' ).read().split('\n')
nvss_data = open('nvss_simple.txt' ).read().split('\n')
mwacs_data = open('mwacs_simple.txt' ).read().split('\n')
fhd_data = open('simple_comp_v5.txt').read().split('\n')
	
##======================================================================================

reject_pile = open('reject_pile.txt','w+')
choice_log = open('choice_log.txt','w+')
sources = []
sources_stats = []


##Input eyeball file
bayes_comp = open(options.input_bayes).read().split('END_GROUP')
del bayes_comp[-1]

for comp in bayes_comp:

	##Get the information into nice usable forms, and get rid of empty/pointless
	##entries
	chunks = comp.split('START_COMP')
	all_info = chunks[0].split('\n')
	
	##FOR SOME REASON CAN'T DO BOTH OF THESE LINES IN THE SAME FOR LOOP?!?!?!
	for entry in all_info:   
		if entry=='': del all_info[all_info.index(entry)]
	for entry in all_info:
		if 'START' in entry: del all_info[all_info.index(entry)]
		
	matches = chunks[1].split('\n')
	
	del matches[0],matches[-1]
	stats = matches[-1]
	del matches[-2:]
	
	match1 = matches[0].split()
	src_g = mkl.get_srcg(match1)
	
	
	##Get some info and find which catalogues are present 
	src_all = mkl.get_allinfo(all_info)
	present_cats = [cat for cat in src_all.cats if cat!='-100000.0']
	
	#if 330.0<=src_all.ras[0]<=360.0 or 0.0<=src_all.ras[0]<=30.0: #if 30.0<=src_all.ras[0]<=90.0: #  #if 330.0<=src_all.ras[0]<=360.0: if 0.0<=src_all.ras[0]<=30.0
		#if -57.0<=src_all.decs[0]<=3.0:
		
	##Currently only have images from vlssr, sumss and nvss. These are saved in a specifc format so can be found
	##if they exist with os.path.exist
	image_num = bayes_comp.index(comp)+1
			
	image_cats = [cat for cat in ['vlssr','sumss','nvss'] if os.path.exists("./extended_fits/%s_%s.fits" %(src_all.names[0],cat))==True]
	#print image_cats
	image_files =["./extended_fits/%s_%s.fits" %(src_all.names[0],cat) for cat in ['vlssr','sumss','nvss'] if os.path.exists("./extended_fits/%s_%s.fits" %(src_all.names[0],cat))==True]

	#Fiddle this number to plot certain information
	i=0
	
	#if src_all.names[0] in HACK_sources:
		#THIS_HACK = HACKS[HACK_sources.index(src_all.names[0])]
		#print "Already exists as %s EXT%s%s " %(THIS_HACK.hack_name,THIS_HACK.cat_name,THIS_HACK.cat_index)
		##Do the image plots
		#do_plot_image(all_info,image_cats,image_files,present_cats,src_all.names,src_g,matches)
		#if src_all.names[0] in mwacs_names: index = mwacs_names.index(src_all.names[0])
		#elif src_all.names[0] in mrc_names: index = mrc_names.index(src_all.names[0])
			
		##print image_num
		##print "(EXT%s%s) %s %s '%s' '%s' ID %s %.5f Jy" %(src_all.cats[i],index,src_all.cats[i],src_all.names[i],deg_to_hour(src_all.ras[i],'info'),deg_to_degmins(src_all.decs[i],'info'),src_all.IDs[i], src_all.fluxs[i])
			
		##for i in xrange(len(src_all.cats)):
			##print "matched ['%s', %.5f, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f]" %(src_all.names[i],src_all.ras[i],src_all.rerrs[i],src_all.decs[i],src_all.derrs[i], src_all.freqs[i][0], src_all.fluxs[i][0],src_all.ferrs[i][0])
		###if type(comb_source)!=str:
			###for i in xrange(len(comb_source.cats)):
				###print "combined '%s' %.5f %.5f %.5f %.5f %.5f %.5f %.5f" %(comb_source.names[i],comb_source.ras[i],comb_source.rerrs[i],comb_source.decs[i],comb_source.derrs[i], comb_source.freqs[i][0],comb_source.fluxs[i][0],comb_source.ferrs[i][0])
		####THIS CAUSES WINDOW TO AUTOMATICALLY MAXIMISE
		#mng = plt.get_current_fig_manager()
		#mng.resize(*mng.window.maxsize())
		###print plt.get_backend()
		

		#plt.show()
	#else:
	
	#if src_all.names[0] in mwacs_names: index = mwacs_names.index(src_all.names[0])
	#elif src_all.names[0] in mrc_names: index = mrc_names.index(src_all.names[0])
	
	#exisiting_hacks.append("EXT%s%s" %(src_all.cats[i],index))
	
	#oot.write("\nEXT%s%s" %(src_all.cats[i],index))
	#for image in image_files:
		#oot.write(" %s" %image.split("/")[-1])
	
	
	
	##if "EXT%s%s" %(src_all.cats[i],index) in hack_names:
		##pass
	##else:
	
	
	
	
	
		##Do the image plots
	
	##if src_all.names[0] in mwacs_names: index = mwacs_names.index(src_all.names[0])
	##elif src_all.names[0] in mrc_names: index = mrc_names.index(src_all.names[0])
		
	#print image_num
	#print "(EXT%s%s) %s %s '%s' '%s' ID %s %.5f Jy" %(src_all.cats[i],index,src_all.cats[i],src_all.names[i],deg_to_hour(src_all.ras[i],'info'),deg_to_degmins(src_all.decs[i],'info'),src_all.IDs[i], src_all.fluxs[i])
		
	#for i in xrange(len(src_all.cats)):
		#print "matched ['%s', %.5f, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f]" %(src_all.names[i],src_all.ras[i],src_all.rerrs[i],src_all.decs[i],src_all.derrs[i], src_all.freqs[i][0], src_all.fluxs[i][0],src_all.ferrs[i][0])
	###if type(comb_source)!=str:
		###for i in xrange(len(comb_source.cats)):
			###print "combined '%s' %.5f %.5f %.5f %.5f %.5f %.5f %.5f" %(comb_source.names[i],comb_source.ras[i],comb_source.rerrs[i],comb_source.decs[i],comb_source.derrs[i], comb_source.freqs[i][0],comb_source.fluxs[i][0],comb_source.ferrs[i][0])
	
	##print plt.get_backend()
	
	####THIS CAUSES WINDOW TO AUTOMATICALLY MAXIMISE
	#mng = plt.get_current_fig_manager()
	#mng.resize(*mng.window.maxsize())
	#plt.show()
	
	
	flags = src_all.flags[0].split('+')
	#fig = do_plot_image(comp,all_info,image_cats,image_files,present_cats,src_all.names,src_g,matches)
	
	#if flags[1]!='0.0':
		#fig,fig2 = do_plot_image(comp,all_info,image_cats,image_files,present_cats,src_all.names,src_g,matches)
		#fig.savefig('%s_KGS.png' %src_all.names[0],bbox_inches='tight',dpi=200)
		#print src_all.names[0]
	#else:
		#fig.savefig('%s_KGS+PUMA.png' %src_all.names[0],bbox_inches='tight',dpi=200)
		#print src_all.names[0]
		#print "--------------------------------------------------------------------------------------------------------------------------------------"
	#else:
		#pass
	#plt.close()
	
	fig,fig2,comb_sources = do_plot_image(comp,all_info,image_cats,image_files,present_cats,src_all.names,src_g,matches)
	
	plt.show(block=False)
	answer = raw_input('Enter either number of accepted combination, letter "c" to combine sources, or "r" to reject sources: ')
	print "--------------------------------"
	#answer='r'
	choice_log.write("%s %s\n" %(src_all.names[0],str(answer)))
	
	meh,num_matches,accept_matches,accepted_inds,accept_type,stage = stats.split()
	g_stats = mkl.group_stats()
	g_stats.num_matches = num_matches
	g_stats.retained_matches = accept_matches
	
	num_matches = len(matches)
	if answer == 'r':
		plt.close()
		plt.close()
		comp = comp[1:]
		reject_pile.write(comp)
		reject_pile.write('END_GROUP\n')
	elif answer == 'c':
		plt.close()
		plt.close()
		sources.append(comb_sources[0])
		g_stats.accept_type = 'user-combine'
		sources_stats.append(g_stats)
	elif answer == 'quit':
		plt.close()
		plt.close()
		break
	else:
		accepted_match = matches[int(answer)-1].split()
		#print accepted_match
		jstat_resids,params,bses,chi_resids = mkl.calculate_resids([accepted_match])
		source = mkl.get_srcg(accepted_match)
		source.SI = params[0][0]
		source.intercept = params[0][1]
		source.SI_err = bses[0][0]
		source.intercept_err = bses[0][1]
		if chi_resids[0]<=2:
			source.low_resids = 0
		else:
			source.low_resids = 1
		sources.append(source)
		g_stats.accept_type = 'user-single'
		sources_stats.append(g_stats)
		plt.close()
		plt.close()
		
updated_ras = []
updated_decs = []
updated_rerrs = []
updated_derrs = []

for source in sources:
	##this is for choosing the most 'reliable' position as defined by
	##the input preference cats
	for p_cat in pref_cats:
		if p_cat in source.cats:
			pos_ind = source.cats.index(p_cat)
			break
	else:
		print 'no catalogue match in positional preferences - OMG'
	updated_ras.append(source.ras[pos_ind])
	updated_decs.append(source.decs[pos_ind])
	updated_rerrs.append(source.rerrs[pos_ind])
	updated_derrs.append(source.derrs[pos_ind])

##Make an 'KGS' name based on position
def make_name(ra,dec,ext):
	ra_str = mkl.deg_to_hour(ra)
	dec_str = mkl.deg_to_degmins(dec)
	if ext=='':
		return 'KGS'+ra_str+dec_str
	else:
		return 'KGS'+ra_str+dec_str+ext


original_ras = [source.ras[0] for source in sources]
original_decs = [source.decs[0] for source in sources]
type_matches = [source.accept_type for source in sources_stats]

##Name the sources based on position. If the sources were the result of a split,
##name by the position of the base source, and add the letter that was defined by the
##split letter in the stat report object
names = []
for i in xrange(len(type_matches)):
	if 'split' in type_matches[i]:
		ext = type_matches[i][-1]
		name = make_name(original_ras[i],original_decs[i],ext)
	else:
		name = make_name(updated_ras[i],updated_decs[i],'')
	names.append(name)

t=atpy.Table(masked=True)

t.add_column('Name',names,description='Name based on position of combined source')

##Add lots of data
prim_names = [source.names[0] for source in sources]
t.add_column('%s_name' %matched_cats[0],prim_names,description='Name of %s component' %matched_cats[0])
t.add_column('updated_RA_J2000',np.array(updated_ras),description='Updated Right ascension of source',unit='deg')
t.add_column('updated_DEC_J2000',np.array(updated_decs),description='Updated Declination of source',unit='deg')
t.add_column('updated_RA_err',np.array(updated_rerrs),description='Error on Updated Right ascension of source',unit='deg')
t.add_column('updated_DEC_err',np.array(updated_derrs),description='Error on Updated Declination of source',unit='deg')
t.add_column('original_RA_J2000',np.array(original_ras),description='Original Right ascension of source',unit='deg')
t.add_column('original_DEC_J2000',np.array(original_decs),description='Original Declination of source',unit='deg')

##For every catalogue in the match
for cat in xrange(len(num_freqs)):
	##See how many frequencies that source has
	num_freq = num_freqs[cat]
	##For every frequency, make a column of fluxes and flux errors, masking every value with -100000.0
	for freq in xrange(num_freq):
		fluxs = np.array([src.fluxs[cat][freq] for src in sources])
		ferrs = np.array([src.ferrs[cat][freq] for src in sources])
		t.add_column('S_%.1f' %float(cat_freqs[cat][freq]),fluxs,description='Flux at %.1fMHz' %float(cat_freqs[cat][freq]),mask=fluxs==-100000.0, fill=None,unit='Jy')
		t.add_column('e_S_%.1f' %float(cat_freqs[cat][freq]),ferrs,description='Flux error at %.1fMHz' %float(cat_freqs[cat][freq]),mask=ferrs==-100000.0, fill=None,unit='Jy')

##Extrapolate the base catalogue frequency flux density, using the fitted values for late comparison
extrap_freq = cat_freqs[0][0]
		
def extrap(freq,SI,intercept):
	return np.exp((np.log(freq)*SI)+intercept)
	
def extrap_error(freq,ext_flux,int_err,SI_err):
	denom = ext_flux**2
	numer = int_err**2 + ((np.log(freq)**2)*(SI_err)**2)
	return np.sqrt(numer/denom)
	
SIs = np.array([source.SI for source in sources])
SI_errs = np.array([source.SI_err for source in sources])
intercepts = np.array([source.intercept for source in sources])
intercept_errs = np.array([source.intercept_err for source in sources])

##create a mask for all of the errors reported for matched with only two cats - can't have an error on
##a fit to just two data points, so statsmodel spits out nonsenseCatalogue & Frequency (MHz) & Number of sources & Sky Coverage
SI_err_mask = [err==float('inf') or err==float('nan') for err in SI_errs]

t.add_column('SI',SIs,description='Spectral Index of Fit')
t.add_column('e_SI',SI_errs,description='Std error on Spectral Index of Fit',mask=SI_err_mask,fill=None)
t.add_column('Intercept',intercepts,description='Intercept of Fit')
t.add_column('e_Intercept',intercept_errs,description='Std error on Intercept of Fit',mask=SI_err_mask,fill=None)

extrap_base = np.array([extrap(extrap_freq,source.SI,source.intercept) for source in sources])
extrap_base_err = np.array([extrap_error(extrap_freq,extrap_flux,source.intercept_err,source.SI_err) for extrap_flux, source in zip(extrap_base,sources)])

t.add_column('S_%.1f_ext' %extrap_freq,extrap_base,unit='Jy',description='Flux at extrapolted to base catalogue frequency using fitted values')
t.add_column('e_S_%.1f_ext' %extrap_freq,extrap_base_err,unit='Jy',description='Error on flux extrapolted to base frequency using error on fitted values',mask=SI_err_mask,fill=None)

##See how many unique catalogues there are
num_cats = [len(set([cat for cat in source.cats if cat!='-100000.0'])) for source in sources]
t.add_column('Number_cats',np.array(num_cats),description='Number of matched catalogues')

numb_matches = [source.num_matches for source in sources_stats]
t.add_column('Number_matches',np.array(numb_matches),description='Number of possible combinations in the group match')

retained_matches = [source.retained_matches for source in sources_stats]
t.add_column('Retained_matches',np.array(retained_matches),description='Number of retained (by position) combinations in the group match')

type_matches = [source.accept_type for source in sources_stats]
t.add_column('Type_matches',np.array(type_matches),description='At which stage the match was accepted')

low_resids = [source.low_resids for source in sources]
t.add_column('Low_resids',np.array(low_resids),description='Whether the fitted data had residuals above the threshold of chi_reduced<=2. 1 is above threshold, 0 otherwise')

flags = [source.flags[0] for source in sources]
t.add_column('Flags',np.array(flags),description='Any flags propagated through')

t.write("%s.vot" %options.output_name,votype='ascii',overwrite=True)
