import numpy as np
import optparse
import subprocess
import os

parser = optparse.OptionParser()

parser.add_option('-i', '--input_bayes', 
	help='Enter name of eyeball bayes file')

options, args = parser.parse_args()


def dec_to_deg(time): 
	'''converts dd:mm:ss.ss in to degrees, must input as a string
	   returns a float in units of degrees'''
	negtest=time[0]
	time=time.split(':')
	degr=float(time[0])
	mins=float(time[1])*(1.0/60.0)
	secs=float(time[2])*(1.0/(3600.0))
	if negtest=='-':
		deg=degr-mins-secs
	if negtest!='-':
		deg=degr+mins+secs
	return deg

##USED TO CREATE NEW SOURCE NAMES
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

class source_group:
	def __init__(self):
		self.cats = []
		self.names = []
		self.ras = []
		self.rerrs = []
		self.decs = []
		self.derrs = []
		self.freqs = []
		self.fluxs = []
		self.ferrs = []
		self.majors = []
		self.minors = []
		self.PAs = []
		self.flags = []
		self.IDs = []
		self.SI = None
		self.intercept = None
		self.prob = None
		self.num_matches = None
		self.type_match = None
		self.SI_err = None
		self.intercept_err = None
		self.low_resids = None

##INFORMAION GATHERING FUCNTIONS-------------------------------------------------------------------------------------
def get_allinfo(all_info):
	'''Takes a list of strings. Each string is a line containing the information for a single source
	in a matched group from an output file of calculate_bayes.py. Gets all of the information
	from each entry and returns them in a source_group() class'''
	src_all = source_group()
	for entry in all_info:
		info = entry.split()
		src_all.cats.append(info[0])
		src_all.names.append(info[1])
		src_all.ras.append(float(info[2]))
		src_all.rerrs.append(float(info[3]))
		src_all.decs.append(float(info[4]))
		src_all.derrs.append(float(info[5]))
		src_all.majors.append(float(info[-5]))
		src_all.minors.append(float(info[-4]))
		src_all.PAs.append(float(info[-3]))
		src_all.flags.append(info[-2])
		src_all.IDs.append(info[-1])
		
		##If the source only has one frequency. Append as an array
		##so that all entries to src_all.freqs etc are of the same
		##type. This deals with cats with multiple freqs, otherwise
		##the position, name etc will have to be repeated for each
		##frequency
		if len(info)==14:
			src_all.freqs.append(np.array([float(info[6])]))          
			src_all.fluxs.append(np.array([float(info[7])]))
			src_all.ferrs.append(np.array([float(info[8])]))
			#src_all.freqs.append(float(info[6]))  ##Left here in case ever want to append just the freq, not
			#src_all.fluxs.append(float(info[7]))  ##an array
			#src_all.ferrs.append(float(info[8]))
			
		##If not, work out how many freqs there are and append to lists
		else:
			extra = (len(info)-14) / 3
			freqs = []
			fluxs = []
			ferrs = []
			for i in xrange(extra+1):
				freqs.append(float(info[6+(3*i)]))
				fluxs.append(float(info[7+(3*i)]))
				ferrs.append(float(info[8+(3*i)]))
			src_all.freqs.append(np.array(freqs))
			src_all.fluxs.append(np.array(fluxs))
			src_all.ferrs.append(np.array(ferrs))
				#src_all.freqs.append(float(info[6+(3*i)]))
				#src_all.fluxs.append(float(info[7+(3*i)]))
				#src_all.ferrs.append(float(info[8+(3*i)]))
	return src_all
	
def get_srcg(info):
	'''Takes a string which contains the information for all sources in particular combination
	Uses num_freqs to work out where each piece of information is, and then return the relevant
	information in a source_group() class. Any catalogues with no matches are entered as -10000.0'''
	src_g = source_group()
	##Work out where in the string each different catalogue will start, using num_freqs
	##to work out how many entries each catalogue will have
	indexes = [(14+((i-1)*3)) for i in num_freqs]
	starts = [0]
	for i in xrange(len(indexes)-1): starts.append(sum(indexes[:i+1]))
	#all_freqs = []
	#all_fluxs = []
	#all_ferrs = []
	for j in xrange(len(starts)): 
		num_freq = num_freqs[j]
		ind = starts[j]
		#cat = info[ind]
		freqss = []
		fluxss = []
		ferrss = []
		for k in xrange(num_freq):
			freqss.append(float(info[6+ind+(3*k)]))
			fluxss.append(float(info[7+ind+(3*k)]))
			ferrss.append(float(info[8+ind+(3*k)]))
			#if float(info[6+ind+(3*k)])!=-100000.0: all_freqs.append(float(info[6+ind+(3*k)]))
			#if float(info[7+ind+(3*k)])!=-100000.0: all_fluxs.append(float(info[7+ind+(3*k)]))
			#if float(info[8+ind+(3*k)])!=-100000.0: all_ferrs.append(float(info[8+ind+(3*k)]))
		src_g.freqs.append(freqss)
		src_g.fluxs.append(fluxss)
		src_g.ferrs.append(ferrss)
		src_g.cats.append(info[ind])
		src_g.names.append(info[ind+1])
		src_g.ras.append(float(info[ind+2]))
		src_g.rerrs.append(float(info[ind+3]))
		src_g.decs.append(float(info[ind+4]))
		src_g.derrs.append(float(info[ind+5]))
		src_g.majors.append(float(info[ind+9+((num_freq-1)*3)]))
		src_g.minors.append(float(info[ind+10+((num_freq-1)*3)]))
		src_g.PAs.append(float(info[ind+11+((num_freq-1)*3)]))
		src_g.flags.append(info[ind+12+((num_freq-1)*3)])
		src_g.IDs.append(info[ind+13+((num_freq-1)*3)])
		src_g.prob = float(info[-1])
	return src_g



##Input eyeball file
bayes_comp = open(options.input_bayes).read().split('END_GROUP')
del bayes_comp[-1]

change_names = []

for comp in bayes_comp:
#for comp in bayes_comp[1:2]:
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
	
	##Get some info and find which catalogues are present 
	src_all = get_allinfo(all_info)
	present_cats = [cat for cat in src_all.cats if cat!='-100000.0']
	
	flags = src_all.flags[0].split('+')
	
	if flags[1]!='0.0':
		print flags
	
	##if 330.0<=src_all.ras[0]<=360.0 or 0.0<=src_all.ras[0]<=90.0: #30.0<=src_all.ras[0]<=90.0 or 
		##if -57.0<=src_all.decs[0]<=3.0:
	
		##Currently only want images from vlssr, sumss and nvss 
		image_cats = [cat for cat in ['vlssr','sumss','nvss'] if cat in present_cats]
		
		##Test to see if images have already been downloaded
		image_locs = "/home/jline/Documents/cataloguing/FHD/carroll/paper/complete/extended_fits"
		image_files =[1 for cat in ['vlssr','sumss','nvss'] if os.path.exists("%s/%s_%s.fits" %(image_locs,src_all.names[0],cat))==True]
		
		if len(image_files)==0:
			position = "position='%s, %s'" %(deg_to_hour(src_all.ras[0],'info'),deg_to_degmins(src_all.decs[0],'info'))
			size = 'size=0.3'
			#size = 'pixels=1000'
			coords = 'coordinates=J2000'
			proj = 'projection=Sin'
			outname = 'output=%s.fits' %src_all.names[0]
			surveys = 'survey='
			for name in image_cats: surveys+=name+','
			surveys = surveys[:-1]
			#scale='scale=0.00002777'
			
			sky_command = 'java -jar /usr/local/SKYVIEW/skyview.jar %s %s %s %s %s %s' %(position,coords,proj,outname,surveys,size) 
			print sky_command
			subprocess.call(sky_command,shell=True)
			
			##Sky view names are stoopid, you can't specify every name - it just whacks on an integer
			##Let's change that for easierness later on
			if len(image_cats)==1:
				change_names.append('mv %s.fits %s_%s.fits' %(src_all.names[0],src_all.names[0],image_cats[0]))
			else:
				for i in xrange(len(image_cats)):
					print 'mv %s_%d.fits %s_%s.fits' %(src_all.names[0],i+1,src_all.names[0],image_cats[i])
					change_names.append('mv %s_%d.fits %s_%s.fits' %(src_all.names[0],i+1,src_all.names[0],image_cats[i]))
		else:
			pass

for name in change_names: subprocess.call(name,shell=True)
