"""
data_plot.py

SuperClass module for the YProfile and Table Classes
It contains numerous plot function for the YProfile and Table Classes

If one in the future wants their class to inherit this superclasses methods
this is what is required:
A. Place 'from data_table import *' at the top of the module
B. If the class is defined like 'class MyClass:', change that to
   'class MyClass(DataTable):'
C. To properly use DataTable's methods properly one will need these methods:
	a get(atri) that returns a numpy array of Data, or a
	list of numpy arrays of data.  The arguments of this function would need to be
	atri which is the name of the data one is looking for.
"""
from numpy import *
from math import *
import matplotlib.pylab as pyl
import matplotlib.pyplot as pl
from matplotlib import colors,cm
from matplotlib.patches import Rectangle, Arrow
from matplotlib.collections import PatchCollection
from matplotlib.offsetbox import AnchoredOffsetbox, TextArea
from matplotlib.lines import Line2D
from matplotlib.ticker import *
import os
import threading
import time
import sys

def padding_model_number(number,max_num):
    ''' this method returns a zero-front padded string

    It makes out of str(45) -> '0045' if 999 < max_num < 10000. This
    is meant to work for reasonable integers (maybe less than 10^6).
    number   number that the string should represent
    max_num  max number of cycle list, implies how many 0s have
	     be padded
    '''

    cnum = str(number)
    clen = len(cnum)

    cmax = int(log10(max_num)) + 1

    return (cmax - clen)*'0' + cnum


class DataPlot():

	def classTest(self):
		'''
		Determines what the type of class instance the subclass is, so
		we can dynamically determine the behaviour of methods.

		This method NEEDS to be modified if any names of files or classes
		are changed
		'''

		c=str(self.__class__)
		Tmp=''
		if 'ppm.yprofile' == c:
			tmp='YProfile'
		elif 'ascii_table.ascii_table' == c:
			tmp='AsciiTable'
		elif 'nugridse.se' == c:
			tmp='se'
		elif 'mesa.mesa_profile' == c:
			tmp='mesa_profile'
		elif 'mesa.star_log' == c or 'mesa.history_data' == c:
			tmp='mesa.star_log'
		elif 'ppn.xtime' == c:
			tmp='xtime'
		elif 'ppn.abu_vector' == c:
			tmp='PPN'

		return tmp

	def which(self, program):
	    '''
	    Mimics which in the unix shell
	    '''
	    def is_exe(fpath):
		return os.path.exists(fpath) and os.access(fpath, os.X_OK)

	    fpath, fname = os.path.split(program)
	    if fpath:
		if is_exe(program):
		    return program
	    else:
		for path in os.environ["PATH"].split(os.pathsep):
		    exe_file = os.path.join(path, program)
		    if is_exe(exe_file):
			return exe_file

	    return None

	def logarithm(self,tmpX,tmpY,logX,logY,base):
		logXER=False
		logYER=False
		for i in range(len(tmpX)):
				if tmpX[i]<=0. and logX:
					#print 'We can not log a number less than or equal to zero'
					#print 'Attempting to remove incompatible values from X'
					logXER=True
				if tmpY[i]<=0. and logY:
					#print 'We can not log a number less than or equal to zero'
					#print 'Attempting to remove incompatible values from Y'
					logYER=True
		tmX=[]
		tmY=[]

		if logXER:
			for i in range(len(tmpX)):
				if tmpX[i]>0.:
					tmX.append( tmpX[i])
					tmY.append(tmpY[i])
			tmpX=tmX
			tmpY=tmY
		elif logYER:
			for i in range(len(tmpY)):
				if tmpY[i]>0.:
					tmX.append( tmpX[i])
					tmY.append(tmpY[i])
			tmpX=tmX
			tmpY=tmY

		tmX=tmpX
		tmY=tmpY

		if logX:
			tmX=tmpX
			try:
				for i in range(len(tmpX)):
					tmX[i]=log(tmpX[i],base)
			except ValueError:
				#print 'We can not log a number less than or equal to zero'
				#print 'Attempting to remove incompatible values from X'
				logXER=True
		if logY:
			tmY=tmpY
			try:
				for i in range(len(tmpY)):
					tmY[i]=log(tmpY[i],base)
			except ValueError:
				#print 'We can not log a number less than or equal to zero'
				#print 'Attempting to remove incompatible values from Y'
				logYER=True

		if logX:
			tmpX=tmX
		if logY:
			tmpY=tmY

		return tmpX,tmpY

	def sparse(self,x,y,sparse):
			"""
			Method that removes every non sparse th element.  For example
			if this argument was 5, This method would plot the 0th, 5th, 10th
			... elements.
			Input:
			x: list of x values, of lenthe j
			y: list of y values, of lenthe j
			sparse: Argument that skips every so many data points
			"""
			tmpX=[]
			tmpY=[]

			for i in range(len(x)):
				if sparse == 1:
					return x,y
				if (i%sparse)==0:
					tmpX.append(x[i])
					tmpY.append(y[i])
			return tmpX, tmpY

	def plotMulti(self,atrix,atriy, cyclist,title,path='/',legend=None,labelx=None, labely=None,logx=False, logy=False, \
			      base=10,sparse=1,pdf=False,limits=None,):
		'''
		Method for plotting multiple plots and saving it to multiple pngs
		or PDFs
		Input:
		atrix: The name of the attribute you want on the x axis
		atriy: The name of the attribute you want on the Y axis
		cyclist: List of cycles that you would like plotted
		title: The title of the graph and the name of the file.
		Legend: A list of legends for each of your cycles, or one legend
			for all of the cycles
		pdf: A boolean of if the image should be saved to a pdf file.
			xMin,xMax, yMin, YMax:  plot coordinates.
		logx: A boolean of whether the user wants the x axis logarithmically
		logy: A boolean of whether the user wants the Y axis logarithmically
		base: The base of the logarithm. Default = 10
		sparse: Argument that skips every so many data points. For
			example if this argument was 5, This method would plot
			the 0th, 5th, 10th ... elements.
		limits: The length four list of the x and y limits. The order of
			the list is xmin,xmax,ymin,ymax
		'''
		if str(legend.__class__)!="<type 'list'>":# Determines the legend is a list
			legendList=False
		else:
			legendList=True

		if legendList and len(cyclist) !=len(legend): #if it is a list, make sure there is an entry for each cycle
			print 'Please input a proper legend, with correct length, aborting plot'
			return None
		for i in xrange(len(cyclist)):
			if legendList:
				self.plot(atrix,atriy,cyclist[i],'ndump',legend[i],labelx,labely,base=base,sparse=sparse, \
						  logx=logx,logy=logy,show=False,limits=limits)
			else:
				self.plot(atrix,atriy,cyclist[i],'ndump',legend,labelx,labely,base=base,sparse=sparse, \
						  logx=logx,logy=logy,show=False,limits=limits)

			pl.title(title)
			if not pdf:
				currentDir = os.getcwd()
				os.chdir(path)
				pl.savefig(title+str(cyclist[i])+'.png', dpi=400)
				os.chdir(currentDir)
			else:
				currentDir = os.getcwd()
				os.chdir(path)
				pl.savefig(title+cyclist[i]+'.pdf', dpi=400)
				os.chdir(currentDir)
			pl.clf()
		return None

	def plot(self,atrix,atriy, fname=None,numtype='ndump',legend=None,labelx=None, labely=None ,\
		indexx=None, indexy=None, title=None, shape='.',logx=False, logy=False, base=10,sparse=1, \
                show=True,limits=None,markevery=None,linewidth=1):
		"""
		Simple function that plots atriy as a function of atrix
		This method will automatically find and plot the requested data.
		Input:
		atrix, The name of the attribute you want on the x axis
		atriy, The name of the attribute you want on the Y axis
		Fname: Be the filename, Ndump or time, or cycle,  If fname is a
		       list, this method will then save a png for each cycle in
		       the list. Warning, this must be a list of cycles and not
		       a list of filenames
		numtype: designates how this function acts and how it interprets
			 fname. Defaults to file
		if numtype is 'file', this function will get the desird
		attribute from that file
		if numtype is 'NDump' function will look at the cycle with that
		nDump
		if numtype is 't' or 'time' function will find the _cycle with
		the closest time stamp
		labelx: The label on the X axis
		labely: The label on the Y axis
		indexx: Depreciated: If the get method returns a list of lists,
			indexx would be the list at the index indexx in the list.
		indexy: Depreciated: If the get method returns a list of lists,
			indexy would be the list at the index indexx in the list.
		shape: What shape and colour the user would like their plot in.
		       Please see
		       http://matplotlib.sourceforge.net/api/pyplot_api.html#matplotlib.pyplot.plot
		       for all possible choices
		title: The Title of the Graph
		logx: A boolean of weather the user wants the x axis logarithmically
		logy: A boolean of weather the user wants the Y axis logarithmically
		base: The base of the logarithm. Default = 10
		sparse: Argument that skips every so many data points. For
			example if this argument was 5, This method would plot
			the 0th, 5th, 10th ... elements.
		show: A boolean of if the plot should be displayed useful with
			the multiPlot method
		WARNING: Unstable if get returns a list with only one element (x=[0])
		limits: The length four list of the x and y limits. The order of
			the list is xmin,xmax,ymin,ymax

                markevery Set the markevery property to subsample the
                          plot when using markers. Eg if markevery=5,
                          every 5-th marker will be plotted. every can
                          be
                          None           very point will be plotted
                          an integer N   Every N-th marker will be plotted 
                                         starting with marker 0
                          every=(start, N) will start at point start and plot 
                linewidth set linewidth, default=1
		"""
		t1=time.time()
		#Setting the axis labels

		if labelx== None :
			labelx=atrix
		if labely== None :
			labely=atriy

		if title!=None:
			title=title
		else:
			title=labely+' vs '+labelx

		if str(fname.__class__)=="<type 'list'>":
			self.plotMulti(atrix,atriy,fname,title,legend ,labelx,labely,logx, logy, 10,1,limits=limits)
			return
		tmpX=[]
		tmpY=[]
		singleX=False
		singleY=False
		#Getting data
		plotType=self.classTest()
		if plotType=='YProfile':
			if fname==None:
				fname=self.cycles[-1]

			listY=self.get(atriy,fname, numtype,resolution='a')
			listX=self.get(atrix,fname, numtype,resolution='a')
		elif plotType=='se':
			if fname==None:
				listY=self.get( atriy,sparse=sparse)
				listX=self.get(atrix,sparse=sparse)
			else:
				listY=self.get(fname, atriy,sparse=sparse)
				listX=self.get(fname, atrix,sparse=sparse)

			t2= time.time()
			print t2 -t1
		elif plotType=='PPN' :
			if fname==None and atrix not in self.cattrs and atriy not in self.cattrs:
				fname=len(self.files)-1
			if numtype=='ndump':
				numtype='cycNum'
			listY=self.get(atriy,fname,numtype)
			listX=self.get(atrix,fname,numtype)
		elif plotType=='xtime' or plotType=='mesa_profile' or plotType=='AsciiTable' or plotType=='mesa.star_log':
			listY=self.get(atriy)
			listX=self.get(atrix)
		else:
			listY=self.get(atriy)
			listX=self.get(atrix)
		tmpX=[]
		tmpY=[]
		if isinstance(listX[0], basestring) or isinstance(listY[0], basestring):
			for i in range(len(listX)):
				if '*****' == listX[i] or '*****' == listY[i]:
					print 'There seems to be a string of * in the lists'
					print 'Cutting out elements in both the lists that have an index equal to or greater than the index of the location of the string of *'
					break
				tmpX.append(float(listX[i]))
				tmpY.append(float(listY[i]))

			listX=tmpX
			listY=tmpY




		#Determining if listX is a list or a list of lists
		try:
			j=listX[0][0]
		except:
			singleX = True

		if len(listX) == 1:  # If it is a list of lists with one element.
			tmpX=listX[0]
		elif singleX == True:# If it is a plain list of values.
			tmpX=listX
		elif indexx==None and len(listX)>1: # If it is a list of lists of values.
						    # take the largest
			tmpX=listX[0]
			for i in range(len(listX)):
				if len(tmpX)<len(listX[i]):
					tmpX=listX[i]
		elif indexx<len(listX): # If an index is specified, use that index
			tmpX=listX[indexx]
		else:
			print 'Sorry that indexx does not exist, returning None'
			return None

		#Determining if listY is a list or a list of lists
		try:
			j=listY[0][0]
		except:
			singleY = True

		if len(listY) == 1: # If it is a list of lists with one element.
			#print 'hello'
			tmpY=listY[0]
		elif singleY == True: # If it is a plain list of values.
			#print 'world'
			tmpY=listY
		elif indexy==None and len(listY)>1:# If it is a list of lists of values.
						    # take the largest
			#print 'fourth'
			tmpY=listY[0]
			for i in range(len(listY)):
				if len(tmpY)<len(listY[i]):
					tmpY=listY[i]
		elif indexy<len(listY): # If an index is specified, use that index
			#print 'sixth'
			tmpY=listY[indexy]
		else:
			print 'Sorry that indexy does not exist, returning None'
			return None
		'''
		elif indexy==None and len(listY)==1:
			#print 'fifth'
			tmpY=listY
		'''




		#Here, if we end up with different sized lists to plot, it
		#searches for a list that is of an equal length
		if len(tmpY)!=len(tmpX):
			found=False
			print "It seems like the lists are not of equal length"
			print "Now attempting to find a compatible list for ListX"
			for i in range(len(listY)):
				if not singleY and len(tmpX)==len(listY[i]):
					tmpY=listY[i]
					found=True

			if not found:
				print "Now attempting to find a compatible list for ListY"
				for i in range(len(listX)):

					if not singleX and len(tmpY)==len(listX[i]):
						tmpX=listX[i]
						found=True

			if found:
				print "Suitable list found"
			else:

				print "There is no suitalble list, returning None"
				return None
		if len(tmpY)!=len(tmpX) and single == True:
			print 'It seems that the selected lists are of different\nsize, now returning none'
			return None
		# Sparse stuff
		if plotType!='se':
			tmpX,tmpY=self.sparse(tmpX,tmpY, sparse)

		# Logarithm stuff
		if logy or logx:
			tmpX,tmpY=self.logarithm(tmpX,tmpY,logx,logy,base)

		# Here it ensures that if we are plotting ncycle no values of '*' will be plotted
		tmX=[]
		tmY=[]
		for i in range(len(tmpX)):
			tmX.append(str(tmpX[i]))
			tmY.append(str(tmpY[i]))

		tmpX=[]
		tmpY=[]
		for i in range(len(tmX)):
			if '*' in tmX[i] or '*' in tmY[i]:
				print 'There seems to be a string of * in the lists'
				print 'Cutting out elements in both the lists that have an index equal to or greater than the index of the location of the string of *'
				break
			tmpX.append(float(tmX[i]))
			tmpY.append(float(tmY[i]))
		listX=tmpX
		listY=tmpY

		#Setting the axis labels

		if logx:
			labelx='log '+labelx
		if logy:
			labely='log '+labely

		if legend!=None:
			legend=legend
		else:
			legend=labely+' vs '+labelx



		pl.plot(listX,listY,shape,label=legend,markevery=markevery,linewidth=linewidth)
		pl.legend()
		pl.title(title)
		pl.xlabel(labelx)
		pl.ylabel(labely)
		if show:
			pl.show()

		if limits != None and len(limits)==4:

			pl.xlim(limits[0],limits[1])
			pl.ylim(limits[2],limits[3])

	def plot_ratios(self,misosx,misosy,solsysx=None,solsysy=None,graindata=None,m_co=None,C_star_only=False,misosxname=None,misosyname=None,deltax=True,deltay=True,logx=False,logy=False,title=None,legend=True,extra_legend='',errbar=True,iniabufile='../../frames/mppnp/USEEPP/iniab2.0E-02GN93.ppn',plt_show=True,modlegend=None,plt_symb='o',plt_col='b',plt_modms=10.,plt_grms=7.,plt_modlw=3.,plt_sparse=0,plt_mrng=False,calling_routine='all',private_legend=None):
		'''
		Method for plotting ratio data from model output as well as grain data.
		Important: You have to give some input to the routine!
		RT, October 2012
		graindata:	presolar grain data -> private is for a private.txt database file, same structure as other files required!
		misosx:		model x data, set None for grain plot only
		misosy:		model y data
      solsysx: 	solar system ratio of x-axis - only necessary if deltax=True, the two solsysx,y variable are necessary to avoid importing utils into DataPlot class. If you import it, mesa.py does not work anymore
      solsysy: 	solar system ratio of y-axis - only necessary if deltay=True
		m_co:			model C/O ratio
		deltax:		Delta values on x-axis
		deltay:		Delta values on y-axis
		logx:			logarithmic x-axis
		logy:			logarithmic y-axis
		title:		Title of plot
		legend:		True or False. Use legend(loc=?) command to move legend around after plotting.
		iniabufile:	initial abundance file - necessary
		plt_show:	Show the plot or not, default is of course
		modlegend:	legend for model data
		plt_symb:	symbold for plotting model data
		plt_modms:	markersize for models
      plt_grms:   markersize for grains
		plt_modlw:	linewidht for models
		plt_col:		color for plotting model data
		plt_sparse:	sparse function for model data
		plt_mrng:   Plot mass range for massive stars (used by plot4iso_exp routine in nugridse)
		calling_routine:	to identify where it comes from for some special treatment
		private_legend:	If you want your own legend for the private file, then you're right here
		'''
		# compatibility
		if misosxname != None and len(misosxname) == 2:
			dumb = []
			dumb = misosxname[0].split('-')
			dumb.append(misosxname[1].split('-')[0])
			dumb.append(misosxname[1].split('-')[1])
			misosxname = dumb
		if misosyname != None and len(misosyname) == 2:
			dumb = []
			dumb = misosyname[0].split('-')
			dumb.append(misosyname[1].split('-')[0])
			dumb.append(misosyname[1].split('-')[1])
			misosyname = dumb
		# style
		# Size of font etc.
		params = {'axes.labelsize':  20,
          'text.fontsize':   14,
          'legend.fontsize': 14,
          'xtick.labelsize': 14,
          'ytick.labelsize': 14}
		pl.rcParams.update(params)

		# prepare model data (if necessary) - here PPN / Nugridse difference if included at some point!
		if m_co != None:
			# find co_ratio, where it gets > 1
			index = -1
			for i in range(len(m_co)):
				if m_co[i] >= 1.:
					index = i
					break
			if index == -1:   # star doesn't get C rich in TPs
				mxdata_orich = misosx
				mydata_orich = misosy
				mxdata_crich = []
				mydata_crich = []
			else:
				mxdata_orich = misosx[0:index+1]
				mydata_orich = misosy[0:index+1]
				mxdata_crich = misosx[index:len(misosx)]
				mydata_crich = misosy[index:len(misosy)]
		# plotting ifs, depending on available data
		if graindata==None:
			if m_co==None:
				#markersize, and style stuff
				# plot
				if modlegend == None:
					pl.plot(misosx,misosy,'o--')
				else:
					if calling_routine == 'general':
						pl.plot(misosx,misosy,'o--',label=modlegend)
					elif calling_routine == '4iso_exp':
						plt_symb = plt_symb + '-'
						for it in range(len(misosx)):
							if it == 0:
								pl.plot(misosx[it],misosy[it],plt_symb,color=plt_col,markevery=plt_sparse,markersize=plt_modms,lw=plt_modlw,label=modlegend)
							else:
								pl.plot(misosx[it],misosy[it],plt_symb,color=plt_col,markevery=plt_sparse,markersize=plt_modms,lw=plt_modlw)
						# annotate mass range
						if plt_mrng != False:
							for mrng_ind in range(len(plt_mrng)):
								pl.text(plt_mrng[mrng_ind][0],plt_mrng[mrng_ind][1],str(round(plt_mrng[mrng_ind][2],2)),ha='right',va='bottom',color=plt_col)

				# axis
				if logx and logy:
					pl.loglog()
				elif logx:
					pl.semilogx()
				elif logy:
					pl.semilogy()
				# labels
				if misosxname != None:
					if deltax:
						pl.xlabel('$\delta$($^{' + str(misosxname[1]) + '}$' + misosxname[0] + '/$^{' + str(misosxname[3]) + '}$' + misosxname[2] + ')')
					else:
						pl.xlabel('$^{' + str(misosxname[1]) + '}$' + misosxname[0] + '/$^{' + str(misosxname[3]) + '}$' + misosxname[2])
				if misosyname != None:
					if deltay:
						pl.ylabel('$\delta$($^{' + str(misosyname[1]) + '}$' + misosyname[0] + '/$^{' + str(misosyname[3]) + '}$' + misosyname[2] + ')')
					else:
						pl.ylabel('$^{' + str(misosyname[1]) + '}$' + misosyname[0] + '/$^{' + str(misosyname[3]) + '}$' + misosyname[2])
				if title != None:
					pl.title(title)
				if legend != None or modlegend != None:
					pl.legend(loc=5)
			else:
				# plot
				if C_star_only ==False:
					pl.plot(mxdata_orich,mydata_orich,linestyle='--',marker=plt_symb,c=plt_col,label='C/O$<$1'+extra_legend)
				pl.plot(mxdata_crich,mydata_crich,'o-',marker=plt_symb,c=plt_col,label='C/O$>$1'+extra_legend,lw=plt_modlw,markersize=plt_modms,markeredgecolor='k',markerfacecolor=plt_col)
				# axis
				if logx and logy:
					pl.loglog()
				elif logx:
					pl.semilogx()
				elif logy:
					pl.semilogy()
				# labels
				if misosxname != None:
					if deltax:
						pl.xlabel('$\delta$($^{' + str(misosxname[1]) + '}$' + misosxname[0] + '/$^{' + str(misosxname[3]) + '}$' + misosxname[2] + ')')
					else:
						pl.xlabel('$^{' + str(misosxname[1]) + '}$' + misosxname[0] + '/$^{' + str(misosxname[3]) + '}$' + misosxname[2])
				if misosyname != None:
					if deltay:
						pl.ylabel('$\delta$($^{' + str(misosyname[1]) + '}$' + misosyname[0] + '/$^{' + str(misosyname[3]) + '}$' + misosyname[2] + ')')
					else:
						pl.ylabel('$^{' + str(misosyname[1]) + '}$' + misosyname[0] + '/$^{' + str(misosyname[3]) + '}$' + misosyname[2])
				if title != None:
					pl.title(title)
				if legend:
					pl.legend(loc=5)
		else:
			# transform data
			gtypelist = graindata[0]
			gdatax    = graindata[1]
			gdataxerr = graindata[2]
			gdatay    = graindata[3]
			gdatayerr = graindata[4]

			### PLOTS ###

			# grains
			for i in range(len(gtypelist)):
				# determine plot symbols and color
				if gtypelist[i][0].lower() == 'sic':
					msymb = 'o'
					if gtypelist[i][1] == 'M':
						mcol = '0.4'
					elif gtypelist[i][1] == 'X':
						mcol = 'b'
					elif gtypelist[i][1] == 'Y':
						mcol = 'g'
					elif gtypelist[i][1] == 'Z':
						mcol = 'r'
					elif gtypelist[i][1] == 'AB':
						mcol = 'c'
					elif gtypelist[i][1] == 'C' or gtypelist[i][1] == 'U/C':
						mcol = 'y'
					elif gtypelist[i][1] == 'N':
						mcol = 'm'
					else:
						mcol = '0.7'
				elif gtypelist[i][0].lower() == 'oxides':
					msymb = '^'
					if gtypelist[i][1] == '1':
						mcol = '0.4'
					elif gtypelist[i][1] == '2':
						mcol = 'b'
					elif gtypelist[i][1] == '3':
						mcol = 'g'
					elif gtypelist[i][1] == '4':
						mcol = 'r'
					else:
						mcol = '0.7'
				elif gtypelist[i][0].lower() == 'silicates':
					msymb = 'v'
					if gtypelist[i][1] == '1':
						mcol = '0.4'
					elif gtypelist[i][1] == '2':
						mcol = 'b'
					elif gtypelist[i][1] == '3':
						mcol = 'g'
					elif gtypelist[i][1] == '4':
						mcol = 'r'
					else:
						mcol = '0.7'
				elif gtypelist[i][0].lower() == 'graphites':
					msymb = 's'
					if gtypelist[i][1] == 'HD':
						mcol = '0.4'
					elif gtypelist[i][1] == 'LD':
						mcol = 'b'
					else:
						mcol = '0.7'
				elif gtypelist[i][0].lower() == 'idp':
					msymb = '>'
					mcol = '0.4'
				elif gtypelist[i][0].lower() == 'private':
					msymb = 'o'
					if gtypelist[i][1] == 'M':
						mcol = '0.4'
					elif gtypelist[i][1] == 'X':
						mcol = 'b'
					elif gtypelist[i][1] == 'Y':
						mcol = 'g'
					elif gtypelist[i][1] == 'Z':
						mcol = 'r'
					elif gtypelist[i][1] == 'AB':
						mcol = 'c'
					elif gtypelist[i][1] == 'C' or gtypelist[i][1] == 'U/C':
						mcol = 'y'
					elif gtypelist[i][1] == 'N':
						mcol = 'm'
					else:
						mcol = '0.7'
				else:
					msymb = '+'
					mcol = '0.4'
				# make nice labels in gtypelist now
				for jt in range(len(gtypelist[i])):
					if gtypelist[i][jt] == 'sic':
						gtypelist[i][jt] = 'SiC'
					elif gtypelist[i][jt] == 'oxides':
						gtypelist[i][jt] = 'Oxides'
					elif gtypelist[i][jt] == 'silicates':
						gtypelist[i][jt] = 'Silicated'
					elif gtypelist[i][jt] == 'graphites':
						gtypelist[i][jt] = 'Graphites'
					elif gtypelist[i][jt] == 'N':
						gtypelist[i][jt] = 'Nova'
					elif gtypelist[i][jt] == 'M':
						gtypelist[i][jt] = 'Mainstream'
					elif gtypelist[i][jt] == 'U':
						gtypelist[i][jt] = 'Unclassified'
					elif gtypelist[i][jt] == 'LD':
						gtypelist[i][jt] = 'low density'
					elif gtypelist[i][jt] == 'HD':
						gtypelist[i][jt] = 'high density'
					elif gtypelist[i][jt] == 'idp':
						gtypelist[i][jt] = 'IDP'
					elif gtypelist[i][jt] == 'private':
						if private_legend==None:
							gtypelist[i][jt] = 'Private'
						else:
							gtypelist[i][jt]=private_legend
				# now plot the grain data!
				# IDP grain labels properly
				if gtypelist[i][0] == 'IDP':
					grainlabelclean = 'IDP'
				else:
					grainlabelclean = gtypelist[i][0] + ' ' + gtypelist[i][1]
				if errbar:
					pl.errorbar(gdatax[i],gdatay[i],xerr=gdataxerr[i],yerr=gdatayerr[i],marker=msymb,c=mcol,linestyle='')
					pl.plot(gdatax[i],gdatay[i],msymb,c=mcol,label=grainlabelclean,markersize=plt_grms)
				else:
					pl.plot(gdatax[i],gdatay[i],msymb,c=mcol,label=grainlabelclean,markersize=plt_grms)
			# plot model:
			if m_co != None and misosx!=None:
				# actual plotting
				if C_star_only==False:
					pl.plot(mxdata_orich,mydata_orich,'--',marker=plt_symb,c=plt_col,label='C/O$<$1'+extra_legend,lw=plt_modlw)
				pl.plot(mxdata_crich,mydata_crich,'*-',marker=plt_symb,c=plt_col,label='C/O$>$1'+extra_legend,lw=plt_modlw,markersize=plt_modms,markeredgecolor='k',markerfacecolor=plt_col)
			elif misosx != None:
				if modlegend != None:
					if calling_routine == 'general':
						pl.plot(misosx,misosy,'o--',label=modlegend,markersize=plt_modms,lw=plt_modlw)
					elif calling_routine == '4iso_exp':
						plt_symb = plt_symb + '-'
						for it in range(len(misosx)):
							if it == 0:
								pl.plot(misosx[it],misosy[it],plt_symb,color=plt_col,markevery=plt_sparse,markersize=plt_modms,lw=plt_modlw,label=modlegend)
							else:
								pl.plot(misosx[it],misosy[it],plt_symb,color=plt_col,markevery=plt_sparse,markersize=plt_modms,lw=plt_modlw)
						# annotate mass range
						if plt_mrng != False:
							for mrng_ind in range(len(plt_mrng)):
								pl.text(plt_mrng[mrng_ind][0],plt_mrng[mrng_ind][1],str(round(plt_mrng[mrng_ind][2],2)),ha='right',va='bottom',color=plt_col)
				else:
					pl.plot(misosx,misosy,'o--',markersize=plt_modms,lw=plt_modlw)
			# axis
			if logx and logy:
				pl.loglog()
			elif logx:
				pl.semilogx()
			elif logy:
				pl.semilogy()
			# labels
			if misosxname != None:
				if deltax:
					pl.xlabel('$\delta$($^{' + str(misosxname[1]) + '}$' + misosxname[0] + '/$^{' + str(misosxname[3]) + '}$' + misosxname[2] + ')')
				else:
					pl.xlabel('$^{' + str(misosxname[1]) + '}$' + misosxname[0] + '/$^{' + str(misosxname[3]) + '}$' + misosxname[2])
			if misosyname != None:
				if deltay:
					pl.ylabel('$\delta$($^{' + str(misosyname[1]) + '}$' + misosyname[0] + '/$^{' + str(misosyname[3]) + '}$' + misosyname[2] + ')')
				else:
					pl.ylabel('$^{' + str(misosyname[1]) + '}$' + misosyname[0] + '/$^{' + str(misosyname[3]) + '}$' + misosyname[2])
			if title != None:
				pl.title(title)
			if legend:
				pl.legend(loc=5)

		# plot horizontal and vertical lines
		if deltax:
			pl.axhline(0,color='k')
		else:
			pl.axhline(solsysy,color='k')
		if deltay:
			pl.axvline(0,color='k')
		else:
			pl.axvline(solsysx,color='k')

		if plt_show:
			pl.show()



	def clear(self, title=True, xlabel=True, ylabel=True):
		'''
		Method for removing the title and/or xlabel and/or Ylabel
		input:
		Title -  boolean of if title will be cleared
		xlabel - boolean of if xlabel will be cleared
		ylabel - boolean of if ylabel will be cleared
		'''
		if title:
			pyl.title('')
		if xlabel:
			pyl.xlabel('')
		if ylabel:
			pyl.ylabel('')
	# From mesa.py
	def xlimrev(self):
		''' reverse xrange'''
		xmax,xmin=pyl.xlim()
		pyl.xlim(xmin,xmax)

	def abu_chartMulti(self,cyclist, mass_range=None ,ilabel = True,imlabel = True,\
			   imagic =  False,boxstable=True,lbound=20,plotaxis=[0,0,0,0],\
			   color_map='jet',pdf=False,title=None):
		'''
		Method that plots abundence chart and saves those
		figures to a .png file (by default). Plots a figure
		for each cycle in the argument cycle

		input:
		cyclist: The list of cycles we are plotting
		ilabel: elemental labels off/on [0/1]
		imlabel: label for isotopic masses off/on [0/1]
		imagic:  turn lines for magic numbers off/on [0/1]
		plotaxis: Set axis limit: If default [0,0,0,0] the complete
			  range in (N,Z) will be plotted
			  format: What format will this be saved in ['pdf'/'png']
		title: The title of the plots and the saved images
		'''

		if self.which('dvipng')==None:
			print "This method may need the third party program dvipng to operate"
			print 'It is located at http://sourceforge.net/projects/dvipng/'

		max_num = max(cyclist)
		for i in xrange(len(cyclist)):
			self.abu_chart( cyclist[i], mass_range ,ilabel,imlabel,imagic,\
					boxstable,lbound,plotaxis,False,color_map)
			if title !=None:
				pl.title(title)
			else:
				name='AbuChart'
			number_str=padding_model_number(cyclist[i],max_num)
			if not pdf:
				pl.savefig(name+number_str+'.png', dpi=100)
			else:
				pl.savefig(name+number_str+'.pdf', dpi=200)
			pl.close()

		return None
	#from mppnp.se
	def abu_chart(self, cycle, mass_range=None ,ilabel = True,imlabel = True,imagic =  False,
		boxstable=True,lbound=20,plotaxis=[0,0,0,0],show=True,color_map='jet',ifig=None):
		'''
		Plots an abundance chart
		input:
		cycle: The cycle we are looking in. It it is a list of cycles,
			this method will then do a plot for each of these cycles
			and save them all to a file
		ilabel: elemental labels off/on [False/True] defaults to True
		imlabel: label for isotopic masses off/on [False/True], defaults to True
		imagic:  turn lines for magic numbers off/on [False/True] defaults to False
		boxstable: plot the black boxes around the stable elements,
			   defaults to true
		lbound: The lower bound of the colour spectrum plotted. Defaults
			to 20
		plotaxis: Set axis limit: If default [0,0,0,0] the complete
			range in (N,Z) will be plotted. It equates to
			[xMin,xMax,Ymin,Ymax]
		mass_range - a 1x2 array containing the lower and upper mass range.
		    		 If this is an instance of abu_vector this will
		    		 only plot isotopes that have an atomic mass
		    		 within this range. This will throw an error if
		    		 this range does not make sence ie [45,2]
			 	if None, it will plot over the entire range
				Defaults to None
		show:  boolean of if the plot should be displayed useful with
		       saving multiple plots using abu_chartMulti
                color_map  color map according to choices in matplotlib
                           (e.g. www.scipy.org/Cookbook/Matplotlib/Show_colormaps)
                ifig  figure number, defaults to cycle number
		'''
                
                if ifig == None:
                    ifig=cycle

		if str(cycle.__class__)=="<type 'list'>":
			self.abu_chartMulti(cycle, mass_range,ilabel,imlabel,imagic,boxstable,\
					    lbound,plotaxis,color_map)
			return
		plotType=self.classTest()

		if mass_range!=None and mass_range[0]>mass_range[1]:
			print 'Please input a proper mass range'
			print 'Returning None'
			return None

		if plotType=='se':
			cycle=self.se.findCycle(cycle)
			nin=zeros(len(self.se.A))
			zin=zeros(len(self.se.Z))
			for i in xrange(len(nin)):
				nin[i]=self.se.A[i]
				zin[i]=self.se.Z[i]
			for i in xrange(len(nin)):
				nin[i]=nin[i]-zin[i]
			yin=self.get(cycle, 'iso_massf')
			isom=self.se.isomeric_states

			masses = self.se.get(cycle,'mass')
			if mass_range != None:
				masses = self.se.get(cycle,'mass')
				masses.sort()

			if mass_range != None:
				tmpyps=[]
				masses = self.se.get(cycle,'mass')
				masses = self.se.get(cycle,'mass')
				masses.sort()
				for i in xrange(len(masses)):
					if (masses[i] >mass_range[0] and masses[i]<mass_range[1]) or\
                                                (masses[i]==mass_range[0] or masses[i]==mass_range[1]):
						tmpyps.append(yin[i])
				yin=tmpyps


			tmp=zeros(len(yin[0]))
			for i in xrange(len(yin)):
				for j in xrange(len(yin[i])):
					tmp[j]+=yin[i][j]

			tmp=tmp/len(yin)

			yin=tmp

		elif plotType=='PPN':

			ain=self.get('A',cycle)
			zin=self.get('Z',cycle)
                        nin=ain-zin
			yin=self.get('ABUNDANCE_MF',cycle)
			isom=self.get('ISOM',cycle)

			if mass_range != None:
				tmpA=[]
				tmpZ=[]
				tmpIsom=[]
				tmpyps=[]
				for i in xrange(len(nin)):
					if (ain[i] >mass_range[0] and ain[i]<mass_range[1])\
					or (ain[i]==mass_range[0] or ain[i]==mass_range[1]):
						tmpA.append(nin[i])
						tmpZ.append(zin[i])
						tmpIsom.append(isom[i])
						tmpyps.append(yin[i])
				zin=tmpZ
				nin=tmpA
				yin=tmpyps
				isom=tmpIsom

		else:
			print 'This method, abu_chart, is not supported by this class'
			print 'Returning None'
			return None
		# in case we call from ipython -pylab, turn interactive on at end again
		turnoff=False
		if not show:
			try:
				ioff()
				turnoff=True
			except NameError:
				turnoff=False

		nnmax = int(max(nin))+1
		nzmax = int(max(zin))+1
		nzycheck = zeros([nnmax,nzmax,3])
		for i in range(len(nin)):
			if isom[i]==1:
				ni = int(nin[i])
				zi = int(zin[i])

				nzycheck[ni,zi,0] = 1
				nzycheck[ni,zi,1] = yin[i]




		#######################################################################
		# elemental names: elname(i) is the name of element with Z=i

		elname=self.elements_names

		#### create plot
		## define axis and plot style (colormap, size, fontsize etc.)
		if plotaxis==[0,0,0,0]:
		  xdim=10
		  ydim=6
		else:
		  dx = plotaxis[1]-plotaxis[0]
		  dy = plotaxis[3]-plotaxis[2]
		  ydim = 6
		  xdim = ydim*dx/dy


		params = {'axes.labelsize':  15,
			  'text.fontsize':   12,
			  'legend.fontsize': 15,
			  'xtick.labelsize': 15,
			  'ytick.labelsize': 15,
			  'text.usetex': True}
		#pl.rcParams.update(params) #May cause Error, someting to do with tex
		fig=pl.figure(ifig,figsize=(xdim,ydim),dpi=100)
		axx = 0.10
		axy = 0.10
		axw = 0.85
		axh = 0.8
		ax=pl.axes([axx,axy,axw,axh])
		# Tick marks
		xminorlocator = MultipleLocator(1)
		xmajorlocator = MultipleLocator(5)
		ax.xaxis.set_major_locator(xmajorlocator)
		ax.xaxis.set_minor_locator(xminorlocator)
		yminorlocator = MultipleLocator(1)
		ymajorlocator = MultipleLocator(5)
		ax.yaxis.set_major_locator(ymajorlocator)
		ax.yaxis.set_minor_locator(yminorlocator)

		# color map choice for abundances

                cmapa = cm.get_cmap(name=color_map)
		# color map choice for arrows
		cmapr = cm.autumn
		# if a value is below the lower limit its set to white
		cmapa.set_under(color='w')
		cmapr.set_under(color='w')
		# set value range for abundance colors (log10(Y))
		norma = colors.Normalize(vmin=-lbound,vmax=0)
		# set x- and y-axis scale aspect ratio to 1
		ax.set_aspect('equal')
		#print time,temp and density on top
		temp = ' '#'%8.3e' %ff['temp']
		time = ' '#'%8.3e' %ff['time']
		dens = ' '#'%8.3e' %ff['dens']

		#May cause Error, someting to do with tex
		'''
		#box1 = TextArea("t : " + time + " s~~/~~T$_{9}$ : " + temp + "~~/~~$\\rho_{b}$ : " \
		#	      + dens + ' g/cm$^{3}$', textprops=dict(color="k"))
		anchored_box = AnchoredOffsetbox(loc=3,
				child=box1, pad=0.,
				frameon=False,
				bbox_to_anchor=(0., 1.02),
				bbox_transform=ax.transAxes,
				borderpad=0.,
				)
		ax.add_artist(anchored_box)
		'''
		## Colour bar plotted

		patches = []
		color = []

		for i in range(nzmax):
		    for j in range(nnmax):
		      if nzycheck[j,i,0]==1:
			xy = j-0.5,i-0.5

			rect = Rectangle(xy,1,1,)

			# abundance
			yab = nzycheck[j,i,1]
			if yab == 0:

				yab=1e-99


			col =log10(yab)

			patches.append(rect)
			color.append(col)


		p = PatchCollection(patches, cmap=cmapa, norm=norma)
		p.set_array(array(color))
		p.set_zorder(1)
		ax.add_collection(p)
		cb = pl.colorbar(p)

		# colorbar label
		cb.set_label('log$_{10}$(X)',fontsize='x-large')

		# plot file name
		graphname = 'abundance-chart'+str(cycle)

		# Add black frames for stable isotopes
		if boxstable:
			for i in xrange(len(self.stable_el)):
				if i == 0:
					continue


				tmp = self.stable_el[i]
				try:
					zz= self.elements_names.index(tmp[0]) #charge
				except:
					continue

				for j in xrange(len(tmp)):
					if j == 0:
						continue

					nn = int(tmp[j]) #atomic mass
					nn=nn-zz

					xy = nn-0.5,zz-0.5
					rect = Rectangle(xy,1,1,ec='k',fc='None',fill='False',lw=3.)
					rect.set_zorder(2)
					ax.add_patch(rect)




		# decide which array to take for label positions
		iarr = 0

		# plot element labels
		if ilabel:
		  for z in range(nzmax):
		    try:
		      nmin = min(argwhere(nzycheck[:,z,iarr]))[0]-1
		      ax.text(nmin,z,elname[z],horizontalalignment='center',verticalalignment='center',\
			      fontsize='x-small',clip_on=True)
		    except ValueError:
		      continue

		# plot mass numbers
		if imlabel:
		  for z in range(nzmax):
		     for n in range(nnmax):
			a = z+n
			if nzycheck[n,z,iarr]==1:
			  ax.text(n,z,a,horizontalalignment='center',verticalalignment='center',\
				  fontsize='xx-small',clip_on=True)

		# plot lines at magic numbers
		if imagic:
		  ixymagic=[2, 8, 20, 28, 50, 82, 126]
		  nmagic = len(ixymagic)
		  for magic in ixymagic:
		    if magic<=nzmax:
		      try:
			xnmin = min(argwhere(nzycheck[:,magic,iarr]))[0]
			xnmax = max(argwhere(nzycheck[:,magic,iarr]))[0]
			line = ax.plot([xnmin,xnmax],[magic,magic],lw=3.,color='r',ls='-')
		      except ValueError:
			dummy=0
		    if magic<=nnmax:
		      try:
			yzmin = min(argwhere(nzycheck[magic,:,iarr]))[0]
			yzmax = max(argwhere(nzycheck[magic,:,iarr]))[0]
			line = ax.plot([magic,magic],[yzmin,yzmax],lw=3.,color='r',ls='-')
		      except ValueError:
			dummy=0

		# set axis limits
		if plotaxis==[0,0,0,0]:

		  xmax=max(nin)
		  ymax=max(zin)
		  ax.axis([-0.5,xmax+0.5,-0.5,ymax+0.5])
		else:
		  ax.axis(plotaxis)

		# set x- and y-axis label
		ax.set_xlabel('neutron number (A-Z)',fontsize=14)
		ax.set_ylabel('proton number Z',fontsize=14)
		pl.title('Isotopic Chart for cycle '+str(int(cycle)))
		#fig.savefig(graphname)
		print graphname,'is done'
		if show:
			pl.show()
		if turnoff:
			ion()
		return

	def abu_flux_chart(self, cycle,ilabel = True,imlabel = True,imagic =  False,
		boxstable=True,lbound=15,plotaxis=[0,0,0,0],which_flux=None,prange=None,profile='charged',show=True):
		'''
		Plots an abundance and flux chart
		input:
		cycle: The cycle we are looking in. It it is a list of cycles,
			this method will then do a plot for each of theses cycles
			and save them all to a file
		ilabel: elemental labels off/on [False/True] defaults to True
		imlabel: label for isotopic masses off/on [False/True], defaults to True
		imagic:  turn lines for magic numbers off/on [False/True] defaults to False
		boxstable: plot the black boxes around the stable elements,
			   defaults to true
		lbound: The lower bound of the colour spectrum ploted. Defaults
			to 20
		plotaxis: Set axis limit: If default [0,0,0,0] the complete
			range in (N,Z) will be plotted. It equates to
			[xMin,xMax,Ymin,Ymax]
		which_flux = 0 is for nucleosynthesis flux plot, which_flux = 1 is for energy flux plot, None is default, that is option 0.
		prange is the range of fluxes to be considered.
		format: What format will this be saved in ['pdf'/'png']
				Defaults to None
		profile: 'charged' is ideal setting to show charged particle reactions flow.
		    'neutron' is ideal setting for neutron captures flows.
		show:  boolean of if the plot should be displayed useful with
		       saving multiple plots using abu_chartMulti
		'''
		#######################################################################
		#### plot options
		# Set axis limit: If default [0,0,0,0] the complete range in (N,Z) will
		# be plotted, i.e. all isotopes, else specify the limits in
		# plotaxis = [xmin,xmax,ymin,ymax]

		#######################################################################

		# read data file
		#inpfile = cycle
		#ff = fdic.ff(inpfile)
		# with the flux implementation I am not using mass range for now.
		# It may be introduced eventually.
		mass_range = None
		if str(cycle.__class__)=="<type 'list'>":
			self.abu_chartMulti(cycle, mass_range,ilabel,imlabel,imagic,boxstable,\
					    lbound,plotaxis)
			return
		plotType=self.classTest()

		#if mass_range!=None and mass_range[0]>mass_range[1]:
			#print 'Please input a proper mass range'
			#print 'Returning None'
			#return None

		if plotType=='se':
			cycle=self.se.findCycle(cycle)
			nin=zeros(len(self.se.A))
			zin=zeros(len(self.se.Z))
			for i in xrange(len(nin)):
				nin[i]=self.se.A[i]
				zin[i]=self.se.Z[i]
			for i in xrange(len(nin)):
				nin[i]=nin[i]-zin[i]
			yin=self.get(cycle, 'iso_massf')
			isom=self.se.isomeric_states

			masses = self.se.get(cycle,'mass')
			if mass_range != None:
				masses = self.se.get(cycle,'mass')
				masses.sort()

			if mass_range != None:
				tmpyps=[]
				masses = self.se.get(cycle,'mass')
				masses = self.se.get(cycle,'mass')
				masses.sort()
				for i in xrange(len(masses)):
					if (masses[i] >mass_range[0] and masses[i]<mass_range[1]) or\
                                                (masses[i]==mass_range[0] or masses[i]==mass_range[1]):
						tmpyps.append(yin[i])
				yin=tmpyps


			tmp=zeros(len(yin[0]))
			for i in xrange(len(yin)):
				for j in xrange(len(yin[i])):
					tmp[j]+=yin[i][j]

			tmp=tmp/len(yin)

			yin=tmp

		elif plotType=='PPN':

			ain=self.get('A',cycle)
			zin=self.get('Z',cycle)
                        nin=ain-zin
			yin=self.get('ABUNDANCE_MF',cycle)
			isom=self.get('ISOM',cycle)

			if mass_range != None:
				tmpA=[]
				tmpZ=[]
				tmpIsom=[]
				tmpyps=[]
				for i in xrange(len(nin)):
					if (ain[i] >mass_range[0] and ain[i]<mass_range[1])\
					or (ain[i]==mass_range[0] or ain[i]==mass_range[1]):
						tmpA.append(nin[i])
						tmpZ.append(zin[i])
						tmpIsom.append(isom[i])
						tmpyps.append(yin[i])
				zin=tmpZ
				nin=tmpA
				yin=tmpyps
				isom=tmpIsom

		else:
			print 'This method, abu_chart, is not supported by this class'
			print 'Returning None'
			return None
		# in case we call from ipython -pylab, turn interactive on at end again
		turnoff=False
		if not show:
			try:
				ioff()
				turnoff=True
			except NameError:
				turnoff=False

		nnmax = int(max(nin))+1
		nzmax = int(max(zin))+1
		nnmax_plot = nnmax
		nzmax_plot = nzmax
		nzycheck = zeros([nnmax,nzmax,3])
		nzycheck_plot = zeros([nnmax,nzmax,3])
		for i in range(len(nin)):
			if isom[i]==1:
				ni = int(nin[i])
				zi = int(zin[i])

				nzycheck[ni,zi,0] = 1
				nzycheck[ni,zi,1] = yin[i]
				nzycheck_plot[ni,zi,0] = 1



		#######################################################################
		# elemental names: elname(i) is the name of element with Z=i

		elname=self.elements_names

		#### create plot
		## define axis and plot style (colormap, size, fontsize etc.)
		if plotaxis==[0,0,0,0]:
		  xdim=10
		  ydim=6
		else:
		  dx = plotaxis[1]-plotaxis[0]
		  dy = plotaxis[3]-plotaxis[2]
		  ydim = 6
		  xdim = ydim*dx/dy


		params = {'axes.labelsize':  15,
			  'text.fontsize':   12,
			  'legend.fontsize': 15,
			  'xtick.labelsize': 15,
			  'ytick.labelsize': 15,
			  'text.usetex': True}
		#pl.rcParams.update(params) #May cause Error, someting to do with tex
		#fig=pl.figure(figsize=(xdim,ydim),dpi=100)
		fig=pl.figure()
		if profile == 'charged':
		    ax1 = fig.add_subplot(1, 2, 1)
		elif profile == 'neutron':
		    ax1 = fig.add_subplot(2, 1, 1)
		#axx = 0.10
		#axy = 0.10
		#axw = 0.85
		#axh = 0.8
		#ax1=pl.axes([axx,axy,axw,axh])
		# Tick marks
		xminorlocator = MultipleLocator(1)
		xmajorlocator = MultipleLocator(5)
		ax1.xaxis.set_major_locator(xmajorlocator)
		ax1.xaxis.set_minor_locator(xminorlocator)
		yminorlocator = MultipleLocator(1)
		ymajorlocator = MultipleLocator(5)
		ax1.yaxis.set_major_locator(ymajorlocator)
		ax1.yaxis.set_minor_locator(yminorlocator)

		# color map choice for abundances
		#cmapa = cm.jet
		cmapa = cm.summer
		# color map choice for arrows
		cmapr = cm.summer
		# if a value is below the lower limit its set to white
		cmapa.set_under(color='w')
		cmapr.set_under(color='w')
		# set value range for abundance colors (log10(Y))
		norma = colors.Normalize(vmin=-lbound,vmax=0)
		# set x- and y-axis scale aspect ratio to 1
		#ax1.set_aspect('equal')
		#print time,temp and density on top
		temp = ' '#'%8.3e' %ff['temp']
		time = ' '#'%8.3e' %ff['time']
		dens = ' '#'%8.3e' %ff['dens']

		#May cause Error, someting to do with tex
		'''
		#box1 = TextArea("t : " + time + " s~~/~~T$_{9}$ : " + temp + "~~/~~$\\rho_{b}$ : " \
		#	      + dens + ' g/cm$^{3}$', textprops=dict(color="k"))
		anchored_box = AnchoredOffsetbox(loc=3,
				child=box1, pad=0.,
				frameon=False,
				bbox_to_anchor=(0., 1.02),
				bbox_transform=ax.transAxes,
				borderpad=0.,
				)
		ax.add_artist(anchored_box)
		'''
		## Colour bar plotted

		patches = []
		color = []

		for i in range(nzmax):
		    for j in range(nnmax):
		      if nzycheck[j,i,0]==1:
			xy = j-0.5,i-0.5

			rect = Rectangle(xy,1,1,)

			# abundance
			yab = nzycheck[j,i,1]
			if yab == 0:

				yab=1e-99


			col =log10(yab)

			patches.append(rect)
			color.append(col)


		p = PatchCollection(patches, cmap=cmapa, norm=norma)
		p.set_array(array(color))
		p.set_zorder(1)
		ax1.add_collection(p)
		cb = pl.colorbar(p)

		# colorbar label
		if profile == 'neutron':
		    cb.set_label('log$_{10}$(X)',fontsize='x-large')

		# plot file name
		graphname = 'abundance-flux-chart'+str(cycle)

		# Add black frames for stable isotopes
		if boxstable:
			for i in xrange(len(self.stable_el)):
				if i == 0:
					continue


				tmp = self.stable_el[i]
				try:
					zz= self.elements_names.index(tmp[0]) #charge
				except:
					continue

				for j in xrange(len(tmp)):
					if j == 0:
						continue

					nn = int(tmp[j]) #atomic mass
					nn=nn-zz

					xy = nn-0.5,zz-0.5
					rect = Rectangle(xy,1,1,ec='k',fc='None',fill='False',lw=4.)
					rect.set_zorder(2)
					ax1.add_patch(rect)




		# decide which array to take for label positions
		iarr = 0

		# plot element labels
		if ilabel:
		  for z in range(nzmax):
		    try:
		      nmin = min(argwhere(nzycheck[:,z,iarr]))[0]-1
		      nmax = max(argwhere(nzycheck[:,z,iarr]))[0]+1
		      ax1.text(nmin,z,elname[z],horizontalalignment='center',verticalalignment='center',\
			      fontsize='medium',clip_on=True)
		      ax1.text(nmax,z,elname[z],horizontalalignment='center',verticalalignment='center',\
			      fontsize='medium',clip_on=True)
		    except ValueError:
		      continue

		# plot mass numbers
		if imlabel:
		  for z in range(nzmax):
		     for n in range(nnmax):
			a = z+n
			if nzycheck[n,z,iarr]==1:
			  ax1.text(n,z,a,horizontalalignment='center',verticalalignment='center',\
				  fontsize='xx-small',clip_on=True)


		# plot lines at magic numbers
		if imagic:
		  ixymagic=[2, 8, 20, 28, 50, 82, 126]
		  nmagic = len(ixymagic)
		  for magic in ixymagic:
		    if magic<=nzmax:
		      try:
			xnmin = min(argwhere(nzycheck[:,magic,iarr]))[0]
			xnmax = max(argwhere(nzycheck[:,magic,iarr]))[0]
			line = ax1.plot([xnmin,xnmax],[magic,magic],lw=3.,color='r',ls='-')
		      except ValueError:
			dummy=0
		    if magic<=nnmax:
		      try:
			yzmin = min(argwhere(nzycheck[magic,:,iarr]))[0]
			yzmax = max(argwhere(nzycheck[magic,:,iarr]))[0]
			line = ax1.plot([magic,magic],[yzmin,yzmax],lw=3.,color='r',ls='-')
		      except ValueError:
			dummy=0

		# set axis limits
		if plotaxis==[0,0,0,0]:

		  xmax=max(nin)
		  ymax=max(zin)
		  ax1.axis([-0.5,xmax+0.5,-0.5,ymax+0.5])
		else:
		  ax1.axis(plotaxis)

		# set x- and y-axis label
		ax1.set_ylabel('Proton number',fontsize='xx-large')
		if profile == 'charged':
		    ax1.set_xlabel('Neutron number',fontsize='xx-large')
		#pl.title('Isotopic Chart for cycle '+str(int(cycle)))

		#
		# here below I read data from the flux_*****.DAT file.
		#
		file_name = 'flux_'+str(cycle).zfill(5)+'.DAT'
		print file_name
		f = open(file_name)
		lines = f.readline()
		lines = f.readlines()
		f.close()

		print_max_flux_in_plot =  False
		# color map choice for fluxes
		#cmapa = cm.jet
		cmapa = cm.autumn
		# color map choice for arrows
		cmapr = cm.autumn
		# starting point of arrow
		coord_x_1 = []
		coord_y_1 = []
		# ending point of arrow (option 1)
		coord_x_2 = []
		coord_y_2 = []
		# ending point of arrow (option 2)
		coord_x_3 = []
		coord_y_3 = []
		# fluxes
		flux_read = []
		flux_log10 = []

		if which_flux == None or which_flux == 0:
			print 'chart for nucleosynthesis fluxes [dYi/dt]'
			line_to_read = 9
		elif which_flux == 1:
			print 'chart for energy fluxes'
			line_to_read = 10
		elif which_flux > 1:
			print "you have only option 0 or 1, not larger than 1"

		single_line = []
		for i in range(len(lines)):
			single_line.append(lines[i].split())
			coord_y_1.append(float(single_line[i][1]))
			coord_x_1.append(float(single_line[i][2])-coord_y_1[i])
			coord_y_2.append(float(single_line[i][5]))
			coord_x_2.append(float(single_line[i][6])-coord_y_2[i])
			coord_y_3.append(float(single_line[i][7]))
			coord_x_3.append(float(single_line[i][8])-coord_y_3[i])
			try:
				flux_read.append(float(single_line[i][line_to_read]))
			except ValueError: # this is done to avoid format issues like 3.13725-181...
				flux_read.append(1.0E-99)
			flux_log10.append(log10(flux_read[i]+1.0e-99))

		print file_name,' read!'

		# I need to select smaller sample, with only fluxes inside plotaxis.
		if plotaxis!=[0,0,0,0]:
			coord_y_1_small=[]
			coord_x_1_small=[]
			coord_y_2_small=[]
			coord_x_2_small=[]
			coord_y_3_small=[]
			coord_x_3_small=[]
			flux_log10_small = []
			for i in range(len(flux_log10)):
				I_am_in = 0
				if coord_y_1[i] > plotaxis[2] and coord_y_1[i] < plotaxis[3] and coord_x_1[i] > plotaxis[0] and coord_x_1[i] < plotaxis[1]:
					I_am_in = 1
					coord_y_1_small.append(coord_y_1[i])
					coord_x_1_small.append(coord_x_1[i])
					coord_y_2_small.append(coord_y_2[i])
					coord_x_2_small.append(coord_x_2[i])
					coord_y_3_small.append(coord_y_3[i])
					coord_x_3_small.append(coord_x_3[i])
					flux_log10_small.append(flux_log10[i])
				if coord_y_3[i] > plotaxis[2] and coord_y_3[i] < plotaxis[3] and coord_x_3[i] > plotaxis[0] and coord_x_3[i] < plotaxis[1] and I_am_in == 0:
					I_am_in = 1
					coord_y_1_small.append(coord_y_1[i])
					coord_x_1_small.append(coord_x_1[i])
					coord_y_2_small.append(coord_y_2[i])
					coord_x_2_small.append(coord_x_2[i])
					coord_y_3_small.append(coord_y_3[i])
					coord_x_3_small.append(coord_x_3[i])
					flux_log10_small.append(flux_log10[i])



		# elemental labels off/on [0/1]
		ilabel = 1

		# label for isotopic masses off/on [0/1]
		imlabel = 1

		# turn lines for magic numbers off/on [0/1]
		imagic = 0

		# flow is plotted over "prange" dex. If flow < maxflow-prange it is not plotted
		if prange == None:
			print 'plot range given by default'
			prange = 8.

		#############################################
		#print flux_log10_small
		# we should scale prange on plot_axis range, not on max_flux!
		max_flux = max(flux_log10)
		ind_max_flux = flux_log10.index(max_flux)
		if plotaxis!=[0,0,0,0]:
			max_flux_small = max(flux_log10_small)

		if plotaxis==[0,0,0,0]:
			nzmax = int(max(max(coord_y_1),max(coord_y_2),max(coord_y_3)))+1
			nnmax = int(max(max(coord_x_1),max(coord_x_2),max(coord_x_3)))+1
			coord_x_1_small = coord_x_1
			coord_x_2_small = coord_x_2
			coord_x_3_small = coord_x_3
			coord_y_1_small = coord_y_1
			coord_y_2_small = coord_y_2
			coord_y_3_small = coord_y_3
			flux_log10_small= flux_log10
			max_flux_small  = max_flux
		else:
			nzmax = int(max(max(coord_y_1_small),max(coord_y_2_small),max(coord_y_3_small)))+1
			nnmax = int(max(max(coord_x_1_small),max(coord_x_2_small),max(coord_x_3_small)))+1

		for i in range(nzmax):
			for j in range(nnmax):
				if nzycheck[j,i,0]==1:
					xy = j-0.5,i-0.5
					rect = Rectangle(xy,1,1,)
					patches.append(rect)


		nzycheck = zeros([nnmax_plot,nzmax,3])
		coord_x_out = zeros(len(coord_x_2_small))
		coord_y_out = zeros(len(coord_y_2_small))
		for i in range(len(flux_log10_small)):
			nzycheck[coord_x_1_small[i],coord_y_1_small[i],0] = 1
			nzycheck[coord_x_1_small[i],coord_y_1_small[i],1] = flux_log10_small[i]
			if coord_x_2_small[i] >= coord_x_3_small[i]:
				coord_x_out[i] = coord_x_2_small[i]
				coord_y_out[i] = coord_y_2_small[i]
				nzycheck[coord_x_out[i],coord_y_out[i],0] = 1
				nzycheck[coord_x_out[i],coord_y_out[i],1] = flux_log10_small[i]
			elif coord_x_2_small[i] < coord_x_3_small[i]:
				coord_x_out[i] = coord_x_3_small[i]
				coord_y_out[i] = coord_y_3_small[i]
				nzycheck[coord_x_out[i],coord_y_out[i],0] = 1
				nzycheck[coord_x_out[i],coord_y_out[i],1] = flux_log10_small[i]
			if flux_log10_small[i]>max_flux_small-prange:
				nzycheck[coord_x_1_small[i],coord_y_1_small[i],2] = 1
				nzycheck[coord_x_out[i],coord_y_out[i],2] = 1

		#### create plot
		if profile == 'charged':
		    ax2 = fig.add_subplot(1, 2, 2)
		elif profile == 'neutron':
		    ax2 = fig.add_subplot(2, 1, 2)
		# Tick marks
		xminorlocator = MultipleLocator(1)
		xmajorlocator = MultipleLocator(5)
		ax2.xaxis.set_major_locator(xmajorlocator)
		ax2.xaxis.set_minor_locator(xminorlocator)
		yminorlocator = MultipleLocator(1)
		ymajorlocator = MultipleLocator(5)
		ax2.yaxis.set_major_locator(ymajorlocator)
		ax2.yaxis.set_minor_locator(yminorlocator)
		## define axis and plot style (colormap, size, fontsize etc.)
		if plotaxis==[0,0,0,0]:
			xdim=10
			ydim=6
		else:
			dx = plotaxis[1]-plotaxis[0]
			dy = plotaxis[3]-plotaxis[2]
			ydim = 6
			xdim = ydim*dx/dy

		format = 'pdf'
		# set x- and y-axis scale aspect ratio to 1
		#ax2.set_aspect('equal')

		# Add black frames for stable isotopes
		# Add black frames for stable isotopes
		if boxstable:
			for i in xrange(len(self.stable_el)):
				if i == 0:
					continue


				tmp = self.stable_el[i]
				try:
					zz= self.elements_names.index(tmp[0]) #charge
				except:
					continue

				for j in xrange(len(tmp)):
					if j == 0:
						continue

					nn = int(tmp[j]) #atomic mass
					nn=nn-zz

					xy = nn-0.5,zz-0.5
					rect = Rectangle(xy,1,1,ec='k',fc='None',fill='False',lw=4.)
					rect.set_zorder(2)
					ax2.add_patch(rect)


		apatches = []
		acolor = []
		m = 0.8/prange
		vmax=ceil(max(flux_log10_small))
		vmin=max(flux_log10_small)-prange
		b=-vmin*m+0.1
		normr = colors.Normalize(vmin=vmin,vmax=vmax)
		ymax=0.
		xmax=0.

		for i in range(len(flux_log10_small)):
			x = coord_x_1_small[i]
			y = coord_y_1_small[i]
			dx = coord_x_out[i]-coord_x_1_small[i]
			dy = coord_y_out[i]-coord_y_1_small[i]
			if flux_log10_small[i]>=vmin:
				arrowwidth = flux_log10_small[i]*m+b
				arrow = Arrow(x,y,dx,dy, width=arrowwidth)
				if xmax<x:
					xmax=x
				if ymax<y:
					ymax=y
				acol = flux_log10_small[i]
				apatches.append(arrow)
				acolor.append(acol)
			xy = x-0.5,y-0.5
			rect = Rectangle(xy,1,1,ec='k',fc='None',fill='False',lw=1.)
			patches.append(rect)
			xy = x+dx-0.5,y+dy-0.5
			rect = Rectangle(xy,1,1,ec='k',fc='None',fill='False',lw=1.)
			patches.append(rect)


		p = PatchCollection(patches,norm=0,facecolor='w')
		p.set_zorder(1)
		ax2.add_collection(p)




		a = PatchCollection(apatches, cmap=cmapr, norm=normr)
		a.set_array(array(acolor))
		a.set_zorder(3)
		ax2.add_collection(a)
		cb = pl.colorbar(a)

		# colorbar label
		cb.set_label('log$_{10}$($x$)',fontsize='x-large')
		if profile == 'neutron':
		    cb.set_label('log$_{10}$(f)',fontsize='x-large')

		# decide which array to take for label positions
		iarr = 2

		# plot element labels
		for z in range(nzmax):
			try:
				nmin = min(argwhere(nzycheck_plot[:,z,iarr-2]))[0]-1
				nmax = max(argwhere(nzycheck_plot[:,z,iarr-2]))[0]+1
				ax2.text(nmin,z,elname[z],horizontalalignment='center',verticalalignment='center',fontsize='medium',clip_on=True)
				ax2.text(nmax,z,elname[z],horizontalalignment='center',verticalalignment='center',fontsize='medium',clip_on=True)
			except ValueError:
				continue

		# plot mass numbers
		if imlabel:
			for z in range(nzmax):
				for n in range(nnmax_plot):
					a = z+n
					if nzycheck_plot[n,z,iarr-2]==1:
						ax2.text(n,z,a,horizontalalignment='center',verticalalignment='center',fontsize='xx-small',clip_on=True)

		# plot lines at magic numbers
		if imagic==1:
			ixymagic=[2, 8, 20, 28, 50, 82, 126]
			nmagic = len(ixymagic)
			for magic in ixymagic:
				if magic<=nzmax:
					try:
						xnmin = min(argwhere(nzycheck[:,magic,iarr-2]))[0]
						xnmax = max(argwhere(nzycheck[:,magic,iarr-2]))[0]
						line = ax2.plot([xnmin,xnmax],[magic,magic],lw=3.,color='r',ls='-')
					except ValueError:
						dummy=0
				if magic<=nnmax:
					try:
						yzmin = min(argwhere(nzycheck[magic,:,iarr-2]))[0]
						yzmax = max(argwhere(nzycheck[magic,:,iarr-2]))[0]
						line = ax2.plot([magic,magic],[yzmin,yzmax],lw=3.,color='r',ls='-')
					except ValueError:
						dummy=0

		# set axis limits
		if plotaxis==[0,0,0,0]:
			ax2.axis([-0.5,xmax+0.5,-0.5,ymax+0.5])
		else:
			ax2.axis(plotaxis)

		# set x- and y-axis label
		ax2.set_xlabel('Neutron number',fontsize='xx-large')
		if profile == 'neutron':
		    ax2.set_ylabel('Proton number',fontsize='xx-large')
		if which_flux == None or which_flux == 0:
			max_flux_label="max flux = "+str('{0:.4f}'.format(max_flux))
		elif which_flux == 1:
			max_flux_label="max energy flux = "+str('{0:.4f}'.format(max_flux))
		if print_max_flux_in_plot:
			ax2.text(plotaxis[1]-1.8,plotaxis[2]+0.1,max_flux_label,fontsize=10.)

		#fig.savefig(graphname)
		print graphname,'is done'
		if show:
			pl.show()
		if turnoff:
			ion()
		return

	def iso_abundMulti(self,cyclist, stable=False,amass_range=None,mass_range=None,
		ylim=[0,0],ref=-1,decayed=False,include_title=False,title=None,pdf=False,\
                color_plot=True,grid=False,point_set=1):
		'''
		Method that plots figures and saves those figures to a .png file
		(by default). Plots a figure for each cycle in the argument cycle. 
                Can be called via iso_abund method by passing a list to cycle. 
                Documentation there.
                 '''
		max_num = max(cyclist)
		for i in xrange(len(cyclist)):
			self.iso_abund(cyclist[i],stable,amass_range,mass_range,ylim,ref,\
			     decayed=decayed,show=False,color_plot=color_plot,grid=False,\
                             point_set=1,include_title=include_title)
			if title !=None:
				pl.title(title)
			else:
				name='IsoAbund'
			number_str=padding_model_number(cyclist[i],max_num)
			if not pdf:
				pl.savefig(name+number_str+'.png', dpi=200)
			else:
				pl.savefig(name+number_str+'.pdf', dpi=200)
			pl.clf()

		return None

	def iso_abund(self,cycle, stable=False,amass_range=None,mass_range=None,ylim=[0,0],
			ref=-1,show=True,log_logic=True,decayed=False,color_plot=True,grid=False,\
            point_set=1,include_title=False):
		''' plot the abundance of all the chemical species

		    cycle      - a string/integer of the cycle of interest.
                                 If it is a list of cycles, this method will do a plot 
                                 for each cycle and save them to a file
		    stable     - a boolean of whether to filter out the unstables.
		    	         Defaults to False
		    amass_range -a 1x2 array containing the lower and upper Atomic
		    		 mass range. optional. if None plot entire available
                                 atomic mass range
		    mass_range - a 1x2 array containing the lower and upper mass range.
		    		 If this is an instance of abu_vector this will
		    		 only plot isotopes that have an atominc mass
		    		 within this range. This will throw an error if
		    		 this range does not make sense ie [45,2]
			 	 if None, it will plot over the entire range
				 Defaults to None
                    ref  - reference cycle. If it is not -1 (default), this method will
                           plot the abundences of cycle devided by the cycle of the same
                           instance given in the ref variable. If ref is a list it will 
                           be interpreted to have two elements: ref=['dir/of/ref/run',cycle]
                           which uses a refernece cycle from another run. If any
		    	   abundence in the reference cycle is zero, it will replace it with
		    	   1e-99. The default is -1, it will do nothing.
		    ylim - A 1x2 array containing the lower and upper Y limits.
		    	   Defaults to [0,0], in which case ylim will be determined automatically
			log_logic  = True[/False] to plot abundances in log scale or linear.
		    decayed    = [True/]False to plot decayed distributions (True), or life 
                                 distribution
		    color_plot - color dots and lines, True or False
                    grid       - add grid, True or default False
                    point_set  - 0 (default), 1 or 2 to select one of three point sets, useful 
                                 for multiple abundances or ratios in one plot
		'''

		plotType=self.classTest()
		if str(cycle.__class__)=="<type 'list'>":
			self.iso_abundMulti(cycle, stable,amass_range,mass_range,ylim,ref,\
                             decayed,include_title,color_plot=color_plot,grid=False,point_set=point_set)
			return

		if mass_range!=None and mass_range[0]>mass_range[1]:
			print 'Please input a proper mass range'
			print 'Returning None'
			return None
		if amass_range!=None and amass_range[0]>amass_range[1]:
			print 'Please input a proper Atomic mass range'
			print 'Returning None'
			return None
		if plotType=='se':
			if decayed:
				print 'Decay option not yet implemented for mppnp - but it is easy do! Consider investing the time!'
				return None


                        # get things as arrays
			cycle=self.se.findCycle(cycle)
                        a_iso_to_plot   = array(self.se.A)
                        abunds          = self.get(cycle,'iso_massf')
                        isotope_to_plot = array(self.se.isotopes)
                        z_iso_to_plot   = array(self.se.Z)
                        isomers_to_plot = array(self.se.isomeric_states)
			if ref >-1:
				ref=self.se.findCycle(ref)
				abundsRef=self.se.get(ref,'iso_massf')


			masses = self.se.get(cycle,'mass')
			if mass_range == None:
			    print 'Using default mass range'
			    mass_range = [min(masses),max(masses)]
			masses.sort()
			mass_range.sort()

                        if amass_range == None:
                            amass_range=[min(a_iso_to_plot),max(a_iso_to_plot)]

                        # remove neutrons - this could move in the non- se/PPN specific part below
                        if 0 in z_iso_to_plot:
                            ind_neut        = where(z_iso_to_plot==0)[0][0]
                            a_iso_to_plot   = delete(a_iso_to_plot,ind_neut)
                            z_iso_to_plot   = delete(z_iso_to_plot,ind_neut)
                            isomers_to_plot = delete(isomers_to_plot,ind_neut)
                            isotope_to_plot = delete(isotope_to_plot,ind_neut)
                            abunds = delete(abunds,ind_neut,1)
                            if ref >-1:
                                abundsRef = delete(abundsRef,ind_neut,1)

                        # extract amass_range
                        acon=(a_iso_to_plot>=amass_range[0]) & (a_iso_to_plot<=amass_range[1])
                        isomers_to_plot = isomers_to_plot[acon]
                        isotope_to_plot = isotope_to_plot[acon]
                        z_iso_to_plot   = z_iso_to_plot[acon]
                        abunds          = abunds.T[acon].T
                        if ref >-1:
                            abundsRef = abundsRef.T[acon].T
                        a_iso_to_plot   = a_iso_to_plot[acon]
                        el_iso_to_plot  = array([x.split('-')[0] for x in isotope_to_plot.tolist()])
                        # apply mass range
                        if mass_range == None:
                            print 'Using default mass range'
                            mass_range = [min(masses),max(masses)]
                        mass_range.sort()
                        aabs = []
                        if ref >-1:
                            cyc  = [cycle,ref]
                            abus = [abunds,abundsRef]
                        else:
                            cyc  = [cycle]
                            abus = [abunds]
                        for cc,aa in zip(cyc,abus):
                            masses = self.se.get(cc,'mass')
                            masses.sort()
                            dmass  = masses[1:] - masses[:-1]    # I should check the grid definition
                            dmass  = append(dmass,0.)
                            mcon   = (masses>=mass_range[0]) & (masses<=mass_range[1])
                            dmass  = dmass[mcon]
                            aa = aa[mcon]

                            # average over mass range:
                            aa = (aa.T*dmass).T.sum(0)
                            aa = aa / (mass_range[1] - mass_range[0])
                            # abunds has now length of isotope_to_plot
                            aabs.append(aa)
                        if ref >-1:
                            abunds = aabs[0]/(aabs[1]+1.e-99)
                        else:
                            abunds = aabs[0]

			self.a_iso_to_plot=a_iso_to_plot
			self.isotope_to_plot=isotope_to_plot
                        self.z_iso_to_plot=z_iso_to_plot
			self.el_iso_to_plot=el_iso_to_plot
			self.abunds=abunds
			self.isomers_to_plot=isomers_to_plot

#                        self.isotopes = self.se.isotopes

		elif plotType=='PPN':
			print "This method adds the following variables to the instance:"
			print "a_iso_to_plot      mass number of plotted range of species"
			print "isotope_to_plot    corresponding list of isotopes"
			print "z_iso_to_plot      corresponding charge numbers"
			print "el_iso_to_plot     corresponding element names"
			print "abunds             corresponding abundances"
			print "isom               isomers and their abundance"

                        self.get(cycle,decayed=decayed)
			if ref is not -1:
                            if type(ref) is list: # reference cycle from other run
                                import ppn
                                pp=ppn.abu_vector(ref[0])
                                abunds_pp=pp.get(ref[1],decayed=decayed)
                                self.abunds=self.abunds/pp.abunds
                            else:
                                abunds=self.abunds
                                self.get(ref,decayed=decayed)
                                self.abunds=abunds/(self.abunds+1.e-99)
                        if amass_range == None:
                            amass_range=[min(self.a_iso_to_plot),max(self.a_iso_to_plot)]

                        aa=ma.masked_outside(self.a_iso_to_plot,amass_range[0],amass_range[1])
                        isotope_to_plot=ma.array(self.isotope_to_plot,mask=aa.mask).compressed()
                        z_iso_to_plot=ma.array(self.z_iso_to_plot,mask=aa.mask).compressed()
                        el_iso_to_plot=ma.array(self.el_iso_to_plot,mask=aa.mask).compressed()
                        abunds=ma.array(self.abunds,mask=aa.mask).compressed()
                        a_iso_to_plot=aa.compressed()
                        isomers_to_plot=[]
                        for i in xrange(len(self.isom)):
                                if int(self.isom[i][0].split('-')[1])>100:
                                        isomers_to_plot.append(self.isom[i])

			self.a_iso_to_plot=a_iso_to_plot
			self.isotope_to_plot=isotope_to_plot
                        self.z_iso_to_plot=z_iso_to_plot
			self.el_iso_to_plot=el_iso_to_plot
			self.abunds=abunds
			self.isomers_to_plot=isomers_to_plot
		else:
			print 'This method, iso_abund, is not supported by this class'
			print 'Returning None'
			return None

		print 'Using the following conditions:'
		if plotType=='se':
                    print '\tmass_range:', mass_range[0], mass_range[1]
                print '\tAtomic mass_range:', amass_range[0], amass_range[1]
		print '\tcycle:           ',cycle
		print '\tplot only stable:',stable
		print '\tplot decayed:    ',decayed


                if stable: # remove unstables:
                    # For the element that belongs to the isotope at index 5 in isotope_to_plot
                    # (C-12) the following gives the mass numbers of stable elements:
                    # self.stable_el[self.stable_names.index(el_iso_to_plot[5])][1:]
                    ind_delete=[]
                    for i in range(len(isotope_to_plot)):
                        if a_iso_to_plot[i] not in self.stable_el[self.stable_names.index(el_iso_to_plot[i])][1:]:
                            ind_delete.append(i)
                    a_iso_to_plot   = delete(a_iso_to_plot,  ind_delete)
                    z_iso_to_plot   = delete(z_iso_to_plot,  ind_delete)
                    isomers_to_plot = delete(isomers_to_plot,ind_delete)
                    isotope_to_plot = delete(isotope_to_plot,ind_delete)
                    el_iso_to_plot  = delete(el_iso_to_plot, ind_delete)
                    abunds          = delete(abunds,         ind_delete)

                el_list=[] # list of elements in el_iso_to_plot
                for el in self.elements_names:
                    if el in el_iso_to_plot:
                        el_list.append(el)

                abund_plot = []  # extract for each element an abundance and associated
                mass_num   = []  # mass number array, sorted by mass number
                for el in el_list:
                    numbers = a_iso_to_plot[(el_iso_to_plot==el)]
                    abund_plot.append(abunds[(el_iso_to_plot==el)][argsort(numbers)])
                    mass_num.append(sort(numbers))
                # now plot:
		plot_type = ['-','--','-.',':','-']
		pl_index = 0
		points = [['o','^','p','h','*'],['x','+','D','>','s'],['H','v','<','*','3']]
		if color_plot:
			colors = ['g','r','c','m','k']
		else:
			colors = ['k','k','k','k','k']

                ylim1 =  1.e99        
                ylim2 = -1.e99        
		for j in xrange(len(abund_plot)):        #Loop through the elements of interest
		    #    Process the line
		    #print 'processing line'
		    for l in xrange(len(abund_plot[j])):
                        # print mass_num[j][l]
                        # print abund_plot[j][l]
			if abund_plot[j][l] == 0:
			    abund_plot[j][l] = 1e-99

                    a_dum=zeros(len(abund_plot[j]))   # this I (FH) have to do because for some
                    if log_logic == False:            # reason log10(abu_abund[j]) does not work
                        a_dum = abund_plot[j]         # although abu_abund[j] is a numpy array?!?
                    else: 
                        for ii in range(len(abund_plot[j])):
                            a_dum[ii]=log10(abund_plot[j][ii])
                    this_label=str(colors[pl_index]+points[point_set][pl_index]+\
                               plot_type[pl_index])
                    pl.plot(mass_num[j],a_dum,this_label,markersize=10)
                    abu_max = max(a_dum)
                    max_index=where(a_dum==abu_max)[0][0]
                    coordinates=[mass_num[j][max_index],abu_max]
                    pl.text(coordinates[0]+0.1,1.05*coordinates[1],el_list[j]) 

                    pl_index+=1
                    if pl_index > 4:
                        pl_index = 0
                    ylim1=min(ylim1,min(a_dum))
                    ylim2=max(ylim2,max(a_dum))
                # now trimming the ylims    
                if log_logic:
                    dylim=0.05*(ylim2-ylim1)
                    ylim1 = ylim1 -dylim
                    ylim2 = ylim2 +dylim
                    if ref is not -1:
                        ylim2 = min(ylim2,4)
                        ylim1 = max(ylim1,-4)
                    else:
                        ylim2 = min(ylim2,0.2)
                        ylim1 = max(ylim1,-13)
                else:
                    ylim1 = ylim1 *0.8
                    ylim2 = ylim2 *1.1
                if include_title:
                    if plotType=='se':
                        if ref == -1:
                            title = str('Range %4.2f' %mass_range[0]) + str('-%4.2f' %mass_range[1]) +\
                                str(' for cycle %d' %int(cycle))
                        else:
                            title = str('Range %4.2f' %mass_range[0]) + \
                                str('-%4.2f' %mass_range[1]) + str(' for cycle %d' %int(cycle))+\
                                str(' relative to cycle %d'  %int(ref))
                    else:
                        if ref == -1:
                            title = str('Cycle %d' %int(cycle))
                        else:
                            title = str('Cycle %d' %int(cycle))+\
                                    str(' relative to cycle %d'  %int(ref))
                    print "including title: ..."
                    pl.title(title)                           
                if ylim[0] == 0 and ylim[1] == 0:
                    pl.ylim(ylim1,ylim2)
                else:
                    pl.ylim(ylim[0],ylim[1])
                pl.xlim([amass_range[0]-.5,amass_range[1]+.5])
                pl.xlabel('mass number (A)',fontsize=14)
                if ref is not -1:
                        if log_logic:
                                pl.ylabel(r'log abundance ratio',fontsize=14)
                        else:
                                pl.ylabel(r'abundance ratio',fontsize=14)
                else:
                        if log_logic:
                                pl.ylabel(r'log mass fraction ',fontsize=14)
                        else:
                                pl.ylabel(r'mass fraction',fontsize=14)
		if grid:
                    pl.grid()
		if show:
                    pl.show()

                if amass_range != None:
                    minimum_mass = amass_range[0]
                    maximum_mass = amass_range[1]

                elif mass_range != None:
                    minimum_mass = mass_range[0]
                    maximum_mass = mass_range[1]

                else:
                    minimum_mass = 0
                    maximum_mass = 200

                if log_logic == False:  
                    pl.plot([amass_range[0]-.5,amass_range[1]+.5],[1,1],'k-')
                else: 
                    pl.plot([amass_range[0]-.5,amass_range[1]+.5],[0,0],'k-')
                    
		ax=pl.axes()
                labelsx=[]
                if (maximum_mass-minimum_mass) > 100:
                    delta_labelsx = 10
                else:
                    delta_labelsx = 5

                iii = amass_range[0]%delta_labelsx
                if iii == 0:
                    labelsx.append(str(amass_range[0]))
                else:
                    labelsx.append(' ')
                iii = iii+1
                kkk = 0
                for label1 in range(amass_range[1]-amass_range[0]):
                    if iii == 5:
                        kkk = kkk+1
                        labelsx.append(str((iii*kkk)+amass_range[0]-(amass_range[0]%5)))
                        iii = 0
                        iii = iii+1
                    else:
                        labelsx.append(' ')
                        iii = iii+1

                if delta_labelsx == 5:
                    xticks = arange(amass_range[0],amass_range[1],1)
                    pl.xticks(xticks,labelsx)
                else:
                    pl.xticks()

##!!FOR!!###### print 'LEN LABELS= ', len(labelsx)
##DEBUGGING####
####!!!######## for bbb in range (len(labelsx)):
###############     print labelsx[bbb]
		return

	def plotprofMulti(self,ini,end,delta,what_specie,xlim1,xlim2,ylim1,ylim2,symbol=None):

		''' create a movie with mass fractions vs mass coordinate
		between xlim1 and xlim2, ylim1 and ylim2. Only works with instances of se

		ini          - initial model
		end          - final model
		delta        - sparsity factor of the frames
		what_specie  - array with species in the plot
		xlim1, xlim2 - mass coordinate range
		ylim1, ylim2 - mass fraction coordinate range
	        symbol       - array indicating which symbol you want to use. Must be of the same len of what_specie array
		'''
		plotType=self.classTest()
		if plotType=='se':
			for i in range(ini,end+1,delta):
			    step = int(i)
			    #print step
			    if symbol==None:
				symbol_dummy = '-'
			    	for j in range(len(what_specie)):
					self.plot_prof_1(step,what_specie[j],xlim1,xlim2,ylim1,ylim2,symbol_dummy)
			    else:
			    	for j in range(len(what_specie)):
					symbol_dummy = symbol[j]
					self.plot_prof_1(step,what_specie[j],xlim1,xlim2,ylim1,ylim2,symbol_dummy)

			    #
			    filename = str('%03d' % step)+'_test.png'
			    pl.savefig(filename, dpi=400)
			    print 'wrote file ', filename
			    #
			    pl.clf()

		else:
			print 'This method is not supported for '+str(self.__class__)
			return
	# From mesa_profile
    	def plot_prof_1(self,species,keystring,xlim1,xlim2,ylim1,ylim2,symbol=None, show=False):
		''' plot one species for cycle between xlim1 and xlim2
		    Only works with instances of se and mesa _profile

		species      - which species to plot
		keystring    - label that appears in the plot or in the cas of se,
			       A cycle or list of cycles
		xlim1, xlim2 - mass coordinate range
		ylim1, ylim2 - mass fraction coordinate range
		symbol       - indicate which symbol you want to use, if required. '''
		plotType=self.classTest()
		if plotType=='se':
			#tot_mass=self.se.get(keystring,'total_mass')
			tot_mass=self.se.get('mini')
			age=self.se.get(keystring,'age')
			mass=self.se.get(keystring,'mass')
			Xspecies=self.se.get(keystring,'iso_massf',species)

			mod=keystring
		elif plotType=='mesa_profile':
			tot_mass=self.header_attr['star_mass']
			age=self.header_attr['star_age']
			mass=self.get('mass')
			mod=self.header_attr['model_number']
			Xspecies=self.get(species)
		else:
			print 'This method is not supported for '+str(self.__class__)
			return

		if symbol == None:
			symbol = '-'

		x,y=self.logarithm(Xspecies,mass,True,False,10)
		#print x
		pl.plot(y,x,symbol,label=str(species))
		pl.xlim(xlim1,xlim2)
		pl.ylim(ylim1,ylim2)
		pl.legend()

		pl.xlabel('$Mass$ $coordinate$', fontsize=20)
		pl.ylabel('$X_{i}$', fontsize=20)
		pl.title('Mass='+str(tot_mass)+', Time='+str(age)+' years, cycle='+str(mod))
		if show:
			pl.show()

	# From mesa.star_log

def flux_chart(file_name,plotaxis,plot_type,which_flux=None,I_am_the_target=None,prange=None):
	'''
	Plots a chart with fluxes
	input:
	file_name: name of the file of fluxes we are looking at.
        plotaxis: [xmin,xmax,ymin,ymax], where on x axis there is neutron number and on y axis there is Z.
        plot_types: 0 for standard flux plot, 1 if fluxes focused on one specie.
        Note: the script is terribly slow, need to be improved. For now I put here in data_plot:
        [1]: import data_plot
        [2]: data_plot.flux_chart('file_name',[xmin,xmax,ymin,ymax],int,which_flux,I_am_the_target,prange)
	The pdf is created, but an error bumped up and the gui is empty. To avoid this, I had to set 'text.usetex': False. See below.
        Also, for the same reason no label in x axys is written using 'text.usetex': True.
        Note also that the GUI works really slow with this plot. so, we need to optimize from the graphic point of view.
        This need to be included in ppn.py I think, and set in multi option too, in case we want to read more flux files at the same time.
        Finally, you need to have stable.dat to read in to make it work....
	which_flux = 0 is for nucleosynthesis flux plot, which_flux = 1 is for energy flux plot, None is default, that is option 0.
	I_am_the_target is a 2Xarray used only if plot_type=1, and is given by [neutron number,proton number].
	prange is the range of fluxes to be considered.
        '''

	import numpy as np
	import matplotlib.pyplot as plt
	from matplotlib.mpl import colors,cm
	from matplotlib.patches import Rectangle, Arrow
	from matplotlib.collections import PatchCollection
	from matplotlib.offsetbox import AnchoredOffsetbox, TextArea
	import sys

	print_max_flux_in_plot =  True


    	f = open(file_name)
    	lines = f.readline()
        lines = f.readlines()
    	f.close()

        # starting point of arrow
        coord_x_1 = []
        coord_y_1 = []
        # ending point of arrow (option 1)
        coord_x_2 = []
        coord_y_2 = []
        # ending point of arrow (option 2)
        coord_x_3 = []
        coord_y_3 = []
        # fluxes
        flux_read = []
        flux_log10 = []

	if which_flux == None or which_flux == 0:
		print 'chart for nucleosynthesis fluxes [dYi/dt]'
		line_to_read = 9
	elif which_flux == 1:
		print 'chart for energy fluxes'
		line_to_read = 10
	elif which_flux == 2:
		print 'chart for timescales'
		line_to_read = 11
	elif which_flux > 2:
		print "you have only option 0, 1 or 2, not larger than 2"

        single_line = []
        for i in range(len(lines)):
                single_line.append(lines[i].split())
 		coord_y_1.append(float(single_line[i][1]))
                coord_x_1.append(float(single_line[i][2])-coord_y_1[i])
 		coord_y_2.append(float(single_line[i][5]))
                coord_x_2.append(float(single_line[i][6])-coord_y_2[i])
 		coord_y_3.append(float(single_line[i][7]))
                coord_x_3.append(float(single_line[i][8])-coord_y_3[i])
		try:
                	flux_read.append(float(single_line[i][line_to_read]))
		except ValueError: # this is done to avoid format issues like 3.13725-181...
			flux_read.append(1.0E-99)
                flux_log10.append(np.log10(flux_read[i]+1.0e-99))

    	print 'file read!'

	# I need to select smaller sample, with only fluxes inside plotaxis.
	coord_y_1_small=[]
	coord_x_1_small=[]
	coord_y_2_small=[]
	coord_x_2_small=[]
	coord_y_3_small=[]
	coord_x_3_small=[]
	flux_log10_small = []
	for i in range(len(flux_log10)):
		I_am_in = 0
		if coord_y_1[i] > plotaxis[2] and coord_y_1[i] < plotaxis[3] and coord_x_1[i] > plotaxis[0] and coord_x_1[i] < plotaxis[1]:
			I_am_in = 1
			coord_y_1_small.append(coord_y_1[i])
			coord_x_1_small.append(coord_x_1[i])
			coord_y_2_small.append(coord_y_2[i])
			coord_x_2_small.append(coord_x_2[i])
			coord_y_3_small.append(coord_y_3[i])
			coord_x_3_small.append(coord_x_3[i])
			flux_log10_small.append(flux_log10[i])
		if coord_y_3[i] > plotaxis[2] and coord_y_3[i] < plotaxis[3] and coord_x_3[i] > plotaxis[0] and coord_x_3[i] < plotaxis[1] and I_am_in == 0:
			I_am_in = 1
			coord_y_1_small.append(coord_y_1[i])
			coord_x_1_small.append(coord_x_1[i])
			coord_y_2_small.append(coord_y_2[i])
			coord_x_2_small.append(coord_x_2[i])
			coord_y_3_small.append(coord_y_3[i])
			coord_x_3_small.append(coord_x_3[i])
			flux_log10_small.append(flux_log10[i])


        if plot_type == 1:
		print 'I_am_the_target=',I_am_the_target
		#I_am_the_target = [56.-26.,26.]
        # here below need for plotting
	# plotaxis = [xmin,xmax,ymin,ymax]
	#plotaxis=[1,20,1,20]
	#plotaxis=[0,0,0,0]

	# elemental labels off/on [0/1]
	ilabel = 1

	# label for isotopic masses off/on [0/1]
	imlabel = 1

	# turn lines for magic numbers off/on [0/1]
	imagic = 0

	# flow is plotted over "prange" dex. If flow < maxflow-prange it is not plotted
	if prange == None:
		print 'plot range given by default'
		prange = 8.

        #############################################

	# we should scale prange on plot_axis range, not on max_flux!
        max_flux = max(flux_log10)
	ind_max_flux = flux_log10.index(max_flux)
	max_flux_small = max(flux_log10_small)
        min_flux = min(flux_log10)
	ind_min_flux = flux_log10.index(min_flux)
	min_flux_small = min(flux_log10_small)

  	#nzmax = int(max(max(coord_y_1),max(coord_y_2),max(coord_y_3)))+1
  	#nnmax = int(max(max(coord_x_1),max(coord_x_2),max(coord_x_3)))+1
	nzmax = int(max(max(coord_y_1_small),max(coord_y_2_small),max(coord_y_3_small)))+1
  	nnmax = int(max(max(coord_x_1_small),max(coord_x_2_small),max(coord_x_3_small)))+1


        nzycheck = np.zeros([nnmax,nzmax,3])
        #coord_x_out = np.zeros(len(coord_x_2))
        #coord_y_out = np.zeros(len(coord_y_2))
        #for i in range(len(flux_log10)):
  	#	nzycheck[coord_x_1[i],coord_y_1[i],0] = 1
  	#	nzycheck[coord_x_1[i],coord_y_1[i],1] = flux_log10[i]
  	#	if coord_x_2[i] >= coord_x_3[i]:
        #                coord_x_out[i] = coord_x_2[i]
        #                coord_y_out[i] = coord_y_2[i]
    	#		nzycheck[coord_x_out[i],coord_y_out[i],0] = 1
    	#		nzycheck[coord_x_out[i],coord_y_out[i],1] = flux_log10[i]
        #        elif coord_x_2[i] < coord_x_3[i]:
        #                coord_x_out[i] = coord_x_3[i]
        #                coord_y_out[i] = coord_y_3[i]
    	#		nzycheck[coord_x_out[i],coord_y_out[i],0] = 1
    	#		nzycheck[coord_x_out[i],coord_y_out[i],1] = flux_log10[i]
    	#	if flux_log10[i]>max_flux-prange:
      	#		nzycheck[coord_x_1[i],coord_y_1[i],2] = 1
      	#		nzycheck[coord_x_out[i],coord_y_out[i],2] = 1
	coord_x_out = np.zeros(len(coord_x_2_small))
        coord_y_out = np.zeros(len(coord_y_2_small))
        for i in range(len(flux_log10_small)):
  		nzycheck[coord_x_1_small[i],coord_y_1_small[i],0] = 1
  		nzycheck[coord_x_1_small[i],coord_y_1_small[i],1] = flux_log10_small[i]
  		if coord_x_2_small[i] >= coord_x_3_small[i]:
                        coord_x_out[i] = coord_x_2_small[i]
                        coord_y_out[i] = coord_y_2_small[i]
    			nzycheck[coord_x_out[i],coord_y_out[i],0] = 1
    			nzycheck[coord_x_out[i],coord_y_out[i],1] = flux_log10_small[i]
                elif coord_x_2_small[i] < coord_x_3_small[i]:
                        coord_x_out[i] = coord_x_3_small[i]
                        coord_y_out[i] = coord_y_3_small[i]
    			nzycheck[coord_x_out[i],coord_y_out[i],0] = 1
    			nzycheck[coord_x_out[i],coord_y_out[i],1] = flux_log10_small[i]
    		if which_flux == None or which_flux < 2 and flux_log10_small[i]>max_flux_small-prange:
      			nzycheck[coord_x_1_small[i],coord_y_1_small[i],2] = 1
      			nzycheck[coord_x_out[i],coord_y_out[i],2] = 1
      		elif which_flux == 2 and flux_log10_small[i]<min_flux_small+prange:
      			nzycheck[coord_x_1_small[i],coord_y_1_small[i],2] = 1
      			nzycheck[coord_x_out[i],coord_y_out[i],2] = 1

	#######################################################################
	# elemental names: elname(i) is the name of element with Z=i
	elname=	('none','H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al','Si','P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr','Mn','Fe',
        'Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb',
        'Te', 'I','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','W','Re','Os',
        'Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th','Pa','U','Np','Pu')


	#### create plot

	## define axis and plot style (colormap, size, fontsize etc.)
	if plotaxis==[0,0,0,0]:
  		xdim=10
  		ydim=6
	else:
  		dx = plotaxis[1]-plotaxis[0]
  		dy = plotaxis[3]-plotaxis[2]
  		ydim = 6
  		xdim = ydim*dx/dy

	format = 'pdf'
        # note that I had to set 'text.usetex': False, to avoid Exception in Tkinter callback.
        # and to make the GUI work properly. Why? some missing package?
	params = {'axes.labelsize':  15,
          'text.fontsize':   15,
          'legend.fontsize': 15,
          'xtick.labelsize': 15,
          'ytick.labelsize': 15,
          'text.usetex': False}
	plt.rcParams.update(params)
	fig=plt.figure(figsize=(xdim,ydim),dpi=100)
	axx = 0.10
	axy = 0.10
	axw = 0.85
	axh = 0.8
	ax=plt.axes([axx,axy,axw,axh])

	# color map choice for abundances
	cmapa = cm.jet
	# color map choice for arrows
	if which_flux == None or which_flux < 2:
	    cmapr = cm.autumn
	elif  which_flux == 2:
	    cmapr = cm.autumn_r
	# if a value is below the lower limit its set to white
	cmapa.set_under(color='w')
	cmapr.set_under(color='w')
	# set value range for abundance colors (log10(Y))
	norma = colors.Normalize(vmin=-20,vmax=0)
	# set x- and y-axis scale aspect ratio to 1
	ax.set_aspect('equal')
	#print time,temp and density on top
	#temp = '%8.3e' %ff['temp']
	#time = '%8.3e' %ff['time']
	#dens = '%8.3e' %ff['dens']

	#box1 = TextArea("t : " + time + " s~~/~~T$_{9}$ : " + temp + "~~/~~$\\rho_{b}$ : " \
        #      + dens + ' g/cm$^{3}$', textprops=dict(color="k"))
	#anchored_box = AnchoredOffsetbox(loc=3,
        #        child=box1, pad=0.,
        #        frameon=False,
        #        bbox_to_anchor=(0., 1.02),
        #        bbox_transform=ax.transAxes,
        #        borderpad=0.,
        #        )
	#ax.add_artist(anchored_box)

	# Add black frames for stable isotopes
	f = open('stable.dat')

	head = f.readline()
	stable = []

	for line in f.readlines():
  		tmp = line.split()
  		zz = int(tmp[2])
  		nn = int(tmp[3])
  		xy = nn-0.5,zz-0.5
  		rect = Rectangle(xy,1,1,ec='k',fc='None',fill='False',lw=3.)
  		rect.set_zorder(2)
  		ax.add_patch(rect)

 	apatches = []
  	acolor = []
  	m = 0.8/prange#0.8/prange
  	if which_flux == None or which_flux < 2:
  	    vmax=np.ceil(max(flux_log10_small))
  	    vmin=max(flux_log10_small)-prange
  	    b=-vmin*m+0.1
  	elif which_flux == 2:
  	    vmin=min(flux_log10_small)
  	    vmax=np.ceil(min(flux_log10_small)+prange)
  	    b=vmax*m+0.1
  	if which_flux == None or which_flux < 3:
  	    normr = colors.Normalize(vmin=vmin,vmax=vmax)
  	    print 'vmin and vmax =',vmin,vmax
  	ymax=0.
  	xmax=0.

	for i in range(len(flux_log10_small)):
    		x = coord_x_1_small[i]
    		y = coord_y_1_small[i]
    		dx = coord_x_out[i]-coord_x_1_small[i]
    		dy = coord_y_out[i]-coord_y_1_small[i]
                if plot_type == 0:
                    if which_flux == None or which_flux < 2:
    			    if flux_log10_small[i]>=vmin:
      				    arrowwidth = flux_log10_small[i]*m+b
      				    arrow = Arrow(x,y,dx,dy, width=arrowwidth)
     				    if xmax<x:
        			    	xmax=x
      				    if ymax<y:
        			    	ymax=y
      				    acol = flux_log10_small[i]
      				    apatches.append(arrow)
      				    acolor.append(acol)
                    elif which_flux == 2:
    			    if flux_log10_small[i]<=vmax:
      				    arrowwidth = -flux_log10_small[i]*m+b
      				    arrow = Arrow(x,y,dx,dy, width=arrowwidth)
     				    if xmax<x:
        			    	xmax=x
      				    if ymax<y:
        			    	ymax=y
      				    acol = flux_log10_small[i]
      				    apatches.append(arrow)
      				    acolor.append(acol)
    		elif plot_type == 1 and which_flux != 2:
			if x==I_am_the_target[0] and y==I_am_the_target[1] and flux_log10_small[i]>=vmin:
      				arrowwidth = flux_log10_small[i]*m+b
      				arrow = Arrow(x,y,dx,dy, width=arrowwidth)
     				if xmax<x:
        				xmax=x
      				if ymax<y:
        				ymax=y
      				acol = flux_log10_small[i]
      				apatches.append(arrow)
      				acolor.append(acol)
			if x+dx==I_am_the_target[0] and y+dy==I_am_the_target[1] and flux_log10_small[i]>=vmin:
      				arrowwidth = flux_log10_small[i]*m+b
      				arrow = Arrow(x,y,dx,dy, width=arrowwidth)
     				if xmax<x:
        				xmax=x
      				if ymax<y:
        				ymax=y
      				acol = flux_log10_small[i]
      				apatches.append(arrow)
      				acolor.append(acol)
    		elif plot_type == 1 and which_flux == 2:
			if x==I_am_the_target[0] and y==I_am_the_target[1] and flux_log10_small[i]<=vmax:
      				arrowwidth = -flux_log10_small[i]*m+b
      				arrow = Arrow(x,y,dx,dy, width=arrowwidth)
     				if xmax<x:
        				xmax=x
      				if ymax<y:
        				ymax=y
      				acol = flux_log10_small[i]
      				apatches.append(arrow)
      				acolor.append(acol)
			if x+dx==I_am_the_target[0] and y+dy==I_am_the_target[1] and flux_log10_small[i]<=vmax:
      				arrowwidth = -flux_log10_small[i]*m+b
      				arrow = Arrow(x,y,dx,dy, width=arrowwidth)
     				if xmax<x:
        				xmax=x
      				if ymax<y:
        				ymax=y
      				acol = flux_log10_small[i]
      				apatches.append(arrow)
      				acolor.append(acol)

 	#apatches = []
  	#acolor = []
  	#m = 0.8/prange
  	#vmax=np.ceil(max(flux_log10))
  	#vmin=max(flux_log10)-prange
  	#b=-vmin*m+0.1
  	#normr = colors.Normalize(vmin=vmin,vmax=vmax)
  	#ymax=0.
  	#xmax=0.

 	#for i in range(len(flux_log10)):
    	#	x = coord_x_1[i]
    	#	y = coord_y_1[i]
    	#	dx = coord_x_out[i]-coord_x_1[i]
    	#	dy = coord_y_out[i]-coord_y_1[i]
        #       if plot_type == 0:
    	#		if flux_log10[i]>=vmin:
      	#			arrowwidth = flux_log10[i]*m+b
      	#			arrow = Arrow(x,y,dx,dy, width=arrowwidth)
     	#			if xmax<x:
        #				xmax=x
      	#			if ymax<y:
        #				ymax=y
      	#			acol = flux_log10[i]
      	#			apatches.append(arrow)
      	#			acolor.append(acol)
    	#	elif plot_type == 1:
	#		if x==I_am_the_target[0] and y==I_am_the_target[1] and flux_log10[i]>=vmin:
      	#			arrowwidth = flux_log10[i]*m+b
      	#			arrow = Arrow(x,y,dx,dy, width=arrowwidth)
     	#			if xmax<x:
        #				xmax=x
      	#			if ymax<y:
        #				ymax=y
      	#			acol = flux_log10[i]
      	#			apatches.append(arrow)
      	#			acolor.append(acol)
	#		if x+dx==I_am_the_target[0] and y+dy==I_am_the_target[1] and flux_log10[i]>=vmin:
      	#			arrowwidth = flux_log10[i]*m+b
      	#			arrow = Arrow(x,y,dx,dy, width=arrowwidth)
     	#			if xmax<x:
        #				xmax=x
      	#			if ymax<y:
        #				ymax=y
      	#			acol = flux_log10[i]
      	#			apatches.append(arrow)
      	#			acolor.append(acol)
	#


        	xy = x-0.5,y-0.5
        	rect = Rectangle(xy,1,1,ec='k',fc='None',fill='False',lw=1.)
        	rect.set_zorder(2)
        	ax.add_patch(rect)
        	xy = x+dx-0.5,y+dy-0.5
        	rect = Rectangle(xy,1,1,ec='k',fc='None',fill='False',lw=1.)
        	rect.set_zorder(2)
        	ax.add_patch(rect)


  	a = PatchCollection(apatches, cmap=cmapr, norm=normr)
  	a.set_array(np.array(acolor))
  	a.set_zorder(3)
  	ax.add_collection(a)
  	cb = plt.colorbar(a)

  	# colorbar label
  	if which_flux == None or which_flux == 0:
  	    cb.set_label('log$_{10}$(f)')
  	elif which_flux ==1:
  	    cb.set_label('log$_{10}$(E)')
  	elif which_flux ==2:
  	    cb.set_label('log$_{10}$(timescale)')


  	# plot file name
  	graphname = 'flow-chart.'+format

	# decide which array to take for label positions
	iarr = 2

	# plot element labels
  	for z in range(nzmax):
    		try:
      			nmin = min(np.argwhere(nzycheck[:,z,iarr-2]))[0]-1
      			ax.text(nmin,z,elname[z],horizontalalignment='center',verticalalignment='center',fontsize='medium',clip_on=True)
    		except ValueError:
      			continue

	# plot mass numbers
	if imlabel==1:
  		for z in range(nzmax):
     			for n in range(nnmax):
        			a = z+n
				if nzycheck[n,z,iarr-2]==1:
          				ax.text(n,z,a,horizontalalignment='center',verticalalignment='center',fontsize='small',clip_on=True)

	# plot lines at magic numbers
	if imagic==1:
  		ixymagic=[2, 8, 20, 28, 50, 82, 126]
  		nmagic = len(ixymagic)
  		for magic in ixymagic:
    			if magic<=nzmax:
      				try:
        				xnmin = min(np.argwhere(nzycheck[:,magic,iarr-2]))[0]
        				xnmax = max(np.argwhere(nzycheck[:,magic,iarr-2]))[0]
        				line = ax.plot([xnmin,xnmax],[magic,magic],lw=3.,color='r',ls='-')
      				except ValueError:
        				dummy=0
    			if magic<=nnmax:
      				try:
        				yzmin = min(np.argwhere(nzycheck[magic,:,iarr-2]))[0]
        				yzmax = max(np.argwhere(nzycheck[magic,:,iarr-2]))[0]
        				line = ax.plot([magic,magic],[yzmin,yzmax],lw=3.,color='r',ls='-')
      				except ValueError:
        				dummy=0

	# set axis limits
	if plotaxis==[0,0,0,0]:
		ax.axis([-0.5,xmax+0.5,-0.5,ymax+0.5])
	else:
  		ax.axis(plotaxis)

	# set x- and y-axis label
	ax.set_xlabel('neutron number')
	ax.set_ylabel('proton number')
	if which_flux == None or which_flux == 0:
		max_flux_label="max flux = "+str('{0:.4f}'.format(max_flux))
	elif which_flux == 1:
		max_flux_label="max energy flux = "+str('{0:.4f}'.format(max_flux))
	elif which_flux == 2:
		min_flux_label="min timescale [s] = "+str('{0:.4f}'.format(min_flux))
	if print_max_flux_in_plot:
	    if which_flux == None or which_flux < 2:
		    ax.text(plotaxis[1]-1.8,plotaxis[2]+0.1,max_flux_label,fontsize=10.)
	    elif which_flux == 2:
		    ax.text(plotaxis[1]-1.8,plotaxis[2]+0.1,min_flux_label,fontsize=10.)



	fig.savefig(graphname)
	print graphname,'is done'
	if which_flux == None or which_flux < 2:
	    print max_flux_label,'for reaction =',ind_max_flux+1
	elif which_flux == 2:
	    print min_flux_label,'for reaction =',ind_min_flux+1

	plt.show()

