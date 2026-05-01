# Author: Jo Vazquez
# Date: 1 September 2022
# Gaseous Astrophysics Fall 2022 Pre-Project 1: Kinematics of HI Gas in the Milky Way
############ I am going to begin this notebook by importing crucial software in order to fit the spectra  ###########

############ I will also define universal constants and functions used throughout the code
import numpy as np  # This package allows me to manipulate arrays of data easily and extract data from files
import matplotlib.pyplot as plt # This package allows me to create plots of my spectra
from scipy.optimize import curve_fit  # This package will allow me to fit my spectra to a Gaussian distribution
from scipy.special import erf  #Importing the error function for integral calculation
import os   # This is for navigating the operating system and finding the files
import argparse  # For the command-line options below
from gaussdecomp import spectrum,fitter,utils
# Here I will define a Gaussian of multiplicity N in terms of the coordinate, amplitude at the mean, the mean, 
# standard deviation, and y axis offset respectively. I will define the function with a lambda, but first, let's make the string to execute:


# Input HI4PI spectra are named by sightline, e.g. 'spectrum_l=45_b=-45.txt' for
# galactic longitude 45, latitude -45. Both directories are configurable so the
# script runs from anywhere; the defaults match the layout in this repository.
parser = argparse.ArgumentParser(
	description='Decompose HI21cm brightness-temperature spectra into Gaussian '
	            'components and derive HI column densities along the b=-45 strip.')
parser.add_argument('--data-dir', default=os.path.join('..', 'Data'),
                    help='directory holding the HI4PI spectrum_l=*_b=*.txt files')
parser.add_argument('--image-dir', default=os.path.join('..', 'Images'),
                    help='directory to write figures into')
parser.add_argument('--savefigs', action='store_true',
                    help='write figures to --image-dir instead of only displaying them')
args = parser.parse_args()

data_dir = args.data_dir
image_dir = args.image_dir
savefigs = args.savefigs

if savefigs:
	os.makedirs(image_dir, exist_ok=True)
	print('Saving figures to %s' % image_dir)
else:
	print('Not saving figures (pass --savefigs to write them out)')
print('')



k = 1.82240e18 # atoms K-1 cm-2. This is the constant of integration to find the HI column density from
#                brightness temperature.




def one_gaussian(x,x0,sigma,Amp,y0):
	return Amp * np.exp((-(x-x0)**2)/(2*sigma**2))+y0


EBHI_ls = [-180,60,90,120,150,180]
GASS_ls=  [-150,-120,-90,-30,0,30]

#Making a list of galactic longitudes to iterate through


file_preamble = os.path.join(data_dir, 'spectrum_l=')
file_postamble = '_b=-45.txt'    

#Above is my file naming convention after extracting the HI4PI dataset. The files are names as such: spectrum_l=-180_b=-45.txt
lrange = np.arange(-180,200,30)  #Starts at -210 but does not include -210
lrange = lrange.tolist()         #Changing lrange to a list to make iterating through easier



#Below I begin a list of NHI sums, NHI components, v0 components, and fwhm components for each sightline
NHI_SUMS = []
ALL_NHI = []
V0 = []
FWHM = []

ALL_POPTS = [] #Creating an empty list to input best-fit multi-gaussian fit parameters
for l in lrange:

##########################################################################################################################################################################
	l = int(l)
	# l = -180 degrees: EBHI data is available.

	fname = file_preamble + str(int(l)) + file_postamble


	#Skipping the text within the files to get to that juicy data! 

	if l in [-180,60,180]:
		skip_header,skip_footer = 17,1780
	elif l in [150,120,90,30,0,-150,-120,-90,-60,-30]:
		skip_header,skip_footer = 14,785
	



	vel,T_B = np.genfromtxt(fname,skip_header=skip_header,skip_footer=skip_footer).T  #Extracting the data

	if l == -60:
		vel_mask = vel <= 65. #km s-1
		vel = vel[vel_mask]
		T_B = T_B[vel_mask]

	velrange = np.linspace(vel[0],vel[-1],5000)
	#Here I am extracting data from lines 17 on and ignoring the last 1780 lines as they are NOT EBHI data; .T means a matrix transpose

	#I will define a new multiple gaussian for every sightline, as each sightline requires a different number of components
	if l in [-180,-120,-150,90]:
		N = 2
	elif l in [-30,-60,-90,30,60,150,180]:
		N = 3
	elif l in [0]:
		N = 4
	elif l in [200]:
		N=1
	def multiple_gaussian(x,*popt):   #Here popt must have a length of 3N+1, where N is the number of components being fit. The last element is y0
		length = len(popt) 

		#Beginning the total gaussian sum
		total_gaussian = 0
		for i in range(int(N)):
			#There are three separate parameters to fit for each gaussian: a center, fwhm, and amplitude. 
			index0 = 3 * i
			#The individual gaussian 
			gauss = one_gaussian(x,popt[index0],popt[index0+1],popt[index0+2],y0=0)
			#Adding the individual Gaussian to the total gaussian
			total_gaussian += gauss


		#Adding y offset
		total_gaussian += popt[-1]

		return total_gaussian


	def individual_gaussians(x,popt,extract_popt=False):

		#This function extracts individual gaussians from the popt vector from a multi-Gaussian fit.

		num_components = int((len(popt)-1)/3)
		gauss_list = []
		popt_list = []
		for i in range(num_components):
			index0 = 3 * i
			gauss = one_gaussian(x,popt[index0],popt[index0+1],popt[index0+2],popt[-1])
			gauss_list.append(gauss)

			popt_list.append([popt[index0],popt[index0+1],popt[index0+2],popt[-1]])

		

		if extract_popt:
			return popt_list
		else:
			return gauss_list



	def integral_gaussian(x,a,b,*pars):


		#I am extracting out the integral of a gaussian over a finite range x in [a,b]

		x0,sigma,Amp,y0 = pars

		integral_multiplier = np.sqrt(np.pi) * Amp * sigma/(np.sqrt(2)) #This is the amplitude of the integral of a gaussian of the form Aexp(-(x-x0)^2/(2sigma^2))



		erf_stuff = erf((b-x0)/(sigma*np.sqrt(2))) - erf((a-x0)/(sigma*np.sqrt(2))) #Here I am taking the limits based on the Fundamental Theorem of Calculus

		integral = integral_multiplier * erf_stuff + y0 * (b - a)  # The final form, which includes the integral of the y offset, which is simply y0 * x over the range [a,b]
		return integral





	
	#A specific function for l = -90 to force a fit at v0 = -42 km s-1
	v0 = 10. #km s-1


	#Below I have the fit components guess lists for each longitude. Notice that each list has a length of 3N+1, where N is the # of components fitted

	if l == -180:
		guesses = [2.71277573e+00, 9.05410680e+00, 8.66035061e+00, 7.39925645e+00,2.73167354e+00, 2.05813887e+01, 3.03863478e-03] #Set up in the style[x0_1,sigma_1,Amp_1,x0_2,sigma_2,Amp_2,y0]

	if l == -150:
		guesses = [-16.0,4.0,1.0,-1.5,5.0,2.2,-4.0,14.0,2.0,0.0] # v0_0 = -15. km s-1; sigma_0 = 5 km s-1; Amp_0 = 2.K/(km s-1); v0_1 = 0. km s-1; sigma_1 = 10. km s-1; 
		#     #        #       #      #     # Amp_1 = 3. K/(km s-1); y0 = 0. K/(km s-1)

	if l == -120:
		guesses = [0.,5.,2.,7.,5.,10.,0.]
	if l == -90:
		guesses = [-24.,1.,0.5,10.,5.,2.,-2.5,5.,2.5,0.]

	if l == -60:
		guesses = [10.,5.,2.,-6.5,5.,4.,0.,7.,10.,0.]

	if l == -30:
		guesses = [-7.,5.,7.,1.,7.,5.,0.,5.,4.,0.]

	#A specific function for l = 120 to force a fit at v0 = -42 km s-1
	v0 = -37. #km s-1
	def l_120_gauss(x,*pars):
		gauss = one_gaussian(x,v0,15.,0.7,y0=0) + one_gaussian(x,pars[0],pars[1],pars[2],pars[3])
		return gauss

		#At this longitude, I am forcing the component at v0 = -37 km s-1 to get a good Milky Way fit. 

	if l == 180:
		guesses = [-3.,5.,6.,8.,4.,max(T_B),4.,3.,15.,0.]
	if l == 150:
		guesses = [-17.,4.,2.,-8.,3.,22.,-1.,2.,max(T_B),0.]
	if l == 120:
		multiple_gaussian = l_120_gauss #Renaming the function to the forced gaussian at -37 km s-1.
		guesses = [-1.0,10.,20.,0.]  #in the form [sigma_0,Amp_0,x0_1,sigma_1,Amp_1,y0]

	if l == 90:
		guesses = [-38.,5.,0.75,0.,10.,13.,0.]

	if l == 60:
		guesses = [-7.0,5.0,12.0,1.5,3.0,7.0,9.0,3.0,5.0,0.0]
	if l == 30:
		guesses = [-15.,5.,2.,-0.8,3.,max(T_B),2.,1.0,8.,0.]

	if l == 0:
		guesses = [3.,7.,7.,-15.,4.,10.,-5.,3.,3.,0.,3.,3.,0.]

	popt,pcov = curve_fit(multiple_gaussian,vel,T_B,guesses)  # Extracting the optimized fit parameters from scipy.optimize.curve_fit with guesses

	fit = multiple_gaussian(velrange,*popt)   # Functional form of the best-fit multi-gaussian.

	plt.title(r'$(\ell,b) =(%d \degree,-45 \degree)$' % l,fontsize=20)  # I am setting the title to '(l,b) = (-180 degrees,-45 degrees)', for example
	plt.xlabel(r'$v_{\rm{LSR}} \rm{(km}\;\rm{s}^{-1})$',fontsize=15) #Setting the x-axis to be LSR velocity
	plt.ylabel(r'$T_B$ ($\frac{\rm{K}}{\rm{km}\;\rm{s}^{-1}})$ ',fontsize=15) #Setting the y-axis to be T_B(K/km s-1)

	plt.plot(velrange,fit,'r-',alpha=0.8,linewidth=3.,label='Total Fit')  #Plotting the best-fit multi-gaussian
	plt.step(vel,T_B,'k-',linewidth=3.,label='Data')  #Plotting the data in histogram style.

	gauss_list = individual_gaussians(vel,popt)   #Making a list of the individual gaussians within the total fit
	popt_list = individual_gaussians(vel,popt,extract_popt=True) #Extracting the list of lists of optimized components
	for i in range(len(gauss_list)):
		gauss = gauss_list[i]   #Setting 'gauss' equal to the ith gaussian component
		popt_ind = popt_list[i] #Setting 'popt_ind' equal to the ith list of optimized gaussian components
		plt.plot(vel,gauss,linestyle='dashed',alpha=0.7,linewidth=1.5,label=r'$v_0 = $%3.4f km s$^{-1}$' % popt_ind[0])  
		#Plotting each individual component and labeling each with its center velocity
	plt.xlim([-50,50]) #Setting the velocity limit from -50 km s-1 to +50 km s-1
	plt.legend()
	if savefigs:
		plt.savefig(os.path.join(image_dir, 'HI_MilkyWay_b=-45_l=%d.pdf' % l))       #If the user elects to save a figure, it will save as '../Images/HI_MilkyWay_b=-45_l=-180.pdf', for example
	plt.show()  #Show the figure to the user

	popt_list = individual_gaussians(vel,popt,extract_popt=True)     # List of popt vectors for the N individual Gaussians fit to the spectrum.

	


	NHIs = []        # Beginning a list of all the HI column densities for the individual components

	v0s = []         # List of individual components of center velocity
	fwhms = []       # ...... of FWHMs
	NHI_sum_lin = 0  # Beginning a sum of all NHI components for the spectrum.

	for i in range(len(popt_list)):       #For each component
		NHI_lin = k*integral_gaussian(vel,-50.,50.,*popt_list[i])   # Individual fit NHI
		NHI_sum_lin += NHI_lin       #Adding NHI component to total NHI sum
		NHIs.append(np.log10(NHI_lin))     # Taking the log base 10, as is standard in the field

		v0 = popt_list[i][0]               # Extracting the center velocity, which is the 0th element of the gaussian, as defined in line 44
		v0s.append(v0)                     # Append the list of center velocities with this component's center velocity.

		sigma = popt_list[i][1]            # Extracting the velocity standard deviation, which is the 1st element of the gaussian, as defined in line 44

		fwhm = 2*np.sqrt(2*np.log(2)) * sigma  # Deriving fwhm from standard deviation
		fwhms.append(fwhm)                  #Appending the list of FWHMs

		print('') #blank line
		print('Component %d at l = %d deg:' % (i, l))   #Printing the component number and the longitude
		print('v0 = %2.3f km s-1' % v0)        #Printing (for example) 'v0 = -12.562 km s-1'
		print('FWHM = %2.3f km s-1' % fwhm)    #Printing (for example) 'FWHM = 10.081 km -1'
		print('log(NHI/cm-2) = %2.3f dex' % np.log10(NHI_lin))    #Printing (for example) 'log(NHI/cm-2) = 20.014 dex'
	
	print('_____________________________')   #line to separate the output for different longitudes


	NHI_sum = np.log10(NHI_sum_lin)        # Converting NHI to log space
	NHI_SUMS.append(NHI_sum)               # Appending the list of total NHI
	
	V0.append(v0s)                         # Appending the list of list of center velocities at different longitudes with the list of all center vel components
	
	FWHM.append(fwhms)                     #..... FWHM componnets

	ALL_NHI.append(NHIs)                   #..... NHI components



new_lrange = []      #This will a 1D list of longitudes, which will include repeated longitudes for multiple components
new_NHI = []         #All NHIs, but converted to a 1D array
new_V0 = []          #..... All v0s
new_FWHM = []        #..... All FWHMs
for i in range(len(lrange)):
	for j in range(len(ALL_NHI[i])):

		#Below, I am appending all lists such that they are all 1D lists so I can plot them
		new_lrange.append(lrange[i])    
		new_NHI.append(ALL_NHI[i][j])
		new_V0.append(V0[i][j])
		new_FWHM.append(FWHM[i][j])





##### CENTER VELOCITY VS. LONGITUDE PLOT ############

plt.scatter(new_lrange,new_V0,c=new_NHI,cmap='viridis') #Plot center vel vs. longitude with a colormap dependent on NHI
plt.xlabel(r'Galactic Longitude $(^{\circ})$')          # Make the x label galactic longitude
plt.ylabel(r'$v_0$ (km s$^{-1}$)')						# Make the y label v_0 (km s-1)
plt.title(r'Center velocities and column densities vs. galactic longitude at $b=-45^{\circ}$', fontsize=10.5)   #Make a title w/ fontsize 10.5
for l in lrange:
	plt.vlines(l,ymin=min(new_V0),ymax=max(new_V0),colors='gray',linestyles='dashed',alpha=0.5)    # Plot dashed vertical lines at all longitude vals to make the plot easier to read

plt.colorbar(label=r'$\log_{10}(N_{\rm{HI}}/\rm{cm}^{-3})$')     #Make a colorbar and label with NHI
if savefigs:
	plt.savefig(os.path.join(image_dir, 'v0_vs_long.pdf'))       #If the user so wishes, the figure will save in the '../Images' folder
plt.show()                                        #Showing the plot to the user



##### FWHM VELOCITY VS. LONGITUDE PLOT ############

plt.scatter(new_lrange,new_FWHM,c=new_NHI,cmap='viridis') #Plot FWHM vel vs. longitude with a colormap dependent on NHI
plt.xlabel(r'Galactic Longitude $(^{\circ})$')				# Make the x label galactic longitude
plt.ylabel(r'$v_{\rm{FWHM}}$ (km s$^{-1}$)')				# Make the y label v_FWHM (km s-1)
plt.title(r'FWHM velocities and column densities vs galactic longitude at $b=-45^{\circ}$', fontsize=11.) #Make a title w/ fontsize 11.
for l in lrange:
	plt.vlines(l,ymin=min(new_FWHM),ymax=max(new_FWHM),colors='gray',linestyles='dashed',alpha=0.5) # Plot dashed vertical lines at all longitude vals to make the plot easier to read


plt.colorbar(label=r'$\log_{10}(N_{\rm{HI}}/\rm{cm}^{-3})$')  #Make a colorbar and label with NHI
if savefigs:
	plt.savefig(os.path.join(image_dir, 'fwhm_vs_long.pdf'))     #If the user so wishes, the figure will save in the '../Images' folder
plt.show()                                        #Showing the plot to the user



##### LOG(NHI)  VS. LONGITUDE PLOT ############

plt.plot(lrange,NHI_SUMS,'ko')   #Plot total NHI vs. longitude with a colormap dependent on NHI
plt.xlabel(r'Galactic Longitude $(^{\circ})$')  # Make the x label galactic longitude
plt.ylabel(r'$\log_{10}(N_{\rm{HI},\rm{total}}/\rm{cm}^{-3})$') # Make the y label v_FWHM (km s-1)
plt.title(r'Column Density vs. Galactic longitude at $b=-45^{\circ}$') #Make a title w/ fontsize 11.
for l in lrange:
	plt.vlines(l,ymin=min(NHI_SUMS),ymax=max(NHI_SUMS),colors='gray',linestyles='dashed',alpha=0.5)# Plot dashed vertical lines at all longitude vals to make the plot easier to read

##### NO COLORBAR NEEDED  #####
if savefigs:
	plt.savefig(os.path.join(image_dir, 'NHI_vs_long.pdf'))#If the user so wishes, the figure will save in the '../Images' folder

plt.show()                                        #Showing the plot to the user



