import math
import numpy as np
import matplotlib.pyplot as plt
import sconvert
import PeakCompare
import time


import sys, os
sys.path.append('C:\\Users\\loogo\\Desktop\\FemtoSecond\\Latest')

### Basic Values ###
g = 1.6e-13
c = 3e8

# Raman transition #
omegaRaman=775
dephasingTime=6.6
HWHMRaman = 0.2206/(dephasingTime*1e-12)/3e10/2

#HWHMRaman = 100 in wavenumbers

convolution = 1

## Note that a higher value is more accurate, but also much more taxing                  
## Also note that this will affect the value for g
## Define the fibre characteristics 

FibreLength= 1      ##in meters 
Intervals= 200      ##number of points per meter; this number must be large for the program to work 
plotPoints= 1000    ##when plotting in 3d, it's best to have only 1000 points 
fibreRadius= 250    ##in microns 
Pressure= 3         ##in atmospheres 
u11= 2.405          ##zero of the Bessel function J[0](x) 


dz= FibreLength/Intervals

## Define largest Stokes and Anti-Stokes order you want to include 
S_Orders= 20        ## including the long-wavelength beam 
A_Orders= 20        ## You can't interactively change the scale on a surface plot, so set the number you want to plot 
minorder= -S_Orders ## use -S_Orders for all of them 
maxorder= A_Orders  ## use A_Orders for all of them

## Set the frequency resolution 
## Note that a typical bandwidth is about 80 wavenumbers, so a spacing of 16 is appropriate 
omegaRes= 7         ## In wavenumbers 

## Set the initial intensity of the Pump, and then the energy ratio between it and the Stokes 
Intensity= 6e12     ##in W/cm^2, but I don't think it matters, as g is a parameter to be determined 
ratio= 1.5
  
## Finally, set the smallest intensity (relative to the Pump), that you want to plot 
plotmin= 1e-4       
## Labeling convention: short-wavelength beam (Pump beam) is order# S_Orders+3 
## long-wavelength beam (Stokes beam) is order# S_Orders+2 
## I'm adding two extra orders on either end for programming considerations 
Pump= S_Orders + 3
Stokes= S_Orders + 2

############################################################################
# Determining resolution and constructing the input spectrum from our data # 
############################################################################

omegaIntervals= math.ceil(omegaRaman/omegaRes) ## number of intervals per order

##  I'm going to make sure that there are an odd number of intervals per order.  By forcing 
#  an odd number, I can define the array elements by using the orders as a reference. 
#  Notice how in the diagram below there are two points on either side of each order 
#  when there's 5 intervals per order. 
#  ..|....|....|....|....|.. 
#  Equivalently, the # of side elements = (# of intervals - 1)/2 

if ( omegaIntervals%2 == 0):
	omegaIntervals= omegaIntervals + 1

side= (omegaIntervals - 1)/2 
## This now sets the actual resolution, which will be as good as or better than 

## The specified resolution 
omegaRes= omegaRaman/omegaIntervals


## Read in the Regen Spectrum
RegenSpec = sconvert.readAndNormalize("C:\\Users\\loogo\\Desktop\\FemtoSecond\\AND064.dat")
omegaReg = RegenSpec[0]
IntensReg = RegenSpec[1]

## Determine the frequency spacing 
points = np.size(omegaReg) 
freqRange = omegaReg[0] - omegaReg[points-1] 
dOmegaReg = freqRange/points 

## Determine the desired index spacing from the specified resolution.  This should 
## give a resolution that is as good or better than that specified 
spacing= math.floor(omegaRes/dOmegaReg) 
omegaRes= spacing*dOmegaReg
RegIntervals= math.floor(points/spacing) 

##Calculate area under the two peaks?
Peaks = PeakCompare.compareBasic(omegaRaman,RegenSpec)

Peak1 = Peaks[0][0]
Peak2 = Peaks[1][0]
Indices = Peaks[2]

PumpArea = Peak1[1]
StokesArea = Peak2[1]

pumpindex = Indices[0]
omegaPump = omegaReg[pumpindex]



## Initialize the Frequency and Intensity arrays 
## I'm adding a couple extra orders on either end as it's needed later 

#omega= (np.arange(omegaPump - (omegaIntervals*(S_Orders + 2) + side)*omegaRes, np.arange(omegaRes, omegaPump + (omegaIntervals*(A_Orders + 2) + side)*omegaRes)))

'''
'''
## added params
simulationMaxW = 40000
simulationMinW = 4000
frewWindowHALF = 200

minW = omegaReg[0]
maxW = omegaReg[-1]
rangeW = maxW - minW
wres= np.min(np.abs(omegaReg - np.roll(omegaReg, 1)))

maxCalcW = (simulationMaxW/wres)//1
minCalcW = (simulationMinW/wres)//1

numPoints  = maxCalcW-minCalcW
omega = np.linspace(simulationMinW,simulationMaxW,int (numPoints))

'''
'''
In= np.zeros((np.size(omega)))
Intensindex = int((minW-simulationMinW)/wres//1)


for n in range(np.size(IntensReg)):
    shiftedIndex = Intensindex + n
    if(n > Indices[0]-frewWindowHALF and n < Indices[0]+frewWindowHALF) or (n > Indices[1]-frewWindowHALF and n < Indices[1]+frewWindowHALF):
        In[shiftedIndex] = IntensReg[n]
    '''
    else:
        In[shiftedIndex] = (1/1e10)
'''
RegenSpecTEST = np.zeros(shape=(np.size(omega[5600:6730]),2)).transpose()

RegenSpecTEST[0] = omega[5600:6730]
RegenSpecTEST[1] = In[5600:6730]

'''
plt.plot(omega[5600:6730],In[5600:6730])
plt.show()
'''
'''
'''

## GVD= 0000e-6; %%I've set it to 10000; for Steve it was -1877 
## TOD= 000000e-9; %%I've set it to 1000000; for Steve it was 18719 
## extraphase= GVD/2.*(omega*100*3e8*2*pi).^2 + TOD/6.*(omega*100*3e8*2*pi).^3; 
## lambda= 1e4./omega; %%this gives the wavelength in microns

##########################
#Set the phase dispersion# 
##########################

ngas= np.zeros((1, np.size(omega)),dtype = "complex_")
k= np.zeros((1, np.size(omega)),dtype = "complex_")
# ngas= ngas + 1; %%comment what's below to set no dispersion
# Sellmeier equations for the gas 
ngas= 1 + Pressure*0.7e30/((31.9e15)**2 - (2.*math.pi*c*100.*omega)**2) ## from physrevA68_023812p15_2003 
# sellmeier equation for fused silica 
# nsilica= sqrt(1.2955 + 0.80985*lambda.^2./(lambda.^2 - .0107945) + 0.91714*lambda.^2./(lambda.^2 - 100)); 
nsilica= 1.5
n1= (nsilica**2 + 1)/(2.*math.sqrt(nsilica**2 - 1))
## Axial wavevector for the EH11 hybrid mode 
k= 2.*math.pi*omega*100.*ngas*(1 - 1/2.*(u11/(omega*100)/(2.*math.pi*fibreRadius))**2)*(1 - 1j*n1/(omega*100)/math.pi/fibreRadius)

###########################
#Electric Field Amplitudes# 
###########################
A= np.zeros(np.size(omega),dtype = "complex_")
B= np.zeros(np.size(omega),dtype = "complex_")
C= np.zeros(np.size(omega),dtype = "complex_")
P= np.zeros((1, Intervals + 1),dtype = "complex_")[0]
A= np.sqrt(In) 
B= np.copy(A)
C= np.copy(A)

## I'll use A, B, and C to keep track of things as I go through the program 

## These next two things are what I'll use to graph the pulses
# 
# 3d plot 
print()
surfacegraph = np.zeros(np.shape(omega),dtype = "complex_")
surfacegraph = surfacegraph.transpose()

## Propagate in z 
for z in range(Intervals):    
    ##for x in range(-(math.ceil(convolution/2/omegaRes)), math.ceil(convolution/2/omegaRes)):
    for x in range(1):                  
        Q= 0
        #omegaintervals is 111 so 334 -> 43*111 = 4773
        for n in range(3*omegaIntervals + 1, (Pump + A_Orders)*omegaIntervals):
        #for n in range(0, np.size(omega)):         
            #Omega= omega[n] - omega[n - omegaIntervals + x]             
            Q= Q + B[n]*np.conjugate(B[n - omegaIntervals + x])*np.exp( -1j *(k[n] - k[n - omegaIntervals + x])*dz)#*1/(omegaRaman**2 - Omega**2 - 2*-1j*Omega*HWHMRaman);                 
            P[z + 1] = Q.real

        #omegaintervals is 111 so 223 -> 43*111 = 4773
        #3*omegaIntervals + 1: (Pump + A_Orders)*omegaIntervals
        for n in range(3*omegaIntervals + 1, (Pump + A_Orders)*omegaIntervals): 
        #for n in range(200, np.size(omega)-200):             
            C[n] = B[n] + g*Intensity/ngas[n]*omega[n]/omega[Pump]*(np.conjugate(Q)*B[n + omegaIntervals - x]*np.exp(1j*(k[n] - k[n + omegaIntervals - x])*dz)- Q*B[n - omegaIntervals + x]*np.exp(1j*(k[n] - k[n - omegaIntervals + x])*dz))*dz
            ## ow to define the values for a 3d plot in frequency and distance  
    '''
    for n in range(2*omegaIntervals + 1, (Pump + A_Orders)*omegaIntervals):
        print(surfacegraph)
        print(np.shape(surfacegraph))
        surfacegraph[z + 1, n] = B[n]*np.conjugate(B[n])*Intensity
        '''
    if(z%10 == 0):
        print("Running: ", z , " / ",Intervals)       
    B = C
    
#######
#Plots#
#######
## Now I'm changing the convention - the Pump is of order 0, and the Stokes is of order -1 
print(-(S_Orders + 2) - side/omegaIntervals)
print( (A_Orders + 2) + side/omegaIntervals)
print( 1/omegaIntervals)
#X= range(-(S_Orders + 2) - side/omegaIntervals, (A_Orders + 2) + side/omegaIntervals, 1/omegaIntervals)
print()

Y= ((B*np.conjugate(B))*Intensity)

print(Y[5600:6730],"")
L= np.log10((B*np.conjugate(B) + plotmin)*Intensity)

plt.figure(1) 

plt.plot(omega, Y) 

##print(np.shape(L))

#omega[5600:6730]
# L[5600:6730]
#np.linspace(5600*wres,6730*wres,6730-5600)
#plt.plot(omega[1000:14000], Y[1000:14000])
#plt.show() 
#plt.xlim([-0.5, 7])

SX= np.arange(0,FibreLength/Intervals, FibreLength) 
SY= (minorder - np.arange(side/omegaIntervals, 1/omegaIntervals, maxorder + side/omegaIntervals)) 

## See diagram below 
##    ..|....|....|....|....|....|....|.. 
##               min   P   max 
## The 3d plot works okay for 500 intervals; above and below that it gets messy 
## plotinterval= floor(Intervals/500)
## Delete elements up to just before the side of the minorder (use the diagram)  
'''
surfacegraph[:, 1:(Pump + minorder - 1)*omegaIntervals - 1] = []     
'''
## Delete elements beyond the side of the maxorder, note that I've already deleted     
## some elements in the line above (the 2 is because of that extra order I added earlier)  
'''
surfacegraph[:, (-minorder + maxorder + 1)*omegaIntervals:(-minorder + A_Orders + 3)*omegaIntervals - 1]= []     
'''
## Delete elements in z to have only 1000 elements - this makes the best plot     
plotInterval= math.floor((Intervals)/plotPoints)

for n in range(1,plotPoints - 1):        
    '''
    surfacegraph[(n + 1) : (n + plotInterval - 1), :] = []         
    '''
    SX[(n + 1): (n + plotInterval - 1)] = []     
    
'''
if (np.size(surfacegraph) > plotPoints):       
	surfacegraph[plotPoints + 1 : np.size(surfacegraph), :] = []     
    
if (np.size(SX) > plotPoints):     
	SX[(plotPoints + 1):np.size(SX)] = []     

SZ= surfacegraph
'''
## Delete some elements so that there's at most 500 points along z 
##   if (Intervals> 500) 
##       for n= Intervals + 1: -1: 1 
##           if (mod(n, plotinterval)~= 1) 
##               SZ(n, :)= [] 
##           end 
##       end 
##   end
## SL= log10(plotmin.*Intensity + surfacegraph)  
## SL = loggraph
##   SL(:, 1: (Pump + minorder - 1)*omegaIntervals - 1) = [] 
##   SL(:, (-minorder + maxorder + 1)*omegaIntervals
##       (-minorder + A_Orders + 3)*omegaIntervals) = [] 
##   if (Intervals> 500) 
##       for n= Intervals + 1: -1: 1 
##           if (mod(n, plotinterval)~= 1) 
##               SL(n, :)= [] 
##               SX(n)= [] 
##           end 
##       end 
##   end
'''
plt.figure(2)
plt.mesh(SY, SX, SZ)
'''
## mesh(SY, SX, SL) 
## shading interp 
 
