import math
from typing import final
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from scipy.fft import fft, ifft
import sconvert
import PeakCompare
import time

### Basic Values ###
g = 1.6e-13
c = 3e8

# File Locations #
filePath = "C:\\Users\\loogo\\Desktop\\FemtoSecond\\AND064.dat"

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
Intervals= 50       ##number of points per meter; this number must be large for the program to work 
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
Intensity = 1
#Intensity= 6e12     ##in W/cm^2, but I don't think it matters, as g is a parameter to be determined 
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

## Read in the Regen Spectrum
#RegenSpec = sconvert.readAndNormalizeOverridePoints("C:\\Users\\loogo\\Desktop\\FemtoSecond\\AND064.dat",150)
## DEFAULT numPoints is 1130

'''
############
'''

# TODO: #RegenSpec = sconvert.manualInterpolateRegenSpec(RegenSpec, 1000)
# ^ Add in support for manual number of points
# TODO: Add support for manual specify Stokes/Pump Frequencies
#RegenSpec = sconvert.readAndNormalize(filePath)
numOrders=30
ramanSep=587
pump=18939
def buildData():
    """
    Import the regen spectrum from the specified file path.
    Determine resolution, frequency spacing, index spacing.
    """
    #RegenSpec = sconvert.readAndNormalize(filePath)
    numOrders=30
    ramanSep=587
    RegenSpec = sconvert.defineBasicTWO(pump,ramanSep,numOrders,0,Intensity)
    omegaReg = RegenSpec[0]
    IntensReg = RegenSpec[1]

    ## The specified resolution 
    omegaRes= omegaRaman/omegaIntervals

    ## Determine the frequency spacing 
    points = np.size(omegaReg) 
    freqRange = omegaReg[0] - omegaReg[points-1] 
    dOmegaReg = freqRange/points 

    ## Determine the desired index spacing from the specified resolution.  This should 
    ## give a resolution that is as good or better than that specified 
    spacing= math.floor(omegaRes/dOmegaReg) 
    omegaRes= spacing*dOmegaReg
    RegIntervals= math.floor(points/spacing) 
    return (RegenSpec,omegaReg,IntensReg,omegaRes,points,freqRange,dOmegaReg,spacing,RegIntervals)

RegenSpec,omegaReg,IntensReg,omegaRes,points,freqRange,dOmegaReg,spacing,RegIntervals = buildData()


# TODO: Rework PeakCompare -- doesnt feel right to me
# ^ perhaps split input in half to begin with then compare
# ^ afterwards to avoid issues

def pumpAndStokesInformation():
    """
    Determine the two peaks from the input RegenSpectrum.
    Returns peaks and indices
    """
    Peaks = PeakCompare.compareNormalized(omegaRaman,RegenSpec)

    PeakData1 = np.ravel(Peaks[0])
    PeakData2 = np.ravel(Peaks[1])
    Peak1 = PeakData1[0]
    Peak2 = PeakData2[0]
    Indices = Peaks[2]

    PumpArea = Peak1
    StokesArea = Peak2

    pumpindex = Indices[0]
    omegaPump = omegaReg[pumpindex]
    return(Peaks,Peak1,Peak2,Indices,PumpArea,StokesArea,pumpindex,omegaPump)
    
Peaks,Peak1,Peak2,Indices,PumpArea,StokesArea,pumpindex,omegaPump = pumpAndStokesInformation()

# TODO: Rework this, slightly inaccurate due to using linspace
# BUG: ^ has been fixed, watch for bugs
def initializeFrequencies():
    """
    Takes a max and min frequency. Find the lowest & highest
    frequencies from the input. Finds the range and resolution.
    Finds the max and min frequencies rectified with resolution.
    Finds number of points between the two ends, uses
    linspace to create a range from rectified min and max
    using numPoints points.
    """
    simulationMaxW = 70000
    simulationMinW = 0
    frewWindowHALF = 100

    minW = omegaReg[0]
    maxW = omegaReg[-1]
    rangeW = maxW - minW
    wres= np.min(np.abs(omegaReg - np.roll(omegaReg, 1)))

    maxCalcW = (simulationMaxW/wres)//1
    minCalcW = (simulationMinW/wres)//1

    numPoints  = maxCalcW-minCalcW
    #omega = np.linspace(simulationMinW,simulationMaxW,int (numPoints))
    omega = np.arange(simulationMinW,simulationMaxW,wres)
    return simulationMaxW,simulationMinW,frewWindowHALF,minW,maxW,wres,maxCalcW,minCalcW,numPoints,omega

simulationMaxW,simulationMinW,frewWindowHALF,minW,maxW,wres,maxCalcW,minCalcW,numPoints,omega = initializeFrequencies()

# BUG: v Fixed; watch for bugs
# TODO: Find more efficient way of copying the values in
def initializeIntensities():
    """
    Creates a blank intensity array of the same size as 
    frequency array, and then plugs all values within 
    frewWindowHALF of the two peaks into th new array
    leaving the rest as 0
    """
    In= np.zeros((np.size(omega)))
    Intensindex = int((minW-simulationMinW)/wres//1)
    In[Intensindex:Intensindex+np.size(IntensReg)-1] = IntensReg[0:-1]
    return In,Intensindex
    
In,Intensindex = initializeIntensities()


## GVD= 0000e-6 %%I've set it to 10000 for Steve it was -1877 
## TOD= 000000e-9 %%I've set it to 1000000 for Steve it was 18719 
## extraphase= GVD/2.*(omega*100*3e8*2*pi).^2 + TOD/6.*(omega*100*3e8*2*pi).^3 
## lambda= 1e4./omega %%this gives the wavelength in microns

##########################
#Set the phase dispersion# 
##########################

ngas= np.zeros((1, np.size(omega)))
k = np.zeros((1, np.size(omega)), dtype = "complex_")
ngas= (ngas + 1)[0]  ##comment what's below to set no dispersion

k = (2 * math.pi) * (omega*100)
# Sellmeier equations for the gas 
# ngas= 1 + Pressure*0.7e30/((31.9e15)**2 - (2.*math.pi*c*100.*omega)**2) ## from physrevA68_023812p15_2003 
# sellmeier equation for fused silica 
# nsilica= sqrt(1.2955 + 0.80985*lambda.^2./(lambda.^2 - .0107945) + 0.91714*lambda.^2./(lambda.^2 - 100)) 
# nsilica= 1.5
# n1= (nsilica**2 + 1)/(2.*math.sqrt(nsilica**2 - 1))
## Axial wavevector for the EH11 hybrid mode 
#k= 2.*math.pi*omega*100.*ngas*(1 - 1/2.*(u11/(omega*100)/(2.*math.pi*fibreRadius))**2)*(1 - 1j*n1/(omega*100)/math.pi/fibreRadius)


###########################
#Electric Field Amplitudes# 
###########################
A= np.zeros(np.size(omega),dtype = "complex")
B= np.zeros(np.size(omega),dtype = "complex")
C= np.zeros(np.size(omega),dtype = "complex")
P= np.zeros((1, Intervals + 1),dtype = "complex")[0]
A= np.sqrt(In)
A = A.astype("complex")
B= np.copy(A)
C= np.copy(A)

## Propagate in z
Q = 0 + 0j
#dz = 10
## COPY OF VARIABLES UP TOP FOR EASY TESTING
FibreLength= 1
Intervals= 3000
dz= FibreLength/Intervals
#g = 1.6e-13
#plotmin = 1e-13
plotmin = 1e-3
g = 0.14
##g= 1.6e-13
#Intervals= 1000
#g = 4e-15
Index1 = Indices[0] + Intensindex
Index2 = Indices[1] + Intensindex

IndexSeperation = Index1 - Index2
Amplitudes = np.zeros(((FibreLength * Intervals)+1,np.size(omega)),dtype = "complex_")
Amplitudes[0] = In
time0 = time.time()
time1 = time0
numTotal = FibreLength * Intervals
print("Start Loop")
for z in range(numTotal): 

    for n in range(int (numPoints) - IndexSeperation):
        Q = Q + B[n + IndexSeperation]*np.conjugate(B[n])*np.exp( -1j *(k[n + IndexSeperation] - k[n])*dz)

    #Down
    for n in range(int (numPoints) - IndexSeperation):
        C[n] = B[n] + g*Intensity/ngas[n]*omega[n]/omega[Index1]*(np.conjugate(Q)*B[n+IndexSeperation]*np.exp(1j*(k[n] - k[n+IndexSeperation])*dz))*dz
    #Up
    for n in range(int (numPoints) - IndexSeperation):
        C[n + IndexSeperation] = C[n + IndexSeperation] + g*Intensity/ngas[n + IndexSeperation]*omega[n + IndexSeperation]/omega[Index1]*(- Q*B[n]*np.exp(1j*(k[n + IndexSeperation] - k[n])*dz))*dz

    #input = ((C*np.conjugate(C))*Intensity)
    np.copyto(Amplitudes[z+1],C)
    np.copyto(B,C)
    #Progress
    if(z%10 == 0):
        time1 = time.time()
        deltaTime = time1 - time0
        timeRemaining = (((numTotal - z + 1)/10)*deltaTime)//1
        print("Running: ", z, " / ",numTotal-1,", ", ((z)/numTotal)*100 , "%, around", timeRemaining, " seconds left...") 
        time0= time.time()
        time1= time.time()
#######
#Plots#
#######
## Now I'm changing the convention - the Pump is of order 0, and the Stokes is of order -1 
#X= range(-(S_Orders + 2) - side/omegaIntervals, (A_Orders + 2) + side/omegaIntervals, 1/omegaIntervals)

SurfaceGraph = Amplitudes*np.conjugate(Amplitudes)
#print(C.imag)
tempindex = (numTotal-1)

##
finalIntens = C*np.conjugate(C)
print(np.max(finalIntens))
finalIntens = finalIntens/np.max(finalIntens)
finalIntens = np.log10((finalIntens.real) + plotmin)
components = (omega - pump)/ramanSep
plt.subplot(211)
plt.xlim(-30, 30)
plt.plot(components,finalIntens)
##

'''
plt.subplot(211)
end = int (pumpindex +((omegaRaman//wres) * 20))
tmp = np.log10((SurfaceGraph.real[tempindex]*np.conjugate(SurfaceGraph.real[tempindex]) + plotmin))
plt.ylim=(0,np.max(tmp))
val = int (pump + ramanSep*numOrders)
print(val)
plt.xlim(0.0, val)
plt.plot(omega,tmp)
'''

eField = fft(Amplitudes[tempindex])
eField = np.fft.ifftshift(eField)

plt.subplot(212)
#[100:200]
plt.plot(omega,eField.real)
plt.show()


for i in range(np.shape(SurfaceGraph)[0]):
    tmp = SurfaceGraph[i]
    tmp = np.log10((tmp*np.conjugate(tmp) + plotmin)*Intensity)
    SurfaceGraph[i] = tmp


###

