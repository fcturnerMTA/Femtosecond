import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage.filters import maximum_filter
import sys
import PeakCompare
import sconvert

def compareBasic(omegaRaman,RegenSpec):
    ##
    ##Find the row with the highest intensity, this will just be the maximum
    ## intensity value and gives peak 1
    highesIntensity_ROW =  np.where(RegenSpec[1] == np.amax(RegenSpec[1]))
    RegenSpec = RegenSpec.transpose()

    ## peak 2 is to the left of peak 1, find half raman and look left of that
    peakOmega = RegenSpec[highesIntensity_ROW][0][0]
    splitOmega = int ( peakOmega//1 - (omegaRaman//2))

    RegenSpec = RegenSpec.transpose()
    difference = np.abs(RegenSpec - splitOmega)
    index = np.argmin(difference)
    RegenSpec = RegenSpec.transpose()

    #print(RegenSpec[index])

    lowerHalf = RegenSpec[index :,:]

    ## same as above but for left half
    lowerHalf = lowerHalf.transpose()
    highesIntensity_ROW2 =  np.where(lowerHalf[1] == np.amax(lowerHalf[1]))
    lowerHalf = lowerHalf.transpose()

    indices = [highesIntensity_ROW[0][0],highesIntensity_ROW2[0][0]]

    return(RegenSpec[highesIntensity_ROW],lowerHalf[highesIntensity_ROW2], indices)

def wip(omegaRaman,RegenSpec):
    local_max = maximum_filter(RegenSpec[1], omegaRaman)
    plt.plot(RegenSpec[0],local_max,RegenSpec[0],RegenSpec[1])
    plt.show()

def compareNormalized(omegaRaman,RegenSpec):
        ##
    ##Find the row with the highest intensity, this will just be the maximum
    ## intensity value and gives peak 1
    highesIntensity_ROW =  np.where(RegenSpec[1] == np.amax(RegenSpec[1]))
    if(np.shape(highesIntensity_ROW)[1] == 2):
        indices = [highesIntensity_ROW[0][1],highesIntensity_ROW[0][0]]
        RegenSpec = RegenSpec.transpose()
        
        return(RegenSpec[highesIntensity_ROW[0][0]],RegenSpec[highesIntensity_ROW[0][1]], indices) 
    else:
        ## peak 2 is to the left of peak 1, find half raman and look left of that
        RegenSpec = RegenSpec.transpose()
        peakOmega = RegenSpec[highesIntensity_ROW[0][0]][0]
        splitOmega = int ( peakOmega//1 - (omegaRaman//2))
        RegenSpec = RegenSpec.transpose()
        difference = np.abs(RegenSpec - splitOmega)
        index = np.argmin(difference)
        RegenSpec = RegenSpec.transpose()

        #print(RegenSpec[index])

        lowerHalf = RegenSpec[0:index,:]
        ## same as above but for left half
        lowerHalf = lowerHalf.transpose()
        highesIntensity_ROW2 =  np.where(lowerHalf[1] == np.amax(lowerHalf[1]))
        lowerHalf = lowerHalf.transpose()

        RegenSpec = RegenSpec.transpose()

        indices = [highesIntensity_ROW[0][0],highesIntensity_ROW2[0][0]]
        RegenSpec = RegenSpec.transpose()

        return(RegenSpec[highesIntensity_ROW],lowerHalf[highesIntensity_ROW2], indices) 

def compareNormalizedAndPlot(omegaRaman,RegenSpec):
    ##
    ##Find the row with the highest intensity, this will just be the maximum
    ## intensity value and gives peak 1
    highesIntensity_ROW =  np.where(RegenSpec[1] == np.amax(RegenSpec[1]))
    if(np.shape(highesIntensity_ROW)[1] == 2):
        indices = [highesIntensity_ROW[0][0],highesIntensity_ROW[0][1]]
        RegenSpec = RegenSpec.transpose()
        plt.plot(RegenSpec[highesIntensity_ROW[0][0]][0],RegenSpec[highesIntensity_ROW[0][0]][1], 'ro')
        plt.plot(RegenSpec[highesIntensity_ROW[0][1]][0],RegenSpec[highesIntensity_ROW[0][1]][1], 'ro')
        
        RegenSpec = RegenSpec.transpose()
        plt.plot(RegenSpec[0],RegenSpec[1])
        plt.show()
        RegenSpec = RegenSpec.transpose()
        return(RegenSpec[highesIntensity_ROW[0][0]],RegenSpec[highesIntensity_ROW[0][1]], indices) 
    else:
        ## peak 2 is to the left of peak 1, find half raman and look left of that
        RegenSpec = RegenSpec.transpose()
        peakOmega = RegenSpec[highesIntensity_ROW[0][0]][0]
        splitOmega = int ( peakOmega//1 - (omegaRaman//2))
        RegenSpec = RegenSpec.transpose()
        difference = np.abs(RegenSpec - splitOmega)
        index = np.argmin(difference)
        RegenSpec = RegenSpec.transpose()

        #print(RegenSpec[index])

        lowerHalf = RegenSpec[0:index,:]
        ## same as above but for left half
        lowerHalf = lowerHalf.transpose()
        highesIntensity_ROW2 =  np.where(lowerHalf[1] == np.amax(lowerHalf[1]))
        lowerHalf = lowerHalf.transpose()

        plt.plot(RegenSpec[highesIntensity_ROW][0][0],RegenSpec[highesIntensity_ROW][0][1], 'ro')
        plt.plot(lowerHalf[highesIntensity_ROW2][0][0],lowerHalf[highesIntensity_ROW2][0][1], 'ro')
        
        RegenSpec = RegenSpec.transpose()
        plt.plot(RegenSpec[0],RegenSpec[1])
        plt.show()

        indices = [highesIntensity_ROW[0][0],highesIntensity_ROW2[0][0]]
        RegenSpec = RegenSpec.transpose()

        return(RegenSpec[highesIntensity_ROW],lowerHalf[highesIntensity_ROW2], indices) 

#RegenSpec = sconvert.readAndNormalize("C:\\Users\\loogo\\Desktop\\FemtoSecond\\AND064.dat")
#RegenSpec = sconvert.defineBasicTWO(12000,775,20,1,7)
#print(compareNormalizedAndPlot(775,RegenSpec))
