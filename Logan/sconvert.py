import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
import math

def readBasic(fName):
    dataFile = open(fName,'r')
    DataValues = np.loadtxt(dataFile, dtype=float,comments='#')
    lambdaReg = np.array(DataValues[:,0])
    IntensReg = np.array(DataValues[:,1])
    omegaReg= 1e7/lambdaReg; 

    RegenSpec = np.zeros(shape=(np.size(lambdaReg),2)).transpose()

    RegenSpec[0] = omegaReg
    RegenSpec[1] = IntensReg
    
    return RegenSpec

def defineBasic(pumpF, ramanSeperation, intensity):
    lambdaReg = np.arange(0, pumpF+ramanSeperation+100)
    IntensReg = np.zeros(np.size(lambdaReg))

    IntensReg[pumpF] = intensity
    IntensReg[pumpF - ramanSeperation] = intensity
    
    RegenSpec = np.zeros(shape=(np.size(lambdaReg),2)).transpose()
    
    RegenSpec[0] = lambdaReg
    RegenSpec[1] = IntensReg

    return RegenSpec

def defineBasicTWO(pumpF, ramanSeperation,numOrder,numPointsBetweenPeaks,intensity):
    ptsLeft = int ((math.ceil(pumpF/ramanSeperation)*numPointsBetweenPeaks)+math.ceil(pumpF/(ramanSeperation)))
    if(pumpF - ramanSeperation < 0 ):
        raise ValueError('Pump - Seperation resulated in a negative raman Frequency')
        
    ptsRight = int ((numOrder * (numPointsBetweenPeaks + 1)))
    ptsTotal = ptsLeft + ptsRight + 1
    finalF = pumpF + (numOrder * ramanSeperation)
    
    indexPump = ptsLeft
    indexRaman = ptsLeft - (numPointsBetweenPeaks+1)

    lambdaReg = np.linspace(0, finalF, ptsTotal)
    IntensReg = np.zeros(np.size(lambdaReg))
    IntensReg[indexPump] = intensity
    IntensReg[indexRaman] = intensity

    RegenSpec = np.zeros(shape=(np.size(lambdaReg),2)).transpose()
    
    RegenSpec[0] = lambdaReg
    RegenSpec[1] = IntensReg

    return RegenSpec

def readAndNormalize(fName):
    dataFile = open(fName,'r')
    DataValues = np.loadtxt(dataFile, dtype=float,comments='#')
    lambdaReg = np.matrix(DataValues[:,0])
    IntensReg = np.matrix(DataValues[:,1])

    flatLambda = np.array(lambdaReg)[0, :]
    flatIntens= np.array(IntensReg)[0, :]

    l= flatLambda
    Il= flatIntens

    w= 1e7/l
    Iw=1e7*Il/w**2

    wres= np.min(np.abs(w - np.roll(w, 1)))
    D= np.max(w) - np.min(w)
    Nfloat= np.ceil(D/wres)
    N= Nfloat.astype(int) 
   
    we= np.linspace(np.min(w), np.max(w), N)
    lin= interpolate.interp1d(w, Iw, kind='linear')
    Ie= lin(we)


    dw= D/Nfloat
    dl= np.min(np.abs(l - np.roll(l,1)))

    Ite= np.trapz(Ie, dx=dw)
    Itl= np.trapz(Il, dx=dl)
    
    RegenSpec = np.zeros(shape=(np.size(we),2)).transpose()


    RegenSpec[0] = we
    RegenSpec[1] = Ie
    
    return RegenSpec

def readAndNormalizeOverridePoints(fName,numPoints):
    dataFile = open(fName,'r')
    DataValues = np.loadtxt(dataFile, dtype=float,comments='#')
    lambdaReg = np.matrix(DataValues[:,0])
    IntensReg = np.matrix(DataValues[:,1])

    flatLambda = np.array(lambdaReg)[0, :]
    flatIntens= np.array(IntensReg)[0, :]

    l= flatLambda
    Il= flatIntens

    w= 1e7/l
    Iw=1e7*Il/w**2

    wres= np.min(np.abs(w - np.roll(w, 1)))
    D= np.max(w) - np.min(w)
    Nfloat= np.ceil(D/wres)
    N= numPoints

    we= np.linspace(np.min(w), np.max(w), N)
    lin= interpolate.interp1d(w, Iw, kind='linear')
    Ie= lin(we)


    dw= D/Nfloat
    dl= np.min(np.abs(l - np.roll(l,1)))

    Ite= np.trapz(Ie, dx=dw)
    Itl= np.trapz(Il, dx=dl)
    
    RegenSpec = np.zeros(shape=(np.size(we),2)).transpose()


    RegenSpec[0] = we
    RegenSpec[1] = Ie
    
    return RegenSpec

def manualInterpolateRegenSpec(regenSpec, numPoints):
    regenSpec[0] = 1e7/regenSpec[0]
    DataValues = regenSpec.transpose()

    lambdaReg = np.matrix(DataValues[:,0])
    IntensReg = np.matrix(DataValues[:,1])

    flatLambda = np.array(lambdaReg)[0, :]
    flatIntens= np.array(IntensReg)[0, :]

    l= flatLambda
    Il= flatIntens

    w= 1e7/l
    Iw=1e7*Il/w**2

    wres = np.min(np.abs(w - np.roll(w, 1)))
    D = np.max(w) - np.min(w)
    Nfloat= np.ceil(D/wres)
    N= numPoints 

    we= np.linspace(np.min(w), np.max(w), N)
    lin= interpolate.interp1d(w, Iw, kind='linear')
    Ie= lin(we)
    
    RegenSpec = np.zeros(shape=(np.size(we),2)).transpose()
    RegenSpec[0] = we
    RegenSpec[1] = Ie

    return RegenSpec

def readAndNormalizeAndDebug(fName):
    dataFile = open(fName,'r')
    DataValues = np.loadtxt(dataFile, dtype=float,comments='#')
    lambdaReg = np.matrix(DataValues[:,0])
    IntensReg = np.matrix(DataValues[:,1])

    flatLambda = np.array(lambdaReg)[0, :]
    flatIntens= np.array(IntensReg)[0, :]

    l= flatLambda
    Il= flatIntens

    w= 1e7/l
    Iw=1e7*Il/w**2

    wres= np.min(np.abs(w - np.roll(w, 1)))
    D= np.max(w) - np.min(w)
    Nfloat= np.ceil(D/wres)
    N= Nfloat.astype(int) 

    we= np.linspace(np.min(w), np.max(w), N)
    lin= interpolate.interp1d(w, Iw, kind='linear')
    Ie= lin(we)


    dw= D/Nfloat
    dl= np.min(np.abs(l - np.roll(l,1)))

    Ite= np.trapz(Ie, dx=dw)
    Itl= np.trapz(Il, dx=dl)

    plt.subplot(211)
    plt.plot(we,Ie)
    plt.subplot(212)
    plt.plot(l,Il)

    #print("Approximately", np.round(100/Itl*(Ite - Itl),1), "% Error between Wavenumber and Wavelength Integrals")
    
    
    RegenSpec = np.zeros(shape=(np.size(we),2)).transpose()
    #print(RegenSpec)

    RegenSpec[0] = we
    RegenSpec[1] = Ie
    
    return RegenSpec

#manualInterpolate(readAndNormalize("C:\\Users\\loogo\\Desktop\\FemtoSecond\\AND064.dat"),100)
#results = defineBasic(12000,800,1)

#pumpF, ramanSeperation,numOrder,numPointsBetweenPeaks,intensity):
#results = defineBasicTWO(12000,775,20,1,7)
#plt.plot(results[0],results[1])