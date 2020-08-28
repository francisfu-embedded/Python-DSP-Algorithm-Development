#import modules
from matplotlib import pyplot as plt
from scipy import signal
from matplotlib import style as st
import numpy as np
import math


#################################signal handle class##########################
class signalHandle:
    
    def __init__(self, name, signalArray):
        self.__name = name
        self.__signalArray = signalArray
        self.__convoSig = []
        self.__rsumSig = []
        self.__fdiffSig = []
        self.__outsideSig = []
        self.__dftSigRe = []
        self.__dftSigIm = []
        self.__dftSigMag = []
        self.__filteredSig = []

    #Signal Processing Methods
    
    def getSignal(self):
        return self.__signalArray
    
    def getMean(self):
        mean = 0;
        
        #accumulate the values in sigalArray
        for x in self.__signalArray:
            mean += x
            
        #mean calculation
        mean = mean/len(self.__signalArray)
        return mean
    
    def getName(self):
        return self.__name
    
    def getPlot(self):
        plt.plot(self.__signalArray)
        plt.show()
        
    def getVariance(self):
        #get mean first
        variance = 0.0
        mean = self.getMean()
        
        #variance = summing (x-mean)^2/item_number
        for x in self.__signalArray:
            variance = variance + (x-mean)**2
        variance /= len(self.__signalArray)
        return variance
    
    def getStDeviation(self):
        #standart deviation = sqrt(variance)
        variance = self.getVariance()
        return variance**(0.5)
    
    def getConvolutionWith (self, outsideSig):
        self.__outsideSig = outsideSig.getSignal()
        #self.__convoSig = signal.convolve(self.__outsideSig,self.__signalArray)
        for x in range(len(self.__signalArray) + len(self.__outsideSig)):
            self.__convoSig.append(0)
        
        for x in range(len(self.__signalArray)):
            for y in range(len(self.__outsideSig)):
                self.__convoSig[x+y] += (self.__signalArray[x]*self.__outsideSig[y])

        return self.__convoSig

    def getRunningSum (self):
        for x in range(len(self.__signalArray)):
            self.__rsumSig.append(0)

        for x in range(len(self.__signalArray)):
            self.__rsumSig[x] = self.__rsumSig[x-1] + self.__signalArray[x]
        
        return self.__rsumSig

    def getFirstDifference(self):
        for x in range(len(self.__signalArray)):
            self.__fdiffSig.append(0)

        for x in range(len(self.__signalArray)):
            self.__fdiffSig[x] = self.__signalArray[x] - self.__signalArray[x-1]

        return self.__fdiffSig
    
    #get real DFT
    def getDFT(self):
        length = len(self.__signalArray)
        halfLength = int(length / 2)
        for i in range (halfLength):
            self.__dftSigRe.append(0)
            self.__dftSigIm.append(0)
            self.__dftSigMag.append(0)

        for j in range (halfLength):
            for k in range (length):
                 self.__dftSigRe[j] += self.__signalArray[k] * math.cos(2*math.pi/length*k*j)
                 self.__dftSigIm[j] -= self.__signalArray[k] * math.sin(2*math.pi/length*k*j)
            self.__dftSigRe[j] = self.__dftSigRe[j] * 2 / length
            self.__dftSigIm[j] = self.__dftSigIm[j] * 2 / length
            
        for index in range (halfLength):
            self.__dftSigMag[index] = math.sqrt((self.__dftSigRe[index])**2 + (self.__dftSigIm[index])**2)
           
        return self.__dftSigRe, self.__dftSigIm, self.__dftSigMag
        

    #freq is the frequency to be rejected in the frequency spectrum,
    #a is the constant to time with for the poles at frequency freq
    def processWithNotchFilter(self, freq, a, sampRate):
        
        sigLen = len(self.__signalArray)
        
        for i in range (sigLen):
            self.__filteredSig.append(0)

        #discrete time angular frequency omega = 2*pi*f*Ts
        omega  = 2*math.pi*freq/sampRate

        for j in range (sigLen):
            self.__filteredSig[j] = (
                                    self.__signalArray[j]-
                                    2*math.cos(omega)*self.__signalArray[j-1]+
                                    self.__signalArray[j-2]-
                                    a*a*self.__filteredSig[j-2]+
                                    2*a*math.cos(omega)*self.__filteredSig[j-1]
                                    )
        
        return self.__filteredSig

    # a comb filter has zeros at 2*pi*k/N0, damping the periodic interference
    #freq is the fundamental frequency of frequencies to be rejected 
    #a is the constant to time with for the poles at frequency freq
    #which has a range between 0 and 1
    def processWithCombFilter(self, freq, a, sampRate):
        sigLen = len(self.__signalArray)

        N0 = int(sampRate/freq)
        
        for i in range (sigLen):
            self.__filteredSig.append(0)

        for j in range (sigLen):
            self.__filteredSig[j] = (self.__signalArray[j]
                                     - self.__signalArray[j-N0]
                                     + (a**N0) * self.__filteredSig[j - N0])
        return self.__filteredSig
                 
    #plotting methods
    def plotConvolution(self):
        f, pltArr = plt.subplots(3,sharex = True)
        f.suptitle('Two Signals Convolved')
        pltArr[0].plot(self.__signalArray , color = 'blue')
        pltArr[0].set_title('Object Signal')
        pltArr[1].plot(self.__outsideSig , color = 'red')
        pltArr[1].set_title('Outside Signal')
        pltArr[2].plot(self.__convoSig, color = 'green')
        pltArr[2].set_title('Convolved Signal')
        plt.show()

    def plotRunningSum(self):
        f, pltArr = plt.subplots(2,sharex = True)
        f.suptitle('Object Signal and Its Running Sum')
        pltArr[0].plot(self.__signalArray , color = 'blue')
        pltArr[0].set_title('Object Signal')
        pltArr[1].plot(self.__rsumSig , color = 'red')
        pltArr[1].set_title('Running Sum')
        plt.show()

    def plotFirstDifference(self):
        f, pltArr = plt.subplots(2,sharex = True)
        f.suptitle('Object Signal and Its First Difference')
        pltArr[0].plot(self.__signalArray , color = 'blue')
        pltArr[0].set_title('Object Signal')
        pltArr[1].plot(self.__fdiffSig , color = 'red')
        pltArr[1].set_title('First Difference')
        plt.show()

    def plotDFT(self):
        f, pltArr = plt.subplots(4,sharex = True)
        f.suptitle('Object Signal and Its Running Sum')
        pltArr[0].plot(self.__signalArray , color = 'blue')
        pltArr[0].set_title('Object Signal')
        pltArr[1].plot(self.__dftSigRe , color = 'red')
        pltArr[1].set_title('DFT Real Domain')
        pltArr[2].plot(self.__dftSigIm , color = 'red')
        pltArr[2].set_title('DFT Imaginary Domain')
        pltArr[3].plot(self.__dftSigMag, color = 'red')
        pltArr[3].set_title('DFT Magnitude')
        plt.show()

    def plotFilteredSig(self):
        f, pltArr = plt.subplots(2,sharex = True)
        f.suptitle('Object Signal and Its Filtered Signal')
        pltArr[0].plot(self.__signalArray , color = 'blue')
        pltArr[0].set_title('Object Signal')
        pltArr[1].plot(self.__filteredSig, color = 'red')
        pltArr[1].set_title('Filtered Signal')
        plt.show()

##############################costumary functions################################

def getIDFT (realSig, imgSig):
    
    if(len(realSig) != len(imgSig)):
        print('length of both signal sequence not identical')

    reconSig = []
    kLength = len(realSig)
    nLength = 2 * kLength
    
    for i in range(nLength):
        reconSig.append(0)

    for j in range(nLength):
        for k in range(kLength):
            reconSig[j] = (reconSig[j] +
                          realSig[k] * math.cos(2 * math.pi * k * j / nLength)+
                          imgSig[k] * math.sin(2 * math.pi * k * j / nLength))
    return reconSig
            


    
