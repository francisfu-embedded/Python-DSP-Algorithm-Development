import mysignals as sigs
import sigproc as sp
from matplotlib import style as st
from matplotlib import pyplot as plt
import numpy as np
import signal
import math


#arr1 = []
#arr2 = []
#arr3 = []

#arr1,arr2, arr3= a.getDFT()

#a.plotDFT()
#signal = a.processWithCombFilter(1000,0.5, 48000)

#a.plotFilteredSig()


windSinc1 = sp.getWindowedSinc('Hamming',0.2,len(sigs.InputSignal_1kHz_15kHz))
windSinc2 = sp.getWindowedSinc('Blackman',0.2,len(sigs.InputSignal_1kHz_15kHz))
windSinc3 = sp.getWindowedSinc('Rectangular',0.2,len(sigs.InputSignal_1kHz_15kHz))
a=sp.signalHandle('sig1',sigs.InputSignal_1kHz_15kHz)
c = a.processWithFIRSincLpFilter(13000, 48000)
a.plotFilteredSig()


