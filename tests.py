import mysignals as sigs
import sigproc as pro
from matplotlib import style as st
from matplotlib import pyplot as plt
a=pro.signalHandle('sig1',sigs.InputSignal_1kHz_15kHz)

#arr1 = []
#arr2 = []
#arr3 = []

#arr1,arr2, arr3= a.getDFT()

#a.plotDFT()
signal = a.processWithCombFilter(1000,0.5, 48000)

a.plotFilteredSig()



