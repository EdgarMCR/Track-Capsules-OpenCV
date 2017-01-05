# -*- coding: utf-8 -*-
"""
Multiple runs in pinching setup

Created on Sun Apr 10 15:36:21 2016

@author: mbbxkeh2
"""
from __future__ import absolute_import, division, print_function
#import cv2
import matplotlib.pyplot as plt
plt.close('all')
#import os
import sys
sys.path.append('C:\\Users\\Edgar\\Dropbox\\PhD\\Python\\OpenCV\\capsule_tracking') #Machine Schuster G.21
sys.path.append('C:\\Users\\mbbxkeh2\\Dropbox\\PhD\\Python\\OpenCV\\capsule_tracking') #Aland Turing 2.105
sys.path.append('/home/magda/Dropbox/PhD/Python/OpenCV/capsule_tracking') #Aland Turing 2.105

import general as gen
import analysis as ana
import semi_sphere as ss
import pinching as pin
import several_runs_pinching as srp

##=============================================================================
## GelBeads20160324-1 #1 04-04-2016
##=============================================================================
directory = 'M:\\EdgarHaener\\Capsules\\GelBead20160324-1\\Pinching-Setup\\GB160324-1_#1_4Apr_Pinching\\'

RCPin_GelBeads20160324d1_nr1_04042016 = pin.ResultsClassPinching(directory)
RCPin_GelBeads20160324d1_nr1_04042016._findAverages()
#=============================================================================

#=============================================================================
# Batch20160223-3 #7 08-04-2016
#=============================================================================
directory = 'M:\\EdgarHaener\\Capsules\\Batch20160223-3\\Pinching\\#7_16-04-08\\'

RCPin_Batch20160223d3_nr7_08042016 = pin.ResultsClassPinching(directory)
RCPin_Batch20160223d3_nr7_08042016._findAverages()
#=============================================================================

#=============================================================================
# Batch20160412-1 #3 15-04-2016
#=============================================================================
directory = 'M:\\EdgarHaener\\Capsules\\Batch20160412-1\\Pinching\\Pin_15Apr_#3\\'

RCPin_Batch20160412d1_nr3_15042016 = pin.ResultsClassPinching(directory)
RCPin_Batch20160412d1_nr3_15042016._findAverages()
#=============================================================================



#=============================================================================
# Collection into Lists
#=============================================================================

instances = [RCPin_GelBeads20160324d1_nr1_04042016, RCPin_Batch20160223d3_nr7_08042016, RCPin_Batch20160412d1_nr3_15042016]
savePath = 'M:\\EdgarHaener\\Capsules\\PlotScripts\\Pinching\\'
nameList=['GelBead20160324-1 #1', 'Batch20160223-3 #7', 'Batch20160412-1 #3']

#instances = [RCPin_Batch20160223d3_nr7_08042016, RCPin_Batch20160412d1_nr3_15042016]
#savePath = 'M:\\EdgarHaener\\Capsules\\PlotScripts\\Pinching\\'
#nameList=['Batch20160223-3 #7', 'Batch20160412-1 #3']

#=============================================================================
# Plotting
#=============================================================================

srp.plotQVsFinalOffset(instances, savePath, forPub=False, names=nameList)
srp.plotQVsFinalOffsetAve(instances, savePath, forPub=False, names=nameList)
srp.plotQVsGapSpaceing(instances, savePath, forPub=False, names=nameList)
