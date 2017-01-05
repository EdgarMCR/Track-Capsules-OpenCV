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
import pinching_mk2 as pin2
import several_runs_pinching as srp2


#=============================================================================
# Batch160510-6 #4
#=============================================================================
directory = 'M:\\EdgarHaener\\Capsules\\Batch20160510-6\\Pinching_mk2\\PinMk2_16-5-12_2Way\\'

RCPin_B160510d6_nr4_120516 = pin2.ResultsClassPinching(directory)
RCPin_B160510d6_nr4_120516._findAverages()
#=============================================================================

#=============================================================================
# Batch20160517-1 #2
#=============================================================================
directory = 'M:\\EdgarHaener\\Capsules\\Batch20160517-1\\Pinching_MkII\\PinMkII_27-5-16\\'

RCPin_B160517d1_nr2_270516 = pin2.ResultsClassPinching(directory)
RCPin_B160517d1_nr2_270516._findAverages()
#=============================================================================
#=============================================================================
# Batch20160531-4 #8
#=============================================================================
directory = 'M:\\EdgarHaener\\Capsules\\Batch20160531-4\\Pinching_MkII\\PinMkII_02-06-2016\\'
RCPin_B160531d4_nr8_020616 = pin2.ResultsClassPinching(directory)
RCPin_B160531d4_nr8_020616._findAverages()
#=============================================================================



#=============================================================================
# Collection into Lists
#=============================================================================

instances = [RCPin_B160510d6_nr4_120516, RCPin_B160517d1_nr2_270516, RCPin_B160531d4_nr8_020616]
savePath = 'M:\\EdgarHaener\\Capsules\\PlotScripts\\Pinching\\'
nameList=['Batch20160510-6 #4 d=3.9mm', 'Batch20160517-1 #2 d=3.8mm', 'Batch20160531-4 #8 d=3.9mm']
force= [3.8, 5.8, 10.01]

#instances = [RCPin_Batch20160223d3_nr7_08042016, RCPin_Batch20160412d1_nr3_15042016]
#savePath = 'M:\\EdgarHaener\\Capsules\\PlotScripts\\Pinching\\'
#nameList=['Batch20160223-3 #7', 'Batch20160412-1 #3']

svn='PinchingMode-FixedBaseFlow'

leng = 0
lengPinFB = 0
for inst in instances:
    leng += len(inst.finalDistCentre)    
    lengPinFB +=len(inst.finalDistCentre)
    
print("lengPinFB = %d" %(lengPinFB))
#=============================================================================
# Plotting
#=============================================================================

#srp2.plotQVsFinalOffset(instances, savePath, forPub=False, names=nameList, savename=svn )
#srp2.plotQVsFinalOffsetAve(instances, savePath, forPub=True, names=nameList, savename=svn)
#srp2.plotQVsFinalOffsetAve(instances, savePath, forPub=False, names=nameList, savename=svn)
#
#srp2.plotCaVsFinalOffsetAve(instances, force, savePath, forPub=True, names=nameList, savename=svn)
#
#srp2. plotAbsFinalVsInitialOffset(instances, savePath, forPub=False)
#print("Initial Centring:")
#for inst in instances:
#    print(inst.aveOffset)
    
##srp2.plotQVsGapAndD12(instances, savePath, forPub=False, names=nameList, savename=svn, TJ=False)
#srp2.plotGapD12VsFinalOffsetAve(instances, savePath, forPub=True, names=nameList, savename=svn)
#srp2.plotCaVsGapAndQvsD12(instances, force, savePath, forPub=True, names=nameList, savename=svn)
#srp2.plotQVsGapAndQvsD12(instances, savePath, forPub=True, names=nameList, savename=svn)
#    
#srp2.plotQVsGapAndQvsD12(instances, savePath, forPub=False, names=nameList, savename=svn, TJ=False, top=True)
#srp2.plotQVsGapAndQvsD12(instances, savePath, forPub=True, names=nameList, savename=svn, TJ=False, top=True)
#=============================================================================
#=============================================================================
# Pinching - Fixed Ratio Mode
#=============================================================================
#=============================================================================

#=============================================================================
# Batch20160505-2 #9
#=============================================================================
directory = 'M:\\EdgarHaener\\Capsules\\Batch20160505-2\\Pinching_MkII\\PinMkII_17Jun\\'
RCPin_B160505d2_nr9_170616 = pin2.ResultsClassPinching(directory)
RCPin_B160505d2_nr9_170616._findAverages()
#=============================================================================

#=============================================================================
# Batch20160622-2 #8
#=============================================================================
directory = 'M:\\EdgarHaener\\Capsules\\Batch20160622-2\\Pinching_MkII\\Pin_MkII_28Jun\\'
RCPin_B160622d2_nr8_280616 = pin2.ResultsClassPinching(directory)
RCPin_B160622d2_nr8_280616._findAverages()
#=============================================================================

#=============================================================================
# Batch20160622-2 #8
#=============================================================================
directory = 'M:\\EdgarHaener\\Capsules\\Batch20160622-2\\Pinching_MkII\\PinMkII_1Jul_FixedRatio\\'
RCPin_B160622d2_nr11_010716 = pin2.ResultsClassPinching(directory)
RCPin_B160622d2_nr11_010716._findAverages()
#=============================================================================


#=============================================================================
# Collection into Lists - Fixed Ratio Mode
#=============================================================================

instances = [RCPin_B160505d2_nr9_170616, RCPin_B160622d2_nr8_280616, RCPin_B160622d2_nr11_010716]
savePath = 'M:\\EdgarHaener\\Capsules\\PlotScripts\\Pinching\\PinchingMkII-FixedRatio\\'
nameList=['Batch20160505-2 #9', 'Batch20160622-2 #8', 'Batch20160622-2 #11']
force= [1.4, 2.3, 14.9]

#instances = [RCPin_Batch20160223d3_nr7_08042016, RCPin_Batch20160412d1_nr3_15042016]
#savePath = 'M:\\EdgarHaener\\Capsules\\PlotScripts\\Pinching\\'
#nameList=['Batch20160223-3 #7', 'Batch20160412-1 #3']

svn='PinchingMode-FixedRatio'

lengPinFR=0
for inst in instances:
    leng += len(inst.finalDistCentre)
    lengPinFR +=len(inst.finalDistCentre)
    
print("lengPinFR = %d" %(lengPinFR))
#=============================================================================
# Plotting
#=============================================================================

#srp2.plotQVsFinalOffset(instances, savePath, forPub=False, names=nameList, savename=svn )
#srp2.plotQVsFinalOffsetAve(instances, savePath, forPub=True, names=nameList, savename=svn)
#srp2.plotQVsFinalOffsetAve(instances, savePath, forPub=False, names=nameList, savename=svn)
#
#srp2.plotCaVsFinalOffsetAve(instances, force, savePath, forPub=True, names=nameList, savename=svn)
#
#srp2.plotQVsGapAndD12(instances, savePath, forPub=False, names=nameList, savename=svn, TJ=False)
#srp2.plotGapD12VsFinalOffsetAve(instances, savePath, forPub=True, names=nameList, savename=svn, pnb=3)
#srp2.plotQVsGapAndQvsD12(instances, savePath, forPub=False, names=nameList, savename=svn)
#srp2.plotQVsGapAndQvsD12(instances, savePath, forPub=False, names=nameList, savename=svn, TJ=False, top=True)
#srp2.plotQVsGapAndQvsD12(instances, savePath, forPub=True, names=nameList, savename=svn, TJ=False, top=True)
#=============================================================================
#=============================================================================
# T-Junction Mode
#=============================================================================
#=============================================================================

#=============================================================================
# Batch20160505-2 #9
#=============================================================================
directory = 'M:\\EdgarHaener\\Capsules\\Batch20160505-2\\Pinching_MkII\\PinMkII_T-J_Mode_17Jun\\'
RCPin_B160505d2_nr9_170616 = pin2.ResultsClassPinching(directory)
RCPin_B160505d2_nr9_170616._findAverages()
#=============================================================================

#=============================================================================
# Batch20160531-4 #8
#=============================================================================
directory = 'M:\\EdgarHaener\\Capsules\\Batch20160531-4\\Pinching_MkII\\T-J_PinMkII_02-06-2016\\'
RCPin_B160531d4_nr8_020616 = pin2.ResultsClassPinching(directory)
RCPin_B160531d4_nr8_020616._findAverages()
#=============================================================================

#=============================================================================
# Batch20160517-1 #2
#=============================================================================
directory = 'M:\\EdgarHaener\\Capsules\\Batch20160517-1\\Pinching_MkII\\T-PinMkII_27-5-16\\'
RCPin_B160517d1_nr2_270516 = pin2.ResultsClassPinching(directory)
RCPin_B160517d1_nr2_270516._findAverages()
#=============================================================================


#=============================================================================
# Batch20160517-1 #2
#=============================================================================
directory = 'M:\\EdgarHaener\\Capsules\\Batch20160622-2\\Pinching_MkII\\PinMkII_1Jul_T-J_Mode\\'
RCPin_B160622d2_nr11_2010716 = pin2.ResultsClassPinching(directory)
RCPin_B160622d2_nr11_2010716._findAverages()
#=============================================================================



#=============================================================================
# Collection into Lists
#=============================================================================

#instances = [RCPin_B160505d2_nr9_170616, RCPin_B160531d4_nr8_020616, RCPin_B160517d1_nr2_270516, RCPin_B160622d2_nr11_2010716]
#savePath = 'M:\\EdgarHaener\\Capsules\\PlotScripts\\Pinching\\T-Junction_Mode\\'
#nameList=['Batch20160505-2 #9  ~1.3 mN', 'Batch20160531-4 #8 ~10 mN', 'Batch20160517-1 #2 ~5.8 mN', 'Batch201600622-2 #11 ~14.9 mN']
#force= [1.35, 10.0, 5.8, 14.9]
instances = [RCPin_B160505d2_nr9_170616, RCPin_B160531d4_nr8_020616, RCPin_B160622d2_nr11_2010716]
savePath = 'M:\\EdgarHaener\\Capsules\\PlotScripts\\Pinching\\T-Junction_Mode\\'
nameList=['Batch20160505-2 #9  ~1.3 mN', 'Batch20160531-4 #8 ~10 mN',  'Batch201600622-2 #11 ~14.9 mN']
force= [1.35, 10.0, 14.9]

#=============================================================================
# Plotting
#=============================================================================
svn = 'T-Junction_Mode'

print("instances[0].aveGap_d12: ", end='')
print(instances[0].aveGap_d12)


for inst in instances:
    leng += len(inst.finalDistCentre)
print("Total runs = %d" %(leng))
##srp2.plotQVsFinalOffset(instances, savePath, forPub=False, names=nameList, savename=svn, TJ=True)
#srp2.plotQVsFinalOffsetAve(instances, savePath, forPub=True, names=nameList, savename=svn, TJ=True)
##srp2.plotQVsFinalOffsetAve(instances, savePath, forPub=False, names=nameList, savename=svn, TJ=True)
##srp2.plotQVsFinalOffsetAve(instances= [RCPin_B160517d1_nr2_270516], savePath=savePath, forPub=False, names=['Batch20160517-1 #2 ~5.8 mN'], savename=svn, TJ=True)
#srp2.plotCaVsFinalOffsetAve(instances, force, savePath, forPub=True, names=nameList, savename=svn, TJ=True)
#
#srp2.fit_exponentiol_LfoVsCa(instances, force, savePath, forPub=True, names=nameList, savename=svn, TJ=True)

#srp2.plotQVsGapSpaceing(instances, savePath, forPub=False, names=nameList, TJ=True)
#srp2.plotGapD12VsFinalOffsetAve(instances, savePath, forPub=False, names=nameList, savename=svn)

#srp2.plotQVsGapAndQvsD12(instances, savePath, forPub=False, names=nameList, savename=svn, TJ=True, top=True)
#srp2.plotQVsGapAndQvsD12(instances, savePath, forPub=False, names=nameList, savename=svn, TJ=True, top=True, force=force)
#srp2.plotQVsGapAndQvsD12(instances, savePath, forPub=True, names=nameList, savename=svn, TJ=True, top=True)
#srp2.plotQVsGapAndQvsD12(instances, savePath, forPub=True, names=nameList, savename=svn, TJ=True, top=True, force=force)

srp2. plotTopGapVsGapD12Ave(instances, savePath, forPub=True, names=nameList, savename=svn)

