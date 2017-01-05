# -*- coding: utf-8 -*-
"""
Created on Mon Feb 08 15:31:09 2016

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

##=============================================================================
## GelBeads20160324-1 #1 04-04-2016
##=============================================================================
#FPS=120
#
#PPMM=15.8
#D0=3.8
#top = 148; bottom = 210; end = 570; leftCorner = 604; rightCorner = 665;
#
#directory = 'M:\\EdgarHaener\\Capsules\\GelBead20160324-1\\Pinching-Setup\\GB160324-1_#1_4Apr_Pinching\\'
#folder =  'GelBeads20160324-1_#1_PS_4Apr_%dFPS_10mlPmin_60mlPmin_1\\' %FPS #\Selection\\
#
#path = directory + folder #Test2\\'
#
#
#backPath = directory   +'Background%dFPS.png' %FPS
##=============================================================================

##=============================================================================
## Batch20160223-3 #7 08-04-2016
##=============================================================================
#FPS=20
#
#PPMM=16.0
#D0=3.9
#top = 148; bottom = 210; end = 590; leftCorner = 627; rightCorner = 688;
#
#directory = 'M:\\EdgarHaener\\Capsules\\Batch20160223-3\\Pinching\\#7_16-04-08\\'
#folder =  'Batch20160223-3_#7_8Apr_Pin_%dFPS_10mlPmin_0mlPmin_2\\' %FPS #\Selection\\
#
##directory = '/mnt/MCND/EdgarHaener/Capsules/Batch20160223-3/Pinching/#7_16-04-08/'
##folder =  'Batch20160223-3_#7_8Apr_Pin_%dFPS_10mlPmin_25mlPmin_6/' %FPS #\Selection\\
#
#path = directory + folder #Test2\\'
#
#
#backPath = directory   +'Background%dFPS.png' %FPS
##backPath = directory   + 'Background80FPS_10mlPmin_1.png'
##backPath = directory   + 'Background60FPS_5mlPmin_13.png'
##backPath = directory   + 'Background60FPS_5mlPmin_15.png'
##backPath = directory   + 'Background60FPS_5mlPmin_16.png'
##=============================================================================

##=============================================================================
## Batch160412-1 #3
##=============================================================================
#FPS=220
#
#PPMM=16.0
#D0=3.9
#top = 138; bottom = 202; end = 579; leftCorner = 612; rightCorner = 676;
#
#directory = 'M:\\EdgarHaener\\Capsules\\Batch20160412-1\\Pinching\\Pin_15Apr_#3\\'
#folder =  'B160412-1_#3_15Apr_Pin_%dFPS_10mlPmin_60mlPmin_5\\' %FPS #\Selection\\
#
##directory = '/mnt/MCND/EdgarHaener/Capsules/Batch20160412-1/Pinching/Pin_15Apr_#3/'
##folder =  'B160412-1_#3_15Apr_Pin_%dFPS_10mlPmin_60mlPmin_5/' %FPS #\Selection\\
#
#path = directory + folder #Test2\\'
#
#
#backPath = directory   +'Background%dFPS.png' %FPS
##=============================================================================


#=============================================================================
# Batch160531-4 #9  07-06-2016
#=============================================================================
FPS=5

PPMM=16.8
D0=3.9
top = 150; bottom = 214; end = 592; leftCorner = 626; rightCorner = 690;

directory = 'C:\\Users\\Edgar\\Desktop\\Experiments\\2016-06-07\\Pinching_7Jun\\'
folder =  'B160531-4_#9_Pin_T-JMode_7Jun_%dFPS_-2p5mlPmin_5mlPmin_3\\' %FPS #\Selection\\

#directory = '/mnt/MCND/EdgarHaener/Capsules/Batch20160412-1/Pinching/Pin_15Apr_#3/'
#folder =  'B160412-1_#3_15Apr_Pin_%dFPS_10mlPmin_60mlPmin_5/' %FPS #\Selection\\

path = directory + folder #Test2\\'


backPath = directory   +'Background%dFPS.png' %FPS
#=============================================================================

gen.setUsingOpenCV3(True)

"""Collect all information need for Pinching Setup. Assumes inflow from right,
towards expansion on the left.


ChannelTop         Pixel y-position of top of channel
ChannelBottom      Pixel y-position of top of channel   
EndChannel         Pixel x-position of end of channel, start of
                   diffuser
LeftCornerInflow   Pixel x-position of left corner of pinching inflow
RightCornerInflow  Pixel x-position of right corner of pinching inflow
                      
"""
PSS = pin.ParametersPinching(listingImageFunc=pin.sortPhotosPCO, coverImgFunc=pin.coverSidePinching,
    baseDirectory= directory, folder= folder, d0= D0, pPmm= PPMM, rotation= -0.3, 
    FPS= FPS, backgroundImage= backPath, ChannelTop = top, ChannelBottom=bottom, 
    LeftCornerInflow =leftCorner, RightCornerInflow = rightCorner,
    EndChannel=end, readParametersFromFile=None)

#PSS = ss.ParametersSemiSphere(listingImageFunc=ss.sortPhotosPCO, coverImgFunc=ss.coverSideSemiSphere,
#    baseDirectory= None, folder= None, d0= None, pPmm= None, rotation= None, 
#    FPS= None, backgroundImage= None, ChannelTop = None, ChannelBottom=None, 
#    EndChannel=None, readParametersFromFile=path)

PSS.plot = False
PSS.printDebugInfo=False


#fileList, leng = ss.sortPhotosPCO(path)
#fileList, leng = PSS.listingImageFunc(path)
#print('leng = %d ' %leng)
#testImg = cv2.imread(backPath,0)
#print(ss.coverSideSemiSphere(testImg, PSS))
##
gen.track_capsules(PSS)
#gen.runOneFPS(PSS)
#gen.wholeRun(PSS)


#
DRP = pin.DataRunPinching(path); DRP.extractMeasures()
#ana.wholeRun(pin.runFunction, directory)
#ana.runOneFPS(pin.runFunction, directory, FPS)

#pin.runFunction(path)

#print(ana._pathExperimentRsltFile(path, isDirectory=False))
#print(ana._pathExperimentRsltFile(directory, isDirectory=True))

#RCPin = pin.ResultsClassPinching(directory); RCPin.plotRslts()
#print(RCPin.batchName)

#for ii in range(9,23):
#    folder = 'Batch20160306-2_#6_080316_%dFPS_5mlPmin_%d\\' %(FPS, ii)
#    path = directory + folder
#    PSS = ss.ParametersSemiSphere(listingImageFunc=ss.sortPhotosPCO, coverImgFunc=ss.coverSideSemiSphere,
#                                  baseDirectory= directory, folder= folder, d0= D0, pPmm= PPMM, rotation= -0.3, 
#                                  FPS= FPS, backgroundImage= backPath, ChannelTop = top, ChannelBottom=bottom, 
#                                  EndChannel=end, readParametersFromFile=None)
#    
#    ss.runFunction(path)
    
#for ii in range(11, 21):
#    PSS.folder = 'Batch20160306-2_#6_080316_%dFPS_30mlPmin_%d\\' %(FPS, ii)
#    PSS.updatePath();gen.track_capsules(PSS)









