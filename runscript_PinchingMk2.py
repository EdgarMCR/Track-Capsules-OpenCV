# -*- coding: utf-8 -*-
"""
Created on Mon Feb 08 15:31:09 2016

@author: mbbxkeh2
"""

from __future__ import absolute_import, division, print_function
#import cv2
import matplotlib.pyplot as plt
plt.close('all')
import numpy as np
#import os
import sys
sys.path.append('C:\\Users\\Edgar\\Dropbox\\PhD\\Python\\OpenCV\\capsule_tracking') #Machine Schuster G.21
sys.path.append('C:\\Users\\mbbxkeh2\\Dropbox\\PhD\\Python\\OpenCV\\capsule_tracking') #Aland Turing 2.105
sys.path.append('/home/magda/Dropbox/PhD/Python/OpenCV/capsule_tracking') #Aland Turing 2.105

import general as gen
import analysis as ana
import semi_sphere as ss
import pinching_mk2 as pin2



def directories_TJ():
    dir_list=[]
    
    #=============================================================================
    # Batch20160517-1 #2 - T-Junction Mode
    #=============================================================================
    directory = 'M:\\EdgarHaener\\Capsules\\Batch20160517-1\\Pinching_MkII\\T-PinMkII_27-5-16\\'
    #=============================================================================
    dir_list.append(directory)
    
    #=============================================================================
    # Batch20160531-4 #8   T-Junction Mode
    #=============================================================================
    directory = 'M:\\EdgarHaener\\Capsules\\Batch20160531-4\\Pinching_MkII\\T-J_PinMkII_02-06-2016\\'
    #=============================================================================
    dir_list.append(directory)

    #=============================================================================
    # Batch20160505-2 #9 17-06-2016 T-Junction
    #=============================================================================
    directory = 'M:\\EdgarHaener\\Capsules\\Batch20160505-2\\Pinching_MkII\\PinMkII_T-J_Mode_17Jun\\'
    #=============================================================================
    dir_list.append(directory)

    #=============================================================================
    # Batch20160622-2 #11 T-Junction Mode
    #=============================================================================
    directory = 'M:\\EdgarHaener\\Capsules\\Batch20160622-2\\Pinching_MkII\\PinMkII_1Jul_T-J_Mode\\'
    #=============================================================================
    dir_list.append(directory)
    
    return dir_list
    
def directories_PIN():
    dir_list=[]

    #=============================================================================
    # Batch160510-6 #4
    #=============================================================================   
    directory = 'M:\\EdgarHaener\\Capsules\\Batch20160510-6\\Pinching_mk2\\PinMk2_16-5-12_2Way\\'
    #=============================================================================
    dir_list.append(directory)
    
    #=============================================================================
    # Batch20160517-1 #2
    #=============================================================================  
    directory = 'M:\\EdgarHaener\\Capsules\\Batch20160517-1\\Pinching_MkII\\PinMkII_27-5-16\\'
    #=============================================================================
    dir_list.append(directory)

    #=============================================================================
    # Batch20160531-4 #8
    #=============================================================================directories_TJ()
    directory = 'M:\\EdgarHaener\\Capsules\\Batch20160531-4\\Pinching_MkII\\PinMkII_02-06-2016\\'
    #=============================================================================
    dir_list.append(directory)

    #=============================================================================
    # Batch20160505-2 #9 17-06-2016
    #=============================================================================
    directory = 'M:\\EdgarHaener\\Capsules\\Batch20160505-2\\Pinching_MkII\\PinMkII_17Jun\\'
    #=============================================================================
    dir_list.append(directory)

    #=============================================================================
    # Batch20160622-3 #6 
    #=============================================================================
    directory = 'M:\\EdgarHaener\\Capsules\\Batch20160622-2\\Pinching_MkII\\Pin_MkII_28Jun\\'
    #=============================================================================
    dir_list.append(directory)

    #=============================================================================
    # Batch20160622-2 #11 
    #=============================================================================    
    directory = 'M:\\EdgarHaener\\Capsules\\Batch20160622-2\\Pinching_MkII\\PinMkII_1Jul_FixedRatio\\'
    #=============================================================================
    dir_list.append(directory)
    
    return dir_list
    
    
def dir_list_TJ2():
    ''' Those used for writing up'''
    dir_list=[]
    #=============================================================================
    # Batch20160505-2 #9
    #=============================================================================
    directory = 'M:\\EdgarHaener\\Capsules\\Batch20160505-2\\Pinching_MkII\\PinMkII_T-J_Mode_17Jun\\'
    dir_list.append(directory)
    #=============================================================================
    
    #=============================================================================
    # Batch20160531-4 #8
    #=============================================================================
    directory = 'M:\\EdgarHaener\\Capsules\\Batch20160531-4\\Pinching_MkII\\T-J_PinMkII_02-06-2016\\'
    dir_list.append(directory)
    #=============================================================================
    
    #=============================================================================
    # Batch20160517-1 #2
    #=============================================================================
    directory = 'M:\\EdgarHaener\\Capsules\\Batch20160517-1\\Pinching_MkII\\T-PinMkII_27-5-16\\'
    dir_list.append(directory)
    #=============================================================================
    
    
    #=============================================================================
    # Batch20160517-1 #2
    #=============================================================================
    directory = 'M:\\EdgarHaener\\Capsules\\Batch20160622-2\\Pinching_MkII\\PinMkII_1Jul_T-J_Mode\\'
    dir_list.append(directory)
    #=============================================================================
    
    return dir_list
    
#=============================================================================
# Batch160510-6 #4
#=============================================================================
FPS=40
 
PPMM=16.6
D0=3.9
top = 244; bottom = 308; end = 338; leftCorner = 258; rightCorner = 320;
 
directory = 'M:\\EdgarHaener\\Capsules\\Batch20160510-6\\Pinching_mk2\\PinMk2_16-5-12_2Way\\'
folder =  'B160510-6_#4_%dFPS_2mlPmin_20mlPmin_2\\' %FPS #\Selection\\
 
#directory = '/mnt/MCND/EdgarHaener/Capsules/Batch20160412-1/Pinching/Pin_15Apr_#3/'
#folder =  'B160412-1_#3_15Apr_Pin_%dFPS_10mlPmin_60mlPmin_5/' %FPS #\Selection\\
 
path = directory + folder #Test2\\'
 
 
backPath = directory   +'Background%dFPS.png' %FPS
#=============================================================================
 
##=============================================================================
## Batch20160517-1 #2
##=============================================================================
#FPS=40
#
#PPMM=16.3
#D0=3.8
#top = 234; bottom = 296; end = 320; leftCorner = 240; rightCorner = 306;
#
#directory = 'M:\\EdgarHaener\\Capsules\\Batch20160517-1\\Pinching_MkII\\PinMkII_27-5-16\\'
#folder =  'B160517-1_#2_PinMk2_27-5_%dFPS_2mlPmin_20mlPmin_2\\' %FPS #\Selection\\
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
# Batch20160517-1 #2 - T-Junction Mode
#=============================================================================
FPS=20
 
PPMM=16.3
D0=3.8
top = 234; bottom = 296; end = 320; leftCorner = 240; rightCorner = 306;
 
directory = 'M:\\EdgarHaener\\Capsules\\Batch20160517-1\\Pinching_MkII\\T-PinMkII_27-5-16\\'
folder =  'B160517-1_#2_T-PinMk2_27-5_%dFPS_-4p9mlPmin_10mlPmin_2\\' %FPS #\Selection\\
 
path = directory + folder #Test2\\'
 
 
backPath = directory   +'Background%dFPS.png' %FPS
#=============================================================================
 
 
##=============================================================================
## Batch20160531-4 #8
##=============================================================================
#FPS=25
#
#PPMM=16.0
#D0=3.9
#top = 244; bottom = 307; end = 316; leftCorner = 236; rightCorner = 302;
#
#directory = 'M:\\EdgarHaener\\Capsules\\Batch20160531-4\\Pinching_MkII\\PinMkII_02-06-2016\\'
#folder =  'B160531-4_#8_2-6-16_PinMkII_%dFPS_2mlPmin_15mlPmin_4\\' %FPS #\Selection\\
#
#path = directory + folder #Test2\\'
#
#
#backPath = directory   +'Background%dFPS.png' %FPS
##=============================================================================
#
#=============================================================================
# Batch20160531-4 #8   T-Junction Mode
#=============================================================================
FPS=20

PPMM=16.0
D0=3.9
top = 245; bottom = 308; end = 315; leftCorner = 233; rightCorner = 299;

directory = 'M:\\EdgarHaener\\Capsules\\Batch20160531-4\\Pinching_MkII\\T-J_PinMkII_02-06-2016\\'
folder =  'B160531-4_#8_2-6-16_T-J_PinMkII_%dFPS_-5mlPmin_20mlPmin_5\\' %FPS #\Selection\\

path = directory + folder #Test2\\'


backPath = directory   +'Background%dFPS_2.png' %FPS
#=============================================================================
#
##=============================================================================
## Batch20160531-4 #8   T-Junction Mode Pressure
##=============================================================================
#FPS=20
#
#PPMM=16.0
#D0=3.9
#top = 245; bottom = 308; end = 315; leftCorner = 233; rightCorner = 299;
#
#directory = 'M:\\EdgarHaener\\Capsules\\Batch20160531-4\\Pinching_MkII\\T-J_Presure_PinMkII_02-06-2016\\'
#folder =  'B160531-4_#8_2-6-16_T-J_P_PinMkII_%dFPS_20mlPmin_1\\' %FPS #\Selection\\
#
#path = directory + folder #Test2\\'
#
#
#backPath = directory   +'Background%dFPS.png' %FPS
##=============================================================================
#
##=============================================================================
## Batch20160505-2 #9 17-06-2016
##=============================================================================
#FPS=4
#
#PPMM=16.25
#D0=3.9
#top = 318; bottom = 383; end = 260; leftCorner = 178; rightCorner = 241;
#
#directory = 'M:\\EdgarHaener\\Capsules\\Batch20160505-2\\Pinching_MkII\\PinMkII_17Jun\\'
#folder =  'B160505-2_#9_17Jun_PinMkII_%dFPS_0p5mlPmin_1mlPmin_5\\' %FPS #\Selection\\
#
#path = directory + folder #Test2\\'
#
#
#backPath = directory   +'Background%dFPS.png' %FPS
##=============================================================================
#
##=============================================================================
## Batch20160505-2 #9 17-06-2016 T-Junction
##=============================================================================
#FPS=5
#
#PPMM=16
#D0=3.9
#top = 318; bottom = 383; end = 178; leftCorner = 98; rightCorner = 162;
#
#directory = 'M:\\EdgarHaener\\Capsules\\Batch20160505-2\\Pinching_MkII\\PinMkII_T-J_Mode_17Jun\\'
#folder =  'B160505-2_#9_17Jun_PinMkII_TJ-Mode_%dFPS_-2p48mlPmin_5mlPmin_1\\' %FPS #\Selection\\
#
#path = directory + folder #Test2\\'
#
#
#backPath = directory   +'Background%dFPS.png' %FPS
##=============================================================================
#
#=============================================================================
# Batch20160622-3 #6 
#=============================================================================
FPS=80

PPMM=16
D0=3.9
top = 234; bottom=298; end=256; leftCorner = 173; rightCorner = 237;

directory = 'M:\\EdgarHaener\\Capsules\\Batch20160622-2\\Pinching_MkII\\Pin_MkII_28Jun\\'
folder =  'B160622-2_#8_PinMkII_28Jun_%dFPS_1mlPmin_2mlPmin_1\\' %FPS #\Selection\\

#directory = '/mnt/MCND/EdgarHaener/Capsules/Batch20160205-1/SemiSphere-2016-02-12-#5/'
#folder = 'Batch20160205-1_20160212_#4_%dFPS_5mlPmin_1/' %FPS 
path = directory + folder #Test2\\'
#=============================================================================
 
 
##=============================================================================
## Batch20160622-2 #11 
##=============================================================================
#PPMM=16
#D0=3.9
#top = 311; bottom=374; end=365; leftCorner = 285; rightCorner = 348;
# 
#directory = 'M:\\EdgarHaener\\Capsules\\Batch20160622-2\\Pinching_MkII\\PinMkII_1Jul_FixedRatio\\'
#FPS=10; folder =  'B160622-2_#11_PinMkII_FR_1Jul_%dFPS_-2p5mlPmin_5mlPmin_9\\' %FPS #\Selection\\
#FPS=210; folder =  'B160622-2_#11_PinMkII_FR_1Jul_%dFPS_-35mlPmin_70mlPmin_6\\' %FPS #\Selection\\
# 
##directory = '/mnt/MCND/EdgarHaener/Capsules/Batch20160205-1/SemiSphere-2016-02-12-#5/'
##folder = 'Batch20160205-1_20160212_#4_%dFPS_5mlPmin_1/' %FPS 
#path = directory + folder #Test2\\'
##=============================================================================
# 
##=============================================================================
## Batch20160622-2 #11 T-Junction Mode
##=============================================================================
#PPMM=16
#D0=3.9
#top = 309; bottom=373; end=184; leftCorner = 104; rightCorner = 168;
# 
#directory = 'M:\\EdgarHaener\\Capsules\\Batch20160622-2\\Pinching_MkII\\PinMkII_1Jul_T-J_Mode\\'
#FPS=10; folder =  'B160622-2_#11_PinMkII_T-JM_1Jul_%dFPS_-2p47mlPmin_5mlPmin_6\\' %FPS #\Selection\\
##FPS=60; folder =  'B160622-2_#11_PinMkII_T-JM_1Jul_%dFPS_-15p08mlPmin_30mlPmin_5\\' %60 #\Selection\\
# 
# 
##directory = '/mnt/MCND/EdgarHaener/Capsules/Batch20160205-1/SemiSphere-2016-02-12-#5/'
##folder = 'Batch20160205-1_20160212_#4_%dFPS_5mlPmin_1/' %FPS 
#path = directory + folder #Test2\\'
##=============================================================================

backPath = directory +'Background%dFPS.png' %FPS


gen.setUsingOpenCV3(True)

"""Collect all information need for Pinching Setup. Assumes inflow from right,
towards expansion on the left.


ChannelTop         Pixel y-position of top of channel
ChannelBottom      Pixel y-position of top of channel   
EndChannel         Pixel x-position of end of channel, start of
                   diffuser
LeftCornerInflow   Pixel x-position of left corner of pinching inflow
RightCornerInflow  Pixel x-position of right corner of pinching inflow
                      
#"""
#PSS = pin2.ParametersPinchingMk2(listingImageFunc=pin2.sortPhotosPCO, coverImgFunc=pin2.coverSidePinching,
#    baseDirectory= directory, folder= folder, d0= D0, pPmm= PPMM, rotation= -0.3, 
#    FPS= FPS, backgroundImage= backPath, ChannelTop = top, ChannelBottom=bottom, 
#    LeftCornerInflow =leftCorner, RightCornerInflow = rightCorner,
#    EndChannel=end, readParametersFromFile=None)
##
###PSS = ss.ParametersSemiSphere(listingImageFunc=ss.sortPhotosPCO, coverImgFunc=ss.coverSideSemiSphere,
###    baseDirectory= None, folder= None, d0= None, pPmm= None, rotation= None, 
###    FPS= None, backgroundImage= None, ChannelTop = None, ChannelBottom=None, 
###    EndChannel=None, readParametersFromFile=path)
#
#PSS.plot = False
#PSS.printDebugInfo= False


#fileList, leng = ss.sortPhotosPCO(path)
#fileList, leng = PSS.listingImageFunc(path)
#print('leng = %d ' %leng)
#testImg = cv2.imread(backPath,0)
#print(ss.coverSideSemiSphere(testImg, PSS))
##
#gen.track_capsules(PSS)
#gen.runOneFPS(PSS)
##gen.wholeRun(PSS)


#
#sp='M:\\EdgarHaener\\Capsules\\PlotScripts\\Pinching\\PinchingMkII-FixedRatio\\'
#DRP = pin2.DataRunPinching(path); 
#DRP.extractMeasures(); 
#DRP._getFinalPosition(TJ_Mode=False, forPub=True, forPubpath=sp);
#DRP._speedVsXpos(savepath=sp)
#DRP._getFinalPosition(TJ_Mode=True, forPub=False)
ana.wholeRun(pin2.runFunction, directory)

#dir_pin = directories_PIN()
#for dire in dir_pin:
#    ana.wholeRun(pin2.runFunction, dire)
#ana.runOneFPS(pin2.runFunction, directory, FPS)

sp='M:\\EdgarHaener\\Capsules\\PlotScripts\\Pinching\\T-Junction_Mode\\'
#DRP = pin2.DataRunPinching(path); DRP.extractMeasures_TJ_Mode()
#DRP.plot=True
#DRP._findGap(TJ_Mode=True)
#DRP._getFinalPosition(TJ_Mode=True, forPub=True, forPubpath=sp);
#DRP._speedVsXpos(sp)

#dir_TJ2 = dir_list_TJ2()
#for dire in dir_TJ2:
#    ana.wholeRun(pin2.runFunction_TJMode, dire)    
#ana.wholeRun(pin2.runFunction_TJMode, directory)
#pin.runFunction(path)

#print(ana._pathExperimentRsltFile(path, isDirectory=False))
#print(ana._pathExperimentRsltFile(directory, isDirectory=True))

#
#RCPin = pin2.ResultsClassPinching(directory); #RCPin.plotRslts()
#RCPin._plotQVsGapSpace()
#RCPin._plotGapD12VsFinalOffset()
#plt.show()

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

#sp = 'M:\\EdgarHaener\\Capsules\\PlotScripts\\Pinching\\' #T-Junction_Mode\\'
##dir_TJ2 = dir_list_TJ2()
##for dire in dir_TJ2:
##    pin2.averageTrajectories(directory=dire, savePath=sp, forPub=False, savename='', TJ=True)
#
##pin2.averageTrajectories(directory=directory, savePath=sp, forPub=True, savename='', TJ=True)
#
#pin2.averageDeformationAndHeight(directory=directory, savePath=sp, forPub=True, savename='', plot_trajectores=True, plot_deformation=True)
#import os
#DRP = pin2.DataRunPinching(os.path.join(directory, 'B160622-2_#11_PinMkII_T-JM_1Jul_80FPS_-20mlPmin_40mlPmin_2'))
#            
##remove entries with no values from centroids
#cent_x, cent_y = DRP._centroidsWithoutErrors()

#cent_y2, cent_x2 = zip(*sorted(zip(cent_y, cent_x)))
#
#cent_y2 = np.array(cent_y2); cent_x2 = np.array(cent_x2)
#
#x, y = ana._remove_outliers(cent_x,cent_y, window=10, cutoff=2.0)
#x, y = ana._remove_outliers(x,y, window=10, cutoff=2.00)
#x2, y2 = ana._remove_outliers(cent_x[20:-20],cent_y[20:-20], window=10, cutoff=2.5)
#x3, y3 = ana._remove_outliers(cent_x,cent_y, window=15, cutoff=1.5)
#
#plt.figure()
#plt.plot(cent_x,cent_y, '-ro')
#plt.plot(cent_x2,cent_y2, '-go')
#plt.plot(x, y, 'sb')
#plt.plot(cent_x[20:-20],cent_y[20:-20], 'hg')
#plt.plot(x3, y3, 'py')
#
#plt.figure()
#plt.plot(cent_x,cent_y, 'ro')
#plt.plot(cent_x2,cent_y2, '-go')
#plt.plot(x+10, y+10, 'sb')
#plt.plot(cent_x[20:-20]+20,cent_y[20:-20]+20, 'hg')
#plt.plot(x2+30, y2+30, '<c')
#plt.plot(x3+40, y3+40, 'py')
#
#plt.figure()
#plt.plot(np.gradient(cent_x), np.gradient(cent_y), 'ro')


#x2, y2 = ana._1d_remove_outliers(cent_x,cent_y, window=30, cutoff=1.5, check=True)
#y2,x2 = ana._1d_remove_outliers(y2,x2, window=30, cutoff=1.5, check=True)
#
#x, y = ana._1d_remove_outliers(cent_x,cent_y, window=30, cutoff=1.8, check=True)
#y,x = ana._1d_remove_outliers(y,x, window=30, cutoff=1.8, check=True)
#
#
#x3, y3 = ana._1d_remove_outliers(cent_x,cent_y, window=20, cutoff=1.8, check=True)
#y3,x3 = ana._1d_remove_outliers(y3,x3, window=20, cutoff=1.8, check=True)
#
#x10, y10 = ana._1d_remove_outliers(cent_x,cent_y, window=10, cutoff=1.8, check=True)
#y10,x10 = ana._1d_remove_outliers(y10,x10, window=10, cutoff=1.8, check=True)
#
#plt.figure()
#plt.plot(cent_x,cent_y, 'ro')
#plt.plot(x+5,y+5, 'bs')
#plt.plot(x2+10,y2+10, 'yp')
#plt.plot(x2+15,y2+15, 'c<')
#plt.plot(x10,y10, 'hg')

#plt.figure()
#plt.plot(np.arange(len(cent_x)), cent_x, 'ro')

#dir_PIN= directories_PIN()
##dir_TJ= directories_TJ()
##dir_TJ2 = dir_list_TJ2()
##
#for dire in dir_PIN:
#    plt.close('all')
#    ana.wholeRun(pin2.runFunction, dire)
##
#for dire in dir_TJ2:
#    plt.close('all')
#    ana.wholeRun(pin2.runFunction_TJMode, dire)    
#    
#plt.close('all')
##for dire in dir_PIN:
##    RCPin = pin2.ResultsClassPinching(dire)
##    RCPin._plotQVsGapSpace()  
##    RCPin._plotGapD12VsFinalOffset()
#for dire in dir_TJ2:
#    RCPin = pin2.ResultsClassPinching(dire)
##    RCPin._plotQVsGapSpace()  
#    RCPin._plotGapD12VsFinalOffset()




