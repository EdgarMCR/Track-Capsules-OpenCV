# -*- coding: utf-8 -*-
"""
Created on Mon Feb 08 15:31:09 2016

@author: mbbxkeh2
"""

from __future__ import absolute_import, division, print_function
#import cv2
import matplotlib.pyplot as plt
plt.close('all')
import os
import sys
sys.path.append('C:\\Users\\Edgar\\Dropbox\\PhD\\Python\\OpenCV\\capsule_tracking') #Machine Schuster G.21
sys.path.append('C:\\Users\\mbbxkeh2\\Dropbox\\PhD\\Python\\OpenCV\\capsule_tracking') #Aland Turing 2.105
sys.path.append('/home/magda/Dropbox/PhD/Python/OpenCV/capsule_tracking') #Aland Turing 2.105

import general as gen
import analysis as ana
import semi_sphere as ss

##=============================================================================
## Batch20160205-1 #4
##=============================================================================
#FPS=5
#
#PPMM=15.7
#D0=3.9
#top = 120; bottom=375; end=552
#
#directory = 'M:\EdgarHaener\Capsules\Batch20160205-1\SemiSphere-2016-02-06\\'
#folder =  'Batch20160205-1-#4-5mlPmin_%dFPS_6\\' %FPS #\Selection\\
#path = directory + folder #Test2\\'
#
#
#backPath = directory +'Background%dFPS.png' %FPS
##=============================================================================


#=============================================================================
# Batch20160205-1 #5
#=============================================================================
FPS=50

PPMM=15.7
D0=3.9
top = 127; bottom=375; end=561

directory = 'M:\\EdgarHaener\\Capsules\\Batch20160205-1\\SemiSphere-2016-02-12-#5\\'
folder =  'Batch20160205-1_20160212_#4_%dFPS_50mlPmin_7\\' %FPS #\Selection\\
#directory = '/mnt/MCND/EdgarHaener/Capsules/Batch20160205-1/SemiSphere-2016-02-12-#5/'
#folder = 'Batch20160205-1_20160212_#4_%dFPS_5mlPmin_1/' %FPS 
path = directory + folder #Test2\\'


backPath = directory +'Background%dFPS.png' %FPS
#=============================================================================

##=============================================================================
## Batch20160205-1 #5 2016-02-15
##=============================================================================
#FPS=20
#
#PPMM=15.7
#D0=3.9
#top = 127; bottom=377; end=550
#
#directory = 'M:\EdgarHaener\Capsules\Batch20160205-1\Batch20160215-1_#5_Run20160215\\'
#folder =  'Batch20160205-1_#5_%dFPS_20mlPmin_1\\' %FPS #\Selection\\
##directory = '/mnt/MCND/EdgarHaener/Capsules/Batch20160205-1/Batch20160215-1_#5_Run20160215/'
##folder = 'Batch20160205-1_#5_%dFPS_20mlPmin_1/' %FPS 
#
#path = directory + folder #Test2\\'
#
#
#backPath = directory +'Background%dFPS.png' %FPS
##=============================================================================


##=============================================================================
## Batch20160223-3 #7
##=============================================================================
#FPS=5
#
#PPMM=15.63
#D0=3.9
#top = 127; bottom=375; end=551
#
##directory = os.path.join('/mnt', 'MCND', 'EdgarHaener', 'Capsules', 'Batch20160223-3', 'Semi-Sphere', 'SemiSphere_26022016_#7')
##folder =  'Batch20160223-3_#7_60mlPmin_%dFPS_1' %FPS #\Selection\\
#
#directory = 'M:\\EdgarHaener\\Capsules\\Batch20160223-3\\Semi-Sphere\\SemiSphere_26022016_#7\\'
#folder = 'Batch20160223-3_#7_5mlPmin_%dFPS_4\\' %FPS 
##path = os.path.join(directory, folder) #Test2\\'
#
#
#backPath = os.path.join(directory, 'Background%dFPS.png' %FPS)
##=============================================================================


##=============================================================================
## Batch20160223-3 #4
##=============================================================================
#FPS=25
#
#PPMM=15.7
#D0=3.9
#top = 128; bottom=376; end=568
#
#directory = 'M:\\EdgarHaener\\Capsules\\Batch20160223-3\\Semi-Sphere\\SemiSphere_28022016_#4\\'
#folder =  'Batch20160223-3_#4_15mlPmin_%dFPS_6\\' %FPS #\Selection\\
##directory = '/mnt/MCND/EdgarHaener/Capsules/Batch20160205-1/SemiSphere-2016-02-12-#5/'
##folder = 'Batch20160205-1_20160212_#4_%dFPS_5mlPmin_1/' %FPS 
#path = directory + folder #Test2\\'
#
#
#backPath = directory +'Background%dFPS.png' %FPS
##=============================================================================

##=============================================================================
## Batch20160223-3 #4 29-02-2016
##=============================================================================
#FPS=25
#
#PPMM=15.7
#D0=3.9
#top = 127; bottom=376; end=567
#
#directory = 'M:\\EdgarHaener\\Capsules\\Batch20160223-3\\Semi-Sphere\\Semi-Sphere_29022016_#4\\'
#folder =  'Batch20160223-3_#4_290216_25mlPmin_%dFPS_4\\' %FPS #\Selection\\
##directory = '/mnt/MCND/EdgarHaener/Capsules/Batch20160205-1/SemiSphere-2016-02-12-#5/'
##folder = 'Batch20160205-1_20160212_#4_%dFPS_5mlPmin_1/' %FPS 
#path = directory + folder #Test2\\'
#
#
#backPath = directory +'Background%dFPS.png' %FPS
##=============================================================================

##=============================================================================
## Batch20160306-2 #6 08-03-2016
##=============================================================================
#FPS=5
#
#PPMM=15.9
#D0=3.9
##top = 128; bottom=378; end=501
##top = 128; bottom=380; end=501 #30 ml/min Runs 1-8
##top = 125; bottom=375; end=500 #30 ml/min Runs 9-20
##top = 130; bottom=380; end=501 #20 ml/min Runs 1-8
##top = 125; bottom=375; end=498 #20 ml/min Runs 9-36
##top = 129; bottom=380; end=502 #10 ml/min Runs 1-8
##top = 126; bottom=375; end=500 #10 ml/min Runs 9-20
##top = 128; bottom=380; end=501 #5 ml/min Runs 1-8
##top = 123; bottom=374; end=500 #5 ml/min Runs 9-22
#
##directory = 'C:\\Users\\Edgar\\Desktop\\Experiments\\2016-03-08\\SemiSphere_08032016\\'
##folder =  'Batch20160306-2_#6_080316_%dFPS_30mlPmin_20\\' %FPS #\Selection\\
#directory = 'M:\\EdgarHaener\\Capsules\\Batch20160306-2\\SemiSphere\\SemiSphere_08032016\\'
#
###for ii in range(12):
##folder = 'Batch20160306-2_#6_080316_%dFPS_5mlPmin_%d\\' %(FPS, 1 )
##path = directory + folder #Test2\\'
##backPath = directory   +'Background%dFPS_2.png' %FPS
##ss.runFunction(path)
#    
##=============================================================================

##=============================================================================
## Batch20160306-2 #6 12-03-2016
##=============================================================================
#FPS=2
#
#PPMM=15.9
#D0=3.9
#top = 128; bottom=378; end=501
#
#directory = 'M:\\EdgarHaener\\Capsules\\Batch20160306-2\\SemiSphere\\SS_2016-03-12\\'
#folder =  'Batch20160306-2_#6_120316_%dFPS_2mlPmin_1\\' %FPS #\Selection\\
#
#path = directory + folder #Test2\\'
#
#
#backPath = directory   +'Background%dFPS.png' %FPS
##=============================================================================

##=============================================================================
## GelBeads20160324-3 #3 29-03-2016
##=============================================================================
#FPS=30
#
#PPMM=15.8
#D0=3.9
#top = 135; bottom=385; end=524
#
#directory = 'M:\\EdgarHaener\\Capsules\\GelBead20160324-1\\Semi-Sphere\\29032016\\'
#folder =  'GelBead20160324-1_#3_%dFPs_30mlPmin_1\\' %FPS #\Selection\\
#
#path = directory + folder #Test2\\'
#
#
#backPath = directory   +'Background%dFPS.png' %FPS
##=============================================================================


#=============================================================================
# Batch20160505-2 #16
#=============================================================================
FPS=50

PPMM=15.8
D0=3.9
top = 118; bottom=368; end=580

directory = 'M:\\EdgarHaener\\Capsules\\Batch20160505-2\\SemiSphere\\SS_15June\\'
folder =  'Batch20160505-2_#16_15Jun_%dFPS_50mlPmin_6\\' %FPS #\Selection\\
#directory = '/mnt/MCND/EdgarHaener/Capsules/Batch20160205-1/SemiSphere-2016-02-12-#5/'
#folder = 'Batch20160205-1_20160212_#4_%dFPS_5mlPmin_1/' %FPS 
path = directory + folder #Test2\\'

backPath = directory +'Background%dFPS.png' %FPS
#=============================================================================

##=============================================================================
##Batch20160531-4 #9 05-06-2016
##=============================================================================
##FPS=
#
#PPMM=15.8
#D0=3.9
#top = 136; bottom= 388; end= 587  # 5ml.min
#top = 129 ; bottom= 379; end= 572 #other flow speeds
#directory = 'M:\\EdgarHaener\\Capsules\\Batch20160531-4\\SemiSphere\\SS_5-6-16\\'
##folder =  
##=============================================================================


##=============================================================================
## Batch20160622-3 #6 
##=============================================================================
#FPS=24
#
#PPMM=15.75
#D0=3.9
#top = 127; bottom=377; end=566
#
#directory ='M:\\EdgarHaener\\Capsules\\Batch20160622-2\\Semi-Sphere\\SS_24Jun_#6\\'
#folder =  'B160622-3_#6_SS_24Jun_%dFPS_20mlPmin_8\\' %FPS #\Selection\\
#
##directory = '/mnt/MCND/EdgarHaener/Capsules/Batch20160205-1/SemiSphere-2016-02-12-#5/'
##folder = 'Batch20160205-1_20160212_#4_%dFPS_5mlPmin_1/' %FPS 
#path = directory + folder #Test2\\'
#
#
#backPath = directory +'Background%dFPS.png' %FPS
##=============================================================================

##=============================================================================
## Batch20160622-3 #69
##=============================================================================
#FPS=70
#
#PPMM=15.75
#D0=3.9
#top = 128; bottom=380; end= 568
#
#directory = 'M:\\EdgarHaener\\Capsules\\Batch20160622-2\\Semi-Sphere\\SS_24Jun_#9\\'
#folder =  'B160622-3_#6_SS_24Jun_%dFPS_60mlPmin_3\\' %FPS #\Selection\\
#
#path = directory + folder #Test2\\'
#
#
#backPath = directory +'Background%dFPS_2.png' %FPS
##=============================================================================

##=============================================================================
## Batch20160610-3 #1 
##=============================================================================
#FPS=60
#
#PPMM=15.6
#D0=3.0
#top = 119; bottom=369; end=528
#
#directory = 'M:\\EdgarHaener\\Capsules\\Batch20160610-3\\Semi-Sphere\\SS_30Jun\\'
#folder =  'B160610-3_#1_SS_30Jun_%dFPS_30mlPmin_4\\' %FPS #\Selection\\
#
#path = directory + folder #Test2\\'
#
#
#backPath = directory +'Background%dFPS_2.png' %FPS
##=============================================================================

##=============================================================================
## Batch20160610-5 #7 
##=============================================================================
#FPS=20
#
#PPMM=15.9
#D0=3.57
#top = 129; bottom=380; end=548
#
##directory = 'C:\\Users\\Edgar\\Desktop\\Experiments\\2016-07-15\\SS_15Jul\\'
#directory = 'M:\\EdgarHaener\\Capsules\\Batch20160610-5\\Semi-Sphere\\SS_15Jul\\'
#folder =  'B160622-5_#7_SS_%dFPS_10mlPmin_2\\' %FPS #\Selection\\
#
#path = directory + folder #Test2\\'
#
#
#backPath = directory +'Background%dFPS.png' %FPS
##=============================================================================

#gen.setUsingOpenCV3(True)

#gen.setENHANCE_CONTRAST(True)

#inputs arguments
#listingImageFunc,coverImgFunc, baseDirectory= None, folder= None, d0= None, pPmm= None, 
#rotation= None, FPS= None, backgroundImage= None, ChannelTop = None, ChannelBottom=None, 
#EndChannel=None, readParametersFromFile=None)

PSS = ss.ParametersSemiSphere(listingImageFunc=ss.sortPhotosPCO, coverImgFunc=ss.coverSideSemiSphere,
    baseDirectory= directory, folder= folder, d0= D0, pPmm= PPMM, rotation= -0.3, 
    FPS= FPS, backgroundImage= backPath, ChannelTop = top, ChannelBottom=bottom, 
    EndChannel=end, readParametersFromFile=None)

##PSS = ss.ParametersSemiSphere(listingImageFunc=ss.sortPhotosPCO, coverImgFunc=ss.coverSideSemiSphere,
##    baseDirectory= None, folder= None, d0= None, pPmm= None, rotation= None, 
##    FPS= None, backgroundImage= None, ChannelTop = None, ChannelBottom=None, 
##    EndChannel=None, readParametersFromFile=path)
#
PSS.plot = False
PSS.printDebugInfo=False
#
#
##fileList, leng = ss.sortPhotosPCO(path)
##fileList, leng = PSS.listingImageFunc(path)
##print('leng = %d ' %leng)
##testImg = cv2.imread(backPath,0)
##print(ss.coverSideSemiSphere(testImg, PSS))
##

#gen.track_capsules(PSS)
sf= os.path.join(PSS.path, 'coloured')
gen.coverImages(PSS, savefolder=sf, color=(0,255,255))










