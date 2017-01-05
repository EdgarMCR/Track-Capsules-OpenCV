# -*- coding: utf-8 -*-
"""
Created on Sun Apr 10 16:09:16 2016

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
import several_runs_semi_sphere as srss

#=============================================================================
# Batch20160205-1 #5
#=============================================================================
directory = 'M:\EdgarHaener\Capsules\Batch20160205-1\SemiSphere-2016-02-12-#5\\'

RCSS_Batch20160205d1_nr5_12022016 = ss.ResultsClassSemiSphere(directory)
RCSS_Batch20160205d1_nr5_12022016._findAverages()
#=============================================================================

#=============================================================================
# Batch20160205-1 #5 2016-02-15
#=============================================================================
directory = 'M:\EdgarHaener\Capsules\Batch20160205-1\Batch20160215-1_#5_Run20160215\\'

RCSS_Batch20160205d1_nr4_15022016 = ss.ResultsClassSemiSphere(directory)
RCSS_Batch20160205d1_nr4_15022016._findAverages()
##=============================================================================


#=============================================================================
# Batch20160223-3 #4
#=============================================================================
directory = 'M:\\EdgarHaener\\Capsules\\Batch20160223-3\\Semi-Sphere\\SemiSphere_28022016_#4\\'

RCSS_Batch20160223d3_nr4_28022016 = ss.ResultsClassSemiSphere(directory)
RCSS_Batch20160223d3_nr4_28022016._findAverages()
#=============================================================================

#=============================================================================
# Batch20160223-3 #4 29-02-2016
#=============================================================================
directory = 'M:\\EdgarHaener\\Capsules\\Batch20160223-3\\Semi-Sphere\\Semi-Sphere_29022016_#4\\'

RCSS_Batch20160223d3_nr4_29022016 = ss.ResultsClassSemiSphere(directory)
RCSS_Batch20160223d3_nr4_29022016._findAverages()
#=============================================================================

#=============================================================================
# Batch20160223-3 #7
#=============================================================================
directory = 'M:\\EdgarHaener\\Capsules\\Batch20160223-3\\Semi-Sphere\\SemiSphere_26022016_#7\\'

RCSS_Batch20160223d3_nr7_26022016 = ss.ResultsClassSemiSphere(directory)
RCSS_Batch20160223d3_nr7_26022016._findAverages()
#=============================================================================

#=============================================================================
# Batch20160306-2 #6 08-03-2016
#=============================================================================
directory = 'M:\\EdgarHaener\\Capsules\\Batch20160306-2\\SemiSphere\\SemiSphere_08032016\\'

RCSS_Batch20160306d2_nr6_08032016 = ss.ResultsClassSemiSphere(directory)
RCSS_Batch20160306d2_nr6_08032016._findAverages()
#=============================================================================

#=============================================================================
# Batch20160306-2 #6 08-03-2016
#=============================================================================
directory = 'M:\\EdgarHaener\\Capsules\\Batch20160306-2\\SemiSphere\\SS_2016-03-12\\'

RCSS_Batch20160306d2_nr6_12032016 = ss.ResultsClassSemiSphere(directory)
RCSS_Batch20160306d2_nr6_12032016._findAverages()

#=============================================================================

#=============================================================================
# GelBeads20160324-3 #3 29-03-2016
#=============================================================================
directory = 'M:\\EdgarHaener\\Capsules\\GelBead20160324-1\\Semi-Sphere\\29032016\\'

RCSS_GelBeads20160324d1_nr1_29032016 = ss.ResultsClassSemiSphere(directory)
RCSS_GelBeads20160324d1_nr1_29032016._findAverages()
#=============================================================================


#=============================================================================
#Batch20160531-4 #9 05-06-2016
#=============================================================================
directory = 'M:\\EdgarHaener\\Capsules\\Batch20160531-4\\SemiSphere\\SS_5-6-16\\'

RCSS_B20160531d4_nr9_5616 = ss.ResultsClassSemiSphere(directory)
RCSS_B20160531d4_nr9_5616._findAverages()
#=============================================================================


#=============================================================================
#Batch20160505-2 #16 15-06-2016
#=============================================================================
directory = 'M:\\EdgarHaener\\Capsules\\Batch20160505-2\\SemiSphere\\SS_15June\\'

RCSS_B20160505d2_nr16_15616 = ss.ResultsClassSemiSphere(directory)
RCSS_B20160505d2_nr16_15616._findAverages()
#=============================================================================

#=============================================================================
#Batch20160622-3 #9 24-06-2016
#=============================================================================
directory = 'M:\\EdgarHaener\\Capsules\\Batch20160622-2\\Semi-Sphere\\SS_24Jun_#9\\'

RCSS_B20160622d2_nr9_24616 = ss.ResultsClassSemiSphere(directory)
RCSS_B20160622d2_nr9_24616._findAverages()
#=============================================================================

#=============================================================================
#Batch20160622-3 #6 24-06-2016
#=============================================================================
directory ='M:\\EdgarHaener\\Capsules\\Batch20160622-2\\Semi-Sphere\\SS_24Jun_#6\\'
RCSS_B20160622d2_nr6_24616 = ss.ResultsClassSemiSphere(directory)
RCSS_B20160622d2_nr6_24616._findAverages()
#=============================================================================


#=============================================================================
#Batch20160610-3 #1 30-06-2016
#=============================================================================
directory = 'M:\\EdgarHaener\\Capsules\\Batch20160610-3\\Semi-Sphere\\SS_30Jun\\'
RCSS_B20160610d3_nr1_300616 = ss.ResultsClassSemiSphere(directory)
RCSS_B20160610d3_nr1_300616._findAverages()
#=============================================================================


#=============================================================================
#Batch20160610-5 #7 15-07-2016
#=============================================================================
directory = 'M:\\EdgarHaener\\Capsules\\Batch20160610-5\\Semi-Sphere\\SS_15Jul\\'
RCSS_B20160610d5_nr7_150716 = ss.ResultsClassSemiSphere(directory)
RCSS_B20160610d5_nr7_150716._findAverages()
#=============================================================================




#=============================================================================
# Plotting
#=============================================================================
#instances = [RCSS_Batch20160205d1_nr5_12022016, RCSS_Batch20160205d1_nr4_15022016, RCSS_Batch20160223d3_nr4_28022016, RCSS_Batch20160223d3_nr4_29022016, RCSS_Batch20160306d2_nr6_08032016, RCSS_Batch20160306d2_nr6_12032016, RCSS_GelBeads20160324d1_nr1_29032016]
#savePath = 'M:\\EdgarHaener\\Capsules\\PlotScripts\\Semi-Sphere\\'
#nameList = ['Batch20160205-1 #5 12/02/16', 'Batch20160205-1 #5 15/02/16','Batch20160223-3 #4 28/02/2016', 'Batch20160223-3 #4 29/02/2016','Batch20160306-2 #6 08/03/16', 'Batch20160306-2 #6 12/03/16', 'GelBeads20160324-1 #1 29/03/16']
#
#
#srss.plotQVsFinalOffset(instances, savePath, forPub=False, names=nameList)
#srss._plotQVsAngle(instances, savePath, forPub=False, names=nameList)
#srss.plotQVsFinalOffsetAve(instances, savePath, forPub=False, names=nameList)
#srss.plotQVsAngleAve(instances, savePath, forPub=False, names=nameList)
#srss.plotFinalVsInitialOffset(instances, savePath, forPub=False, names=[])
#
#
#instances = [RCSS_Batch20160205d1_nr5_12022016, RCSS_Batch20160205d1_nr4_15022016]
#savePath = 'M:\\EdgarHaener\\Capsules\\PlotScripts\\Semi-Sphere\\'
#nameList = ['Batch20160205-1 #5 12/02/16', 'Batch20160205-1 #5 15/02/16']
#snm='Batch20160205d1_nr4'
#
#srss.plotQVsFinalOffset(instances, savePath, forPub=False, names=nameList, savename=snm)
#srss._plotQVsAngle(instances, savePath, forPub=False, names=nameList, savename=snm)
#srss.plotQVsFinalOffsetAve(instances, savePath, forPub=False, names=nameList, savename=snm)
#srss.plotQVsAngleAve(instances, savePath, forPub=False, names=nameList, savename=snm)
#srss.plotFinalVsInitialOffset(instances, savePath, forPub=False, names=[], savename=snm)

#instances = [ RCSS_Batch20160223d3_nr4_28022016, RCSS_Batch20160223d3_nr4_29022016, RCSS_Batch20160306d2_nr6_08032016, RCSS_Batch20160306d2_nr6_12032016, RCSS_GelBeads20160324d1_nr1_29032016]
#savePath = 'M:\\EdgarHaener\\Capsules\\PlotScripts\\Semi-Sphere\\'
#nameList = ['Batch20160223-3 #4 28/02/2016', 'Batch20160223-3 #4 29/02/2016','Batch20160306-2 #6 08/03/16', 'Batch20160306-2 #6 12/03/16', 'GelBeads20160324-1 #1 29/03/16']
#snm='B0223_B0306_GB0324'
#
##srss.plotQVsFinalOffset(instances, savePath, forPub=False, names=nameList, savename=snm)
#srss._plotQVsAngle(instances, savePath, forPub=False, names=nameList, savename=snm)
##srss.plotQVsFinalOffsetAve(instances, savePath, forPub=False, names=nameList, savename=snm)
#srss.plotQVsAngleAve(instances, savePath, forPub=False, names=nameList, savename=snm)
##srss.plotFinalVsInitialOffset(instances, savePath, forPub=Fadirectory = 'M:\\EdgarHaener\\Capsules\\Batch20160610-5\\Semi-Sphere\\SS_15Jul\\'lse, names=[], savename=snm)

instances = [RCSS_Batch20160205d1_nr5_12022016, RCSS_Batch20160223d3_nr4_28022016, RCSS_Batch20160223d3_nr7_26022016,  RCSS_Batch20160306d2_nr6_08032016, RCSS_B20160531d4_nr9_5616, RCSS_GelBeads20160324d1_nr1_29032016]
savePath = 'M:\\EdgarHaener\\Capsules\\PlotScripts\\Semi-Sphere\\'
nameList = ['Batch20160205-1 #5 12/02/16', 'Batch20160223-3 #4 28/02/2016', 'Batch20160223-3 #7 26/02/2016', 'Batch20160306-2 #6 8/3/2016','Batch20160531-4 #9 05/06/2016',  'GelBeads20160324-1 #1 29/03/16']
force = [8.6, 5.6, 2.4, 4.2, 9.6, 100.0] #old measurments
snm='B0205_B0223_B0306_B0531_GB0324'

#srss.plotQVsFinalOffset(instances, savePath, forPub=False, names=nameList, savename=snm)
#srss._plotQVsAngle(instances, savePath, forPub=False, names=nameList, savename=snm)

#srss.plotQVsFinalOffsetAve(instances, savePath, forPub=False, names=nameList, savename=snm)
#srss.plotQVsAngleAve(instances, savePath, forPub=Fadirectory = 'M:\\EdgarHaener\\Capsules\\Batch20160610-3\\Semi-Sphere\\SS_30Jun\\'lse, names=nameList, savename=snm)
#
#srss.plotCaVsFinalOffsetAve(instances, force, savePath, forPub=False, names=nameList, savename=snm)
#srss.plotCaVsAngleAve(instances, force, savePath, forPub=False, names=nameList, savename=snm)
#
#srss.plotFinalVsInitialOffset(instances, savePath, forPub=False, names=nameList, savename=snm)

#instances = [RCSS_Batch20160205d1_nr5_12022016, RCSS_Batch20160223d3_nr4_29022016, RCSS_B20160531d4_nr9_5616, RCSS_GelBeads20160324d1_nr1_29032016]
#savePath = 'M:\\EdgarHaener\\Capsules\\PlotScripts\\Semi-Sphere\\'
#nameList = ['Batch20160205-1 #5 12'Batch20160505-2 #9 $F_{60%} = 1.3 mN$'/02/16', 'Batch20160223-3 #4 29/02/2016',  'Batch20160531-4 #9 05/06/2016', 'GelBeads20160324-1 #1 29/03/16']
#snm='B0205_B0223_B0531_GB0324'
#
#srss.plotQVsFinalOffset(instances, savePath, forPub=False, names=nameList, savename=snm)
#srss._plotQVsAngle(instances, savePath, forPub=False, names=nameList, savename=snm)
#srss.plotQVsFinalOffsetAve(instances, savePath, forPub=True, names=nameList, savename=snm)
#srss.plotQVsAngleAve(instances, savePath, forPub=True, names=nameList, savename=snm)
#srss.plotFinalVsInitialOffset(instances, savePath, forPub=False, names=[], savename=snm)

#for plotting
instances = [RCSS_Batch20160205d1_nr5_12022016,  RCSS_B20160531d4_nr9_5616, RCSS_B20160505d2_nr16_15616, RCSS_B20160622d2_nr6_24616, RCSS_GelBeads20160324d1_nr1_29032016]
savePath = 'M:\\EdgarHaener\\Capsules\\PlotScripts\\Semi-Sphere\\'
nameList = ['Batch20160205-1 #5 12/02/16','Batch20160531-4 #9 05/06/2016', 'Batch20160505-2 #16 15/06/2016',  'Batch20160622-3 #6 24/06/16', 'GelBeads20160324-1 #1 29/03/16']
force = [8.6, 9.6, 2.8, 35.0, 100.0] #old measurments
#force = [7.7, 6.2, 4.6, 100.0, 30.0]
snm='B0205_B0223_B0505_B0531_GB0324'

srss.plotQVsFinalOffsetAve(instances, savePath, forPub=True, names=nameList, savename=snm)
#srss.plotVelocitiesVsFinalOffsetAve(instances, savePath, forPub=False, names=nameList, savename=snm)
#srss._plotQVsAngle(instances, savePath, forPub=False, names=nameList, savename=snm)
#srss.plotQVsFinalOffsetAve(instances, savePath, forPub=False, names=nameList, savename=snm)
srss.plotQVsAngleAve(instances, savePath, forPub=True, names=nameList, savename=snm)
#srss.plotVelocityVsAngleAve(instances, savePath, forPub=False, names=nameList, savename=snm)
##
srss.plotCaVsFinalOffsetAve(instances, force, savePath, forPub=True, names=nameList, savename=snm)

srss.plotCaVsFinalOffsetAveSelected(instances, force, savePath, initialCentringLimit=0.1, forPub=False,  names=nameList, savename=snm)
srss.plotCaVsFinalOffsetAveSelected(instances, force, savePath, initialCentringLimit=0.2, forPub=False,  names=nameList, savename=snm)
srss.plotCaVsFinalOffsetAveSelected(instances, force, savePath, initialCentringLimit=0.3, forPub=False,  names=nameList, savename=snm)
srss.plotCaVsFinalOffsetAveSelected(instances, force, savePath, initialCentringLimit=3, forPub=False,  names=nameList, savename=snm)
srss.plotCaVsAngleAve(instances, force, savePath, forPub=True, names=nameList, savename=snm)
#
#srss.plotFinalVsInitialOffset(instances, savePath, forPub=False, names=nameList, savename=snm)
#
#srss.plotAveAngleVsInitialOffset(instances, savePath, forPub=False, names=nameList, savename=snm)
inst = [RCSS_B20160505d2_nr16_15616, RCSS_Batch20160205d1_nr5_12022016,  RCSS_B20160531d4_nr9_5616,  RCSS_B20160622d2_nr6_24616, RCSS_GelBeads20160324d1_nr1_29032016]
nl = ['HC3','HC1', 'HC4',  'HC7', 'HC9.']
srss.plotQVsFinalOffset_specificQ(inst, savePath, VolFlux=10.0, forPub=True, names=nl, savename=snm)
#srss.plotNormedFinalVsInitialOffset(instances, savePath, forPub=False, names=nameList, savename=snm)
#srss.plotQvsMeanHeightAtSS(instances, savePath, forPub=True, names=nameList, savename=snm)

#instances = [RCSS_B20160531d4_nr9_5616, RCSS_B20160505d2_nr16_15616,  RCSS_B20160622d2_nr6_24616, RCSS_GelBeads20160324d1_nr1_29032016]
#savePath = 'M:\\EdgarHaener\\Capsules\\PlotScripts\\Semi-Sphere\\'
#nameList = ['Batch20160531-4 #9 05/06/2016', 'Batch20160505-2 #16 15/06/2016', 'GelBeads20160324-1 #1 29/03/16', 'Batch20160622-3 #6 24/06/16']
#force = [9.6, 2.8,  35.0, 100.0,] #old measurments
##force = [7.7, 6.2, 4.6, 100.0, 30.0]
#snm='B0223_B0505_B0531_GB0324'
#srss.plotAveInitialOffsetVsQ(instances, savePath, forPub=False, names=nameList, savename=snm)
srss.plotFinalNormedVsInitialOffset(instances, force, savePath, forPub=True, names=nameList, savename=snm)

#==============================================================================
# All Batches
#==============================================================================
#RCSS_Batch20160223d3_nr7_26022016.pixelsPmm.fill(15.4)
instances = [RCSS_Batch20160205d1_nr5_12022016,  RCSS_B20160531d4_nr9_5616, RCSS_B20160505d2_nr16_15616, RCSS_GelBeads20160324d1_nr1_29032016, RCSS_Batch20160223d3_nr4_28022016, RCSS_Batch20160223d3_nr7_26022016,  RCSS_Batch20160306d2_nr6_08032016, RCSS_B20160622d2_nr6_24616]
savePath = 'M:\\EdgarHaener\\Capsules\\PlotScripts\\Semi-Sphere\\'
nameList = ['Batch20160205-1 #4 12/02/16','Batch20160531-4 #9 05/06/2016',  'Batch20160505-2 #15 15/06/2016', 'GelBeads20160324-1 #1 29/03/16',  'Batch20160223-3 #4 28/02/2016', 'Batch20160223-3 #7 26/02/2016', 'Batch20160306-2 #6 8/3/2016', 'Batch20160622-3 #6 24/06/16']
#force = [8.6, 9.6, 2.8, 100.0, 5.6, 2.4, 4.2, 35.0] #old measurments
force = [7.7, 6.2, 4.6, 100.0, 6.2, 7.6, 4.7, 30.0]
snm='All'

totalRuns=0
totalImages=0
for hh in range(len(instances)):
    print('%s: %d runs' %(nameList[hh], len(instances[hh].volumeFlux)))
    totalRuns += len(instances[hh].volumeFlux)
#    totalImages += len
print("Total number of runs = %d" %totalRuns)

#srss.plotQVsFinalOffset(instances, savePath, forPub=False, names=nameList, savename=snm)
###srss._plotQVsAngle(instances, savePath, forPub=False, names=nameList, savename=snm)
##srss.plotQVsFinalOffsetAve(instances, savePath, forPub=False, names=nameList, savename=snm)
###srss.plotQVsAngleAve(instances, savePath, forPub=False, names=nameList, savename=snm)
###
#srss.plotCaVsFinalOffsetAve(instances, force, savePath, forPub=False, names=nameList, savename=snm)
##srss.plotCaVsFinalOffsetAveSelected(instances, force, savePath, initialCentringLimit=0.2, forPub=False, names=nameList, savename=snm)
###srss.plotCaVsAngleAve(instances, force, savePath, forPub=False, names=nameList, savename=snm)
##
##srss.plotInitialOffsetVsQ(instances, savePath, forPub=False, names=nameList, savename=snm)
##srss.plotAveInitialOffsetVsQ(instances, savePath, forPub=False, names=nameList, savename=snm)
##srss.plotFinalVsInitialOffset(instances, savePath, forPub=False, names=nameList, savename=snm)
##srss.plotFinalNormedVsInitialOffset(instances, savePath, forPub=False, names=nameList, savename=snm)




#==============================================================================
#  Effect of Size in Semi-Cylinder
#==============================================================================
instances = [RCSS_B20160505d2_nr16_15616, RCSS_B20160531d4_nr9_5616,  RCSS_B20160622d2_nr9_24616, RCSS_B20160610d3_nr1_300616, RCSS_B20160610d5_nr7_150716]
savePath = 'M:\\EdgarHaener\\Capsules\\PlotScripts\\Semi-Sphere\\SizeEffect\\'
nameList = ['Batch20160505-2 #9 15/6/16 3.9mm 2.8mN','Batch20160531-1 #9 5/6/16 3.9mm 6.2 mN', 'Batch20160622-2 #9 24/06/16 3.8mm 20.5 mN', 'Batch20160610-3 #1 30/06/16 3mm 14.9 mN', 'Batch20160610-5 #7 15/07/2016 3.6mm 4.8 mN']
force = [2.8,6.2, 20.5, 14.9, 4.8] #old measurments
#force = [7.7, 20.9, 6.6]
snm='size_SS'



#print('len(nameList) = %d' %len(nameList))
#print('len(instances) = %d' %len(instances))

#srss.plotQVsFinalOffsetAve(instances, savePath, forPub=False, names=nameList, savename=snm)
##srss.plotQVsAngleAve(instances, savePath, forPub=False, names=nameList, savename=snm)
##
#srss.plotCaVsFinalOffsetAve(instances, force, savePath, forPub=False, names=nameList, savename=snm)
##srss.plotCaVsAngleAve(instances, force, savePath, forPub=False, names=nameList, savename=snm)
#
#srss.plotAveInitialOffsetVsQ(instances, savePath, forPub=False, names=nameList, savename=snm)


#==============================================================================
# Graphs for Effect of Centring in Half-Cylinder
#==============================================================================
instances = [  RCSS_Batch20160306d2_nr6_08032016, RCSS_Batch20160205d1_nr5_12022016,  RCSS_B20160531d4_nr9_5616, RCSS_B20160505d2_nr16_15616]
savePath = 'M:\\EdgarHaener\\Capsules\\PlotScripts\\Semi-Sphere\\Centring\\'
nameList = ['Batch20160306-2 #6 8/3/2016', 'Batch20160205-1 #4 12/02/16','Batch20160531-4 #9 05/06/2016',  'Batch20160505-2 #15 15/06/2016' ]
#force = [8.6, 9.6, 2.8, 100.0, 5.6, 2.4, 4.2, 35.0] #old measurments
force = [4.7, 7.7, 6.2, 4.6]
snm='Centring'

#srss.plotNormedFinalVsInitialOffset(instances, savePath, forPub=False, names=nameList, savename=snm)
#srss.plotNormedFinalVsInitialOffset(instances, savePath, forPub=True, names=nameList, savename=snm)
instances = [RCSS_Batch20160205d1_nr5_12022016,  RCSS_B20160531d4_nr9_5616, RCSS_B20160505d2_nr16_15616, RCSS_GelBeads20160324d1_nr1_29032016,  RCSS_Batch20160306d2_nr6_08032016, RCSS_B20160622d2_nr6_24616,  RCSS_B20160622d2_nr9_24616, RCSS_B20160610d3_nr1_300616]
totalRuns=0
for hh in range(len(instances)):
    print('%s: %d runs' %(instances[hh].name[0], len(instances[hh].volumeFlux)))
    totalRuns += len(instances[hh].volumeFlux)
print("Total number of runs = %d" %totalRuns)
