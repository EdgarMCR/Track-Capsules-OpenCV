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
import t_junction as tj
import ReadOutputFile as TJ_ROF
import ReadResultsFile as TJ_RRF
import track_capsule_TJ as TR
import CompareRuns as CP_TJ


def runfunc():
    #    directory = 'M:\\EdgarHaener\\Capsules\\Batch170615-002\\T-Junction\\2015-06-22\\Batch170615-002_#2\\'
    #    folder = 'Batch170615-002-#2-%dFPS-35mlPmin-2\\' %FPS
    
    #    directory = 'M:\\EdgarHaener\\Capsules\\Batch170615-002\\T-Junction\\2015-06-20\\Batch170615-002_#5\\'    
    #    folder = 'Batch170615-002_#5_%dFPS_5mlPmin-1\\' %FPS
        
    #    directory = 'M:\\EdgarHaener\\Capsules\\Batch270715-001\\T-Junction\\2015-08-04\\Batch270715-001-#5\\'
    #    folder =  'Batch270715-001-#5-%dFPS-15mlPmin-4\\' %FPS
        
    #    directory = 'M:\\EdgarHaener\\Capsules\\Batch040615-002\\T-Junction\\'
    
    #    directory = 'M:\\EdgarHaener\\Capsules\\GelBeads150715-1\\T-Junction\\2015-07-22\\GelBead150715-1-#4\\'
    #    folder =  'GelBead150715-1-#4-%dFPS-35mlPmin-8\\' %FPS
        
    #    directory = 'M:\\EdgarHaener\\Capsules\\Batch120615-004\\T-Junction\\'
    #    folder =  'Batch120615-001-#4-%dFPS-5mlPmin-4\\' %FPS
        
    #    
    #    directory = 'M:\\EdgarHaener\\Capsules\\Batch260615-001\\T-Junction\\2015-07-01\\Batch260615-001-#17\\'
    #    folder =  'Batch260615-001-#17-%dFPS-70mlPmin-1\\' %FPS
        
    #    directory = 'M:\\EdgarHaener\\Capsules\\GelBeads150730-1\\T-Junction\\2015-08-04\\GelBead150730-1-#1\\'
    #    folder =  'GelBead150730-1-#1-%dFPS-35mlPmin-1\\' %FPS
        
    #    directory = 'M:\\EdgarHaener\\Capsules\\Batch010615-001\\T-Junction\\Capsule#1\\'
    #    folder = 'Batch010615-001_#1-030615-5kcSt-1S-%dFPS-70mlPmin-5\\' %FPS
    
     
        
        
    
    #    centerline, width, pPmm, geometryTJ = None, None, None, None
    #    directory = 'M:\\EdgarHaener\\Capsules\\GelBeads150715-1\\T-Junction\\2015-07-22\\GelBead150715-1-#4\\'; geometryTJ=[35, 121, 535, 709] #GelBeads150715-1 
    #    centerline, width, pPmm = 622, 174, 22.3; directory = 'M:\\EdgarHaener\\Capsules\\GelBeads150715-1\\T-Junction\\2015-07-22\\GelBead150715-1-#4\\'; geometryTJ=[37, 121, 535, 709]; #GelBeads150715-1
    #    folder = 'GelBead150715-1-#4-%dFPS-70mlPmin-7\\' %FPS
    #    path=directory+folder
    #    directory = 'M:\\EdgarHaener\\Capsules\\GelBeads150730-1\\T-Junction\\2015-08-04\\GelBead150730-1-#1\\'; centerline, width, pPmm = 624.5, 175, 22.4; geometryTJ=[32, 119, 537, 712] #GelBeads150730-1 #11
    #    directory = '/mnt/MCND/EdgarHaener/Capsules/GelBeads150730-1/T-Junction/2015-08-04/GelBead150730-1-#1/'; 
    #    path = '/mnt/MCND/EdgarHaener/Capsules/GelBeads150730-1/T-Junction/2015-08-04/GelBead150730-1-#1/GelBead150730-1-#1-10FPS-5mlPmin-1' #GelBeads150730-1 #1
    #    centerline, width, pPmm = 646, 174, 22.3; geometryTJ=[28, 114, 559,  733] #Batch260615-001 #17
    
    #    centerline, width, pPmm = 637, 175, 22.4; geometryTJ=[72, 159, 550,  725] #Batch120615-004 #4 15ml/min
    #    centerline, width, pPmm = 637, 175, 22.4; geometryTJ=[35, 118, 549,  724] #Batch120615-004 #4 Other Runs
    #    centerline, width, pPmm = 631.5, 175, 22.4; geometryTJ=[73, 158, 544,  719] #Batch170615-002 #5 5ml/min
    #    centerline, width, pPmm = 632, 176, 22.4; geometryTJ=[34, 120, 543,  719] #Batch170615-002 #5 Other
    #    geometryTJ=[33, 118, 541, 717] #Batch170615-002 #2
    #    centerline, width, pPmm = 631, 174, 22.3; geometryTJ=[39, 125, 544, 718] #Batch010615-001 #1
    #    centerline, width, pPmm = 623, 176, 22.5; geometryTJ=[120, 535, 711] #Batch270715-001 #5
    
    #    centerline, width, pPmm = 628, 176, 22.6; geometryTJ=[119, 540, 716] #Batch100815-001-#8
    #    centerline, width, pPmm = 633, 166, 20.2; geometryTJ=[105, 550, 716] #Batch100815-001 #6 & #7
    #    centerline, width, pPmm = 634.5, 165, 20.2; geometryTJ=[105, 552, 717] #Batch100815-001 #3 & #4
    
    
     #geometryTJ = [Top Daugther Channel, Bottom Daugther Channel, left side Main Channel, right side Main Channel]

    #=============================================================================
    # GelBead150715-1
    #=============================================================================
    name4 = 'GelBead150715-1#1'
        
    centerline, width, pPmm = 622, 174, 22.3;  geometryTJ=[37, 121, 535, 709]; #GelBeads150715-1
    directory = 'M:\\EdgarHaener\\Capsules\\GelBeads150715-1\\T-Junction\\2015-07-22\\GelBead150715-1-#4\\';  
    FPS=40; folder =  'Batch040615-002-#1-1S-5kcSt-%dFPS-70mlPmin-6\\' %FPS
        
    path = directory + folder #Test2\\'
    pathRESLT4 = directory + '2015-07-22-GelBead150715-1-#4_Results.txt'
    #=============================================================================
    
    #=============================================================================
    # Batch040615-002
    #=============================================================================
    name3 = 'Batch040615-002#1'
    
     
    centerline, width, pPmm3 = 635.5, 173, 22.2; geometryTJ=[36, 122, 549,  722] #
     
    directory = 'M:\\EdgarHaener\\Capsules\\Batch040615-002\\T-Junction\\Capsule#1\\'
    FPS3=10; folder =  'Batch040615-002-#1-1S-5kcSt-%dFPS-5mlPmin-4\\' %FPS3
        
    path3 = directory + folder #Test2\\'
    pathRESLT3 = directory + 'T-Junction-Capsule#1_Results.txt'
    #=============================================================================
    
    #=============================================================================
    # Batch260615-001
    #=============================================================================
    name2 = 'Batch260615-001#17'
    
    centerline, width, pPmm2 = 646, 174, 22.3; geometryTJ2=[28, 114, 559,  733] #Batch260615-001 #17
     
    directory = 'M:\\EdgarHaener\\Capsules\\Batch260615-001\\T-Junction\\2015-07-01\\Batch260615-001-#17\\'
    FPS2=140; folder =  'Batch260615-001-#17-%dFPS-70mlPmin-1\\' %FPS2
        
    path2 = directory + folder #Test2\\'
    pathRESLT2 = directory + '2015-07-01-Batch260615-001-#17_Results.txt'
    #=============================================================================
    
    #=============================================================================
    # Batch120615-004
    #=============================================================================
    name1 = 'Batch120615-004#4'
    
    directory = 'M:\\EdgarHaener\\Capsules\\Batch120615-004\\T-Junction\\'
    FPS = 10; folder =  'Batch120615-001-#4-%dFPS-5mlPmin-4\\' %FPS
    FPS1 = 100; folder =  'Batch120615-001-#4-%dFPS-50mlPmin-1\\' %FPS1
    
    if FPS1 ==  30:
        centerline, width, pPmm1 = 637, 175, 22.4; geometryTJ1=[72, 159, 550,  725] #Batch120615-004 #4 15ml/min
    else:
        centerline, width, pPmm1 = 637, 175, 22.4; geometryTJ1=[35, 118, 549,  724] #Batch120615-004 #4 Other Runs
        
    path1 = directory + folder #Test2\\'
    pathRESLT1 = directory + 'Batch120615-004-T-Junction_Results.txt'
    #=============================================================================
    
    
    
    
    #=============================================================================
    # Run stuff
    #=============================================================================
    
    sp = "M:\\EdgarHaener\\Capsules\\PlotScripts\\T-Junction-RESLT\\"
    
    #tj.tryout(directory)
    #tj.averageTrajectories(directory = directory,
    #                       savePath = sp, 
    #                       pPmm = pPmm, 
    #                       batchName = name, 
    #                       forPub=False)
    
    
    # Plot Trajectory
#    tj.plotTrajecotry(path = path, pPmm=pPmm, geometryTJ = geometryTJ, savepath = sp, forPub=True)
#    commandsToGeneratePlotsForThesis()
#    areaAndPerimeter(sp)
         
    TJ1 = TJ_RRF.ResultsClass(pathRESLT1)
    TJ2 = TJ_RRF.ResultsClass(pathRESLT2)
    TJ3 = TJ_RRF.ResultsClass(pathRESLT3)
    TJ4 = TJ_RRF.ResultsClass(pathRESLT4)
    
    inst = [TJ1, TJ2, TJ3, TJ4]
    names = [name1, name2, name3, name4]
    dia = [3.77, 3.77, 3.87, 4.00]
    pPmm = [22.4, 22.3, 22.2, 22.3]
    force = [3.1, 6.2, 9.2, 180]

    CP_TJ.plotMeanMaxAcceleration(listOfInstances =inst, 
                              listOfNames = names, pPmm=pPmm,
                              savepath=sp, title  = '', savename ='',
                              show=True, forPub=True)    
#    TandNDTvsQ(inst, names, sp)
#    plotGauss()
#    plotTimeAndRelaxation(inst, names,dia, pPmm, force, sp)
#    migrationVelocity(sp)
#    velocity(sp)
#    acceleration(sp)
    
def velocity(sp):
    #=============================================================================
    # Batch040615-002
    #=============================================================================
    centerline, width, pPmm3 = 635.5, 173, 22.2; geometryTJ=[36, 122, 549,  722] #   
    directory = 'M:\\EdgarHaener\\Capsules\\Batch040615-002\\T-Junction\\Capsule#1\\'
    FPS3=10; folder =  'Batch040615-002-#1-1S-5kcSt-%dFPS-5mlPmin-4\\' %FPS3        
    path3 = directory + folder #Test2\\'
    #=============================================================================
    tj.plotVelocity(path3, pPmm3, FPS3,  savepath = sp, 
                    forPub=True, TJlim=[5.8,10.1], lim=[0, 16, 0, 4])
    FPS3=140; folder =  'Batch040615-002-#1-1S-5kcSt-%dFPS-70mlPmin-4\\' %FPS3        
    path3 = directory + folder #Test2\\'
    tj.plotVelocity(path3, pPmm3, FPS3,  savepath = sp, 
                    forPub=True, TJlim=[0.7, 1.1], lim=[0.3, 1.6, 0, 60])
                    
                    
    #=============================================================================
    # Batch120615-004
    #=============================================================================
    directory = 'M:\\EdgarHaener\\Capsules\\Batch120615-004\\T-Junction\\'
    FPS1 = 10; folder =  'Batch120615-001-#4-%dFPS-5mlPmin-10\\' %FPS1
    if FPS1 ==  30:
        centerline, width, pPmm1 = 637, 175, 22.4; geometryTJ1=[72, 159, 550,  725] #Batch120615-004 #4 15ml/min
    else:
        centerline, width, pPmm1 = 637, 175, 22.4; geometryTJ1=[35, 118, 549,  724] #Batch120615-004 #4 Other Runs
    path1 = directory + folder #Test2\\'
    tj.plotVelocity(path1, pPmm1, FPS1,  savepath = sp, 
                    forPub=True, TJlim=[7.5, 12.5], lim=[2, 21, 0, 4.5], cutoff=21)
    FPS1 = 10; folder =  'Batch120615-001-#4-%dFPS-5mlPmin-12\\' %FPS1
    path1 = directory + folder #Test2\\'
    tj.plotVelocity(path1, pPmm1, FPS1,  savepath = sp, 
                    forPub=True, TJlim=[7.5, 15], lim=[2, 21, 0, 4.5])

def acceleration(sp):
    #=============================================================================
    # Batch120615-004
    #=============================================================================
    directory = 'M:\\EdgarHaener\\Capsules\\Batch120615-004\\T-Junction\\'
    FPS1 = 10; folder =  'Batch120615-001-#4-%dFPS-5mlPmin-10\\' %FPS1
    if FPS1 ==  30:
        centerline, width, pPmm1 = 637, 175, 22.4; geometryTJ1=[72, 159, 550,  725] #Batch120615-004 #4 15ml/min
    else:
        centerline, width, pPmm1 = 637, 175, 22.4; geometryTJ1=[35, 118, 549,  724] #Batch120615-004 #4 Other Runs
    path1 = directory + folder #Test2\\'
    tj.plotAcceleration(path1, pPmm1, FPS1,  savepath = sp, 
                    forPub=True, TJlim=[7.5, 12.5], lim=[2, 21, -4.5, 4.5], cutoff=21)
    FPS1 = 100; folder =  'Batch120615-001-#4-%dFPS-50mlPmin-1\\' %FPS1
    path1 = directory + folder #Test2\\'
    tj.plotAcceleration(path1, pPmm1, FPS1,  savepath = sp, 
                    forPub=True, TJlim=[0.9, 1.4], lim=[0.0, 2, -200, 200.0])
    
def migrationVelocity(sp):    
    #=============================================================================
    # Batch260615-001
    #=============================================================================   
    centerline, width, pPmm2 = 646, 174, 22.3; geometryTJ2=[28, 114, 559,  733] #Batch260615-001 #17
    directory = 'M:\\EdgarHaener\\Capsules\\Batch260615-001\\T-Junction\\2015-07-01\\Batch260615-001-#17\\'
    FPS2=140; folder =  'Batch260615-001-#17-%dFPS-70mlPmin-1\\' %FPS2
    path2 = directory + folder 
    
    #=============================================================================
    # Batch120615-004
    #=============================================================================
    directory = 'M:\\EdgarHaener\\Capsules\\Batch120615-004\\T-Junction\\'
    FPS1 = 100; folder =  'Batch120615-001-#4-%dFPS-50mlPmin-1\\' %FPS1
    if FPS1 ==  30:
        centerline, width, pPmm1 = 637, 175, 22.4; geometryTJ1=[72, 159, 550,  725] #Batch120615-004 #4 15ml/min
    else:
        centerline, width, pPmm1 = 637, 175, 22.4; geometryTJ1=[35, 118, 549,  724] #Batch120615-004 #4 Other Runs
    path1 = directory + folder #Test2\\'
    
    #Plots of Migration Velocity
    anPos=[0.69, -3.5, 0.8, -6.0, 1.10, 0.2, 0.8, -2.0,]
    tj.plotMigrationVelocity(path2, pPmm2, FPS2, geometryTJ2, savepath = sp, 
                             forPub=True, lim=[0.6, 1.25, -10, 3.5], anPos= anPos)
    #Batch120615-001-nr4-100FPS-50mlPmin-1
    anPos=[1.29, -2, 1.4, -6, 1.92, 0.3, 1.4, -2,]
    tj.plotMigrationVelocity(path1, pPmm1, FPS1, geometryTJ1, savepath = sp, 
                             forPub=True, lim=[1.1, 2.1, -10, 3.5], anPos= anPos)

    
def plotTimeAndRelaxation(inst, names, dia, pPmm, force, sp):

    CP_TJ.plotTimeVsOffcentre(listOfClassInstances =inst, 
                              listOfNames = names, 
                              listOfD0 = dia, pPmm = pPmm,
                              nrCap=0, 
                              savepath=sp, 
                              show=True, forPub=False,
                              absoluteV=False, NDT=True)

    CP_TJ.plotTimeVsOffcentre(listOfClassInstances =inst, 
                              listOfNames = names, 
                              listOfD0 = dia, pPmm = pPmm,
                              nrCap=0, 
                              savepath=sp, 
                              show=True, forPub=True,
                              absoluteV=False, NDT=True)

    CP_TJ.plotTimeVsOffcentre(listOfClassInstances =inst, 
                              listOfNames = names, 
                              listOfD0 = dia, pPmm = pPmm,
                              nrCap=0, 
                              savepath=sp, 
                              show=True, forPub=False, 
                              absoluteV=True, NDT=True)
                              
    CP_TJ.plotTimeVsOffcentre(listOfClassInstances =inst, 
                              listOfNames = names, 
                              listOfD0 = dia, pPmm = pPmm,
                              nrCap=0, 
                              savepath=sp, 
                              show=True, forPub=True, 
                              absoluteV=True, NDT=True)

    CP_TJ.plotMinSpeed(listOfInstances =inst, listOfNames = names, 
                      listOfd0 = dia, savepath=sp, title  = '', 
                      savename ='', show=True, nonDim=True, forPub=False)
    
    CP_TJ.plotRelaxationDistanceVsCa(listOfInstances =inst, listOfNames = names, 
                                     pPmm = pPmm, force = force, 
                                     listOfD0 = dia, savepath=sp, title  = '', 
                                     savename ='', 
                                     show=True, forPub=False, ND=True)
    CP_TJ.plotRelaxationDistanceVsCa(listOfInstances =inst, listOfNames = names, 
                                     pPmm = pPmm, force = force, 
                                     listOfD0 = dia, savepath=sp, title  = '', 
                                     savename ='', 
                                     show=True, forPub=True, ND=False)
                                     
    CP_TJ.plotRelaxationTimeVsCa(listOfInstances =inst, listOfNames = names, 
                                 force = force, savepath=sp, title  = '', 
                                 savename ='', 
                                 show=True, forPub=False, ND=True)
                                     
    CP_TJ.plotRelaxationTime(listOfInstances =inst, listOfNames = names, 
                             savepath=sp, title  = '', savename ='', 
                             show=True, forPub=True, ND=False)                                

    CP_TJ.plotRelaxationTimeVsCa(listOfInstances =inst, listOfNames = names, 
                                 force = force, savepath=sp, title  = '', 
                                 savename ='', 
                                 show=True, forPub=True, ND=False)     

    CP_TJ.plotRelaxationTime(listOfInstances =inst, listOfNames = names, 
                             savepath=sp, title  = '', savename ='', 
                             show=True, forPub=True, ND=True)                                

    CP_TJ.plotRelaxationTimeVsCa(listOfInstances =inst, listOfNames = names, 
                                 force = force, savepath=sp, title  = '', 
                                 savename ='', 
                                 show=True, forPub=True, ND=True)                                   

    CP_TJ.plotRelaxationSpeed(listOfInstances =inst, listOfNames = names, 
                              pPmm=pPmm, savepath=sp, title  = '', 
                              savename ='', show=True, forPub=True)

    
def gaus(x,a,x0,sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))
    
def plotGauss():
    a=1.0
    x0=0.0
    x=np.arange(-1.0, 1.0, 0.001)
    plt.figure()
    
    sigma=1.0
    plt.plot(x, gaus(x,a,x0,sigma), label='1.0')
    sigma=0.1
    plt.plot(x, gaus(x,a,x0,sigma), label='0.1')
    sigma=10.0
    plt.plot(x, gaus(x,a,x0,sigma), label='10.0')
    plt.legend()
    
def TandNDTvsQ(inst, names, sp):
    CP_TJ.plotTvsQ(listOfInstances = inst, 
                   listOfNames = names, 
                   savepath=sp, 
                   title  = '', 
                   savename ='', 
                   show=True, 
                   forPub=True, 
                   nondimensionalise=False)

    CP_TJ.plotTvsQ(listOfInstances = inst, 
                   listOfNames = names, 
                   savepath=sp, 
                   title  = '', 
                   savename ='', 
                   show=True, 
                   forPub=True, 
                   nondimensionalise=True)

def averageOffset(pathRESLT1, pathRESLT2, pathRESLT3,pathRESLT4):
    TJ1 = TJ_RRF.ResultsClass(pathRESLT1)
    print(TJ1.offCentre)
    Lic = TJ1.offCentre/22.4
    TJ2 = TJ_RRF.ResultsClass(pathRESLT2)
    TJ3 = TJ_RRF.ResultsClass(pathRESLT3)
    TJ4 = TJ_RRF.ResultsClass(pathRESLT4)
    Lic = np.concatenate((Lic, TJ2.offCentre/22.3))
    Lic = np.concatenate((Lic, TJ3.offCentre/22.2))
    Lic = np.concatenate((Lic, TJ4.offCentre/22.3))
    
    print("Mean Lic = %f" %(np.mean(Lic)))
    
def commandsToGeneratePlotsForThesis():
    sp = "M:\\EdgarHaener\\Capsules\\PlotScripts\\T-Junction-RESLT\\"
    trajectoryPlots(sp)

def trajectoryPlots(savepath):
    # Batch120615-004
    directory = 'M:\\EdgarHaener\\Capsules\\Batch120615-004\\T-Junction\\'
    FPS = 10; folder =  'Batch120615-001-#4-%dFPS-5mlPmin-4\\' %FPS
    path1 = directory + folder #Test2\\'
    FPS = 80; folder =  'Batch120615-001-#4-%dFPS-40mlPmin-6\\' %FPS
    path2 = directory + folder #Test2\\'
    centerline, width, pPmm = 637, 175, 22.4; geometryTJ=[35, 118, 549,  724] #Batch120615-004 #4 Other Runs
    
    tj.plotTrajecotry(path = path1, pPmm=pPmm, geometryTJ = geometryTJ, savepath = savepath, forPub=True)
    tj.plotTrajecotry(path = path2, pPmm=pPmm, geometryTJ = geometryTJ, savepath = savepath, forPub=True)
    
def areaAndPerimeter(savepath):
    # Batch120615-004
    directory = 'M:\\EdgarHaener\\Capsules\\Batch120615-004\\T-Junction\\'
    FPS1 = 20; folder =  'Batch120615-001-#4-%dFPS-10mlPmin-4\\' %FPS1 ;    path1 = directory + folder #Test2\\'
    FPS2 = 100; folder =  'Batch120615-001-#4-%dFPS-50mlPmin-1\\' %FPS2 ;    path2 = directory + folder #Test2\\'
    centerline, width, pPmm = 637, 175, 22.4; geometryTJ=[35, 118, 549,  724] #Batch120615-004 #4 Other Runs
    
    TJ_ROF.plotAreaAndPerimeter(path1, FPS1, pPmm, 3.0*FPS1, 5.75*FPS1 , xmax=10.0, savepath=savepath)
    TJ_ROF.plotAreaAndPerimeter(path2, FPS2, pPmm, 0.9*FPS2, 1.4*FPS2, xmax=2.2, savepath=savepath)
    

    # Batch040615-002
    centerline, width, pPmm = 635.5, 173, 22.2; geometryTJ=[36, 122, 549,  722] #
     
    directory = 'M:\\EdgarHaener\\Capsules\\Batch040615-002\\T-Junction\\Capsule#1\\'
    FPS1=10; folder =  'Batch040615-002-#1-1S-5kcSt-%dFPS-5mlPmin-2\\' %FPS1 ; path1 = directory + folder       
    FPS2=140; folder =  'Batch040615-002-#1-1S-5kcSt-%dFPS-70mlPmin-6\\' %FPS2 ; path2 = directory + folder       
    
    TJ_ROF.plotAreaAndPerimeter(path1, FPS1, pPmm, 5.0*FPS1, 9*FPS1 , xmax=16.0, savepath=savepath)
    TJ_ROF.plotAreaAndPerimeter(path2, FPS2, pPmm, 0.575*FPS2, 0.875*FPS2, xmax=1.5, savepath=savepath)


if __name__ == "__main__":
#    print("Starting 'funfunc'")
    runfunc()