# -*- coding: utf-8 -*-
"""
Semi-Sphere
Functions to track capsules in top-view images using OpenCV.

Created on 08.02.2016

@author: Edgar Haener
edgar.haner@gmail.com

"""
from __future__ import absolute_import, division, print_function
import numpy as np
import matplotlib.pyplot as plt
import cv2
import os
import time
#import scipy.stats
from scipy import stats


import general as gen
from general import filenameClass

import analysis as ana
#from analysis import DataRun

ERROR_CONST=-1

class ParametersPinching(gen.Parameters):
    """Collect all information need for Pinching Setup. Assumes inflow from right,
    towards expansion on the left.
    

    ChannelTop         Pixel y-position of top of channel
    ChannelBottom      Pixel y-position of top of channel   
    EndChannel         Pixel x-position of end of channel, start of
                       diffuser
    LeftCornerInflow   Pixel x-position of left corner of pinching inflow
    RightCornerInflow  Pixel x-position of right corner of pinching inflow
                          
    """
    
    def __init__(self, listingImageFunc,coverImgFunc, baseDirectory= None, folder= None, d0= None, pPmm= None, 
                 rotation= None, FPS= None, backgroundImage= None, 
                 ChannelTop = None, ChannelBottom=None, EndChannel=None, 
                 LeftCornerInflow = None, RightCornerInflow = None,
                 readParametersFromFile=None):
        
        if baseDirectory != None:
            gen.Parameters.__init__(self, baseDirectory, folder, d0, pPmm, 
                 rotation, FPS, backgroundImage)
            self.ChannelTop = ChannelTop 
            self.ChannelBottom  =ChannelBottom 
            self.EndChannel = EndChannel
            self.LeftCornerInflow = LeftCornerInflow
            self.RightCornerInflow = RightCornerInflow
            self._writeToFile()
            
        elif readParametersFromFile != None:
            self._readFromFile(readParametersFromFile)
        
        self.setFunc(listingImageFunc,coverImgFunc)
        
        
        
        
    def _writeToFile(self):
        '''         Write a parameters file in path         '''
        pathParameters = self.path + 'ImageAnalysis_Parameters_'+time.strftime("%Y-%m-%d")+'.txt' 
        print('os.path.isdir(self.path):', end='');print(os.path.isdir(self.path))
        fileParameters = open(pathParameters, 'w')
        fileParameters.write('# Parameters for the analysis of the images. Pinching Setup. \n')
        fileParameters.write('# Lines starting with # are ignored \n')
        fileParameters.write('# The first position holds the parameter, the second a description, tab separated  \n')
        
        fileParameters.write('%s \t self.baseDirectory \n' %(self.baseDirectory))
        fileParameters.write('%s \t self.folder \n' %(self.folder))
        fileParameters.write('%s \t self.path \n' %(self.path))
        
        fileParameters.write('%f \t self.d0 \n' %(self.d0))
        fileParameters.write('%f \t self.pPmm \n' %(self.pPmm))
        fileParameters.write('%f \t self.rotation \n' %(self.rotation))
        
        fileParameters.write('%d \t self.FPS \n' %(self.FPS))
        fileParameters.write('%s \t self.backgroundImage \n' %(self.backgroundImage))
        
        
        fileParameters.write('%d \t self.ChannelTop \n' %(self.ChannelTop))
        fileParameters.write('%d \t self.ChannelBottom \n' %(self.ChannelBottom))
        fileParameters.write('%d \t self.EndChannel \n' %(self.EndChannel))
        fileParameters.write('%d \t self.LeftCornerInflow \n' %(self.LeftCornerInflow))
        fileParameters.write('%d \t self.RightCornerInflow \n' %(self.RightCornerInflow))

        fileParameters.close()
        
    def _printParameter(self):
        '''        Print Parameters         '''

        print('%s \t self.baseDirectory \n' %(self.baseDirectory))
        print('%s \t self.folder \n' %(self.folder))
        print('%s \t self.path \n' %(self.path))
        
        print('%f \t self.d0 \n' %(self.d0))
        print('%f \t self.pPmm \n' %(self.pPmm))
        print('%f \t self.rotation \n' %(self.rotation))
        
        print('%d \t self.FPS \n' %(self.FPS))
        print('%s \t self.backgroundImage \n' %(self.backgroundImage))
        
        
        print('%d \t self.ChannelTop \n' %(self.ChannelTop))
        print('%d \t self.ChannelBottom \n' %(self.ChannelBottom))
        print('%d \t self.EndChannel \n' %(self.EndChannel))
        print('%d \t self.LeftCornerInflow \n' %(self.LeftCornerInflow))
        print('%d \t self.RightCornerInflow \n' %(self.RightCornerInflow))
        
        
    def _readFromFile(self, path):
        '''Find and read parameter file in path '''
        
        
        listOfPamameterfiles=[]        
        for fileName in os.listdir(path):
            if fileName.endswith(".txt"):
                if fileName.find("ImageAnalysis_Parameters_") != -1:
                    listOfPamameterfiles.append(fileName)
#                    print(fileName)
#        print('befor sorting')
#        print(listOfPamameterfiles)
#        print('after sorting')
#        print(sorted(listOfPamameterfiles))
#        print(listOfPamameterfiles[:-1])
        
        fileParameters = open(os.path.join(path,listOfPamameterfiles[-1]), 'r')
        lines = fileParameters.readlines()
        fileParameters.close()       
#        print(lines)
        
        orderedLines=[]
        for line in lines:
            entriesLine=line.split('\t')
            if entriesLine[0][0].strip() != '#' and entriesLine[0][0].strip() != '': 
                orderedLines.append(entriesLine[0].strip())
#        print(orderedLines)
        if len(orderedLines) != 13:
            print('Error reading file, incorrect number of arguments!')
        else:
            baseDirectory = orderedLines[0]
            folder = orderedLines[1]
            path = orderedLines[2]
            
            d0 = float(orderedLines[3])
            pPmm = float(orderedLines[4])
            rotation = float(orderedLines[5])
            
            FPS = int(float(orderedLines[6]))
            backgroundImage = orderedLines[7]
            
            gen.Parameters.__init__(self, baseDirectory, folder, d0, pPmm, 
                 rotation, FPS, backgroundImage)
            
            #Pinching specific            
            self.ChannelTop =  int(float(orderedLines[8]))
            self.ChannelBottom  = int(float(orderedLines[9]))
            self.EndChannel =  int(float(orderedLines[10]))
            self.LeftCornerInflow = int(float(orderedLines[11]))
            self.RightCornerInflow = int(float(orderedLines[12]))
            
def sortPhotosPCO(path, fileType = '.png', prefixleng=10):
    """
    PCO camera
    
    Typical filename :     GelBeads20160324-1_#1_PS_4Apr_20FPS_10mlPmin_0mlPmin_3_0002.png
    #TODO: fix this
    """
    sperator='_' #'-'

    dirs = os.listdir(path)
    numberOfJPG=0

#    filenameList = [ filenameClass() for i in range(numberOfJPG)]
    filenameList=[]
    
    i=-1
    for fname in dirs:
        if (fname[-4:len(fname)] == fileType):
            try:
                d=fname.rfind(os.sep, 1, -1)
                if fname[d+1:d+3] == '._':
#                    print(fname[d+1:d+3])
                    continue
                i+=1            
                name=fname[:-4]
                #find first dash (day)
                d1=name.find(sperator,-5,-1)
#                print('number = ' + name[d1+1:])
                
#                print('abs(int(float(name[d1+1:]))) = ', end=''); print(abs(int(float(name[d1+1:]))))
    
                temp=filenameClass(fname, -1, abs(int(float(name[d1+1:]))) )
                filenameList.append(temp)
                numberOfJPG += 1
            except:
                print('\t Not inluding %s' %fname)
    
    #sort by milliseconds
    filenameList.sort(key=lambda filenameClass: filenameClass.number)
#    newlist = sorted(filenameList, key=lambda gen.filenameClass: gen.filenameClass.number) 
#    for obj in newlist:
#        obj.printOut()
    return filenameList, numberOfJPG
            
            
def coverSidePinching(img, offset, ParameterClass):
    '''Take image and delete areas not of intereste'''

    top = ParameterClass.ChannelTop
    bottom = ParameterClass.ChannelBottom 
    end = ParameterClass.EndChannel
    lc = ParameterClass.LeftCornerInflow
    rc = ParameterClass.RightCornerInflow
    
    #Get size of image
    yImg,xImg = img.shape
    
#    cv2.rectangle(img, (x1, y1), (x2, y2), (255,0,0), 2)  
#    x1,y1 ------
#    |          |
#    |          |
#    |          |
#    --------x2,y2
    #cover area above and below main channel
    cv2.rectangle(img, (end+offset, top+offset), (xImg, 0), 255, thickness=-1)   
    
    cv2.rectangle(img, (end+offset, bottom+offset), (lc+offset, yImg),  255, thickness=-1)
    cv2.rectangle(img, (rc+offset, bottom+offset), (xImg, yImg),  255, thickness=-1)
    
    
    #Cover outside
    side=5
    cv2.rectangle(img, (0, offset+side), (xImg, 0), 255, thickness=-1)   
    cv2.rectangle(img, (0, yImg - offset-side), (xImg, yImg), 255, thickness=-1)   
#    cv2.rectangle(img, (end+offset, bottom+offset), (xImg, yImg),  255, thickness=-1)
    
    #through rotation, there is a black border that seems to throw off the auto-
    # matic threshold setting, add a 3 pixel wide margine
#    crop=8
#    cv2.rectangle(img, (0, 0), (crop, yImg), 255, thickness=-1)  
#    cv2.rectangle(img, (0, yImg), (xImg, yImg-crop), 255, thickness=-1)  
    
    #time stampe
#    cv2.rectangle(img, (50, 12), (300, 6), 255, thickness=-1)  
    return img



# =============================================================================
# Extract relevant measurments from found data
# =============================================================================
def runFunction(path):
    '''Runs DataRunPinching on path'''
#    print('Starting runFunction')
    DRSS = DataRunPinching(path)
    DRSS.plot = False
    DRSS.extractMeasures()
#    DRSS = DataRunPinching(path)
#    DRSS._writeToFilePinching()
#    DRSS.writeNameToFile()
#    del DRSS
    
class DataRunPinching(ana.DataRun):
    
    def __init__(self, path):
        #load parameters from image analysis
        self.PP = ParametersPinching(listingImageFunc = sortPhotosPCO,
                                   coverImgFunc = coverSidePinching,
                                   readParametersFromFile=path)
        
        if self.PP.path != path:
            self.PP.path = path
            self.PP.baseDirectory = os.path.dirname(os.path.dirname(path))
            pa , folder = os.path.split(path)
            p , folder = os.path.split(pa)
            self.PP.folder = folder
            
        #use base class to load file
        ana.DataRun.__init__(self, self.PP.path, self.PP.pPmm, self.PP.FPS)
        
        self.plot=True
#        self.plot=False
        self.volumFlux2=-1
        self._getVolumnFluxes()
        
#        self.extractMeasures()
        
    def extractMeasures(self):
        '''Extract measures of interest from the raw tracking data and write to 
        file. '''
        self._findCentering()
        self._getFinalPosition()
        self._maxDeformation()        
#        self._findFinalDiameter()
        self.d=-1; self.errd = -1
        self._findRelaxationSize()
        self._findGap()
        self._writeToFilePinching()
        
    def writeNameToFile(self):
        '''Debug function'''
        directory = os.path.dirname(self.path)
        directory = os.path.dirname(directory)
        fileData = open(os.path.join(directory,'Names.txt'), 'a')
        fileData.write('%s \n' %(self.name))
        fileData.close()
        
    def _getVolumnFluxes(self):
        '''Find volume flux from folder name'''
        name=self.path[:-2]
        d1=name.find('mlPmin')
        d2=name[:d1].rfind('_')
        #I'm inconsistent, sometimes us _ or - as seperator, check for both
        if d2 ==-1 or (d1-d2)>20:
            d2 = name[:d1].rfind('-')
        d3=name[d2:d1].find('m')
        if d3==-1:
            volFlux=name[d2+1:d1]
        else:
            volFlux=name[d2+1:d2+d3]
        
        n2 = name[:d1+1]+name[d1+6:]
        print('n2 = %s' %n2)
        #try to find second volume flux
        d1=n2.find('mlPmin')
        d2=n2[:d1].rfind('_')
        #I'm inconsistent, sometimes us _ or - as seperator, check for both
        if d2 ==-1 or (d1-d2)>20:
            d2 = name[:d1].rfind('-')
        d3=n2[d2:d1].find('m')

        if d3==-1:
            volFlux2=n2[d2+1:d1]
        else:
            volFlux2=n2[d2+1:d2+d3]
            
        print('Volume Fluxes:  %s and %s ' %(volFlux, volFlux2))
            
        self.volumFlux  = float(volFlux.replace("p", "."))
        self.volumFlux2  = float(volFlux2.replace("p", "."))
        
    def _getFinalPosition(self):
        '''Finds final position of capsule after the diffuser and check whether
        capsule went right (up) or left (down) (from the capsules point of view.'''
        #Central line can be evaluated from PP.ChannelTop  and PP.ChannelBottom
        
        centreline = (self.PP.ChannelTop +  self.PP.ChannelBottom)/2.0
        
        # Somewhat arbitrary end position to measure: 6 diameters
        xposAfterDiffuser = self.PP.EndChannel - \
            (self.PP.d0*6.0)*self.PP.pPmm
            
#        xposAfterDiffuser = self.PP.EndChannel - \
#            (self.PP.d0*4.0)*self.PP.pPmm
        
        #3/4 of a diameter from the edge of the image to only capture full outlines
        minXpos = (3.0*self.PP.d0/4.0)*self.PP.pPmm
        #remove entries with no values
        cent_x, cent_y = self._centroidsWithoutErrors()
        ri = int(len(cent_x)/4.0) #relevant indexes are the last 4th
        cent_xs, cent_ys =cent_x[-ri:], cent_y[-ri:] 
#        print("cent_xs = )
        print('xposAfterDiffuser = %d' %xposAfterDiffuser)
        print('minXpos = %d' %minXpos)
        #Find if we have measurments at that position
        indexes = [i for i in range(len(cent_xs)) if cent_xs[i] < xposAfterDiffuser \
                    and cent_xs[i] > minXpos ]
        
#        print('indexes = ', end=''); print(indexes)
        if np.mean(cent_ys[indexes]) < centreline:
            self.wentRight = True
            distFromCentreline = centreline - cent_ys[indexes] 
        else:
            self.wentRight = False
            distFromCentreline = cent_ys[indexes] - centreline
#        print('self.wentRight :  ' +str(self.wentRight))
        
        #remove entries below zero
        distFromCentreline = distFromCentreline[distFromCentreline>0.0]
        
        self.distanceFromCentreline = np.mean(distFromCentreline)
        self.distanceFromCentrelineSTD = np.std(distFromCentreline)
        
#        print('distFromCentreline = ' +str(distFromCentreline))
#        plt.figure()
#        plt.plot(cent_x[indexes], distFromCentreline, 'bs')
#        plt.xlabel('x-position Centroid [pixels]'); plt.ylabel('Distance from centreline [pixels]')       
        
        if self.plot:
#            yfit=[];
#            for x in self.pos_diffuser_x:
#                yfit.append(ana.func(x, self.diffuser_m, self.diffuser_c))
            
            fig = plt.figure(figsize=(8, 6), dpi=200); ax = fig.add_subplot(111)
            plt.plot(cent_x, cent_y, 'bs', label='Centroid position')
            plt.plot(cent_xs[indexes], cent_ys[indexes], 'ro', label='Centroid position used for final offset')
            plt.plot(self.intialCentreingXpos, self.intialCentreingYpos, 'g<', label='Centroid position used for initial offset')
#            plt.plot(self.pos_diffuser_x, self.pos_diffuser_y, 'ch',label='Centroid position for angle measurment' )
#            xmax=self.PP.EndChannel
#            ymin=ana.func(0.0, self.diffuser_m, self.diffuser_c)
#            ymax = ana.func(xmax, self.diffuser_m, self.diffuser_c)
#            plt.plot([0, xmax], [ymin, ymax], 'c-', label='Fit angle = %f' %(self.diffuser_angle))

            plt.plot([0, np.max(cent_x)], [centreline, centreline], 'k-', label='Centreline')
            
            if self.wentRight:
                dc = centreline - self.distanceFromCentreline
            else:
                dc = centreline + self.distanceFromCentreline

            ic = centreline + self.aveOffset
                
            plt.plot([0, np.max(cent_x)], [dc, dc], 'r-', label='Final Offset')
            plt.plot([0, np.max(cent_x)], [ic, ic], 'g-', label='Initial Offset')
            plt.xlabel('Centroid x-position [pixels]', fontsize=12); plt.ylabel('Centroid y-position [pixels]', fontsize=12)

            plt.text(0.1, 0.6, 'Final Offset from Centreline = %.3f mm' %(self.distanceFromCentreline/self.PP.pPmm), fontsize = 16, horizontalalignment='left', verticalalignment='center', transform = ax.transAxes)
            plt.text(0.1, 0.4, 'Initial Centreing = %.3f mm' %(self.aveOffset/self.PP.pPmm), fontsize = 16, horizontalalignment='left', verticalalignment='center', transform = ax.transAxes)
            
            plt.title("Centroid position -  " +self.name, fontsize=12)
            plt.legend(loc='best', fontsize=6, ncol=2)
            fig.tight_layout()

            plt.savefig((self.path+"CentroidPosition_" +self.name+".jpg"), dpi=300)
            
        else:
            print('self.plot = %s' %self.plot)

#        print('self.distanceFromCentreline: ' +str(self.distanceFromCentreline) + ' +/- ' + str(self.distanceFromCentrelineSTD))
    
    def _maxDeformation(self):
        '''Find max deformation, as measured by taylor deformation parameter'''
        
        #remove entries with no values
        cent_x, cent_y = self._centroidsWithoutErrors()
 
        #First find area around semi-sphere, which is at x-pos end-channel        
        startSS = self.PP.RightCornerInflow
        endSS = self.PP.EndChannel - 2*self.PP.d0*self.PP.pPmm
        indexes = [i for i in range(len(cent_x)) if cent_x[i] < startSS and cent_x[i] > endSS]
        
#        centreline = (self.PP.ChannelTop + self.PP.ChannelBottom)/2.0
#        plt.figure()
#        plt.plot(cent_x, cent_y, 'bs')
#        plt.plot(cent_x[indexes], cent_y[indexes], 'ro')
#        plt.plot([0, np.max(cent_x)], [centreline, centreline], 'k-')
#        
        self.maxd12 = np.max(self.d12[indexes])
        
        #ax width
        self.maxWidth = np.max(self.width[indexes])
        #minimum width
        self.minHeight = np.min(ana._removeERRORCONST(self.height[indexes]))
        
        #======================================================================
        #
        # Find the deformation between the end of the channel and the end of
        # inflow and at the inflow
        #
        # Consider the width of the capsule at the end of the channel and average
        # over all capsule outlines that are there.
        #======================================================================
        indexesInflow = [i for i in range(len(self.centroid_x)) if (self.centroid_x[i] - self.width[i]/2.0) < self.PP.LeftCornerInflow and (self.centroid_x[i] + self.width[i]/2.0) > self.PP.LeftCornerInflow]

            
        print('Inflow: # of points = %d' %(len(indexesInflow)))
        if len(indexesInflow) >0 :
            self.maxd12Inflow = np.max(self.d12[indexesInflow])
            self.maxWidthInflow = np.max(self.width[indexesInflow])
            self.minHeightInflow = np.min(self.height[indexesInflow])
            self.meanHeightInflow = np.mean(self.height[indexesInflow])
        else:
            self.maxd12Inflow = -1; self.maxWidthInflow = -1; self.minHeightInflow = -1; self.meanHeightInflow = -1
        
        indexes2mm = [i for i in range(len(self.centroid_x)) if (self.centroid_x[i] - self.width[i]/2.0) < self.PP.EndChannel and (self.centroid_x[i] + self.width[i]/2.0) > self.PP.EndChannel]

            
        print('2mm: # of points = %d' %(len(indexes2mm)))
        if len(indexes2mm) >0 :
            self.maxd122mm = np.max(self.d12[indexes2mm])
            self.maxWidth2mm = np.max(self.width[indexes2mm])
            self.minHeight2mm = np.min(self.height[indexes2mm])
            self.meanHeight2mm = np.mean(self.height[indexes2mm])
        else:
            self.maxd122mm = -1; self.maxWidth2mm = -1; self.minHeight2mm = -1; self.meanHeight2mm = -1
            
            
        self.mean_xd_Inflow, self.mean_yd_Inflow = self._findCentralWidthAndHeight(indexesInflow)
        self.mean_xd_2mm, self.mean_yd_2mm = self._findCentralWidthAndHeight(indexes2mm)
        
        self.mean_ydf_Inflow = self._findHeightAtXpos(indexesInflow, self.PP.LeftCornerInflow)
        self.mean_ydf_2mm = self._findHeightAtXpos(indexes2mm, self.PP.EndChannel)
        
        print('Inflow: xd = %f yd=%f ydf=%f d_12=%f \t indexes: %d - %d' %(self.mean_xd_Inflow, self.mean_yd_Inflow , self.mean_ydf_Inflow, abs(self.mean_xd_Inflow - self.mean_yd_Inflow)/(self.mean_xd_Inflow+self.mean_yd_Inflow), np.min(indexesInflow), np.max(indexesInflow)))
        print('2mm: xd = %f yd=%f  ydf=%f d_12=%f \t indexes: %d - %d'%(self.mean_xd_2mm, self.mean_yd_2mm, self.mean_ydf_2mm, abs(self.mean_xd_2mm - self.mean_yd_2mm)/(self.mean_xd_2mm+self.mean_yd_2mm), np.min(indexes2mm), np.max(indexes2mm)))
        
#        if self.plot:
#            fig = plt.figure(figsize=(8, 6), dpi=200); ax = fig.add_subplot(111)
#            plt.plot(self.centroid_x, self.centroid_y, 'bs', label='Centroid position')
#            
#            plt.plot(self.centroid_x[indexesInflow], self.centroid_y[indexesInflow], 'ro', label='Centroid position used for Inflow')
#            plt.plot(self.centroid_x[indexes2mm], self.centroid_y[indexes2mm], 'g<', label='Centroid position used for 2mm Strech')
#            
#            ymin, ymax = ax.get_ylim()
#            
#            plt.plot([startInflow, startInflow], [ymin, ymax], 'r-')
#            plt.plot([endInflow+0.5, endInflow+0.5], [ymin, ymax], 'r-')
#            plt.plot([start2mm-0.5, start2mm-0.5], [ymin, ymax], 'g-')
#            plt.plot([end2mm, end2mm], [ymin, ymax], 'g-')
#            plt.title("Centroid position -  " +self.name, fontsize=12)
#            plt.legend(loc='best', fontsize=6, ncol=2)
#            fig.tight_layout()
#
#            plt.savefig((self.path+"CentroidPositionForHeight_" +self.name+".jpg"), dpi=300)
        
    def _findCentering(self):
        '''find centering of capsule in channel leading up 2dn inflow'''
        
        #find centerline
        centreline = (self.PP.ChannelTop + self.PP.ChannelBottom)/2.0
        
        # We know the start of the pinching inflow
        endSS = self.PP.RightCornerInflow + (self.PP.d0)*self.PP.pPmm
#        endSS = self.PP.RightCornerInflow + (0.25*self.PP.d0)*self.PP.pPmm
        
        #remove entries with no values
        cent_x, cent_y = self._centroidsWithoutErrors()
        
        #find max x position
        maxX = np.max(cent_x)
        
        # lets not take the first half diameter as there might be artefacts
        # from capsule tracking as the capsule enters the image
        startSS_IgnoreLast = maxX - 0.5*self.PP.d0*self.PP.pPmm
        
        #Take 2 capsule diameters
        startSS_3capsules = endSS + 2.0*self.PP.d0*self.PP.pPmm
        
        if startSS_IgnoreLast > startSS_3capsules:
            startSS =startSS_3capsules
        else:
            startSS = startSS_IgnoreLast
        
        #Find if we have measurments at that position
        indexes = [i for i in range(len(cent_x)) if cent_x[i] > endSS and cent_x[i] < startSS]
        
        #TODO: measure undeformed diameter at the end of the diffuser
        #1. get correct index, pretty certain its wronge
#        self.d, self.errd = self.fitCircle(indexes)
#        print('d = %f =/- %f' %(self.d/self.PP.pPmm, self.errd/self.PP.pPmm))
        
        meanYpos=cent_y[indexes]
        
        offset = meanYpos - centreline
        
#        plt.figure()
#        plt.plot(cent_x, cent_y, 'bs')
#        plt.plot(cent_x[indexes], cent_y[indexes], 'ro')
#        plt.plot([0, np.max(cent_x)], [centreline, centreline], 'k-')
#        plt.xlabel('x-pos [pixels]'); plt.ylabel('y-pos [pixels]')
#        
#        plt.figure()
#        plt.plot(cent_x[indexes]/self.PP.pPmm, offset/self.PP.pPmm, 'bs')
#        plt.xlabel('x-pos [mm]'); plt.ylabel('Offset from Centreline [mm]')
        self.intialCentreingXpos=cent_x[indexes]
        self.intialCentreingYpos=cent_y[indexes]
        
        self.aveOffset = np.mean(offset)
        self.aveOffsetSTD = np.std(offset)
        

    def _findFinalDiameter(self):
        '''find diameter of capsule in diffuser'''        
        # lets not take the last half diameter as there might be artefacts
        # from capsule tracking as the capsule enters the image
        start =  0.5*self.PP.d0*self.PP.pPmm
        
        leng=0; counter = 0
        while leng <1:        
            end = (2.0 + 0.5*counter)*self.PP.d0*self.PP.pPmm
            indexes = [i for i in range(len(self.centroid_x)) if self.centroid_x[i] < end and self.centroid_x[i] > start]
            counter += 1
            leng = len(indexes)
            
#        print('start = %d end = %d  and indexes = ' %(start, end), end=''); print(indexes)
        self.d, self.errd = self.fitCircle(indexes)
        print('Average d = %f =/- %f (based on %d outlines), counter = %d' %(self.d/self.PP.pPmm, self.errd/self.PP.pPmm, len(indexes), counter))    
        
    def _findGap(self):
        '''Find gap between capsule and wall during pinching flow'''
        indexesR = [i for i in range(len(self.centroid_x)) if (self.centroid_x[i] - self.width[i]/4.0) > self.PP.EndChannel and (self.centroid_x[i] - self.width[i]/2.0) < self.PP.RightCornerInflow]
        
        gap=[]        
        for kk in indexesR     :
            gap.append(self.centroid_y[kk] - self.height[kk]/2.0 - self.PP.ChannelTop)
            
        self.gapSpace = np.mean(gap)
        self.gapSpaceSTD = np.std(gap)
            
        if self.plot:
            fig = plt.figure(figsize=(8, 6), dpi=200); ax = fig.add_subplot(111)
            plt.plot(indexesR, np.array(gap)/self.PP.pPmm, 'cs', label='Gap width')
            
            plt.plot([np.min(indexesR), np.max(indexesR)], [self.gapSpace/self.PP.pPmm, self.gapSpace/self.PP.pPmm])

            plt.title("Gap Width -  " +self.name, fontsize=12)
            plt.xlabel('Image Number')
            plt.ylabel('Gap [mm]')
#            plt.legend(loc='best', fontsize=6, ncol=2)
            fig.tight_layout()

            plt.savefig((self.path+"Gap_" +self.name+".jpg"), dpi=300)
            
            fig = plt.figure(figsize=(8, 6), dpi=200); ax = fig.add_subplot(111)
            plt.plot(self.centroid_x, self.centroid_y, 'bs')
            plt.plot(self.centroid_x[indexesR], self.centroid_y[indexesR], 'cs', label='Gap width')
            
            from matplotlib.patches import Rectangle
            currentAxis = plt.gca()
            print('%d, %d' %(self.PP.EndChannel, self.PP.ChannelBottom))
            currentAxis.add_patch(Rectangle((self.PP.EndChannel, self.PP.ChannelBottom), 2.0*self.PP.pPmm, 4.0*self.PP.pPmm, facecolor="grey"))
            currentAxis.add_patch(Rectangle((self.PP.RightCornerInflow, self.PP.ChannelBottom), 4.0*self.PP.pPmm, 4.0*self.PP.pPmm, facecolor="grey"))
            currentAxis.add_patch(Rectangle((self.PP.EndChannel, self.PP.ChannelTop), 8.0*self.PP.pPmm, -4.0*self.PP.pPmm, facecolor="grey"))
            
            plt.title("Gap Width -  " +self.name, fontsize=12)
            plt.xlabel('X-position')
            plt.ylabel('Y-position')
#            plt.legend(loc='best', fontsize=6, ncol=2)
            fig.tight_layout()
            plt.show()
            plt.ylim( (self.PP.ChannelBottom +5, self.PP.ChannelTop - 5) )
            plt.xlim( (self.PP.End - 100, self.PP.RightCornerInflow +100) )
        
            plt.savefig((self.path+"GapCentorid_" +self.name+".jpg"), dpi=300)
        
        
    def _findRelaxationSize(self):
        '''Find dimensions of capsule after exiting channel'''
        indexesRelaxt = [i for i in range(len(self.centroid_x)) if (self.centroid_x[i] - self.width[i]/2.0) < self.PP.EndChannel + self.PP.pPmm*2.0 and (self.centroid_x[i] - self.width[i]/2.0) > self.PP.EndChannel - 8.0*self.PP.pPmm]
        
        wt = self.width[indexesRelaxt]
        ht = self.height[indexesRelaxt]
        xt = self.centroid_x[indexesRelaxt]
        yt = self.centroid_y[indexesRelaxt]

        w=[]; h=[]; x=[]; y=[];ind=[]; vx=[]; vy=[];        
        for ii in range(len(indexesRelaxt)):
            if xt[ii] > 0.0 and yt[ii] > 0.0:
                x.append(xt[ii])
                y.append(yt[ii])
                w.append(wt[ii])
                h.append(ht[ii])
                ind.append(indexesRelaxt[ii])
                vx.append(self.vel_x[indexesRelaxt[ii]])
                vy.append(self.vel_y[indexesRelaxt[ii]])
                
        self.widthRelax = np.array(w)
        self.heightRelax = np.array(h)

        if self.plot:
            fig = plt.figure(figsize=(8, 6), dpi=200); ax = fig.add_subplot(111)
            plt.plot(ind, self.widthRelax/self.PP.pPmm, 'bs', label='Width')
            plt.plot(ind, self.heightRelax/self.PP.pPmm, 'ro', label='Height')
            
            
#            ymin, ymax = ax.get_ylim()            
#            plt.plot([startInflow, startInflow], [ymin, ymax], 'r-')

            plt.title("Relaxation -  " +self.name, fontsize=12)
            plt.xlabel('Image #')
            plt.ylabel('Width / Height [mm]')
            plt.legend(loc='best', fontsize=6, ncol=2)
            fig.tight_layout()

            plt.savefig((self.path+"Relaxation_" +self.name+".jpg"), dpi=300)
            
            fig = plt.figure(figsize=(8, 6), dpi=200); ax = fig.add_subplot(111)
            plt.plot(x,y, 'bs', label='Width')

            from matplotlib.patches import Rectangle
            currentAxis = plt.gca()
            print('%d, %d' %(self.PP.EndChannel, self.PP.ChannelBottom))
            currentAxis.add_patch(Rectangle((self.PP.EndChannel, self.PP.ChannelBottom), 2.0*self.PP.pPmm, 4.0*self.PP.pPmm, facecolor="grey"))
            currentAxis.add_patch(Rectangle((self.PP.RightCornerInflow, self.PP.ChannelBottom), 4.0*self.PP.pPmm, 4.0*self.PP.pPmm, facecolor="grey"))
            currentAxis.add_patch(Rectangle((self.PP.EndChannel, self.PP.ChannelTop), 8.0*self.PP.pPmm, -4.0*self.PP.pPmm, facecolor="grey"))
#            plt.plot(ind, self.heightRelax/self.PP.pPmm, 'ro', label='Height')
                        
#            ymin, ymax = ax.get_ylim()            
#            plt.plot([startInflow, startInflow], [ymin, ymax], 'r-')
            plt.ylim( (self.PP.ChannelBottom +5, self.PP.ChannelTop - 5) )
            plt.title("Relaxation  Centroid-  " +self.name, fontsize=12)
            plt.xlabel('X-Position [pixels]')
            plt.ylabel('Y-Position [pixels]')
#            plt.legend(loc='best', fontsize=6, ncol=2)
            fig.tight_layout()

            plt.savefig((self.path+"RelaxationCentroid_" +self.name+".jpg"), dpi=300)
            
            fig = plt.figure(figsize=(8, 6), dpi=200); ax = fig.add_subplot(111)
            plt.plot(ind, np.array(vx)/self.PP.pPmm, 'cs', label='X-Velocity')
            plt.plot(ind, np.array(vy)/self.PP.pPmm, 'mo', label='Y-Velocity')

            plt.title("Relaxation  Velocity-  " +self.name, fontsize=12)
            plt.xlabel('Image Number')
            plt.ylabel('Velocity [mm/s]')
#            plt.legend(loc='best', fontsize=6, ncol=2)
            fig.tight_layout()

            plt.savefig((self.path+"RelaxationVelocity_" +self.name+".jpg"), dpi=300)
        
        
    def _writeToFilePinching(self):
        '''Write measurments of run to file'''

        self.maxd12Inflow, self.maxWidthInflow, self.minHeightInflow, self.maxd122mm, self.maxWidth2mm, self.minHeight2mm
        headerString = '\tVolume Flux 2\tFinal Distance From Centreline \t STD \tmax d12 \t maxWidth \tmin height \t Centring, \t Centring STD \t Final diameter \t Final d STD \t maxd12 Inflow \tmaxWidth Inflow \tminHeight Inflow \t meanHeight Inflow \tmaxd12 2mm \t maxWidth 2mm \tminHeight 2mm \t meanHeight 2mm \t mean_xd_Inflow \tmean_yd_Inflow \tmean_xd_2mm \tmean_yd_2mm \tmean_ydf_Inflow \tmean_ydf_2mm \tGap Height \tGap STD'
        String = '\t%.1f \t%.6f \t%.6f \t%.6f \t%d \t%d \t%.6f \t%.6f \t%.6f \t%.6f \t%.6f \t%.6f \t%.6f \t%.6f \t%.6f \t%.6f \t%.6f \t%.6f \t%.6f \t%.6f \t%.6f \t%.6f \t%.6f \t%.6f \t%.6f \t%.6f' %(self.volumFlux2, self.distanceFromCentreline, self.distanceFromCentrelineSTD, self.maxd12, self.maxWidth,  self.minHeight, self.aveOffset, self.aveOffsetSTD, self.d, self.errd, self.maxd12Inflow, self.maxWidthInflow, self.minHeightInflow, self.meanHeightInflow, self.maxd122mm, self.maxWidth2mm, self.minHeight2mm, self.meanHeight2mm, self.mean_xd_Inflow, self.mean_yd_Inflow, self.mean_xd_2mm, self.mean_yd_2mm, self.mean_ydf_Inflow, self.mean_ydf_2mm, self.gapSpace, self.gapSpaceSTD)
        
#        headerString = '\tVolume Flux 2\tFinal Distance From Centreline \t STD \tmax d12 \t maxWidth \tmin height \t Centring, \t Centring STD '
#        String = '\t%.1f \t%.6f \t%.6f \t%.6f \t%d \t%d \t%.6f \t%.6f ' %(self.volumFlux2, self.distanceFromCentreline, self.distanceFromCentrelineSTD, self.maxd12, self.maxWidth,  self.minHeight, self.aveOffset, self.aveOffsetSTD)
        
        self._writeToFile(headerString, String)
        
        
        
class ResultsClassPinching(ana.ResultsClass):
    '''Holds experiment wide results for Semi Sphere Setup'''
    
    def __init__(self, directory):
        ana.ResultsClass.__init__(self,directory)
        
        self.volumeFlux2 = np.zeros(len(self.volumeFlux))
        self.finalDistCentre = np.zeros(len(self.volumeFlux))
        self.finalDistCentreSTD = np.zeros(len(self.volumeFlux))
        self.maxD12 = np.zeros(len(self.volumeFlux))
        self.maxWidth = np.zeros(len(self.volumeFlux))
        self.minHeight = np.zeros(len(self.volumeFlux))
        self.Offset = np.zeros(len(self.volumeFlux))
        self.OffsetSTD = np.zeros(len(self.volumeFlux))
        self.Diameter = np.zeros(len(self.volumeFlux))
        self.DiameterSTD = np.zeros(len(self.volumeFlux))
        self.maxD12Inflow = np.zeros(len(self.volumeFlux))
        self.maxWidthInflow = np.zeros(len(self.volumeFlux))
        self.minHeightInflow = np.zeros(len(self.volumeFlux))
        self.meanHeightInflow = np.zeros(len(self.volumeFlux))
        self.maxD122mm = np.zeros(len(self.volumeFlux))
        self.maxWidth2mm = np.zeros(len(self.volumeFlux))
        self.minHeight2mm = np.zeros(len(self.volumeFlux))
        self.meanHeight2mm = np.zeros(len(self.volumeFlux))
        
        self.mean_xd_Inflow  = np.zeros(len(self.volumeFlux))
        self.mean_yd_Inflow  = np.zeros(len(self.volumeFlux))
        self.mean_d12xy_Inflow  = np.zeros(len(self.volumeFlux))
        
        self.mean_xd_2mm  = np.zeros(len(self.volumeFlux)) 
        self.mean_yd_2mm  = np.zeros(len(self.volumeFlux))
        self.mean_d12xy_2mm  = np.zeros(len(self.volumeFlux))
        
        self.mean_ydf_Inflow  = np.zeros(len(self.volumeFlux))
        self.mean_ydf_2mm  = np.zeros(len(self.volumeFlux))
        
        self.gapSpace = np.zeros(len(self.volumeFlux))
        self.gapSpaceSTD = np.zeros(len(self.volumeFlux))
        
        self._readDataPinching()
        
        self.PP = ParametersPinching(listingImageFunc = sortPhotosPCO,
                                   coverImgFunc = coverSidePinching,
                                   readParametersFromFile=os.path.join(self.directory, self.name[0] ))
            
    def _readDataPinching(self):
        '''Read in specific data, see above DataRunPinching._writeToFilePinching'''
        file_handle = open(self.filePath, 'r')        
        lines_list = file_handle.readlines()        # Read in all the lines of your file into a list of lines
        file_handle.close()
        
        ii=0
        linecounter=-1
        for line in lines_list:
            linecounter=linecounter+1
            if linecounter==0 : #ignore first line as it contains the header
                continue
            
            entries=line.split()

            if len(entries)>=12:
                self.volumeFlux2[ii] = (float(entries[4]))
                self.finalDistCentre[ii] = (float(entries[5]))
                self.finalDistCentreSTD[ii] = (float(entries[6]))
                self.maxD12[ii] = (float(entries[7]))
                self.maxWidth[ii] = (float(entries[8]))
                self.minHeight[ii] = (float(entries[9]))
                self.Offset[ii] = (float(entries[10]))
                self.OffsetSTD[ii] = (float(entries[11]))
                ii += 1
            if len(entries) >=28:
                self.Diameter[ii-1] = (float(entries[12]))
                self.DiameterSTD[ii-1] = (float(entries[13]))
                
                self.maxD12Inflow[ii-1] = (float(entries[14]))
                self.maxWidthInflow[ii-1] = (float(entries[15]))
                self.minHeightInflow[ii-1] = (float(entries[16]))
                self.meanHeightInflow[ii-1] = (float(entries[17]))
                self.maxD122mm[ii-1] = (float(entries[18]))
                self.maxWidth2mm[ii-1] = (float(entries[19]))
                self.minHeight2mm[ii-1] = (float(entries[20]))
                self.meanHeight2mm[ii-1] = (float(entries[21]))
                
                self.mean_xd_Inflow[ii-1] = (float(entries[22]))
                self.mean_yd_Inflow[ii-1] = (float(entries[23]))
                self.mean_d12xy_Inflow[ii-1] = abs((float(entries[22])) - (float(entries[23]))) / ((float(entries[22])) + (float(entries[23])))
                
                self.mean_xd_2mm[ii-1] = (float(entries[24]))
                self.mean_yd_2mm[ii-1] = (float(entries[25]))
                self.mean_d12xy_2mm[ii-1] = abs((float(entries[24])) - (float(entries[25]))) / ((float(entries[24])) + (float(entries[25])))
                
                self.mean_ydf_Inflow[ii-1] = (float(entries[26]))
                self.mean_ydf_2mm[ii-1] = (float(entries[27]))
                
            if len(entries)  >= 30:
                self.gapSpace[ii-1] = (float(entries[28]))
                self.gapSpaceSTD[ii-1] = (float(entries[29]))
                


    def averageByQ2(self, xInput):
        '''     Average the values for each volumn flux and return average and std     '''
        if not hasattr(self, "volumeFluxes"):
            #get all volumn fluxes
            self.volumeFluxes2=[]
            for VF in self.volumeFlux2:
                if VF not in self.volumeFluxes2:
                    self.volumeFluxes2.append(VF)
            self.volumeFluxes2 = np.array(self.volumeFluxes2)
        
        lenvfl=len(self.volumeFluxes2)
        x=np.zeros((lenvfl))
        xSTD=np.zeros((lenvfl))
        
        counter=0
        for vfl in self.volumeFluxes2:
            tempx=[]
#            print('Counter =%d' %counter)
            for i in range(len(self.volumeFlux2)):
                if self.volumeFlux2[i] == vfl:
                    tempx.append(xInput[i])
                    
            x[counter]=np.average(tempx)
            xSTD[counter]=np.std(tempx, ddof=1)
            
            counter +=1
            
        return x, xSTD
        
        
    def plotRslts(self):
        '''Top level function plotting results'''
        self._findAverages()
        self._plotQVsFinalOffset()
        self._plotQVsFinalOffsetND()
        self._plotQVsInitialOffset()
#        self._plot1oQVsOffset()
        self._plotQVsD12()
        
#        self._plotAbsFinalVsInitialOffset()
#        self._plotFinalVsInitialOffset()
#        self._plotmCvsO()
#        self._plotcCvsO()
#        self._plotQVsFinalOffsetByStreamline()
        self._plotMinHeight()
        self._plotMeanHeight()
        self._plotAveMeanHeight()
        self._plotHeightAtPos()
        self._plotQVsD12xy()
#        self._plotDiameter()
        self._plotQVsGapSpace()
        
    def _findAverages(self):
        '''Find averages of data loaded from file'''
        self.avePPmm, _ = self.averageByQ2(self.pixelsPmm) 
        self.aveFinalDistCentre, self.aveFinalDistCentreSTD = self.averageByQ2(self.finalDistCentre) 
 
        self.aveMaxD12, self.aveMaxD12STD = self.averageByQ2(self.maxD12)
        self.aveMaxWidth, self.aveMaxWidthSTD = self.averageByQ2(self.maxWidth)
        self.aveMinHeight, self.aveMinHeightSTD = self.averageByQ2(self.minHeight)
        self.aveOffset, self.aveOffsetSTD = self.averageByQ2(self.Offset)
        
        self.aveMeanHeightInflow, self.aveMeanHeightInflowSTD = self.averageByQ2(self.meanHeightInflow)
        self.aveMeanHeight2mm, self.aveMeanHeight2mmSTD = self.averageByQ2(self.meanHeight2mm)
        
        self.aveMean_yd_Inflow, self.aveMean_yd_InflowSTD = self.averageByQ2(self.meanHeightInflow)
        self.aveMean_yd_2mm, self.aveMean_yd_2mmSTD = self.averageByQ2(self.meanHeight2mm)

        self.aveMean_ydf_Inflow, self.aveMean_ydf_InflowSTD = self.averageByQ2(self.mean_ydf_Inflow)
        self.aveMean_ydf_2mm, self.aveMean_ydf_2mmSTD = self.averageByQ2(self.mean_ydf_2mm)
        
        self.aveGap, self.aveGapSTD = self.averageByQ2(self.gapSpace)
       
       
       
    def _plotQVsGapSpace(self):
        '''Plot final distance from centreline versus volume flux Q'''
        m = ana.markerSymbols(); c = ana.coloursForMarker()
        fig = plt.figure(figsize=(8, 6), dpi=200,); ax = fig.add_subplot(111)
        
        plt.errorbar(self.volumeFlux2, self.gapSpace/self.pixelsPmm, 
                     #yerr=self.finalDistCentreSTD/self.pixelsPmm,
                    linestyle='None', marker = m[0], color = c[0])
                    
        plt.errorbar(self.volumeFluxes2,  self.aveGap/self.avePPmm, 
                     yerr=self.aveGapSTD/self.avePPmm, 
                     linestyle='None', marker = m[3], markersize=8, color = c[3])

        plt.title('%s Q vs Gap Spacing' %self.batchName)
        plt.xlabel('Pinching Flow Volume Flux [$ml/min$]');   plt.ylabel('Distance betweeen Capsule and Wall [$mm$]')  
        ax = ana._standardPlotSize(ax)        
        
        ana._saveFig(fig, os.path.join(self.directory, 'RESLT') ,'%s_QvsGapSpacing.jpg' %self.batchName)
        
    def _plotQVsFinalOffset(self):
        '''Plot final distance from centreline versus volume flux Q'''
        m = ana.markerSymbols(); c = ana.coloursForMarker()
        fig = plt.figure(figsize=(8, 6), dpi=200,); ax = fig.add_subplot(111)
        
        plt.errorbar(self.volumeFlux2, self.finalDistCentre/self.pixelsPmm, 
                     #yerr=self.finalDistCentreSTD/self.pixelsPmm,
                    linestyle='None', marker = m[0], color = c[0])
                    
        plt.errorbar(self.volumeFluxes2, self.aveFinalDistCentre/self.avePPmm, 
                     yerr=self.aveFinalDistCentreSTD/self.avePPmm, 
                     linestyle='None', marker = m[3], markersize=8, color = c[3])
                    
        plt.title('%s Q vs Offset' %self.batchName)
        plt.xlabel('Pinching Flow Volume Flux [$ml/min$]');   plt.ylabel('Final Offset from centreline [$mm$]')  
        ax = ana._standardPlotSize(ax)        
        
        ana._saveFig(fig, os.path.join(self.directory, 'RESLT') ,'%s_QvsOffset.jpg' %self.batchName)
        
    def _plotQVsFinalOffsetND(self):
        '''Plot final distance from centreline versus volume flux Q'''
        m = ana.markerSymbols(); c = ana.coloursForMarker()
        fig = plt.figure(figsize=(8, 6), dpi=200,); ax = fig.add_subplot(111)
        
        plt.errorbar(self.volumeFlux2, (self.finalDistCentre/self.pixelsPmm)/self.PP.d0, 
                     #yerr=self.finalDistCentreSTD/self.pixelsPmm,
                    linestyle='None', marker = m[0], color = c[0])
                    
        plt.errorbar(self.volumeFluxes2, (self.aveFinalDistCentre/self.avePPmm)/self.PP.d0, 
                     yerr=(self.aveFinalDistCentreSTD/self.avePPmm)/self.PP.d0, 
                     linestyle='None', marker = m[3], markersize=8, color = c[3])
                    
        plt.title('%s Q vs Offset' %self.batchName)
        plt.xlabel('Pinching Flow Volume Flux [$ml/min$]');   plt.ylabel('Final Offset from Centreline (in diameters) [$d_0$]')  
        ax = ana._standardPlotSize(ax)        
        
        ana._saveFig(fig, os.path.join(self.directory, 'RESLT') ,'%s_QvsOffsetND.jpg' %self.batchName)
        
    def _plotMinHeight(self):
        '''Plot min height of capsules'''
        m = ana.markerSymbols(); c = ana.coloursForMarker()
        fig = plt.figure(figsize=(8, 6), dpi=200,); ax = fig.add_subplot(111)
        
        plt.errorbar(self.volumeFlux2, self.minHeight/self.pixelsPmm, 
                     yerr=1.0/self.pixelsPmm, 
                     label= 'Runs +/- 1 pixel ',
                    linestyle='None', marker = m[0], color = c[0])
                    
        plt.errorbar(self.volumeFluxes2, self.aveMinHeight/self.avePPmm, 
                     yerr=self.aveMinHeightSTD/self.avePPmm, 
                     label= 'Average +/- Standard Deviation',
                     linestyle='None', marker = m[3], markersize=8, color = c[3])

        plt.title('%s Q vs Min Height' %self.batchName)
        plt.xlabel('Pinching Flow Volume Flux [$ml/min$]');   plt.ylabel('Min Height of Capsule [$mm$]')  
        plt.legend(loc='best', fontsize=8)
        ax = ana._standardPlotSize(ax)        
        
        ana._saveFig(fig, os.path.join(self.directory, 'RESLT') ,'%s_QvsMinHeight.jpg' %self.batchName)
        
    def _plotMeanHeight(self):
        '''Plot mean height of capsules'''
        m = ana.markerSymbols(); c = ana.coloursForMarker()
        fig = plt.figure(figsize=(8, 6), dpi=200,); ax = fig.add_subplot(111)

        plt.errorbar(self.volumeFlux2, self.meanHeightInflow/self.pixelsPmm, 
                     yerr=1.0/self.pixelsPmm, 
                     label= 'Mean Height at Inflow ',
                    linestyle='None', marker = m[0], color = c[0])
                    
        plt.errorbar(self.volumeFlux2, self.meanHeight2mm/self.pixelsPmm, 
                     yerr=1.0/self.pixelsPmm, 
                     label= 'Mean Height in 2mm Exit ',
                    linestyle='None', marker = m[3], color = c[3])
                    

#        newax.axis('off')
                    
#        plt.errorbar(self.volumeFluxes2, self.aveMinHeight/self.avePPmm, 
#                     yerr=self.aveMinHeightSTD/self.avePPmm, 
#                     label= 'Average +/- Standard Deviation',
#                     linestyle='None', marker = m[3], markersize=8, color = c[3])

        plt.title('%s Q vs Mean Height' %self.batchName)
        plt.xlabel('Pinching Flow Volume Flux [$ml/min$]');   plt.ylabel('Mean Height of Capsule [$mm$]')  
        plt.legend(loc='best', fontsize=8)
        ax = ana._standardPlotSize(ax)      
        
        
#
#        im = plt.imread('C:\\Users\\mbbxkeh2\\Dropbox\\PhD\\Sketch_Pinching.png')
#        from matplotlib.offsetbox import TextArea, DrawingArea, OffsetImage,  AnnotationBbox
#        xy = (0.5, 0.8)
#        img = OffsetImage(im, zoom=0.3)
#        ab = AnnotationBbox(img, xy,
#                        xybox=xy,
#                        xycoords="axes fraction",
#                        boxcoords="axes fraction",
#                        pad=0.1,
#                        arrowprops=dict(arrowstyle="->"))
#
#        ax.add_artist(ab)
        
        
        ana._saveFig(fig, os.path.join(self.directory, 'RESLT') ,'%s_QvsMeanHeight.jpg' %self.batchName)
        
    def _plotAveMeanHeight(self):
        '''Plot mean height of capsules'''
        m = ana.markerSymbols(); c = ana.coloursForMarker()
        fig = plt.figure(figsize=(8, 6), dpi=200,); ax = fig.add_subplot(111)

        plt.errorbar(self.volumeFluxes2, self.aveMeanHeightInflow/self.avePPmm, 
                     yerr=self.aveMeanHeightInflowSTD/self.avePPmm, 
                     label= 'Mean Height at Inflow',
                    linestyle='None', marker = m[1], color = c[1])
                    
        plt.errorbar(self.volumeFluxes2, self.aveMeanHeight2mm/self.avePPmm, 
                     yerr=self.aveMeanHeight2mmSTD/self.avePPmm, 
                     label= 'Mean Height in 2mm Exit',
                    linestyle='None', marker = m[3], color = c[3])
                    
        plt.errorbar(self.volumeFluxes2, self.aveMean_yd_Inflow/self.avePPmm, 
                     yerr=self.aveMean_yd_InflowSTD/self.avePPmm, 
                     label= 'Mean Height at Inflow from contour',
                    linestyle='None', marker = m[0], color = c[0])
                    
        plt.errorbar(self.volumeFluxes2, self.aveMean_yd_2mm/self.avePPmm, 
                     yerr=self.aveMean_yd_2mmSTD/self.avePPmm, 
                     label= 'Mean Height in 2mm Exit from contour',
                    linestyle='None', marker = m[2], color = c[2])
                    
#        plt.errorbar(self.volumeFluxes2, self.aveMinHeight/self.avePPmm, 
#                     yerr=self.aveMinHeightSTD/self.avePPmm, 
#                     label= 'Average +/- Standard Deviation',
#                     linestyle='None', marker = m[3], markersize=8, color = c[3])

        plt.title('%s Q vs Ave Mean Height' %self.batchName)
        plt.xlabel('Pinching Flow Volume Flux [$ml/min$]');   plt.ylabel('Ave Mean Height of Capsule [$mm$]')  
        plt.legend(loc='best', fontsize=8)
        ax = ana._standardPlotSize(ax)        
        
        ana._saveFig(fig, os.path.join(self.directory, 'RESLT') ,'%s_QvsAveMeanHeight.jpg' %self.batchName)

    def _plotHeightAtPos(self):
        m = ana.markerSymbols(); c = ana.coloursForMarker()
        fig = plt.figure(figsize=(8, 6), dpi=200,); ax = fig.add_subplot(111)

        plt.errorbar(self.volumeFlux2, self.mean_ydf_Inflow/self.pixelsPmm, 
                     yerr=1.0/self.pixelsPmm, 
                     label= 'Mean Height at Inflow ',
                    linestyle='None', marker = m[0], color = c[0])
                    
        plt.errorbar(self.volumeFlux2, self.mean_ydf_2mm/self.pixelsPmm, 
                     yerr=1.0/self.pixelsPmm, 
                     label= 'Mean Height in 2mm Exit ',
                    linestyle='None', marker = m[3], color = c[3])

        plt.errorbar(self.volumeFluxes2, self.aveMean_ydf_Inflow/self.avePPmm, 
                     yerr=self.aveMean_ydf_InflowSTD/self.avePPmm, 
                     label= 'Average Mean Height at Inflow',
                    linestyle='None', marker = m[1], color = c[1])
                    
        plt.errorbar(self.volumeFluxes2, self.aveMean_ydf_2mm/self.avePPmm, 
                     yerr=self.aveMean_ydf_2mmSTD/self.avePPmm, 
                     label= 'AverageMean Height in 2mm Exit',
                    linestyle='None', marker = m[2], color = c[2])
                    
        plt.title('%s Q vs Ave Mean Height at position' %self.batchName)
        plt.xlabel('Pinching Flow Volume Flux [$ml/min$]');   plt.ylabel('Ave Mean Height of Capsule [$mm$]')  
        plt.legend(loc='best', fontsize=8)
        ax = ana._standardPlotSize(ax)        
        
        ana._saveFig(fig, os.path.join(self.directory, 'RESLT') ,'%s_QvsAveMeanHeightPos.jpg' %self.batchName)
        
    
    
    def _plotDiameter(self):
        '''Plot min height of capsules'''
        m = ana.markerSymbols(); c = ana.coloursForMarker()
        fig = plt.figure(figsize=(8, 6), dpi=200,); ax = fig.add_subplot(111)
        
        plt.errorbar(self.volumeFlux2, self.Diameter/self.pixelsPmm, 
                     yerr=self.DiameterSTD/self.pixelsPmm, 
                     label= 'Runs +/- 1 pixel ',
                    linestyle='None', marker = m[0], color = c[0])
                    
        vf2=[]; d=[]; dSTD=[]; tol=0.2;        
        
        for ii in range(len(self.Diameter)):
            td = self.Diameter[ii]
            tdSTD = self.DiameterSTD[ii]
            
            if abs(self.PP.d0*self.pixelsPmm[0] - td)/(self.PP.d0*self.pixelsPmm[0]) < tol:
                vf2.append(self.volumeFlux2[ii])
                d.append(td/self.pixelsPmm[0])
                dSTD.append(tdSTD/self.pixelsPmm[0])
                
        d = np.array(d); dSTD = np.array(dSTD); vf2 = np.array(vf2);
        
        plt.errorbar(vf2, d,yerr=dSTD, 
                     label= 'Runs within %.2f of Mean' %tol,
                    linestyle='None', marker = m[3], color = c[3])
            
#        print(d)
        dMean = np.mean(d)
#        print('dMean = %f' %dMean)
        textstring = 'Average Diameter = %f' %(dMean)
        plt.text(0.1, 0.5, textstring, fontsize = 12, 
                     horizontalalignment='left', verticalalignment='center', 
                     transform = ax.transAxes)
            
        xmin, xmax = ax.get_xlim()
        plt.plot([xmin, xmax], [dMean, dMean], 'r-', label='Average')

        plt.title('%s Q vs Diameter in Diffuser' %self.batchName)
        plt.xlabel('Pinching Flow Volume Flux [$ml/min$]');   plt.ylabel('Diameter Capsule (topview) [$mm$]')  
        plt.legend(loc='best', fontsize=8)
        
        ax = ana._standardPlotSize(ax)        
        
        ana._saveFig(fig, os.path.join(self.directory, 'RESLT') ,'%s_QvsDiameter.jpg' %self.batchName)
        
    def _sortByCentering(self, width):
        '''Select the value within width of centreline'''
        vf=[]; finalDist=[]; 
        for ii in range(len(self.volumeFlux2)):
            if np.abs(self.Offset[ii]/self.pixelsPmm[ii]) < width:
                vf.append(self.volumeFlux2[ii])
                finalDist.append(self.finalDistCentre[ii]/self.pixelsPmm[ii])
        
        return vf, finalDist
        
    def _plotQVsFinalOffsetByStreamline(self):
        '''Plot final distance from centreline versus volume flux Q'''
        m = ana.markerSymbols(); c = ana.coloursForMarker()
        fig = plt.figure(figsize=(8, 6), dpi=200,); ax = fig.add_subplot(111)
        
        vf1, finalDist1 =self._sortByCentering(0.1)
        vf2, finalDist2 =self._sortByCentering(0.2)
        vf3, finalDist3 =self._sortByCentering(0.3)
        plt.errorbar(self.volumeFlux2, self.finalDistCentre/self.pixelsPmm, 
                     #yerr=self.finalDistCentreSTD/self.pixelsPmm,
                    linestyle='None', marker = m[0], color = c[0], 
                    markersize=4, label = 'All Data')
                    
#        plt.errorbar(vf3, finalDist3, 
#                     #yerr=self.finalDistCentreSTD/self.pixelsPmm,
#                    linestyle='None', marker = m[1], color = c[1],
#                    label = 'Within 0.3 mm of centreline' )
                    
        plt.errorbar(vf2, finalDist2, 
                     #yerr=self.finalDistCentreSTD/self.pixelsPmm,
                    linestyle='None', marker = m[2], color = c[2],
                    markersize=6, label = 'Within 0.2 mm of centreline' )
                    
#        plt.errorbar(vf1, finalDist1, 
#                     #yerr=self.finalDistCentreSTD/self.pixelsPmm,
#                    linestyle='None', marker = m[3], color = c[3],
#                    label = 'Within 0.1 mm of centreline' )
                    
        plt.title('%s Q vs Offset' %self.batchName)
        plt.xlabel('Volume Flux [$ml/min$]');   plt.ylabel('Final Offset [$mm$]')  
        ax = ana._standardPlotSize(ax)        
        plt.legend(loc='best', fontsize=8, ncol=1)
        ana._saveFig(fig, os.path.join(self.directory, 'RESLT') ,'%s_QvsOffset_centeringSelected.jpg' %self.batchName)
        
    def _plot1oQVsOffset(self):
        '''Plot final distance from centreline versus volume flux Q'''
        m = ana.markerSymbols(); c = ana.coloursForMarker()
        fig = plt.figure(figsize=(8, 6), dpi=200,); ax = fig.add_subplot(111)
        
        plt.errorbar(1.0/self.volumeFlux2, self.finalDistCentre/self.pixelsPmm, 
                     #yerr=self.finalDistCentreSTD/self.pixelsPmm,
                    linestyle='None', marker = m[0], color = c[0])
                    
        plt.errorbar(1.0/self.volumeFluxes2, self.aveFinalDistCentre/self.avePPmm, 
                     yerr=self.aveFinalDistCentreSTD/self.avePPmm, 
                     linestyle='None', marker = m[3], markersize=8, color = c[3])
                    
        plt.title('%s 1/Q vs Offset' %self.batchName)
        plt.xlabel('1/Volume Flux [$min/ml$]');   plt.ylabel('Final Offset [$mm$]')  
        ax = ana._standardPlotSize(ax)        
        
        ana._saveFig(fig, os.path.join(self.directory, 'RESLT') ,'%s_1oQvsOffset.jpg' %self.batchName)
        
    def _plotQVsD12(self):
        '''Plot max taylor deformation versus volume flux Q'''
        m = ana.markerSymbols(); c = ana.coloursForMarker()
        fig = plt.figure(figsize=(8, 6), dpi=200,); ax = fig.add_subplot(111)
        
        plt.plot(self.volumeFlux2, self.maxD12, 
                    linestyle='None', marker = m[0], color = c[0])
        plt.xlabel('Volume Flux [$ml/min$]');   
        plt.ylabel('Taylor deformation parameter $D_{12}$')  
        
        ax.spines['left'].set_color(c[0])
        ax.tick_params(axis='y', colors=c[0])
        ax.yaxis.label.set_color(c[0])
        
        locs, labels =plt.xticks()        
        ax2 = ax.twinx()
        maxW = (self.maxWidth/self.pixelsPmm) / self.PP.d0
        ax2.plot(self.volumeFlux2, maxW,
                linestyle='None', marker = m[3], color = c[3])
        ax2.set_ylabel('Max Width [$d_0$]')
        
        ax2.spines['left'].set_color(c[3])
        ax2.tick_params(axis='y', colors=c[3])
        ax2.yaxis.label.set_color(c[3])
        
        plt.title('%s Q vs $D_{12}$' %self.batchName)
        
        ax = ana._standardPlotSize(ax)        
        
        ana._saveFig(fig, os.path.join(self.directory, 'RESLT') ,'%s_QvsD12.jpg' %self.batchName)
        
    def _plotQVsD12xy(self):
        '''Plot max taylor deformation versus volume flux Q'''
        m = ana.markerSymbols(); c = ana.coloursForMarker()
        fig = plt.figure(figsize=(8, 6), dpi=200,); ax = fig.add_subplot(111)

        plt.errorbar(self.volumeFlux2, self.mean_d12xy_Inflow, 
#                     yerr=self.aveMeanHeightInflowSTD/self.avePPmm, 
                     label= 'Mean $D_{12}$ at Inflow',
                    linestyle='None', marker = m[1], color = c[1])
                    
        plt.errorbar(self.volumeFlux2, self.mean_d12xy_2mm,
#                     yerr=self.aveMeanHeight2mmSTD/self.avePPmm, 
                     label= 'Mean $D_{12}$ in 2mm Exit',
                    linestyle='None', marker = m[3], color = c[3])
                    
#        plt.errorbar(self.volumeFluxes2, self.aveMinHeight/self.avePPmm, 
#                     yerr=self.aveMinHeightSTD/self.avePPmm, 
#                     label= 'Average +/- Standard Deviation',
#                     linestyle='None', marker = m[3], markersize=8, color = c[3])

        plt.title('%s Q vs $D_{12}$ xy' %self.batchName)
        plt.xlabel('Pinching Flow Volume Flux [$ml/min$]');   plt.ylabel('$D_{12}$ xy of Capsule')  
        plt.legend(loc='best', fontsize=8)
        ax = ana._standardPlotSize(ax)        
        
        ana._saveFig(fig, os.path.join(self.directory, 'RESLT') ,'%s_QvsD12xy.jpg' %self.batchName)
        
    
    
    def _plotAbsFinalVsInitialOffset(self):
        '''Plot initial offset from centreline in the channel leading to the 
        semi-sphere versus volume flux Q'''      
        
        p1, p2, fres = self._linearFit(np.abs(self.Offset/self.pixelsPmm),self.finalDistCentre/self.pixelsPmm)
        
        m = ana.markerSymbols(); c = ana.coloursForMarker()
        fig = plt.figure(figsize=(8, 6), dpi=200,); ax = fig.add_subplot(111)
        
        plt.errorbar(np.abs(self.Offset/self.pixelsPmm), self.finalDistCentre/self.pixelsPmm,
                    #yerr=self.finalDistCentreSTD/self.pixelsPmm, 
                    #xerr=self.OffsetSTD/self.pixelsPmm,
                    linestyle='None', marker = m[0], color = c[0], label='Runs')
                    
        plt.plot([0.0, np.max(np.abs(self.Offset/self.pixelsPmm))], [ana.func(0.0, p1,p2), ana.func(np.max(np.abs(self.Offset/self.pixelsPmm)), p1,p2)], 
                  'k-', label='Linear Best Fit with m= %.4f and c = %.4f (fres = %.4f)' %(p1, p2, fres))
        
        plt.legend(loc='best', fontsize=6, ncol=1)
        plt.title('%s Initial vs Final Offset' %self.batchName)
        plt.ylabel('(Absolute) Final Offset [$mm$]');   plt.xlabel('(Absolute) Initial Centering [$mm$]')  
        
        textstring = 'Number of Runs = %d' %(len(self.Offset))
        plt.text(0.1, 0.1, textstring, fontsize = 12, 
                     horizontalalignment='left', verticalalignment='center', 
                     transform = ax.transAxes)
        ax = ana._standardPlotSize(ax)        
        
        ana._saveFig(fig, os.path.join(self.directory, 'RESLT') ,'%s_AbsInitiaVsFinalOffset.jpg' %self.batchName)
        
        counter=0
        self.mCentreingOffset=[]
        self.cCentreingOffset=[]
        self.merrCentreingOffset = []
        self.cerrCentreingOffset = []
        self.msigmaCentreingOffset = []
        self.csigmaCentreingOffset = []
        for vf in self.volumeFluxes2:
            fig = plt.figure(figsize=(8, 6), dpi=200,); ax = fig.add_subplot(111)      
            tempx=[]; tempy=[]
            for i in range(len(self.volumeFlux2)):
                if self.volumeFlux2[i] == vf:
                    tempx.append(self.Offset[i]/self.pixelsPmm[i])
                    tempy.append(self.finalDistCentre[i]/self.pixelsPmm[i])
            tempx = np.array(tempx); tempy = np.array(tempy)
            
            my, by, ry, errmy, errby   = ana.lsqfity(np.abs(tempx),tempy)
            sigmamy2 = stats.norm.interval(0.68, loc=my, scale=errmy)
            sigmamy = stats.norm.interval(0.68, loc=my, scale=errmy/np.sqrt(len(tempx)))
            sigmaby = stats.norm.interval(0.68, loc=my, scale=errby/np.sqrt(len(tempx)))
            print('confidence interval sample = %f,%f population = %f, %f \t sterrs = %f, %f' %(sigmamy[0],sigmamy[1],sigmamy2[0], sigmamy2[1], errmy, errby))
#            print(sigmamy)
            self.mCentreingOffset.append(my);self.cCentreingOffset.append(by); 
            self.merrCentreingOffset.append(errmy);self.cerrCentreingOffset.append(errby); 
            self.msigmaCentreingOffset.append(sigmamy);self.csigmaCentreingOffset.append(sigmaby); 
            
            plt.errorbar(np.abs(tempx), tempy,
                        #yerr=self.finalDistCentreSTD/self.pixelsPmm, 
                        #xerr=self.OffsetSTD/self.pixelsPmm,
                        linestyle='None', marker = m[0], color = c[0])
                        
            labelString = 'Linear Best Fit with m= %.4f and c = %.4f (RSS = %.4f)' %(my, by, ry)

            plt.plot([0.0, np.max(np.abs(tempx))], [ana.func(0.0, my,by), 
                      ana.func(np.max(np.abs(tempx)), my,by)], 'k-', 
                    label=labelString)
        
            plt.legend(loc='best', fontsize=6, ncol=1)

            textstring = 'Number of Runs = %d' %(len(tempx))
            plt.text(0.1, 0.6, textstring, fontsize = 16, 
                     horizontalalignment='left', verticalalignment='center', 
                     transform = ax.transAxes)
            plt.title('%s Initial vs Final Offset, Q = %.f' %(self.batchName, vf))
            plt.ylabel('(Absolute) Final Offset [$mm$]');   plt.xlabel('(Absolute) Initial Centering [$mm$]')  
            ax = ana._standardPlotSize(ax)        
            
            ana._saveFig(fig, os.path.join(self.directory, 'RESLT') ,'%s_Q=%.f_AbsInitiaVsFinalOffset.jpg' %(self.batchName, vf))
            counter +=1
    
    def _plotmCvsO(self):
        ''' Plot the gradient of the linear fit to the initial centreing against
        final offset'''
        m = ana.markerSymbols(); c = ana.coloursForMarker()
        fig = plt.figure(figsize=(8, 6), dpi=200,); ax = fig.add_subplot(111)
        
#        print(len(self.volumeFluxes2))
#        print(len(self.mCentreingOffset))
#        print(len())
#        upperCI=[]        
#        lowerCI=[]  
#        for sg in self.msigmaCentreingOffset:
#            upperCI.append(sg[0])
#            lowerCI.append(sg[1])
#            
#        plt.errorbar(self.volumeFluxes, self.mCentreingOffset, 
#                     yerr = [upperCI,lowerCI],
#                     linestyle='None', marker = m[0], color = c[0])
#        plt.text(0.1, 0.1, 'Where the errorbars indicate one $\sigma$, i.e. the 68% confidence interval.', fontsize = 12, 
#             horizontalalignment='left', verticalalignment='center', 
#             transform = ax.transAxes)
             
        plt.errorbar(self.volumeFluxes2, self.mCentreingOffset, 
                     yerr = self.merrCentreingOffset,
                     linestyle='None', marker = m[0], color = c[0])
        plt.text(0.1, 0.1, 'Where the errorbars indicate the standard error', fontsize = 12, 
             horizontalalignment='left', verticalalignment='center', 
             transform = ax.transAxes)
             
#        plt.errorbar(self.aveOffset/self.avePPmm, self.aveFinalDistCentre/self.avePPmm,
#                     yerr = self.aveFinalDistCentreSTD/self.avePPmm, 
#                     xerr=self.aveOffsetSTD/self.avePPmm, 
#                     linestyle='None', marker = m[3], markersize=8, color = c[3])
                    
        plt.title('%s Gradient Initial vs Final Offset' %self.batchName)
        plt.xlabel('Volume Flux [$ml/min$]');   plt.ylabel('Gradient')  
        
#        textstring = 'Number of Runs = %d' %(len(self.Offset))
#        plt.text(0.1, 0.1, textstring, fontsize = 12, 
#                     horizontalalignment='left', verticalalignment='center', 
#                     transform = ax.transAxes)

        ax = ana._standardPlotSize(ax)        
        
        ana._saveFig(fig, os.path.join(self.directory, 'RESLT') ,'%s_GradientInitiaVsFinalOffset.jpg' %self.batchName)
        
    def _plotcCvsO(self):
        ''' Plot the gradient of the linear fit to the initial centreing against
        final offset'''
        m = ana.markerSymbols(); c = ana.coloursForMarker()
        fig = plt.figure(figsize=(8, 6), dpi=200,); ax = fig.add_subplot(111)
        
#        print(len(self.volumeFluxes))
#        print(len(self.mCentreingOffset))
#        print(len())
        plt.errorbar(self.volumeFluxes2, self.cCentreingOffset, 
                     yerr = self.cerrCentreingOffset,
                     linestyle='None', marker = m[0], color = c[0])
        plt.text(0.1, 0.1, 'Where the errorbars indicate the standard error', fontsize = 12, 
                     horizontalalignment='left', verticalalignment='center', 
                     transform = ax.transAxes)
                     
#        plt.errorbar(self.aveOffset/self.avePPmm, self.aveFinalDistCentre/self.avePPmm,
#                     yerr = self.aveFinalDistCentreSTD/self.avePPmm, 
#                     xerr=self.aveOffsetSTD/self.avePPmm, 
#                     linestyle='None', marker = m[3], markersize=8, color = c[3])
                    
        plt.title('%s Intercept Initial vs Final Offset' %self.batchName)
        plt.xlabel('Volume Flux [$ml/min$]');   plt.ylabel('Intercept [$mm$]')  
        
#        textstring = 'Number of Runs = %d' %(len(self.Offset))

        ax = ana._standardPlotSize(ax)        
        
        ana._saveFig(fig, os.path.join(self.directory, 'RESLT') ,'%s_InterceptInitiaVsFinalOffset.jpg' %self.batchName)
    
    def _plotFinalVsInitialOffset(self):
        '''Plot initial offset from centreline in the channel leading to the 
        semi-sphere versus volume flux Q'''
        m = ana.markerSymbols(); c = ana.coloursForMarker()
        fig = plt.figure(figsize=(8, 6), dpi=200,); ax = fig.add_subplot(111)
        
        plt.errorbar(self.Offset/self.pixelsPmm, self.finalDistCentre/self.pixelsPmm,
                    #yerr=self.finalDistCentreSTD/self.pixelsPmm, 
                    #xerr=self.OffsetSTD/self.pixelsPmm,
                    linestyle='None', marker = m[0], color = c[0])
                    
        plt.errorbar(self.aveOffset/self.avePPmm, self.aveFinalDistCentre/self.avePPmm,
                     yerr = self.aveFinalDistCentreSTD/self.avePPmm, 
                     xerr=self.aveOffsetSTD/self.avePPmm, 
                     linestyle='None', marker = m[3], markersize=8, color = c[3])
                    
        plt.title('%s Initial vs Final Offset' %self.batchName)
        plt.ylabel('Final Offset [$mm$]');   plt.xlabel('Initial Centering [$mm$]')  
        
        textstring = 'Number of Runs = %d' %(len(self.Offset))
        plt.text(0.1, 0.1, textstring, fontsize = 12, 
                     horizontalalignment='left', verticalalignment='center', 
                     transform = ax.transAxes)
        ax = ana._standardPlotSize(ax)        
        
        ana._saveFig(fig, os.path.join(self.directory, 'RESLT') ,'%s_InitiaVsFinalOffset.jpg' %self.batchName)
        counter=0
        for vf in self.volumeFluxes2:
            fig = plt.figure(figsize=(8, 6), dpi=200,); ax = fig.add_subplot(111)      
            tempx=[]; tempy=[]
            for i in range(len(self.volumeFlux2)):
                if self.volumeFlux2[i] == vf:
                    tempx.append(self.Offset[i]/self.pixelsPmm[i])
                    tempy.append(self.finalDistCentre[i]/self.pixelsPmm[i])
            plt.errorbar(tempx, tempy,
                        #yerr=self.finalDistCentreSTD/self.pixelsPmm, 
                        #xerr=self.OffsetSTD/self.pixelsPmm,
                        linestyle='None', marker = m[0], color = c[0])
            
#            tempx2=[]; tempy2=[]; erry=[]; errx=[]
#            for i in range(len(self.volumeFluxes)):
#                if self.volumeFluxes[i] == vf:
#                    tempx2.append(self.Offset[i]/self.pixelsPmm[i])
#                    tempy2.append(self.finalDistCentre[i]/self.pixelsPmm[i])
                    
            plt.errorbar(self.aveOffset[counter]/self.avePPmm[counter], 
                         self.aveFinalDistCentre[counter]/self.avePPmm[counter],
                         yerr = self.aveFinalDistCentreSTD[counter]/self.avePPmm[counter], 
                         xerr=self.aveOffsetSTD[counter]/self.avePPmm[counter], 
                         linestyle='None', marker = m[3], markersize=8, color = c[3])
            textstring = 'Number of Runs = %d' %(len(tempx))
            plt.text(0.1, 0.6, textstring, fontsize = 16, 
                     horizontalalignment='left', verticalalignment='center', 
                     transform = ax.transAxes)
            plt.title('%s Initial vs Final Offset, Q = %.f' %(self.batchName, vf))
            plt.ylabel('Final Offset [$mm$]');   plt.xlabel('Initial Centering [$mm$]')  
            ax = ana._standardPlotSize(ax)        
            
            ana._saveFig(fig, os.path.join(self.directory, 'RESLT') ,'%s_Q=%.f_InitiaVsFinalOffset.jpg' %(self.batchName, vf))
            counter +=1
            
    def _plotQVsInitialOffset(self):
        '''Plot final versus initial offset'''
        m = ana.markerSymbols(); c = ana.coloursForMarker()
        fig = plt.figure(figsize=(8, 6), dpi=200,); ax = fig.add_subplot(111)
        
        plt.errorbar(self.volumeFlux2, self.Offset/self.pixelsPmm, 
                     #yerr=self.OffsetSTD/self.pixelsPmm,
                    linestyle='None', marker = m[0], color = c[0])
                    
        plt.errorbar(self.volumeFluxes2, self.aveOffset/self.avePPmm, 
                     yerr=self.aveOffsetSTD/self.avePPmm, 
                     linestyle='None', marker = m[3], markersize=8, color = c[3])
                    
        plt.title('%s Q vs Centreing' %self.batchName)
        plt.xlabel('Volume Flux [$ml/min$]');   plt.ylabel('Initial Centering [$mm$]')  
        ax = ana._standardPlotSize(ax)        
        
        ana._saveFig(fig, os.path.join(self.directory, 'RESLT') ,'%s_QvsCentering.jpg' %self.batchName)
        
    

        
        
        
        
        