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
import pinching as pin
#from analysis import DataRun

ERROR_CONST=-1

class ParametersPinchingMk2(gen.Parameters):
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
    use sort function in pinching file
    """
    filenameList, numberOfJPG = pin.sortPhotosPCO(path, fileType, prefixleng)
    return filenameList, numberOfJPG
            
            
def coverSidePinching(img, offset, ParameterClass):
    '''Take image and delete areas not of intereste'''

    top = ParameterClass.ChannelTop
    bottom = ParameterClass.ChannelBottom 
    end = ParameterClass.EndChannel
    lc = ParameterClass.LeftCornerInflow
    rc = ParameterClass.RightCornerInflow
    
    width = (bottom - top)  #4mm
    ydown = bottom + width
    #Get size of image
    yImg,xImg = img.shape
    
#    cv2.rectangle(img, (x1, y1), (x2, y2), (255,0,0), 2)  
#    x1,y1 ------
#    |          |
#    |          |
#    |          |
#    --------x2,y2
    #cover area above and below main channel
    cv2.rectangle(img, (end+offset, top+offset), (0, 0), 255, thickness=-1)   
    
    cv2.rectangle(img, (end+offset, bottom+offset), (rc+offset, ydown+offset),  255, thickness=-1)
    cv2.rectangle(img, (lc+offset, bottom+offset), (0, ydown+offset),  255, thickness=-1)
    
    shapeLeft = np.array([ [lc+offset,ydown+offset], [lc+offset - int(1.5*width) ,ydown+offset+2*width], [lc+offset - int(1.5*width) ,yImg], [0 ,yImg], [0, ydown+offset ] ], np.int32)
    cv2.fillPoly(img, [shapeLeft], (255,255,255))

    shapeR = np.array([[rc+offset,bottom+offset], 
                       [rc+offset,ydown+offset], 
                       [rc+offset - width, ydown+offset+2*width],
                       [rc+offset - width, yImg], 
                       [end+offset ,yImg], 
                       [end+offset, bottom+offset ] ], np.int32)
    cv2.fillPoly(img, [shapeR], (255,255,255))    
    
    #Cover outside
    side=5
    cv2.rectangle(img, (0, offset+side), (xImg, 0), 255, thickness=-1)   
    cv2.rectangle(img, (0, yImg - offset-side), (xImg, yImg), 255, thickness=-1)   
    cv2.rectangle(img, (offset+side, 0), (0, yImg), 255, thickness=-1)   
    cv2.rectangle(img, (xImg - offset-side, 0), (xImg, yImg), 255, thickness=-1)   
        
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
#    DRSS.plot = False
    DRSS.extractMeasures()
#    DRSS = DataRunPinching(path)
#    DRSS._writeToFilePinching()
#    DRSS.writeNameToFile()
#    del DRSS
    
def runFunction_TJMode(path):
    '''Runs DataRunPinching on path'''
#    print('Starting runFunction')
    DRSS = DataRunPinching(path)
#    DRSS.plot = False
    DRSS.extractMeasures_TJ_Mode()
#    DRSS = DataRunPinching(path)
#    DRSS._writeToFilePinching()
#    DRSS.writeNameToFile()
#    del DRSS
    
class DataRunPinching(ana.DataRun):
    
    def __init__(self, path):
        #load parameters from image analysis
        self.PP = ParametersPinchingMk2(listingImageFunc = sortPhotosPCO,
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
        '''Extract measures of interest from the raw tracking data self.distanceFromCentrelineand write to 
        file. '''
        self._findGap()
        self._findCentering()
        self._getFinalPosition()        
        
        self._maxDeformation()        
#
#        self.d=-1; self.errd = -1
#        self._findRelaxationSize()
        self._findGap()
        self._writeToFilePinching()
        
    def extractMeasures_TJ_Mode(self):
        '''Extract measures of interest from the raw tracking data self.distanceFromCentrelineand write to 
        file. '''
        self._findGap()
        self._findCentering_TJ_Mode()
        self._getFinalPosition(TJ_Mode=True)        
        
        self._maxDeformation()        
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
        
    def _getFinalPosition(self, TJ_Mode=False, forPub=False, forPubpath=''):
        '''Finds final position of capsule after the diffuser and check whether
        capsule went right (up) or left (down) (from the capsules point of view.'''
        #Central line can be evaluated from PP.ChannelTop  and PP.ChannelBottom
        
        centreline = (self.PP.ChannelTop +  self.PP.ChannelBottom)/2.0
        
        # Somewhat arbitrary end position to measure: 6 diameters.append(interp1d(x, y))
        
        if TJ_Mode:
#            print('Using t-Junction Mode final position')
            xposAfterDiffuser = self.PP.EndChannel + \
                (self.PP.d0*5.0)*self.PP.pPmm
        else:
            xposAfterDiffuser = self.PP.EndChannel + \
                (self.PP.d0*6.0)*self.PP.pPmm

        #remove entries with no values
        cent_x, cent_y = self._centroidsWithoutErrors()
        ri = int(len(cent_x)/2.0) #relevant indexes are the second half
        cent_xs, cent_ys =cent_x[-ri:], cent_y[-ri:] 
#        print("cent_xs = ", end = '')
#        print(cent_xs)
        print('xposAfterDiffuser = %d' %xposAfterDiffuser)

        #find max x position
        maxX = np.max(cent_x)
        
        #ignore half a diameter from the edge
        edge = maxX - (self.PP.d0*0.5)*self.PP.pPmm
        
        #Find if we have measurments at that position
        indexes = [i for i in range(len(cent_xs)) if cent_xs[i] > xposAfterDiffuser 
                and cent_xs[i] < edge]
        print('indexes of positions that are great than %d and smaller than %d' %(xposAfterDiffuser, edge))
        print('indexes = ', end='')
        print(indexes)

        meany =np.mean(cent_ys[indexes])

        #Find if we have measurments at that position
        
        indexes = [i for i in indexes if abs(cent_ys[i] -meany+0.0)/meany < 0.2 ]
                
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
            if forPub:
                fig = plt.figure(figsize=(6, 6), dpi=200); ax = fig.add_subplot(111); FS=18
                plt.plot(cent_x/self.PP.pPmm, cent_y/self.PP.pPmm, 'bs', markeredgecolor='b', label='Centroid position')
                plt.plot(cent_xs[indexes]/self.PP.pPmm, cent_ys[indexes]/self.PP.pPmm, 'ro', markeredgecolor='r', label='Centroid position used for final offset')
                plt.plot(self.intialCentreingXpos/self.PP.pPmm, self.intialCentreingYpos/self.PP.pPmm, 'g<', markeredgecolor='g', label='Centroid position used for initial offset')
                plt.plot(self.x_gap/self.PP.pPmm, self.y_gap/self.PP.pPmm, 'ch', markeredgecolor='c', label='Centroid position used Gap Width')
            else:
                fig = plt.figure(figsize=(8, 6), dpi=200); ax = fig.add_subplot(111); FS=12
                plt.plot(cent_x, cent_y, 'bs', label='Centroid position')
                plt.plot(cent_xs[indexes], cent_ys[indexes], 'ro', label='Centroid position used for final offset')
                plt.plot(self.intialCentreingXpos, self.intialCentreingYpos, 'g<', label='Centroid position used for initial offset')
                plt.plot(self.x_gap, self.y_gap, 'ch', label='Centroid position used Gap Width')

#            plt.plot(self.pos_diffuser_x, self.pos_diffuser_y, 'ch',label='Centroid position for angle measurment' )
#            xmax=self.PP.EndChannel
#            ymin=ana.func(0.0, self.diffuser_m, self.diffuser_c)
#            ymax = ana.func(xmax, self.diffuser_m, self.diffuser_c)
#            plt.plot([0, xmax], [ymin, ymax], 'c-', label='Fit angle = %f' %(self.diffuser_angle))


            
            if self.wentRight:
                dc = centreline - self.distanceFromCentreline
            else:
                dc = centreline + self.distanceFromCentreline

            ic = centreline + self.aveOffset
            maxX = np.max(cent_x)

            if forPub: 
                maxX = maxX/self.PP.pPmm                    
                plt.plot([0, maxX], [centreline/self.PP.pPmm, centreline/self.PP.pPmm], 'k-', label='Centreline')
                plt.plot([0, maxX], [dc/self.PP.pPmm, dc/self.PP.pPmm], 'r-', label='Final Offset')   
                plt.xlabel('Centroid x-position [mm]', fontsize=FS); plt.ylabel('Centroid y-position [mm]', fontsize=FS)
                
            else:
                plt.plot([0, maxX], [centreline, centreline], 'k-', label='Centreline')
                plt.plot([0, maxX], [dc, dc], 'r-', label='Final Offset')                
                plt.xlabel('Centroid x-position [pixels]', fontsize=FS); plt.ylabel('Centroid y-position [pixels]', fontsize=FS)

            if TJ_Mode:
                cl = (self.PP.LeftCornerInflow + self.PP.RightCornerInflow)/2.0
                ic = cl + self.aveOffset
                maxY = np.max(cent_y)
                if forPub: ic = ic/self.PP.pPmm; maxY = maxY/self.PP.pPmm
                plt.plot([ic, ic], [0, maxY],  'g-', label='Initial Offset')
            else:

                if forPub: ic = ic/self.PP.pPmm; 
                plt.plot([0, maxX], [ic, ic], 'g-', label='Initial Offset')
                
                
            
            
            if not forPub:
                plt.text(0.01, 0.8, 'Final Offset from Centreline = %.3f mm' %(self.distanceFromCentreline/self.PP.pPmm), fontsize = 16, horizontalalignment='left', verticalalignment='center', transform = ax.transAxes)
                plt.text(0.01, 0.4, 'Initial Centreing = %.3f mm' %(self.aveOffset/self.PP.pPmm), fontsize = 16, horizontalalignment='left', verticalalignment='center', transform = ax.transAxes)
                
                plt.text(0.55, 0.3, 'Base Flow = %.1f ml/min' %(self.volumFlux), fontsize = 15, horizontalalignment='left', verticalalignment='center', transform = ax.transAxes)
                plt.text(0.55, 0.2, 'Pinching Flow = %.1f ml/min' %(self.volumFlux2), fontsize = 15, horizontalalignment='left', verticalalignment='center', transform = ax.transAxes)
                        
            print('self.PP.ChannelTop - 4.0*self.PP.pPmm = %f' %(self.PP.ChannelTop - 4.0*self.PP.pPmm))
            print('self.distanceFromCentreline = %f' %self.distanceFromCentreline)
            print('self.PP.pPmm = %f' %(self.PP.pPmm))
            


            
            if forPub:
                minylim= 150.0; maxylim = 450.0
                plt.xlim( (0.0/self.PP.pPmm, 900.0/self.PP.pPmm))
            else:
                maxylim = self.PP.ChannelBottom + 4.0*self.PP.pPmm
                minylim = self.PP.ChannelTop - 4.0*self.PP.pPmm
                if minylim > dc:
                    minylim = dc *0.9
                    
            

            from matplotlib.patches import Rectangle
            currentAxis = plt.gca()
            print('self.PP.EndChannel, self.PP.ChannelBottom = %d, %d' %(self.PP.EndChannel, self.PP.ChannelBottom))
            
            XY=[]
            x11 = self.PP.RightCornerInflow ; y11 = self.PP.ChannelBottom
            x12 = 1.0*self.PP.pPmm ; y12 = maxylim - self.PP.ChannelBottom
            XY += [x11, y11, x12, y12]
            
            x21 = 0.0 ; y21 = self.PP.ChannelBottom
            x22 =  self.PP.LeftCornerInflow; y22 = maxylim - self.PP.ChannelBottom
            XY += [x21, y21, x22, y22]
            
            x31 = 0.0; y31 = self.PP.ChannelTop
            x32 = self.PP.RightCornerInflow ; y32 = -1.0*(self.PP.ChannelTop - minylim)
            XY += [x31, y31, x32, y32] ; XY = np.array(XY)
            
            if forPub: XY = XY/self.PP.pPmm; minylim = minylim/self.PP.pPmm; maxylim = maxylim/self.PP.pPmm

            currentAxis.add_patch(Rectangle((XY[0], XY[1]), XY[2], XY[3], facecolor="grey"))
            currentAxis.add_patch(Rectangle((XY[4], XY[5]), XY[6], XY[7], facecolor="grey"))
            currentAxis.add_patch(Rectangle((XY[8], XY[9]), XY[10], XY[11], facecolor="grey"))
            
            plt.ylim( ( minylim, maxylim))
#            currentAxis.add_patch(Rectangle((self.PP.RightCornerInflow, self.PP.ChannelBottom), 
#                                            1.0*self.PP.pPmm, maxylim - self.PP.ChannelBottom, facecolor="grey"))
#
#            currentAxis.add_patch(Rectangle((0, self.PP.ChannelBottom), 
#                                            self.PP.LeftCornerInflow, maxylim - self.PP.ChannelBottom, facecolor="grey"))
#            
#            dist = self.PP.ChannelTop - minylim
#            currentAxis.add_patch(Rectangle((0, self.PP.ChannelTop), 
#                                            self.PP.RightCornerInflow, -dist, facecolor="grey"))

            if not forPub:
                plt.title("Centroid position -  " +self.name, fontsize=12)
                plt.legend(loc='best', fontsize=FS/2, ncol=2)
#            else:
#                plt.legend(loc='best', fontsize=FS/2)
            fig.tight_layout()
            
            if not forPub:
                plt.savefig((self.path+"CentroidPosition_" +self.name+".jpg"), dpi=300)
            else:
                plt.gca().invert_yaxis()
                fig.set_size_inches(6.5, 6)
                fig.tight_layout()
                sn = ana._save_for_latex("CentroidPosition_" +self.name+"_forpub.jpg")
                if   forPubpath == '':
                    plt.savefig(self.path+sn, dpi=400)
                else:
                    plt.savefig(forPubpath+sn, dpi=400)
            
        else:
            print('self.plot = %s' %self.plot)
            
    def _speedVsXpos(self, savepath=''):
        '''Plot speed against the x-position'''  
        dispCutoff=40 #10 pixel max displacmeent between frames
        FS = 18
        dx= np.gradient(self.centroid_x+0.0)
        dy = np.gradient(self.centroid_y+0.0)
                
        xpos = self.centroid_x
        
        speed = np.sqrt(np.power(dx, 2.0) + np.power(dy, 2.0))
        ImgNr = np.arange(len(speed))
        
        delIndex=[]
        for uu in reversed(range(len(speed))):
            if abs(speed[uu]) > dispCutoff:
                delIndex.append(uu)

        speed = np.delete(speed, delIndex)/self.pixelsPmm  * self.FPS
        ImgNr = np.delete(ImgNr, delIndex)
        xpos = np.delete(xpos, delIndex)
        
        aveWindow = 6
        speed_ave =  gen.runningMean(speed, aveWindow)  
        speed_ave =  gen.runningMean(speed_ave, aveWindow)         
        speed_ave =  gen.runningMean(speed_ave, aveWindow)         
        speed_ave =  gen.runningMean(speed_ave, aveWindow)     
        
        fig = plt.figure(figsize=(6, 6), dpi=200); ax = fig.add_subplot(111)
        plt.plot(ImgNr/self.FPS, speed_ave, 'bs'); plt.xlabel('Time [s]', fontsize=FS); plt.ylabel('Capsule Speed [mm/s]', fontsize=FS)
        sn = ana._save_for_latex("SpeedVsTime_" +self.name+".jpg")
        a = plt.axes([0.4, 0.3, .5, .5], axisbg='w')
        plt.plot(ImgNr/self.FPS, speed_ave, 'bs'); plt.xlabel('Time [s]', fontsize=FS/2); plt.ylabel('Capsule Speed [mm/s]', fontsize=FS/2)
#        plt.xlim(0, 1.2)
        plt.xlim(0, 12)

#        plt.plot(xpos[12:-10]/self.PP.pPmm, speed_ave[12:-10], 'bs-'); plt.xlabel('$x$-position of Centroid [mm]', fontsize=FS); 
#        plt.xlim([0.0, 900.0/self.PP.pPmm]); sn = ana._save_for_latex("SpeedVsXpos_" +self.name+".jpg")      
#        plt.ylabel('Capsule Speed [$mm/s$]', fontsize=FS)

#        ax.invert_xaxis()
#        plt.title("Speed -  " +self.name)
#        plt.legend(loc='best', fontsize=6, ncol=2)
        fig.set_size_inches(6.5, 6)
        fig.tight_layout()
        
        if savepath != '':
            plt.savefig(os.path.join(savepath+sn), dpi=400)
        else:
            plt.savefig((self.path+sn), dpi=400)
#        print('self.distanceFromCentreline: ' +str(self.distanceFromCentreline) + ' +/- ' + str(self.distanceFromCentrelineSTD))
    
    def _maxDeformation(self):
        '''Find max deformation, as measured by taylor deformation parameter'''
        
        #remove entries with no values
        cent_x, cent_y = self._centroidsWithoutErrors()
 
        #First find area around semi-sphere, which is at x-pos end-channel        
        print('self.PP.LeftCornerInflow = %d' %(self.PP.LeftCornerInflow) )
        print('self.PP.EndChannel = %d' %(self.PP.EndChannel) )
        
        startSS = self.PP.LeftCornerInflow - 0.5*self.PP.d0*self.PP.pPmm
        endSS = self.PP.EndChannel + 0.5*self.PP.d0*self.PP.pPmm
        indexes = [i for i in range(len(cent_x)) if cent_x[i] > startSS and cent_x[i] < endSS]
        print('maxDeformation')
        print('startSS = %d' %(startSS))
        print('endSS = %d' %(endSS))
        print('indexes = ')        
        print(indexes)
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
#        print(self.height[indexes])
        try:
            self.minHeight = np.min(ana._removeERRORCONST(self.height[indexes]))
        except:
            self.minHeight = ERROR_CONST
                
     
#        if self.plot:
#            fig = plt.figure(figsize=(8, 6), dpi=200); ax = fig.add_subplot(111)
#            plt.plot(self.centroid_x, self.centroid_y, 'bs', label='Centroid position')
#            
#            plt.plot(self.centroid_x[indexes], self.centroid_y[indexes], 'ro', label='Centroid position used for max deformation')
#            
#            ymin, ymax = ax.get_ylim()
#            
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
        endSS = self.PP.LeftCornerInflow - (self.PP.d0)*self.PP.pPmm
#        endSS = self.PP.RightCornerInflow + (0.25*self.PP.d0)*self.PP.pPmm
        
        #remove entries with no values
        cent_x, cent_y = self._centroidsWithoutErrors()
        
        # lets not take the first half diameter as there might be artefacts
        # from capsule tracking as the capsule enters the image
        startSS =  0.5*self.PP.d0*self.PP.pPmm
        
        #Find if we have measurments at that position
        indexes = [i for i in range(len(cent_x)) if cent_x[i] < endSS and cent_x[i] > startSS and cent_y[i] > self.PP.ChannelTop and cent_y[i] < self.PP.ChannelBottom]
        
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
        

    def _findCentering_TJ_Mode(self):
        '''find centering of capsule in channel leading up 2dn inflow'''        
        #find centerline
        centreline = (self.PP.RightCornerInflow  + self.PP.LeftCornerInflow )/2.0

        #remove entries with no values
        cent_x, cent_y = self._centroidsWithoutErrors()
        
        #Find if we have measurments at that position
        indexes = [i for i in range(len(cent_x)) if cent_x[i] < self.PP.RightCornerInflow and cent_x[i] > self.PP.LeftCornerInflow and cent_y[i] > self.PP.ChannelBottom and cent_y[i] < (self.PP.ChannelBottom + 4.0 * self.PP.pPmm)]

        meanXpos=cent_x[indexes]
        
        offset = meanXpos - centreline
        
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
        
    def _findGap(self, TJ_Mode=False, forPub=False):
        '''Find gap between capsule and wall during pinching flow'''
        indexesR=[]
        for jj in range(len(self.centroid_x)):
            #check whether capsule forward central quater is between the end of the channel and the right hand side of the inflow, i.e. within a 1mm strech
                        
#            if (self.centroid_x[jj] + self.width[jj]/2.0 * 0.8) > self.PP.RightCornerInflow and (self.centroid_x[jj] + self.width[jj]/4.0) < self.PP.EndChannel:
            #Within  1 mm of the right hand side of the channel and in the central part of the channel in the y direction
            #Doesn't work well as the capsule is already outside and expanding
#            if(self.centroid_x[jj] > (self.PP.RightCornerInflow - self.PP.pPmm*1.0) 
#                and self.centroid_x[jj] < (self.PP.RightCornerInflow + self.PP.pPmm*1.0)
#                and self.centroid_y[jj] < self.PP.ChannelBottom - self.PP.pPmm*1.0
#                and self.centroid_y[jj] > self.PP.ChannelTop + self.PP.pPmm*1.0):
            if(self.centroid_x[jj] > (self.PP.RightCornerInflow - self.PP.pPmm*3.5) 
                and self.centroid_x[jj] < (self.PP.RightCornerInflow - self.PP.pPmm*0.5)
                and self.centroid_y[jj] < self.PP.ChannelBottom - self.PP.pPmm*1.0
                and self.centroid_y[jj] > self.PP.ChannelTop + self.PP.pPmm*0.5):
                indexesR.append(jj)
                
#        print("In _findGap(): ", end=''); print(indexesR)
        gap=[]; x_gap=[]; y_gap=[];  gap_cap_h=[]; gap_d12=[];    
        topgap=[]; #Distance between capsule and top side wall of straigth channel
        
        if len(indexesR) ==0:
            gap.append(ERROR_CONST)
            gap_cap_h.append(ERROR_CONST)
            gap_d12.append(ERROR_CONST)
            x_gap.append(ERROR_CONST)
            y_gap.append(ERROR_CONST)
            topgap.append(ERROR_CONST)

            indexesR.append(ERROR_CONST)
        else:
            for kk in indexesR:
                #gap.append(self.centroid_y[kk] - self.height[kk]/2.0 - self.PP.ChannelTop) #ChannelTop is at a large y coordinate than y
                gap.append(self.PP.ChannelBottom -  (self.centroid_y[kk] + self.height[kk]/2.0) )  
                topgap.append((self.centroid_y[kk] - self.height[kk]/2.0) - self.PP.ChannelTop    )
                
    #            print("self.PP.ChannelBottom = %d self.centroid_y[kk] = %d self.height[kk]/2.0 = %d gap = %d"  %(self.PP.ChannelBottom, self.centroid_y[kk], self.height[kk]/2.0, gap[-1]))
                gap_cap_h.append(self.height[kk])
                gap_d12.append(self.d12[kk])
                x_gap.append(self.centroid_x[kk])
                y_gap.append(self.centroid_y[kk])
        x_gap = np.array(x_gap); y_gap = np.array(y_gap)            
        
        self.x_gap = x_gap;self.y_gap = y_gap;
        self.gap_d12 = np.mean(gap_d12)
        self.gap_cap_h = np.mean(gap_cap_h)
        self.gapSpace = np.mean(gap)
        self.gapSpaceSTD = np.std(gap)
        self.topGapSpace = np.mean(topgap)
        self.topGapSpaceSTD = np.std(topgap)        
        
        if self.plot:
            fig = plt.figure(figsize=(8, 6), dpi=200); ax = fig.add_subplot(111)
#            print("indexesR = " , end='')            
#            print(indexesR)
#            print(np.array(gap)/self.PP.pPmm)
            plt.plot(indexesR, np.array(gap)/self.PP.pPmm, 'cs', label='Gap width')
            plt.plot(indexesR, np.array(topgap)/self.PP.pPmm, 'ro', label='Top Gap width')
            plt.plot([np.min(indexesR), np.max(indexesR)], [self.topGapSpace/self.PP.pPmm, self.topGapSpace/self.PP.pPmm])
            plt.plot([np.min(indexesR), np.max(indexesR)], [self.gapSpace/self.PP.pPmm, self.gapSpace/self.PP.pPmm])
            
            
            plt.title("Gap Width -  " +self.name, fontsize=12)
            plt.xlabel('Image Number')
            plt.ylabel('Gap [mm]')
            plt.legend(loc='best', fontsize=6, ncol=1)

            ax2 = ax.twinx()
            plt.plot(indexesR, np.array(gap_d12), 'rh', label='Deformation')            
            ax2.set_ylabel("Deformation Parameter $D_{12}$")

            fig.tight_layout()

            plt.savefig((self.path+"Gap_" +self.name+".jpg"), dpi=300)
        
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
#            plt.legend(loc='best', fontsize=6, ncol=1)
            fig.tight_layout()

            plt.savefig((self.path+"RelaxationVelocity_" +self.name+".jpg"), dpi=300)
        
        
    def _writeToFilePinching(self):
        '''Write measurments of run to file'''

        headerString = '\tVolume Flux 2\tFinal Distance From Centreline \t STD \tmax d12 \t maxWidth \tmin height \t Centring, \t Centring STD \t Gap \t Gap STD \t Gap Capsule Height \tgap_d12 \tTop Gap \t Top Gap SDT'
        String = '\t%.1f \t%.6f \t%.6f \t%.6f \t%d \t%d \t%.6f \t%.6f \t%.6f \t%.8f \t%.6f \t%.6f \t%.6f \t%.6f' %(self.volumFlux2, self.distanceFromCentreline, self.distanceFromCentrelineSTD, self.maxd12, self.maxWidth,  self.minHeight, self.aveOffset, self.aveOffsetSTD, self.gapSpace, self.gapSpaceSTD, self.gap_cap_h, self.gap_d12, self.topGapSpace, self.topGapSpaceSTD,)
        
#        headerString = '\tVolume Flux 2\tFinal Distance From Centreline \t STD \tmax d12 \t maxWidth \tmin height \t Centring, \t Centring STD '
#        String = '\t%.1f \t%.6f \t%.6f \t%.6f \t%d \t%d \t%.6f \t%.6f ' %(self.volumFlux2, self.distanceFromCentreline, self.distanceFromCentrelineSTD, self.maxd12, self.maxWidth,  self.minHeight, self.aveOffset, self.aveOffsetSTD)
        
        self._writeToFile(headerString, String)
        
        
        
#==============================================================================
# Other Functions
#==============================================================================
def averageTrajectories(directory, savePath, forPub=False, savename='', TJ=True):
    ''' Average all trajectories for each volumn flux '''    
    assert TJ #other thing not implemented
    
    #step 1
    # Interpolate them y  is a function of x to get equal sampeling points
    from scipy.interpolate import interp1d
    
    listFolders = ana._get_foler_list(directory)
    minX=0; maxX=10000;
    inter_funcs=[]; Q=[]; Q2=[]


#    c=0
#    for f in listFolders:
#        print('%d \t%s' %(c, f))
#        c+=1
        
    for f in listFolders:
#    for  f in list(listFolders[i] for i in [0, 4, 9, 11, 17, 28]): #5 ML/MIN
#    for  f in list(listFolders[i] for i in [3, 8, 23, 25, 26, 33]): #20 ML/MIN
        try:
            DRP = DataRunPinching(os.path.join(directory, f))

        except:
            print("\nDidn't work for \n%s\n\n" %f)
            continue
            
        #remove entries with no values from centroids
        cent_x, cent_y = DRP._centroidsWithoutErrors()
#        x, y = ana._remove_outliers(cent_x,cent_y, window=15, cutoff=1.5)
        x, y = ana._remove_outliers(cent_x,cent_y, window=15, cutoff=1.5, check=True, sp='%s_filtered.jpg' %f)
#        x2, y2 = ana._1d_remove_outliers(cent_x,cent_y, window=30, cutoff=1.5)
#        y2,x2 = ana._1d_remove_outliers(y2,x2, window=30, cutoff=1.5)

        if np.min(x) < 130:
            Q.append(DRP.volumFlux)
            Q2.append(DRP.volumFlux2)
            inter_funcs.append(interp1d(x, y))
            
            if np.min(x) > minX: minX = np.min(x)
            if np.max(x) < maxX: maxX = np.max(x)
            
    Q2, Q, inter_funcs = zip(*sorted(zip(Q2, Q, inter_funcs)))
    
    print(Q2)
    
    qs = np.unique(Q2)
    
    if minX < 70: minX=70
    x = np.arange(minX, maxX, 1)
    
    ave_y=[];ave_y2=[]    ; all_y=[]    
    ave_ey=[];
    print('qs: ')
    print(qs)
    for q in qs:
        notFound = True
        ys=[]  
        append_count=0
        for ii in range(len(Q2)):
            if q == Q2[ii]:
                ys.append(np.array(inter_funcs[ii](x)))
                append_count +=1
#        print("appended %d runs for Q = %.1f" %(append_count, q))
        ys = np.array(ys)
        all_y.append(ys)
#        print('ys : ')
#        print(ys)
        
        threshold = 0.25; runCounter=0
        while notFound and runCounter < 10:
            y=[]; y2=[]
            ey=[]
            for kk in range(len(ys[0])):
#            for kk in range(90, 115):
    #            print("\nys[:, kk]")            
    #            print(ys[:, kk])
                ytemp = np.zeros(len(ys))
                ytemp = ys[:, kk]
                
                outlier = ana.is_outlier(ytemp, thresh=threshold)
    #            if np.any(outlier) == True:
    #                countTrue =0
    #                for o in outlier: 
    #                    if o: 
    #                        countTrue +=1
    #                print("%d / %d are outliers for kk = %d" %(countTrue, len(ys), kk))
                    
                yt=[]; 
                ii=0
                for o in outlier:
                    if not o:
                        yt.append(ytemp[ii])
                    ii +=1
    
#                print("len(yt)/(ytemp) = %d / %d \t\t np.mean(yt)/(ytemp) = %.4f / %.4f" %(len(yt), len(ytemp), np.mean(yt), np.mean(ytemp)))
                y.append(np.mean(yt))  
                ey.append(np.std(yt))
                y2.append(np.mean(ytemp))
                
            if not np.isnan(y).any():
                notFound = False
            else:
                threshold += 0.25
                print("Increasing threshold to %f" %threshold)
                runCounter +=1
            
        ave_y.append(np.array(y))
        ave_ey.append(np.array(ey))
        ave_y2.append(np.array(y2))
        
#    x = x[90:115]
#    print('ave_y / ave_y2:')
#    print(ave_y)
#    print(ave_y2)
    m = ana.markerSymbols(); c = ana.coloursForMarker(len(ave_y))
    fig = plt.figure(figsize=(8, 6), dpi=200); ax = fig.add_subplot(111)
    for ii in range(len(ave_y)):    
        plt.errorbar(x/DRP.PP.pPmm, ave_y[ii]/DRP.PP.pPmm, yerr=ave_ey[ii]/DRP.PP.pPmm, ecolor='k', markersize=6, marker=m[ii], color=c[ii], label='Q = %.1f ml/min' %(qs[ii]))
#        plt.plot(x/DRP.PP.pPmm, ave_y2[ii]/DRP.PP.pPmm, markersize=3, marker=m[ii+1], color=c[ii+1], label='No filter Q = %.1f' %(qs[ii]))    
#    for jj in range(len(all_y)):
#            for ll in range(len(all_y[jj])):
#                plt.plot(x/DRP.PP.pPmm, all_y[jj][ll][90:115]/DRP.PP.pPmm, marker=m[3], color=c[3], label='Q = %.1f' %(qs[ii]))                
#    if forPub:
#        minylim= 150; maxylim = 450
#        plt.xlim( (0, 900))
#    else:
    maxylim = DRP.PP.ChannelBottom + 4.0*DRP.PP.pPmm
    minylim = DRP.PP.ChannelTop - 8.0*DRP.PP.pPmm

                    
    plt.ylim( (  maxylim/DRP.PP.pPmm, minylim/DRP.PP.pPmm))
    
    from matplotlib.patches import Rectangle
    currentAxis = plt.gca()
    print('EndChannel, ChannelBottom = %d, %d' %(DRP.PP.EndChannel/DRP.PP.pPmm, DRP.PP.ChannelBottom/DRP.PP.pPmm))
    currentAxis.add_patch(Rectangle((DRP.PP.RightCornerInflow/DRP.PP.pPmm, DRP.PP.ChannelBottom/DRP.PP.pPmm), 
                                    1.0*DRP.PP.pPmm/DRP.PP.pPmm, (maxylim - DRP.PP.ChannelBottom)/DRP.PP.pPmm, facecolor="grey"))

    currentAxis.add_patch(Rectangle((0, DRP.PP.ChannelBottom/DRP.PP.pPmm), 
                                    DRP.PP.LeftCornerInflow/DRP.PP.pPmm, (maxylim - DRP.PP.ChannelBottom)/DRP.PP.pPmm, facecolor="grey"))
    
    dist = DRP.PP.ChannelTop - minylim
    currentAxis.add_patch(Rectangle((0, DRP.PP.ChannelTop/DRP.PP.pPmm), 
                                    DRP.PP.RightCornerInflow/DRP.PP.pPmm, -dist/DRP.PP.pPmm, facecolor="grey"))
            
    if not forPub:
        plt.title("Average Trajectories-  " +DRP.bathName, fontsize=12)
    plt.xlabel('$X$ [mm]')
    plt.ylabel('$Y$ [mm]')
    if forPub:
        plt.legend(loc='best', fontsize=12, ncol=1)
    else:
        plt.legend(loc='best', fontsize=6, ncol=1)
    fig.tight_layout()

    
    ax_inset=fig.add_axes([0.37,0.18,0.3,0.3])
    for ii in range(len(ave_y)):    
        ax_inset.errorbar(x/DRP.PP.pPmm, ave_y[ii]/DRP.PP.pPmm, yerr=ave_ey[ii]/DRP.PP.pPmm, ecolor='k', marker=m[ii], color=c[ii], label='Q = %.1f' %(qs[ii]))

    currentAxis = plt.gca()
    currentAxis.add_patch(Rectangle((DRP.PP.RightCornerInflow/DRP.PP.pPmm, DRP.PP.ChannelBottom/DRP.PP.pPmm), 
                                    1.0*DRP.PP.pPmm/DRP.PP.pPmm, (maxylim - DRP.PP.ChannelBottom)/DRP.PP.pPmm, facecolor="grey"))

    currentAxis.add_patch(Rectangle((0, DRP.PP.ChannelBottom/DRP.PP.pPmm), 
                                    DRP.PP.LeftCornerInflow/DRP.PP.pPmm, (maxylim - DRP.PP.ChannelBottom)/DRP.PP.pPmm, facecolor="grey"))
    
    dist = DRP.PP.ChannelTop - minylim
    currentAxis.add_patch(Rectangle((0, DRP.PP.ChannelTop/DRP.PP.pPmm), 
                                    DRP.PP.RightCornerInflow/DRP.PP.pPmm, -dist/DRP.PP.pPmm, facecolor="grey"))

    ax_inset.set_ylim((23.5, 20.5))
    ax_inset.set_xlim((8,10))

#    from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
#    from mpl_toolkits.axes_grid1.inset_locator import mark_inset
#    
#    axins = zoomed_inset_axes(ax, 6, loc=1)  # zoom = 6
#    for ii in range(len(ave_y)):    
#        axins.plot(x/DRP.PP.pPmm, ave_y[ii]/DRP.PP.pPmm, marker=m[ii], color=c[ii], label='Q = %.1f' %(qs[ii]))
#    axins.set_ylim((20,23))
#    axins.set_xlim((6,12))
#  
#    # draw a bbox of the region of the inset axes in the parent axes and
#    # connecting lines between the bbox and the inset axes area
#    mark_inset(ax, axins, loc1=2, loc2=4, fc="none", ec="0.5")
    

    if not savePath == '':
        if forPub:
            plt.savefig((savePath+"AverageTrajectory_forPub_" +DRP.bathName+".jpg"), dpi=300)
        else:
            plt.savefig((savePath+"AverageTrajectory_" +DRP.bathName+".jpg"), dpi=300)
    
        
def _1dFilterOutliers(a, threshold = 0.25, Max_Increases=5):
    notFound = True; runCounter=0
    while notFound and runCounter < Max_Increases:
        outlier = ana.is_outlier(a, thresh=threshold)
                
        yt=[]; 
        ii=0
        for o in outlier:
            if not o:
                yt.append(a[ii])
            ii +=1
            
        if len(yt) > 0:
            notFound = False
        else:
            threshold += 0.25
            print("Increasing threshold to %f" %threshold)
            runCounter +=1
            
    return np.mean(yt), np.std(yt)

def _1dfilterOutliers_ArrayOfArray(a):
        temp=[]; etemp=[]
        for kk in range(len(a[0])):
            tsy, tesy = _1dFilterOutliers(a[:, kk], threshold = 0.25, Max_Increases=5)
            temp.append(tsy)
            etemp.append(tesy)
        return np.array(temp), np.array(etemp)
        
def averageDeformationAndHeight(directory, savePath, forPub=False, savename='', plot_trajectores=True, plot_deformation=True):
    ''' Average deformation for each volumn flux along either x or y position'''        
    
    #step 1
    # Interpolate them a function of x to get equal sampeling points
    from scipy.interpolate import interp1d
    
    listFolders = ana._get_foler_list(directory)
    minX=0; maxX=10000;
    inter_funcs=[]; Q=[]; Q2=[]
    inter_funcs_d12=[]; inter_funcs_d12rect=[]; inter_funcs_h=[];
    
#    c=0
#    for f in listFolders:
#        print('%d \t%s' %(c, f))
#        c+=1
        
#    lf = listFolders.sort()
    for f in listFolders:
#    for  f in list(listFolders[i] for i in [0, 1, 16, 30, 48, 51, 53]):
#    for  f in list(listFolders[i] for i in [0, 4, 9, 11, 17, 28]): #5 ML/MIN T-J
#    for  f in list(listFolders[i] for i in [3, 8, 23, 25, 26, 33]): #20 ML/MIN T-J
        try:
            DRP = DataRunPinching(os.path.join(directory, f))

        except:
            print("\nDidn't work for \n%s\n\n" %f)
            continue
            
        #remove entries with no values from centroids
        cent_x, cent_y = DRP._centroidsWithoutErrors()
        x, y = ana._remove_outliers(cent_x,cent_y, window=15, cutoff=1.5, check=False, sp='%s_filtered.jpg' %f)

        if np.min(x) < 130: # variable to ensure that interpolated data spanns range, should be fine for pinching
            Q.append(DRP.volumFlux)
            Q2.append(DRP.volumFlux2)
            inter_funcs.append(interp1d(x, y))
#            print("f:%s \nlen(DRP.centroid_x) = %d, len(DRP.d12) = %d" %(f, len(DRP.centroid_x), len(DRP.d12)))
            inter_funcs_d12.append(interp1d(DRP.centroid_x.squeeze(), DRP.d12.squeeze()))
            inter_funcs_d12rect.append(interp1d(DRP.centroid_x.squeeze(), DRP.d12rect.squeeze()))
            inter_funcs_h.append(interp1d(DRP.centroid_x.squeeze(), DRP.height.squeeze()))
                
            if np.min(x) > minX: minX = np.min(x)
            if np.max(x) < maxX: maxX = np.max(x)
            
#    Q2, Q, inter_funcs = zip(*sorted(zip(Q2, Q, inter_funcs)))
#    print(Q2)
    # Step 2: average over different runs, excluding outliers
    qs = np.unique(Q2)
#    if minX < 70: minX=70
    x = np.arange(minX, maxX, 1)
    
    ave_y=[];ave_d12=[]; ave_d12rect=[]; ave_h=[]; all_y=[]    
    ave_ey=[]; ave_ed12=[]; ave_ed12rect=[]; ave_eh=[];
    print('qs: ')
    print(qs)
    for q in qs:
        notFound = True
        ys=[] ; d12=[]; d12rect=[]; h=[];  
        append_count=0
        for ii in range(len(Q2)):
            if q == Q2[ii]:
                ys.append(np.array(inter_funcs[ii](x)))
                d12.append(np.array(inter_funcs_d12[ii](x)))
                d12rect.append(np.array(inter_funcs_d12rect[ii](x)))
                h.append(np.array(inter_funcs_h[ii](x)))
                append_count +=1
#        print("appended %d runs for Q = %.1f" %(append_count, q))
        ys = np.array(ys); d12 = np.array(d12); d12rect = np.array(d12rect); h= np.array(h)

        reslt, e_reslt = _1dfilterOutliers_ArrayOfArray(ys)
        ave_y.append(reslt); ave_ey.append(e_reslt)

        reslt, e_reslt = _1dfilterOutliers_ArrayOfArray(d12)
        ave_d12.append(reslt); ave_ed12.append(e_reslt)
        reslt, e_reslt = _1dfilterOutliers_ArrayOfArray(d12rect)
        ave_d12rect.append(reslt); ave_ed12rect.append(e_reslt)

        reslt, e_reslt = _1dfilterOutliers_ArrayOfArray(h)
        ave_h.append(reslt); ave_eh.append(e_reslt)        
        
    # Step 3: Plotting

    if plot_trajectores:
        m = ana.markerSymbols(); c = ana.coloursForMarker(len(ave_y))
        fig = plt.figure(figsize=(8, 6), dpi=200); ax = fig.add_subplot(111)
        for ii in range(len(ave_y)):    
            plt.errorbar(x/DRP.PP.pPmm, ave_y[ii]/DRP.PP.pPmm, yerr=ave_ey[ii]/DRP.PP.pPmm, ecolor='k', markersize=6, marker=m[ii], color=c[ii], label='Q = %.1f ml/min' %(qs[ii]))

        maxylim = DRP.PP.ChannelBottom + 4.0*DRP.PP.pPmm
        minylim = DRP.PP.ChannelTop - 8.0*DRP.PP.pPmm
    
        plt.ylim( (  maxylim/DRP.PP.pPmm, minylim/DRP.PP.pPmm))
        
        from matplotlib.patches import Rectangle
        currentAxis = plt.gca()
        print('EndChannel, ChannelBottom = %d, %d' %(DRP.PP.EndChannel/DRP.PP.pPmm, DRP.PP.ChannelBottom/DRP.PP.pPmm))
        currentAxis.add_patch(Rectangle((DRP.PP.RightCornerInflow/DRP.PP.pPmm, DRP.PP.ChannelBottom/DRP.PP.pPmm), 
                                        1.0*DRP.PP.pPmm/DRP.PP.pPmm, (maxylim - DRP.PP.ChannelBottom)/DRP.PP.pPmm, facecolor="grey"))
    
        currentAxis.add_patch(Rectangle((0, DRP.PP.ChannelBottom/DRP.PP.pPmm), 
                                        DRP.PP.LeftCornerInflow/DRP.PP.pPmm, (maxylim - DRP.PP.ChannelBottom)/DRP.PP.pPmm, facecolor="grey"))
        
        dist = DRP.PP.ChannelTop - minylim
        currentAxis.add_patch(Rectangle((0, DRP.PP.ChannelTop/DRP.PP.pPmm), 
                                        DRP.PP.RightCornerInflow/DRP.PP.pPmm, -dist/DRP.PP.pPmm, facecolor="grey"))
                
        if not forPub:
            plt.title("Average Trajectories-  " +DRP.bathName, fontsize=12)
        plt.xlabel('$X$ [mm]')
        plt.ylabel('$Y$ [mm]')
        if forPub:
            plt.legend(loc='best', fontsize=12, ncol=1)
        else:
            plt.legend(loc='best', fontsize=6, ncol=1)
        fig.tight_layout()
    
        ax_inset=fig.add_axes([0.37,0.18,0.3,0.3])
        for ii in range(len(ave_y)):    
            ax_inset.errorbar(x/DRP.PP.pPmm, ave_y[ii], yerr=ave_ey[ii], ecolor='k', marker=m[ii], color=c[ii], label='Q = %.1f' %(qs[ii]))
    
        currentAxis = plt.gca()
        currentAxis.add_patch(Rectangle((DRP.PP.RightCornerInflow/DRP.PP.pPmm, DRP.PP.ChannelBottom/DRP.PP.pPmm), 
                                        1.0*DRP.PP.pPmm/DRP.PP.pPmm, (maxylim - DRP.PP.ChannelBottom)/DRP.PP.pPmm, facecolor="grey"))
    
        currentAxis.add_patch(Rectangle((0, DRP.PP.ChannelBottom/DRP.PP.pPmm), 
                                        DRP.PP.LeftCornerInflow/DRP.PP.pPmm, (maxylim - DRP.PP.ChannelBottom)/DRP.PP.pPmm, facecolor="grey"))
        
        dist = DRP.PP.ChannelTop - minylim
        currentAxis.add_patch(Rectangle((0, DRP.PP.ChannelTop/DRP.PP.pPmm), 
                                        DRP.PP.RightCornerInflow/DRP.PP.pPmm, -dist/DRP.PP.pPmm, facecolor="grey"))
    
        ax_inset.set_ylim((23.5, 20.5))
        ax_inset.set_xlim((8,10))
            
        if not savePath == '':
            if forPub:
                plt.savefig((savePath+"AverageTrajectory_forPub_" +DRP.bathName+".jpg"), dpi=300)
            else:
                plt.savefig((savePath+"AverageTrajectory_" +DRP.bathName+".jpg"), dpi=300)
   
    if plot_deformation:
        m = ana.markerSymbols(); c = ana.coloursForMarker(len(ave_y))
        fig = plt.figure(figsize=(8, 6), dpi=200); ax = fig.add_subplot(111)
        for ii in range(len(ave_y)):    
            plt.errorbar(x/DRP.PP.pPmm, ave_d12[ii]/DRP.PP.pPmm, yerr=ave_ed12[ii]/DRP.PP.pPmm, ecolor='k', markersize=6, marker=m[ii], color=c[ii], label='Q = %.1f ml/min' %(qs[ii]))

        if not forPub:
            plt.title("Average Trajectories-  " +DRP.bathName, fontsize=12)
        plt.xlabel('$X$ [mm]')
        plt.ylabel('Tylor Deformation Parameter $D_{12}$')
        
        if forPub:
            plt.legend(loc='best', fontsize=12, ncol=1)
        else:
            plt.legend(loc='best', fontsize=6, ncol=1)
        fig.tight_layout()

        if not savePath == '':
            if forPub:
                plt.savefig((savePath+"AverageD12_forPub_" +DRP.bathName+".jpg"), dpi=300)
            else:
                plt.savefig((savePath+"AverageD12_" +DRP.bathName+".jpg"), dpi=300)

        fig = plt.figure(figsize=(8, 6), dpi=200); ax = fig.add_subplot(111)
        if forPub: FS=16
        else: FS=12

        #Show where Pinching channel hits        
        plt.plot([DRP.PP.RightCornerInflow/DRP.PP.pPmm, DRP.PP.RightCornerInflow/DRP.PP.pPmm], 
                 [0, 10], '--k')
        plt.plot([DRP.PP.LeftCornerInflow/DRP.PP.pPmm, DRP.PP.LeftCornerInflow/DRP.PP.pPmm], 
                 [0, 10], '--k')
        plt.axvspan(0.0, DRP.PP.LeftCornerInflow/DRP.PP.pPmm, facecolor='g', alpha=0.2)
        plt.axvspan(DRP.PP.LeftCornerInflow/DRP.PP.pPmm, DRP.PP.RightCornerInflow/DRP.PP.pPmm, facecolor='y', alpha=0.2)
        
        for ii in range(len(ave_y)):    
#            plt.errorbar(x/DRP.PP.pPmm, ave_h[ii]/DRP.PP.pPmm, yerr=ave_eh[ii]/DRP.PP.pPmm, ecolor='k', markersize=4, marker=m[len(ave_y)+ii], color=c[ii], label='Q = %.1f ml/min' %(qs[ii]))
            plt.errorbar(x/DRP.PP.pPmm, 
#                         ave_h[ii]/DRP.PP.pPmm,
#                         gen.smooth(ave_h[ii]/DRP.PP.pPmm,window_len=11,window='hanning'), 
                         gen.savitzky_golay(ave_h[ii]/DRP.PP.pPmm, 15, 3, deriv=0, rate=1),
                         yerr=ave_eh[ii]/DRP.PP.pPmm, 
#                         ecolor='k', 
                         markersize=2, 
                         marker=m[ii], color=c[ii], 
                         label='$Q_{Pin}$ = %.1f ml/min' %(qs[ii]))
        if not forPub:
            plt.title("Average Trajectories-  " +DRP.bathName, fontsize=12)
        plt.xlabel('$X$ [mm]', fontsize=FS)
        plt.ylabel("Smoothed Capsule Lateral ($y$) Size [mm]", fontsize=FS)
        
        if forPub:
            plt.legend(loc='best', fontsize=10, ncol=1)
        else:
            plt.legend(loc='best', fontsize=6, ncol=1)
        fig.tight_layout()
        plt.xlim((minX/DRP.PP.pPmm, maxX/DRP.PP.pPmm))
        plt.ylim((2,5.0))
        if not savePath == '':
            if forPub:
                plt.savefig((savePath+"AverageH_forPub_" +DRP.bathName+".jpg"), dpi=300)
            else:
                plt.savefig((savePath+"AverageH_" +DRP.bathName+".jpg"), dpi=300)
                
        print("x-pos of left corner of inflow: %f" %(DRP.PP.LeftCornerInflow/DRP.PP.pPmm))        
        print("x-pos of right corner of inflow: %f" %(DRP.PP.RightCornerInflow/DRP.PP.pPmm))

            

        
    
        
            

#==============================================================================

#==============================================================================





#==============================================================================
# Capsule Level Visualization
#==============================================================================
        
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
        self.gap = np.zeros(len(self.volumeFlux))
        self.gapSTD = np.zeros(len(self.volumeFlux))
        self.gap_capsule_height = np.zeros(len(self.volumeFlux))
        self.gap_d12 = np.zeros(len(self.volumeFlux))
        self.topGap = np.zeros(len(self.volumeFlux))
        self.topGapSTD = np.zeros(len(self.volumeFlux))        
        
        self._readDataPinching()
        
        self.PP = ParametersPinchingMk2(listingImageFunc = sortPhotosPCO,
                                   coverImgFunc = coverSidePinching,
                                   readParametersFromFile=os.path.join(self.directory, self.name[0] ))
        self._findAverages()
            
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
            if len(entries)>=18:
                self.volumeFlux2[ii] = (float(entries[4]))
                self.finalDistCentre[ii] = (float(entries[5]))
                self.finalDistCentreSTD[ii] = (float(entries[6]))
                self.maxD12[ii] = (float(entries[7]))
                self.maxWidth[ii] = (float(entries[8]))
                self.minHeight[ii] = (float(entries[9]))
                self.Offset[ii] = (float(entries[10]))
                self.OffsetSTD[ii] = (float(entries[11]))
                self.gap[ii] = (float(entries[12]))
                self.gapSTD[ii] = (float(entries[13]))
                self.gap_capsule_height[ii] = (float(entries[14]))
                self.gap_d12[ii] = (float(entries[15])) 
                self.topGap[ii] = (float(entries[16]))
                self.topGapSTD[ii] = (float(entries[17]))                
                ii += 1
                


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
#        self._findAverages()
        self._plotQVsFinalOffset()
        self._plotQVsFinalOffsetND()
        self._plotQVsInitialOffset()
#        self._plot1oQVsOffset()
        self._plotQVsD12()
        
#        self._plotAbsFinalVsInitialOffset()
        self._plotFinalVsInitialOffset()
#        self._plotmCvsO()
#        self._plotcCvsO()
#        self._plotQVsFinalOffsetByStreamline()
        self._plotMinHeight()

        
    def _findAverages(self):
        '''Find averages of data loaded from file'''
        self.avePPmm, _ = self.averageByQ2(self.pixelsPmm) 
        self.aveFinalDistCentre, self.aveFinalDistCentreSTD = self.averageByQ2(self.finalDistCentre) 
 
        self.aveMaxD12, self.aveMaxD12STD = self.averageByQ2(self.maxD12)
        self.aveMaxWidth, self.aveMaxWidthSTD = self.averageByQ2(self.maxWidth)
        self.aveMinHeight, self.aveMinHeightSTD = self.averageByQ2(self.minHeight)
        self.aveOffset, self.aveOffsetSTD = self.averageByQ2(self.Offset)
        
        self.aveGap, self.aveGapSTD = self.averageByQ2(self.gap)  
        self.aveGap_capsule_height, self.aveGap_capsule_heightSTD = self.averageByQ2(self.gap_capsule_height)
        self.aveGap_d12, self.aveGap_d12STD = self.averageByQ2(self.gap_d12)
        self.aveTopGap, self.aveTopGapSTD = self.averageByQ2(self.topGap)  
        
    def _plotQVsGapSpace(self):
        '''Plot final distance from centreline versus volume flux Q'''
        m = ana.markerSymbols(); c = ana.coloursForMarker()
        fig = plt.figure(figsize=(8, 6), dpi=200,); ax = fig.add_subplot(111)
        
        #old self.gapSpace
#        plt.errorbar(self.volumeFlux2, self.gap/self.pixelsPmm, 
#                    yerr=self.gapSTD/self.pixelsPmm,
#                    markersize=2,
#                    linestyle='None', marker = m[0], color = c[0])
                    
        plt.errorbar(self.volumeFluxes2,  self.aveGap/self.avePPmm, 
                     yerr=self.aveGapSTD/self.avePPmm, 
                     linestyle='None', marker = m[3], markersize=8, color = c[3],
                     label="Gap Width [mm]")
                     
        
        ax.errorbar(self.volumeFluxes2,  self.aveGap_capsule_height/self.avePPmm, 
                     yerr=self.aveGap_capsule_heightSTD/self.avePPmm, 
                     linestyle='None', marker = m[2], markersize=8, color = c[2],
                     label="Capsule Height [mm]")


        plt.errorbar(self.volumeFluxes2,  self.aveTopGap/self.avePPmm, 
                     yerr=self.aveTopGapSTD/self.avePPmm, 
                     linestyle='None', marker = m[0], markersize=8, color = c[0],
                     label="Top Gap Width [mm]")
                                          
        plt.title('%s Q vs Gap Spacing' %self.batchName)
        plt.xlabel('Pinching Flow Volume Flux [$ml/min$]');   plt.ylabel('Distance betweeen Capsule and Wall [$mm$]')  
        ax = ana._standardPlotSize(ax)        
        plt.legend(loc="lower center", fontsize=6, ncol=1)
        
        ax2 = ax.twinx()
        ax2.errorbar(self.volumeFluxes2,  self.aveGap_d12, 
                     yerr=self.aveGap_d12STD, 
                     linestyle='None', marker = m[0], markersize=8, color = c[0],
                     label="Deformation")
        ax2.set_ylabel("Deformation Parameter $D_{12}$")
        ax.legend(loc="LowerLeft", fontsize=6, ncol=1)
        ax2.legend(loc='best', fontsize=6, ncol=1)
        ax2 = ana._standardPlotSize(ax2)               
        ana._saveFig(fig, os.path.join(self.directory, 'RESLT') ,'%s_QvsGapSpacing2.jpg' %self.batchName)


    def _plotGapD12VsFinalOffset(self):
        '''Plot final distance from centreline versus volume flux Q'''
        m = ana.markerSymbols(); c = ana.coloursForMarker()
        fig = plt.figure(figsize=(8, 6), dpi=200,); ax = fig.add_subplot(111)
        
        plt.errorbar(self.aveGap_d12,  self.aveFinalDistCentre/self.avePPmm, 
                     yerr=self.aveFinalDistCentreSTD/self.avePPmm, 
                     linestyle='None', marker = m[3], markersize=8, color = c[3])
                     
        
        plt.title('%s Gap $D_{12}$ Vs $L_{fo}$' %self.batchName)
        plt.xlabel('Gap Deformation $D_{12}$');   plt.ylabel('Final Distance from Centreline $L_{fo}$ [$mm$]')  
        ax = ana._standardPlotSize(ax)        
#        plt.legend(loc="lower center", fontsize=6, ncol=1)

        ana._saveFig(fig, os.path.join(self.directory, 'RESLT') ,'%s_GapD12vsFinalOffset.jpg' %self.batchName)
        
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
        
    

        
        
        
        
        