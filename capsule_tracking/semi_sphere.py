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
import matplotlib.patches as mpatches
import cv2
import os
import time
#import scipy.stats
from scipy import stats
from scipy import optimize 

import general as gen
from general import filenameClass

import analysis as ana
#from analysis import DataRun

ERROR_CONST=-1

class ParametersSemiSphere(gen.Parameters):
    """Collect all information need for semi-sphere. Assumes inflow from right,
    towards diffuser on the left.
    

    ChannelTop         Pixel y-position of top of channel
    ChannelBottom      Pixel y-position of top of channel   
    EndChannel         Pixel x-position of end of channel, start of
                       diffuser
                       
    Instead of providing parameters, read them from file. 
    readParametersFromFile  
    
    """
    
    def __init__(self, listingImageFunc,coverImgFunc, baseDirectory= None, folder= None, d0= None, pPmm= None, 
                 rotation= None, FPS= None, backgroundImage= None, 
                 ChannelTop = None, ChannelBottom=None, EndChannel=None, 
                 readParametersFromFile=None):
        
        if baseDirectory != None:
            gen.Parameters.__init__(self, baseDirectory, folder, d0, pPmm, 
                 rotation, FPS, backgroundImage)
            self.ChannelTop = ChannelTop 
            self.ChannelBottom  =ChannelBottom 
            self.EndChannel = EndChannel
            self._writeToFile()
            
        elif readParametersFromFile != None:
            self._readFromFile(readParametersFromFile)
        
        self.setFunc(listingImageFunc,coverImgFunc)
        
        
        
        
    def _writeToFile(self):
        '''         Write a parameters file in path         '''
        pathParameters = self.path + 'ImageAnalysis_Parameters_'+time.strftime("%Y-%m-%d")+'.txt' 
        print('os.path.isdir(self.path):', end='');print(os.path.isdir(self.path))
        fileParameters = open(pathParameters, 'w')
        fileParameters.write('# Parameters for the analysis of the images. \n')
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
        if len(orderedLines) != 11:
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
            
            #semi-sphere specific            
            self.ChannelTop =  int(float(orderedLines[8]))
            self.ChannelBottom  = int(float(orderedLines[9]))
            self.EndChannel =  int(float(orderedLines[10]))
            
def sortPhotosPCO(path, fileType = '.png', prefixleng=10):
    """
    PCO camera
    
    Typical filename :     Batch20160205-1-#4-5mlPmin_5FPS_7_0000 
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
            
            
def coverSideSemiSphere(img, offset, ParameterClass, color=(255,255,255)):
    '''Take image and delete areas not of intereste'''

    top=ParameterClass.ChannelTop
    bottom=ParameterClass.ChannelBottom 
    end=ParameterClass.EndChannel
    
    #Get size of image
    try:
        yImg,xImg = img.shape
    except:
        yImg,xImg, _ = img.shape
    
#    cv2.rectangle(img, (x1, y1), (x2, y2), (255,0,0), 2)  
#    x1,y1 ------
#    |          |
#    |          |
#    |          |
#    --------x2,y2
    #cover area above and below main channel
    cv2.rectangle(img, (end+offset, top+offset), (xImg, 0), color, thickness=-1)   
    cv2.rectangle(img, (end+offset, bottom+offset), (xImg, yImg),  color, thickness=-1)
    
    #Cover outside
    cv2.rectangle(img, (0, offset+2), (xImg, 0), color, thickness=-1)   
#    cv2.rectangle(img, (end+offset, bottom+offset), (xImg, yImg),  255, thickness=-1)
    
    #through rotation, there is a black border that seems to throw off the auto-
    # matic threshold setting, add a 3 pixel wide margine
    crop=16
    cv2.rectangle(img, (0, 0), (crop, yImg), color, thickness=-1)  
    cv2.rectangle(img, (0, yImg), (xImg, yImg-crop), color, thickness=-1)  
    
    #left & right
    crop=4
    cv2.rectangle(img, (crop, 0), (0, yImg), color, thickness=-1)   
    cv2.rectangle(img, (xImg - crop, 0), (xImg, yImg), color, thickness=-1) 
    
    #cover traingle bottom
    shapeR = np.array([[end+offset-5,bottom+offset], 
                       [end+offset-5, yImg], 
                       [end+offset-5 - (yImg-bottom-offset), yImg]], np.int32)
    cv2.fillPoly(img, [shapeR], color)   
    
        
    #cover traingle top
    shapeR = np.array([[end+offset-5,top+offset], 
                       [end+offset-5, 0], 
                       [end+offset-5 - (top+offset), 0]], np.int32)
    cv2.fillPoly(img, [shapeR], color) 
    
    #time stampe
    cv2.rectangle(img, (50, 12), (300, 6), color, thickness=-1)  
    return img



# =============================================================================
# Extract relevant measurments from found data
# =============================================================================
def runFunction(path):
    '''Runs DataRunSemiSphere on path'''
    DataRunSemiSphere(path)
#    DRSS = DataRunSemiSphere(path)
#    DRSS._writeToFileSemiSphere()
#    DRSS.writeNameToFile()
#    del DRSS
    
class DataRunSemiSphere(ana.DataRun):
    
    def __init__(self, path):
        #load parameters from image analysis
        self.PSS = ParametersSemiSphere(listingImageFunc = sortPhotosPCO,
                                   coverImgFunc = coverSideSemiSphere,
                                   readParametersFromFile=path)
        
        if self.PSS.path != path:
            self.PSS.path = path
            self.PSS.baseDirectory = os.path.dirname(os.path.dirname(path))
            pa , folder = os.path.split(path)
            p , folder = os.path.split(pa)
            self.PSS.folder = folder
            
        #use base class to load file
        ana.DataRun.__init__(self, self.PSS.path, self.PSS.pPmm, self.PSS.FPS)
        
        self.plot=True
        
        self.extractMeasures()
        
    def extractMeasures(self):
        '''Extract measures of interest from the raw tracking data and write to 
        file. '''
        self._findAngle()
        self._findCentering()
        self._getFinalPosition()
        self._maxDeformation()  
        self._findFinalDiameter()
        self._writeToFileSemiSphere()
#        self._speedVsXpos()
        
    def writeNameToFile(self):
        '''Debug function'''
        directory = os.path.dirname(self.path)
        directory = os.path.dirname(directory)
        fileData = open(os.path.join(directory,'Names.txt'), 'a')
        fileData.write('%s \n' %(self.name))
        fileData.close()
        
    def _getFinalPosition(self, forPub=False, savePath=''):
        '''Finds final position of capsule after the diffuser and check whether
        capsule went right (up) or left (down) (from the capsules point of view.'''
        #Central line can be evaluated from PSS.ChannelTop  and PSS.ChannelBottom
        
        centreline = (self.PSS.ChannelTop + self.PSS.ChannelBottom)/2.0
        
        # length of the diffuser is 27mm and capsule diameter is PSS.d0 and 
        # Position of semi-sphere in channel is given PSS.endChannel
        # Getting x-position at which the capsule has stoped moving in the 
        # y-direction
#        xposAfterDiffuser = self.PSS.EndChannel - \
#            (self.PSS.d0/2.0 + 27.0)*self.PSS.pPmm
        xposAfterDiffuser = self.PSS.EndChannel - \
            (self.PSS.d0 + 27.0)*self.PSS.pPmm
        
        #half a diameter from the edge of the image to only capture full outlines
        minXpos = (self.PSS.d0/2.0)*self.PSS.pPmm
        #remove entries with no values
        cent_x, cent_y = self._centroidsWithoutErrors()
        ri = int(len(cent_x)/4.0) #relevant indexes are the last 4th
        cent_xs, cent_ys =cent_x[-ri:], cent_y[-ri:] 
        print(cent_xs)
        print('xposAfterDiffuser = %d' %(xposAfterDiffuser))
        print('minXpos = %d' %(minXpos))
        #Find if we have measurments at that position
        indexes = [i for i in range(len(cent_xs)) if cent_xs[i] < xposAfterDiffuser \
                    and cent_xs[i] > minXpos ]
        
        print('indexes = ', end=''); print(indexes)
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
        
        if np.isnan(self.distanceFromCentreline):
            self.distanceFromCentreline=-1
        
        if self.plot:
            yfit=[];
            for x in self.pos_diffuser_x:
                yfit.append(ana.func(x, self.diffuser_m, self.diffuser_c))
                    
            
            
            fig = plt.figure(figsize=(8, 6), dpi=200); ax = fig.add_subplot(111)
            if forPub:
                FS=18
            else:
                FS=12
            
            
            plt.plot(cent_x, cent_y, 'bs', markeredgecolor = 'none', label='Centroid position')
            plt.plot(cent_xs[indexes], cent_ys[indexes], 'ro', markeredgecolor = 'none', label='Centroid position used for final offset')
            plt.plot(self.intialCentreingXpos, self.intialCentreingYpos, 'g<', markeredgecolor = 'none', label='Centroid position used for initial offset')
            plt.plot(self.pos_diffuser_x, self.pos_diffuser_y, 'ch',markeredgecolor = 'none', label='Centroid position for angle measurment' )
            xmax=self.PSS.EndChannel
            ymin=ana.func(0.0, self.diffuser_m, self.diffuser_c)
            ymax = ana.func(xmax, self.diffuser_m, self.diffuser_c)
            if not forPub:
                plt.plot([0, xmax], [ymin, ymax], 'c-', label='Fit angle = %f' %(self.diffuser_angle))                       
            
            plt.plot([0, np.max(cent_x)], [centreline, centreline], 'k--', label='Centreline')
            
            if self.wentRight:
                dc = centreline - self.distanceFromCentreline
            else:
                dc = centreline + self.distanceFromCentreline
            
            ic = centreline - self.aveOffset
            
            plt.xlabel('Centroid x-position [pixels]', fontsize=FS); plt.ylabel('Centroid y-position [pixels]', fontsize=FS)
             
            if not forPub:
                plt.plot([0, np.max(cent_x)], [dc, dc], 'r-', label='Final Offset')
                plt.plot([0, np.max(cent_x)], [ic, ic], 'g-', label='Initial Offset')
               
    
                plt.text(0.1, 0.6, 'Final Offset from Centreline = %.3f mm' %(self.distanceFromCentreline/self.PSS.pPmm), fontsize = 16, horizontalalignment='left', verticalalignment='center', transform = ax.transAxes)
                plt.text(0.1, 0.4, 'Initial Centreing = %.3f mm' %(self.aveOffset/self.PSS.pPmm), fontsize = 16, horizontalalignment='left', verticalalignment='center', transform = ax.transAxes)
                
                plt.title("Centroid position -  " +self.name)
                plt.legend(loc='best', fontsize=6, ncol=2)

            
            ymin, ymax = ax.get_ylim()
            xmin, xmax = ax.get_xlim()
            currentAxis = plt.gca()
            xy = np.array([[self.PSS.EndChannel, self.PSS.ChannelBottom], 
                          [self.PSS.EndChannel - int(self.PSS.pPmm*20), self.PSS.ChannelBottom + int(self.PSS.pPmm*20)],
                          [self.PSS.EndChannel, self.PSS.ChannelBottom + int(self.PSS.pPmm*20)],
                          [self.PSS.EndChannel + int(self.PSS.pPmm*30), self.PSS.ChannelBottom + int(self.PSS.pPmm*20.0)],
                          [self.PSS.EndChannel + int(self.PSS.pPmm*30), self.PSS.ChannelBottom]])
            currentAxis.add_patch(mpatches.Polygon(xy, closed=True, facecolor="grey"))

            xy = np.array([[self.PSS.EndChannel, self.PSS.ChannelTop], 
                          [self.PSS.EndChannel - int(self.PSS.pPmm*20), self.PSS.ChannelTop - int(self.PSS.pPmm*20)],
                          [self.PSS.EndChannel, self.PSS.ChannelTop - int(self.PSS.pPmm*20)],
                          [self.PSS.EndChannel + int(self.PSS.pPmm*30), self.PSS.ChannelTop - int(self.PSS.pPmm*20.0)],
                          [self.PSS.EndChannel + int(self.PSS.pPmm*30), self.PSS.ChannelTop]])
            currentAxis.add_patch(mpatches.Polygon(xy, closed=True, facecolor="grey"))
                                            
            currentAxis.add_patch(mpatches.Circle((self.PSS.EndChannel, centreline ), 
                                            self.PSS.pPmm*4.0, facecolor="grey"))
            #Cover half of the circle
            currentAxis.add_patch(mpatches.Rectangle((self.PSS.EndChannel, int(centreline - 4.05*self.PSS.pPmm)), 
                                            - int(2.85*self.PSS.pPmm), int(8.15*self.PSS.pPmm), facecolor="white", edgecolor="none"))
            currentAxis.add_patch(mpatches.Rectangle((self.PSS.EndChannel- int(2.05*self.PSS.pPmm), int(centreline - 3.05*self.PSS.pPmm)), 
                                            - int(2.05*self.PSS.pPmm), int(6.15*self.PSS.pPmm), facecolor="white", edgecolor="none"))

            # Set Axis Limits
            plt.axes().set_aspect('equal', 'datalim')        
            fig.tight_layout()
#            ax.set_ylim((centreline + self.distanceFromCentreline + 3*self.PSS.pPmm, centreline - self.distanceFromCentreline - 3*self.PSS.pPmm))     
            ax.set_ylim(( 2*centreline -20, 10))
#            plt.ylim(( 2*centreline -20, 10))
            xmin, xmax = plt.xlim()
            plt.xlim(20, xmax)
            plt.show()
            
            if forPub:
                if savePath != '':
                    plt.savefig(ana._save_for_latex(savePath+"CentroidPosition_forpub_" +self.name+".jpg"), dpi=300)
                else:
                    plt.savefig(ana._save_for_latex((self.path+"CentroidPosition_forpub_" +self.name+".jpg")), dpi=300)
            else:
                plt.savefig((self.path+"CentroidPosition_" +self.name+".jpg"), dpi=300)

#            fig = plt.figure( dpi=200); ax = fig.add_subplot(111)
#            plt.plot(cent_x, cent_y, 'bo', markersize=3, label='Centroid position')
#            plt.axes().set_aspect('equal', 'datalim')
#            plt.savefig((self.path+"CentroidPosition_Plain" +self.name+".jpg"), dpi=300)            

    def _speedVsXpos(self, savePath=''):
        '''Plot speed against the x-position'''
        
        dispCutoff=10 #10 pixel max displacmeent between frames
        
        dx= np.gradient(self.centroid_x+0.0)
        dy = np.gradient(self.centroid_y+0.0)
                     
        speed = np.sqrt(np.power(dx, 2.0) + np.power(dy, 2.0))
        ImgNr = np.arange(len(speed))
        xpos = self.centroid_x
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
#        speed_ave =  gen.runningMean(speed_ave, aveWindow)         
#        speed_ave =  gen.runningMean(speed_ave, aveWindow)         
#        speed_ave =  gen.runningMean(speed_ave, aveWindow)  

        fig = plt.figure(figsize=(8, 6), dpi=200); ax = fig.add_subplot(111)
        plt.plot(ImgNr/self.FPS, speed_ave, 'bs')
        
        fig = plt.figure(figsize=(8, 6), dpi=200); ax = fig.add_subplot(111)
        plt.plot(xpos, speed_ave, 'bs', label='Centroid position')
        
        xmin, xmax = plt.xlim()
        plt.xlim(20, xmax)
        plt.xlabel('Centroid x-position [pixels]', fontsize=18); plt.ylabel('Capsule Speed [$mm/s$]', fontsize=18)

#        plt.text(0.1, 0.6, 'Final Offset from Centreline = %.3f mm' %(self.distanceFromCentreline/self.PSS.pPmm), fontsize = 16, horizontalalignment='left', verticalalignment='center', transform = ax.transAxes)
#        plt.text(0.1, 0.4, 'Initial Centreing = %.3f mm' %(self.aveOffset/self.PSS.pPmm), fontsize = 16, horizontalalignment='left', verticalalignment='center', transform = ax.transAxes)

#        ax.invert_xaxis()
#        plt.title("Speed -  " +self.name)
#        plt.legend(loc='best', fontsize=6, ncol=2)
        fig.tight_layout()
        plt.show()        
        if savePath == '':
            plt.savefig(ana._save_for_latex(self.path+"SpeedVsXpos_" +self.name+".jpg"), dpi=300)
        else:
            plt.savefig(ana._save_for_latex(savePath+"SpeedVsXpos_" +self.name+".jpg"), dpi=300)
        
        
    
    def _maxDeformation(self):
        '''Find max deformation, as measured by taylor deformation parameter'''
        
        #remove entries with no values
        cent_x, cent_y = self._centroidsWithoutErrors()
        
        #First find area around semi-sphere, which is at x-pos end-channel        
        startSS = self.PSS.EndChannel + 2*self.PSS.d0*self.PSS.pPmm
        endSS = self.PSS.EndChannel - 2*self.PSS.d0*self.PSS.pPmm
        indexes = [i for i in range(len(cent_x)) if cent_x[i] < startSS and cent_x[i] > endSS]
        
#        centreline = (self.PSS.ChannelTop + self.PSS.ChannelBottom)/2.0
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
        # Find the deformation between the end of semisphere and the beginning of
        # semisphere
        #
        #======================================================================
        startBeginningSS = self.PSS.EndChannel + (5.0 + self.PSS.d0/2.0) * self.PSS.pPmm
        endBeginningSS = self.PSS.EndChannel + 4.0 * self.PSS.pPmm
        indexesBeginningSS = [i for i in range(len(self.centroid_x)) if self.centroid_x[i] > endBeginningSS and self.centroid_x[i] < startBeginningSS]
        
        if len(indexesBeginningSS) >0 :
            self.maxd12BeginningSS = np.max(self.d12[indexesBeginningSS])
            self.maxWidthBeginningSS = np.max(self.width[indexesBeginningSS])
            self.minHeightBeginningSS = np.min(self.height[indexesBeginningSS])
            self.meanHeightBeginningSS = np.mean(self.height[indexesBeginningSS])
            self.meanWidthBeginningSS = np.mean(self.width[indexesBeginningSS])
            self.meand12BeginningSS = np.mean(self.d12[indexesBeginningSS])
            self.meand12rectBeginningSS = np.mean(self.d12rect[indexesBeginningSS])
        else:
            self.maxd12BeginningSS = -1; self.maxWidthBeginningSS = -1; self.minHeightBeginningSS = -1; self.meanHeightBeginningSS = -1; self.meanWidthBeginningSS =-1; self.meand12BeginningSS =-1; self.meand12rectBeginningSS = -1;
        
        startEndSS  = self.PSS.EndChannel + (0.5 ) * self.PSS.pPmm
        endEndSS  = self.PSS.EndChannel + 0.0 * self.PSS.pPmm             
        indexesEndSS  = [i for i in range(len(self.centroid_x)) if self.centroid_x[i] > endEndSS  and self.centroid_x[i] < startEndSS ]
        
        if len(indexesEndSS ) >0 :
            self.maxd12EndSS  = np.max(self.d12[indexesEndSS ])
            self.maxWidthEndSS  = np.max(self.width[indexesEndSS ])
            self.minHeightEndSS  = np.min(self.height[indexesEndSS ])
            self.meanHeightEndSS  = np.mean(self.height[indexesEndSS ])
            self.meanWidthEndSS  = np.mean(self.width[indexesEndSS ])
            self.meand12EndSS  = np.mean(self.d12[indexesEndSS ])
            self.meand12rectEndSS  = np.mean(self.d12rect[indexesEndSS ])
        else:
            self.maxd12EndSS  = -1; self.maxWidthEndSS  = -1; self.minHeightEndSS  = -1; self.meanHeightEndSS  = -1; self.meanWidthEndSS  =-1; self.meand12EndSS  =-1; self.meand12rectEndSS=-1;
        
        self.mean_xd_Beginning, self.mean_yd_Beginning = self._findCentralWidthAndHeight(indexesBeginningSS)
        self.mean_xd_End, self.mean_yd_End = self._findCentralWidthAndHeight(indexesEndSS, plot=True, pos=[self.PSS.ChannelTop, self.PSS.ChannelBottom])
        
        print('Beginning of Semi-Sphere: xd = %f yd=%f (%.2f mm) d_12=%f' %(self.mean_xd_Beginning, self.mean_yd_Beginning, self.mean_yd_Beginning/self.pixelsPmm, abs(self.mean_xd_Beginning - self.mean_yd_Beginning)/(self.mean_xd_Beginning+self.mean_yd_Beginning)))
        print('End of Semi-Sphere: xd = %f yd=%f (%.2f mm) d_12=%f' %(self.mean_xd_End, self.mean_yd_End, self.mean_yd_End/self.pixelsPmm ,abs(self.mean_xd_End - self.mean_yd_End)/(self.mean_xd_End+self.mean_yd_End)))
        
        #======================================================================
        # Plot central distance on images to check them
        #======================================================================
        
    def _findFinalDiameter(self):
        '''find diameter of capsule in diffuser'''        
        # Give capsule time to relaxe
        end = 2.0*self.PSS.d0*self.PSS.pPmm

        # lets not take the last half diameter as there might be artefacts
        # from capsule tracking as the capsule enters the image
        start =  0.5*self.PSS.d0*self.PSS.pPmm
        
        #Find if we have measurments at that position
        indexes = [i for i in range(len(self.centroid_x)) if self.centroid_x[i] < end and self.centroid_x[i] > start]
        
#        print('start = %d end = %d  and indexes = ' %(start, end), end=''); print(indexes)
        self.d, self.dSTD = self.fitCircle(indexes)
        print('Average d = %f =/- %f (based on %d outlines)' %(self.d/self.PSS.pPmm, self.dSTD/self.PSS.pPmm, len(indexes)))
        
            
    def _findCentering(self):
        '''find centering of capsule in channel leading up to semi-sphere'''
        
        #find centerline
        centreline = (self.PSS.ChannelTop + self.PSS.ChannelBottom)/2.0
        
        # We know the end of the semi-sphere. Its starts 4mm ealier and we would 
        # want to give the capsule at least 2 diameter 
        endSS = self.PSS.EndChannel + (4 + 2*self.PSS.d0)*self.PSS.pPmm
        
        #remove entries with no values
        cent_x, cent_y = self._centroidsWithoutErrors()
        
        #find max x position
        maxX = np.max(cent_x)
        
        # lets not take the first half diameter as there might be artefacts
        # from capsule tracking as the capsule enters the image
        startSS_IgnoreLast = maxX - 0.5*self.PSS.d0*self.PSS.pPmm
        
        #Take 3 capsule diameters
        startSS_3capsules = endSS + 3*self.PSS.d0*self.PSS.pPmm
        
        if startSS_IgnoreLast > startSS_3capsules:
            startSS =startSS_3capsules
        else:
            startSS = startSS_IgnoreLast
        
        #Find if we have measurments at that positionhttp://www.contributoria.com/issue/2014-05/5319c4add63a707e780000cd/
        indexes = [i for i in range(len(cent_x)) if cent_x[i] > endSS and cent_x[i] < startSS]
        
        meanYpos=cent_y[indexes]
        
        offset = centreline - meanYpos
        
        # ====================================================================        
        # Also find the velocity in the initial channel
        # ====================================================================
        
        if len(cent_x[indexes]) >1 :
            dispCutoff=25 #pixel max displacmeent between frames
            
            dx= np.gradient(cent_x[indexes]+0.0)
            dy = np.gradient(cent_y[indexes]+0.0)
                    
            speed = np.sqrt(np.power(dx, 2.0) + np.power(dy, 2.0))
            ImgNr = np.arange(len(speed))
            
            delIndex=[]
            for uu in reversed(range(len(speed))):
                if abs(speed[uu]) > dispCutoff:
                    delIndex.append(uu)
    
            speed = np.delete(speed, delIndex)/self.pixelsPmm  * self.FPS
            ImgNr = np.delete(ImgNr, delIndex)
            
            self.initial_speed  = np.average(speed)
        else:
            self.initial_speed = ERROR_CONST
#        plt.figure()
#        plt.plot(cent_x, cent_y, 'bs')
#        plt.plot(cent_x[indexes], cent_y[indexes], 'ro')
#        plt.plot([0, np.max(cent_x)], [centreline, centreline], 'k-')
#        plt.xlabel('x-pos [pixels]'); plt.ylabel('y-pos [pixels]')
#        
#        plt.figure()
#        plt.plot(cent_x[indexes]/self.PSS.pPmm, offset/self.PSS.pPmm, 'bs')
#        plt.xlabel('x-pos [mm]'); plt.ylabel('Offset from Centreline [mm]')
        self.intialCentreingXpos=cent_x[indexes]
        self.intialCentreingYpos=cent_y[indexes]
        
        self.aveOffset = np.mean(offset)
        self.aveOffsetSTD = np.std(offset)
        

    def _findAngle(self):
        '''find angle capsule takes away from obstacle'''
        
        # We know the end of the semi-sphere. Its starts 4mm ealier and we would 
        # want to give the capsule at least 2 diameter 
        start = self.PSS.EndChannel - (2.2*self.PSS.d0)*self.PSS.pPmm
        
        #The diffuser is 27mm long, lets neglect the last two capsule diameters
        end = self.PSS.EndChannel + (2*self.PSS.d0)*self.PSS.pPmm - 27.0*self.PSS.pPmm
        
        print('end = %d start = %d' %(end, start))
        
        #remove entries with no values
        cent_x, cent_y = self._centroidsWithoutErrors()
        
        #Find if we have measurments at that position
        indexes = [i for i in range(len(cent_x)) if cent_x[i] > end and cent_x[i] < start]
        
        pos_diffuser_x = cent_x[indexes]
        pos_diffuser_y = cent_y[indexes]
        
        ave_y = np.mean(pos_diffuser_y)
        
        pdx=[]; pdy=[];
        for kk in range(len(pos_diffuser_y)):
            if abs((pos_diffuser_y[kk]-ave_y)/ave_y) < 0.1:
                pdx.append(pos_diffuser_x[kk])
                pdy.append(pos_diffuser_y[kk])
                
        self.pos_diffuser_x = np.array(pdx)
        self.pos_diffuser_y = np.array(pdy)
        
        #        plt.figure()_plotMeanD12(s
#        plt.plot(cent_x, cent_y, 'bs')
#        plt.plot(self.pos_diffuser_x, self.pos_diffuser_y, 'ro')

        p1, p2, fres = ana.linearFit( self.pos_diffuser_x,self.pos_diffuser_y)
        
        
        self.diffuser_m=p1 
        self.diffuser_c = p2
        self.diffuser_angle = np.arctan(abs(p1))#Take the absolute value of the gradient as the capsule can go up or down
        
        
        
    def _writeToFileSemiSphere(self):
        '''Write measurments of run to file'''
        headerString = '\tFinal Distance From Centreline \t STD \tmax d12 \t maxWidth \tmin height \tAve Offset \t Ave Offset STD \t angle \tFinal Diameter \tDiameter STD \t maxd12BeginningSS \tself.maxWidthBeginningSS \tminHeightBeginningSS \tmeanHeightBeginningSS \tmeanWidthBeginningSS \t meand12BeginningSS \t self.maxd12EndSS  \tmaxWidthEndSS  \tminHeightEndSS  \tmeanHeightEndSS  \tmeanWidthEndSS  \tself.meand12EndSS \tmeand12rectBeginningSS \tmeand12rectEndSS \tmean_xd_Beginning \tmean_yd_Beginning \tmean_xd_End \tmean_yd_End \tself.wentRight \tself.initial_speed'
        String = '\t%.6f \t%.6f \t%.6f \t%d \t%d \t%.6f \t%.6f \t%.6f \t%.6f \t%.6f \t%.6f \t%.6f \t%.6f \t%.6f \t%.6f \t%.6f \t%.6f \t%.6f \t%.6f \t%.6f  \t%.6f \t%.6f \t%.6f \t%.6f \t%.6f \t%.6f \t%.6f \t%.6f  \t%r \t%.6f' %(self.distanceFromCentreline, self.distanceFromCentrelineSTD, self.maxd12, self.maxWidth,  self.minHeight, self.aveOffset, self.aveOffsetSTD, self.diffuser_angle, self.d, self.dSTD, self.maxd12BeginningSS , self.maxWidthBeginningSS , self.minHeightBeginningSS , self.meanHeightBeginningSS , self.meanWidthBeginningSS, self.meand12BeginningSS, self.maxd12EndSS  , self.maxWidthEndSS  , self.minHeightEndSS  , self.meanHeightEndSS  , self.meanWidthEndSS, self.meand12EndSS, self.meand12rectBeginningSS, self.meand12rectEndSS, self.mean_xd_Beginning, self.mean_yd_Beginning, self.mean_xd_End, self.mean_yd_End, self.wentRight, self.initial_speed)        
        self._writeToFile(headerString, String)
        
        
        
class ResultsClassSemiSphere(ana.ResultsClass):
    '''Holds experiment wide results for Semi Sphere Setup'''
    
    def __init__(self, directory):
        ana.ResultsClass.__init__(self,directory)
        self.finalDistCentre = np.zeros(len(self.volumeFlux))
        self.finalDistCentreSTD = np.zeros(len(self.volumeFlux))
        self.maxD12 = np.zeros(len(self.volumeFlux))
        self.maxWidth = np.zeros(len(self.volumeFlux))
        self.minHeight = np.zeros(len(self.volumeFlux))
        self.Offset = np.zeros(len(self.volumeFlux))
        self.OffsetSTD = np.zeros(len(self.volumeFlux))
        self.angle = np.zeros(len(self.volumeFlux))
        
        self.d = np.zeros(len(self.volumeFlux))
        self.dSTD = np.zeros(len(self.volumeFlux))
        
        self.maxd12BeginningSS = np.zeros(len(self.volumeFlux))
        self.maxWidthBeginningSS = np.zeros(len(self.volumeFlux))
        self.minHeightBeginningSS = np.zeros(len(self.volumeFlux))
        self.meanHeightBeginningSS = np.zeros(len(self.volumeFlux))
        self.meanWidthBeginningSS = np.zeros(len(self.volumeFlux))
        self.meand12BeginningSS = np.zeros(len(self.volumeFlux))

        self.maxd12EndSS  = np.zeros(len(self.volumeFlux))
        self.maxWidthEndSS  = np.zeros(len(self.volumeFlux))
        self.minHeightEndSS  = np.zeros(len(self.volumeFlux))
        self.meanHeightEndSS  = np.zeros(len(self.volumeFlux))
        self.meanWidthEndSS  = np.zeros(len(self.volumeFlux))
        self.meand12EndSS  = np.zeros(len(self.volumeFlux))
        
        self.meand12rectBeginningSS  = np.zeros(len(self.volumeFlux))
        self.meand12rectEndSS  = np.zeros(len(self.volumeFlux))
        
        self.mean_xd_BeginningSS  = np.zeros(len(self.volumeFlux))
        self.mean_yd_BeginningSS  = np.zeros(len(self.volumeFlux))
        self.mean_d12xy_BeginningSS  = np.zeros(len(self.volumeFlux))

        self.mean_xd_EndSS  = np.zeros(len(self.volumeFlux))
        self.mean_yd_EndSS  = np.zeros(len(self.volumeFlux))
        self.mean_d12xy_EndSS  = np.zeros(len(self.volumeFlux))
        
        self.wentRight  = np.zeros(len(self.volumeFlux), np.bool)
        
        self.inital_velocities   = np.zeros(len(self.volumeFlux))
        
        self._readDataSemiSphere()
        
        self.PSS = ParametersSemiSphere(listingImageFunc = sortPhotosPCO,
                                   coverImgFunc = coverSideSemiSphere,
                                   readParametersFromFile=os.path.join(self.directory, self.name[0] ))
            
    def _readDataSemiSphere(self):
        '''Read in specific data, see above DataRunSemiSphere._writeToFileSemiSphere'''
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
                self.finalDistCentre[ii] = (float(entries[4]))
                self.finalDistCentreSTD[ii] = (float(entries[5]))
                self.maxD12[ii] = (float(entries[6]))
                self.maxWidth[ii] = (float(entries[7]))
                self.minHeight[ii] = (float(entries[8]))
                self.Offset[ii] = (float(entries[9]))
                self.OffsetSTD[ii] = (float(entries[10]))
                self.angle[ii] = (float(entries[11]))
                ii += 1
            if len(entries)>=32:
                self.d[ii-1] = (float(entries[12]))
                self.dSTD[ii-1] = (float(entries[13]))
        
                self.maxd12BeginningSS[ii-1] = (float(entries[14]))
                self.maxWidthBeginningSS[ii-1] = (float(entries[15]))
                self.minHeightBeginningSS[ii-1] = (float(entries[16]))
                self.meanHeightBeginningSS[ii-1] = (float(entries[17]))
                self.meanWidthBeginningSS[ii-1] = (float(entries[18]))
                self.meand12BeginningSS[ii-1] = (float(entries[19]))
        
                self.maxd12EndSS[ii-1] = (float(entries[20]))
                self.maxWidthEndSS[ii-1] = (float(entries[21]))
                self.minHeightEndSS[ii-1] = (float(entries[22]))
                self.meanHeightEndSS[ii-1] = (float(entries[23]))
                self.meanWidthEndSS[ii-1] = (float(entries[24]))
                self.meand12EndSS[ii-1] = (float(entries[25]))

                self.meand12rectBeginningSS[ii-1] = (float(entries[26]))                
                self.meand12rectEndSS[ii-1] = (float(entries[27]))
                
                self.mean_xd_BeginningSS[ii-1]  = (float(entries[28]))
                self.mean_yd_BeginningSS[ii-1]  = (float(entries[29]))
                self.mean_d12xy_BeginningSS[ii-1]  = abs(float(entries[28]) - float(entries[29]))/(float(entries[28]) + float(entries[29]))

                self.mean_xd_EndSS[ii-1]  = (float(entries[30]))
                self.mean_yd_EndSS[ii-1]  = (float(entries[31]))
                self.mean_d12xy_EndSS[ii-1]  = abs(float(entries[30]) - float(entries[31]))/(float(entries[30]) + float(entries[31]))
            if len(entries) >= 34:
#                self.wentRight[ii-1] = self._boolConver(entries[32])
                #TODO: fix this and re-run analyisis for runs
                self.inital_velocities[ii-1] = (float(entries[33]))
                    
    def _boolConver(self, string):
        if string[0]=='T':
            return True
        elif string[0]=='F':
            return False
        else:
            raise NameError('Bool string neither ture or false')
                
    def plotRslts(self):
        '''Top level function plotting results'''
        self._findAverages()
        self._plotQVsFinalOffset()
        self._plotQVsInitialOffset()
#        self._plot1oQVsOffset()
        self._plotQVsD12()
        
        self._plotAbsFinalVsInitialOffset()
        self._plotNormedFinalVsInitialOffset()
        
#        self._plotFinalVsInitialOffset()
        self._plotmCvsO()
        self._plotcCvsO()
        self._plotQVsFinalOffsetByStreamline()
        self._plotQVsAngle()
        self._plotMeanD12()
        self._plotMeanD12Rect()
        self._plotMeanD12xy()
        self._plotDiameter()

    def plotRslts_small(self):
        '''Top level function plotting results'''
        self._findAverages()
        self._plotQVsFinalOffset()
        self._plotQVsInitialOffset()
#        self._plot1oQVsOffset()
        self._plotQVsD12()
        
#        self._plotAbsFinalVsInitialOffset()
#        self._plotFinalVsInitialOffset()
#        self._plotmCvsO()
#        self._plotcCvsO()
        self._plotQVsFinalOffsetByStreamline()
        self._plotQVsAngle()
#        self._plotMeanD12()
#        self._plotMeanD12Rect()
#        self._plotMeanD12xy()
#        self._plotDiamester()
        
    def _findAverages(self):
        '''Find averages of data loaded from file'''
        self.avePPmm, _ = self.averageByQ(self.pixelsPmm) 
        
        self.aveFinalDistCentre, self.aveFinalDistCentreSTD = self.averageByQ(self.finalDistCentre, grZero=True) 
 
        self.aveMaxD12, self.aveMaxD12STD = self.averageByQ(self.maxD12)
        self.aveMaxWidth, self.aveMaxWidthSTD = self.averageByQ(self.maxWidth)
        self.aveOffset, self.aveOffsetSTD = self.averageByQ(self.Offset)
        self.aveAngle, self.aveAngleSTD = self.averageByQ(self.angle)
#        self.aveIntialCentring, self.aveIntialCentringSTD = self.averageByQ(self.angle)
        
        self.aveMinHeightEndSS, self.aveMinHeightEndSSSTD  = self.averageByQ(self.minHeightEndSS, grZero=True)
        self.aveMeanHeightEndSS, self.aveMeanHeightEndSSSTD  = self.averageByQ(self.meanHeightEndSS, grZero=True)
        
        self.aveMean_yd_EndSS, self.aveMean_yd_EndSSSTD  = self.averageByQ(self.mean_yd_EndSS, grZero=True)
        
        self.aveInitialVelocities, self.aveInitialVelocitiesSTD  = self.averageByQ(self.inital_velocities, grZero=True)

        
    def _plotQVsFinalOffset(self, fit=False, mode='Linear'):
        '''Plot final distance from centreline versus volume flux Q'''
        m = ana.markerSymbols(); c = ana.coloursForMarker()
        fig = plt.figure(figsize=(8, 6), dpi=200,); ax = fig.add_subplot(111)
#        fDC = [kk for kk in self.finalDistCentre if kk > 0.0]
        plt.errorbar(self.volumeFlux, self.finalDistCentre/self.pixelsPmm, 
                     #yerr=self.finalDistCentreSTD/self.pixelsPmm,
                    linestyle='None', marker = m[0], color = c[0])
                    
        plt.errorbar(self.volumeFluxes, self.aveFinalDistCentre/self.avePPmm, 
                     yerr=self.aveFinalDistCentreSTD/self.avePPmm, 
                     linestyle='None', marker = m[3], markersize=8, color = c[3])
                     
        if fit:
            x=[]; y=[]
            for ss in range(len(self.finalDistCentre)):
                if self.finalDistCentre[ss] > 0.0:
                    x.append(self.volumeFlux[ss])
                    y.append(self.finalDistCentre[ss]/self.pixelsPmm[ss])
            x = np.array(x);             y = np.array(y);
            if mode.lower() == 'linear':
                slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
                plt.plot([0.0, np.max(x)], 
                          [slope*0.0+intercept, slope*np.max(x)+intercept], 
                          '--k', 
                          label='Linear Fit with slope = %f and intercept = %f' %(slope, intercept) )                
                print("In Q vs Final Offset plot: linear fit results in slope = %f and intercept = %f" %(slope, intercept) )    
                print("std_err: %f" %std_err)
                print("r-squared: %f" %r_value**2)
                
                my, by, ry, smy, sby   = ana.lsqfity(x, y)
                print("2nd way: linear fit results in slope = %f and intercept = %f" %(my, by) )    
                print("ry = %f, smy = %f, sby = %f" %(ry, smy, sby))
            
            elif mode.lower() == 'exp':
                fitParams, fitCovariance = optimize.curve_fit(ana.fitFuncExp, x, y)
                sigma = [fitCovariance[0,0], \
                        fitCovariance[1,1], \
                        fitCovariance[2,2], \
                        ]

                error = [] 
                for i in range(len(fitCovariance)):
                    error.append(np.absolute(fitCovariance[i][i])**0.5)
                perr_leastsq = np.array(error)
                #from:
                # http://stackoverflow.com/a/21844726/1826893
                
                #get r^2
                residuals = y - ana.fitFuncExp(x, fitParams[0], fitParams[1] , fitParams[2])
                ss_res = np.sum(residuals**2)
                ss_tot = np.sum((y-np.mean(y))**2)
                r_squared = 1 - (ss_res / ss_tot)

                print("_plotQVsFinalOffset - fitting exp: y =  c * np.exp(-b*x) + a")
                print("a = %f pm %f \t b = %f pm %f \t c = %f pm %f" %(fitParams[0], perr_leastsq[0], fitParams[1], perr_leastsq[1], fitParams[2], perr_leastsq[2]))
                print("ss_res = %f r^2 = %f" %(ss_res, r_squared))
                xfit = np.arange(0.0, np.max(x))
                yfit = ana.fitFuncExp(xfit, fitParams[0], fitParams[1] , fitParams[2])
                plt.plot(xfit, yfit, '--k')
                plt.plot([0, np.max(x)], [fitParams[0] +fitParams[2],fitParams[0] +fitParams[2]], '-.b')
                
            elif mode.lower() == 'sqrt':
                fitParams, fitCovariance = optimize.curve_fit(ana.fitFuncSqrt, x, y)
                error = [] 
                for i in range(len(fitCovariance)):
                    error.append(np.absolute(fitCovariance[i][i])**0.5)
                perr_leastsq = np.array(error)
                #from:
                # http://stackoverflow.com/a/21844726/1826893
                
                #get r^2
                residuals = y - ana.fitFuncSqrt(x, fitParams[0], fitParams[1])
                ss_res = np.sum(residuals**2)
                ss_tot = np.sum((y-np.mean(y))**2)
                r_squared = 1 - (ss_res / ss_tot)

                print("_plotQVsFinalOffset - fitting sqrt: y =  b * np.sqrt(x) + a")
                print("a = %f pm %f \t b = %f pm %f " %(fitParams[0], perr_leastsq[0], fitParams[1], perr_leastsq[1]))
                print("ss_res = %f r^2 = %f" %(ss_res, r_squared))
                xfit = np.arange(0.0, np.max(x))
                yfit = ana.fitFuncSqrt(xfit, fitParams[0], fitParams[1])
                plt.plot(xfit, yfit, '--k')
                plt.plot([0, np.max(x)], [fitParams[0] +fitParams[1],fitParams[0] +fitParams[1]], '-.b')
            else:                
                print("mode not defined")
                
        plt.title('%s Q vs Offset' %self.batchName)
        plt.xlabel('Volume Flux [$ml/min$]');   plt.ylabel('Final Offset [$mm$]')  
        x1,x2,y1,y2 = plt.axis()
        plt.axis((x1,x2,0,y2))
        ax = ana._standardPlotSize(ax)        
        
        ana._saveFig(fig, os.path.join(self.directory, 'RESLT') ,'%s_QvsOffset.jpg' %self.batchName)
        

    def _plotQVsVelocity(self):
        '''Plot final distance from centreline versus volume flux Q'''
        m = ana.markerSymbols(); c = ana.coloursForMarker()
        fig = plt.figure(figsize=(8, 6), dpi=200,); ax = fig.add_subplot(111)
#        fDC = [kk for kk in self.finalDistCentre if kk > 0.0]
        plt.errorbar(self.inital_velocities, self.finalDistCentre/self.pixelsPmm, 
                     #yerr=self.finalDistCentreSTD/self.pixelsPmm,
                    linestyle='None', marker = m[0], color = c[0])
                    
        plt.errorbar(self.aveInitialVelocities, self.aveFinalDistCentre/self.avePPmm, 
                     xerr = self.aveInitialVelocitiesSTD,
                     yerr=self.aveFinalDistCentreSTD/self.avePPmm, 
                     linestyle='None', marker = m[3], markersize=8, color = c[3])
                    
        plt.title('%s Velocity vs Offset' %self.batchName)
        plt.xlabel('Velocity in inital Channel [$mm/s$]');   plt.ylabel('Final Offset [$mm$]')  
        x1,x2,y1,y2 = plt.axis()            
        plt.axis((0,x2,0,y2))
        ax = ana._standardPlotSize(ax)        
        
        ana._saveFig(fig, os.path.join(self.directory, 'RESLT') ,'%s_VelocityvsOffset.jpg' %self.batchName)
                
    def _sortByCentering(self, width):
        '''Select the value within width of centreline'''
        vf=[]; finalDist=[]; 
        for ii in range(len(self.volumeFlux)):
            if np.abs(self.Offset[ii]/self.pixelsPmm[ii]) < width:
                vf.append(self.volumeFlux[ii])
                finalDist.append(self.finalDistCentre[ii]/self.pixelsPmm[ii])
        
        return vf, finalDist
        
    def _plotQVsFinalOffsetByStreamline(self):
        '''Plot final distance from centreline versus volume flux Q'''
        m = ana.markerSymbols(); c = ana.coloursForMarker()
        fig = plt.figure(figsize=(8, 6), dpi=200,); ax = fig.add_subplot(111)
        
        vf1, finalDist1 =self._sortByCentering(0.1)
        vf2, finalDist2 =self._sortByCentering(0.2)
        vf3, finalDist3 =self._sortByCentering(0.3)
        plt.errorbar(self.volumeFlux, self.finalDistCentre/self.pixelsPmm, 
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
        
    def _plot1oQVsOffset(self, fit=False):
        '''Plot final distance from centreline versus volume flux Q'''
        m = ana.markerSymbols(); c = ana.coloursForMarker()
        fig = plt.figure(figsize=(8, 6), dpi=200,); ax = fig.add_subplot(111)
        
        plt.errorbar(1.0/self.volumeFlux, self.finalDistCentre/self.pixelsPmm, 
                     #yerr=self.finalDistCentreSTD/self.pixelsPmm,
                    linestyle='None', marker = m[0], color = c[0])
                    
        plt.errorbar(1.0/self.volumeFluxes, self.aveFinalDistCentre/self.avePPmm, 
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
        
        plt.plot(self.volumeFlux, self.maxD12, 
                    linestyle='None', marker = m[0], color = c[0])
        plt.xlabel('Volume Flux [$ml/min$]');   
        plt.ylabel('Taylor deformation parameter $D_{12}$')  
        
        
        locs, labels =plt.xticks()
        ax2 = ax.twinx()
        maxW = (self.maxWidth/self.pixelsPmm) / self.PSS.d0
        ax2.plot(self.volumeFlux, maxW,
                linestyle='None', marker = m[3], color = c[3])
        ax2.set_ylabel('Max Width [$d_0$]')
        
        plt.title('%s Q vs $D_{12}$' %self.batchName)
        
        ax = ana._standardPlotSize(ax)        
        
        ana._saveFig(fig, os.path.join(self.directory, 'RESLT') ,'%s_QvsD12.jpg' %self.batchName)
        
    def _plotMeanD12(self):
        '''Plot mean d12 of capsules'''
        m = ana.markerSymbols(); c = ana.coloursForMarker()
        fig = plt.figure(figsize=(8, 6), dpi=200,); ax = fig.add_subplot(111)

        plt.errorbar(self.volumeFlux, self.meand12BeginningSS/self.pixelsPmm, 
                     yerr=1.0/self.pixelsPmm, 
                     label= 'Mean $D_{12}$ at beginning of Semi-Sphere ',
                    linestyle='None', marker = m[0], color = c[0])
                    
        plt.errorbar(self.volumeFlux, self.meand12EndSS/self.pixelsPmm, 
                     yerr=1.0/self.pixelsPmm, 
                     label= 'Mean $D_{12}$ at end of Semi-Sphere',
                    linestyle='None', marker = m[3], color = c[3])

        plt.title('%s Q vAngle in diffuser versus volume flux $Q$ for a near-rigid gel bead and capsules of varying elasticity.s Mean $D_{12}$' %self.batchName)
        plt.xlabel('Volume Flux [$ml/min$]');   plt.ylabel('Mean $D_{12}$ of Capsule')  
        plt.legend(loc='best', fontsize=8)
        ax = ana._standardPlotSize(ax)      

        
        ana._saveFig(fig, os.path.join(self.directory, 'RESLT') ,'%s_QvsMeanD12.jpg' %self.batchName)
        
    def _plotMeanD12Rect(self):
        '''Plot mean d12 of capsules'''
        m = ana.markerSymbols(); c = ana.coloursForMarker()
        fig = plt.figure(figsize=(8, 6), dpi=200,); ax = fig.add_subplot(111)

        plt.errorbar(self.volumeFlux, self.meand12rectBeginningSS/self.pixelsPmm, 
                     yerr=1.0/self.pixelsPmm, 
                     label= 'Mean $D_{12}$ at beginning of Semi-Sphere ',
                    linestyle='None', marker = m[0], color = c[0])
                    
        plt.errorbar(self.volumeFlux, self.meand12rectEndSS/self.pixelsPmm, 
                     yerr=1.0/self.pixelsPmm, 
                     label= 'Mean $D_{12}$ at end of Semi-Sphere',
                    linestyle='None', marker = m[3], color = c[3])

#        newax.axis('off')
                    
#        plt.errorbar(self.volumeFluxes2, self.aveMinHeight/self.avePPmm, 
#                     yerr=self.aveMinHeightSTD/self.avePPmm, 
#                     label= 'Average +/- Standard Deviation',
#                     linestyle='None', marker = m[3], markersize=8, color = c[3])

        plt.title('%s Q vs Mean $D_{12}$ Rect' %self.batchName)
        plt.xlabel('Volume Flux [$ml/min$]');   plt.ylabel('Mean $D_{12}$ Rect of Capsule')  
        plt.legend(loc='best', fontsize=8)
        ax = ana._standardPlotSize(ax)      
        
        ana._saveFig(fig, os.path.join(self.directory, 'RESLT') ,'%s_QvsMeanD12Rect.jpg' %self.batchName)
                
    def _plotMeanD12xy(self):
        '''Plot mean d12 of capsules'''
        m = ana.markerSymbols(); c = ana.coloursForMarker()
        fig = plt.figure(figsize=(8, 6), dpi=200,); ax = fig.add_subplot(111)
        
        plt.errorbar(self.volumeFlux, self.mean_d12xy_BeginningSS/self.pixelsPmm, 
                     yerr=1.0/self.pixelsPmm, 
                     label= 'Mean $D_{12}$ at beginning of Semi-Sphere ',
                    linestyle='None', marker = m[0], color = c[0])
                    
        plt.errorbar(self.volumeFlux, self.mean_d12xy_EndSS/self.pixelsPmm, 
                     yerr=1.0/self.pixelsPmm, 
                     label= 'Mean $D_{12}$ at end of Semi-Sphere',
                    linestyle='None', marker = m[3], color = c[3])
                    

#        newax.axis('off')
                    
#        plt.errorbar(self.volumeFluxes2, self.aveMinHeight/self.avePPmm, 
#                     yerr=self.aveMinHeightSTD/self.avePPmm, 
#                     label= 'Average +/- Standard Deviation',
#                     linestyle='None', marker = m[3], markersize=8, color = c[3])

        plt.title('%s Q vs Mean $D_{12}$ xy' %self.batchName)
        plt.xlabel('Volume Flux [$ml/min$]');   plt.ylabel('Mean $D_{12}$ xy of Capsule')  
        plt.legend(loc='best', fontsize=8)
        ax = ana._standardPlotSize(ax)      
        
        ana._saveFig(fig, os.path.join(self.directory, 'RESLT') ,'%s_QvsMeanD12xy.jpg' %self.batchName)


    def _plotDiameter(self):
        '''Plot min height of capsules'''
        m = ana.markerSymbols(); c = ana.coloursForMarker()
        fig = plt.figure(figsize=(8, 6), dpi=200,); ax = fig.add_subplot(111)
        
        plt.errorbar(self.volumeFlux, self.d/self.pixelsPmm, 
                     yerr=self.dSTD/self.pixelsPmm, 
                     label= 'Runs +/- 1 pixel ',
                    linestyle='None', marker = m[0], color = c[0])
                    
        vf2=[]; d=[]; dSTD=[]; tol=0.2;        
        
        for ii in range(len(self.d)):
            td = self.d[ii]
            tdSTD = self.dSTD[ii]
            
            if abs(self.PSS.d0*self.pixelsPmm[0] - td)/(self.PSS.d0*self.pixelsPmm[0]) < tol:
                vf2.append(self.volumeFlux[ii])
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
    
    
    def _plotAbsFinalVsInitialOffset(self, forPub=False):
        '''Plot initial offset from centreline in the channel leading to the 
        semi-sphere versus volume flux Q'''      
        
        p1, p2, fres = self._linearFit(np.abs(self.Offset/self.pixelsPmm),self.finalDistCentre/self.pixelsPmm)
        
        m = ana.markerSymbols(); c = ana.coloursForMarker()
        
        if forPub:
            FS = 18; MS= 10;
        else:
            FS=12; MS=8;
            
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
        for vf in self.volumeFluxes:
            fig = plt.figure(figsize=(8, 6), dpi=200,); ax = fig.add_subplot(111)      
            tempx=[]; tempy=[];             err_tempx=[]; err_tempy=[]
            for i in range(len(self.volumeFlux)):
                if self.volumeFlux[i] == vf:
                    tempx.append(self.Offset[i]/self.pixelsPmm[i])
                    tempy.append(self.finalDistCentre[i]/self.pixelsPmm[i])
                    err_tempx.append(self.OffsetSTD[i]/self.pixelsPmm[i])
                    err_tempy.append(self.finalDistCentreSTD[i]/self.pixelsPmm[i])

            tempx = np.array(tempx); tempy = np.array(tempy)
            err_tempx = np.array(tempx); err_tempy = np.array(tempy)
            
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
#                        yerr=err_tempy, xerr=err_tempx,
                        markersize = MS,
                        linestyle='None', marker = m[0], color = c[0])
                        
            labelString = 'Linear Best Fit with m= %.4f and c = %.4f (RSS = %.4f)' %(my, by, ry)

            plt.plot([0.0, np.max(np.abs(tempx))], [ana.func(0.0, my,by), 
                      ana.func(np.max(np.abs(tempx)), my,by)], 'k-', 
                    label=labelString)
            if not forPub:
                plt.legend(loc='best', fontsize=6, ncol=1)

                textstring = 'Number of Runs = %d' %(len(tempx))
                plt.text(0.05, 0.4, textstring, fontsize = 16, 
                     horizontalalignment='left', verticalalignment='center', 
                     transform = ax.transAxes)
                plt.title('%s Initial vs Final Offset, $Q = %.f \ ml/min$' %(self.batchName, vf))
            ax = ana._standardPlotSize(ax)        
            plt.ylabel('Final Offset $|L_{fo}|$ [mm]', fontsize=FS);   plt.xlabel('Initial Centring $|L_{ic}|$ [mm]', fontsize=FS)  
            
            ana._saveFig(fig, os.path.join(self.directory, 'RESLT') ,'%s_Q=%.f_AbsInitiaVsFinalOffset.jpg' %(self.batchName, vf))
            
            #I want to plot the actual distance, not its absolute value so I need
            # to use the variable wentRight to determine that. 
            
            fig = plt.figure(figsize=(8, 6), dpi=200,); ax = fig.add_subplot(111)      
            tempx=[]; tempy=[]
            for i in range(len(self.volumeFlux)):
                if self.volumeFlux[i] == vf:
                    if self.wentRight[i]: mp=-1.0
                    else: mp =1.0
                    tempx.append(self.Offset[i]/self.pixelsPmm[i])
                    tempy.append(mp*self.finalDistCentre[i]/self.pixelsPmm[i])
            tempx = np.array(tempx); tempy = np.array(tempy)
                      
            plt.errorbar(np.abs(tempx), tempy,
                        #yerr=self.finalDistCentreSTD/self.pixelsPmm, 
                        #xerr=self.OffsetSTD/self.pixelsPmm,
                        linestyle='None', marker = m[0], color = c[0])

            textstring = 'Number of Runs = %d' %(len(tempx))
            plt.text(0.05, 0.4, textstring, fontsize = 16, 
                     horizontalalignment='left', verticalalignment='center', 
                     transform = ax.transAxes)
            plt.title('%s Initial vs Final Offset, $Q = %.f \ ml/min$' %(self.batchName, vf))
            plt.ylabel('Final Offset [$mm$]');   plt.xlabel('Initial Centering [$mm$]')  
            ax = ana._standardPlotSize(ax)        
            
            ana._saveFig(fig, os.path.join(self.directory, 'RESLT') ,'%s_Q=%.f_InitiaVsFinalOffset.jpg' %(self.batchName, vf))
            counter +=1
            
    def _plotmCvsO(self):
        ''' Plot the gradient of the linear fit to the initial centreing against
        final offset'''
        m = ana.markerSymbols(); c = ana.coloursForMarker()
        fig = plt.figure(figsize=(8, 6), dpi=200,); ax = fig.add_subplot(111)
             
        plt.errorbar(self.volumeFluxes, self.mCentreingOffset, 
                     yerr = self.merrCentreingOffset,
                     linestyle='None', marker = m[0], color = c[0])
        plt.text(0.1, 0.1, 'Where the errorbars indicate the standard error', fontsize = 12, 
             horizontalalignment='left', verticalalignment='center', 
             transform = ax.transAxes)

        plt.title('%s Gradient Initial vs Final Offset' %self.batchName)
        plt.xlabel('Volume Flux [$ml/min$]');   plt.ylabel('Gradient')  


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
        plt.errorbar(self.volumeFluxes, self.cCentreingOffset, 
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
        for vf in self.volumeFluxes:
            fig = plt.figure(figsize=(8, 6), dpi=200,); ax = fig.add_subplot(111)      
            tempx=[]; tempy=[]
            for i in range(len(self.volumeFlux)):
                if self.volumeFlux[i] == vf:
                    tempx.append(self.Offset[i]/self.pixelsPmm[i])
                    tempy.append(self.finalDistCentre[i]/self.pixelsPmm[i])
            plt.errorbar(tempx, tempy,
                        #yerr=self.finalDistCentreSTD/self.pixelsPmm, entries[30]
                        #xerr=self.OffsetSTD/self.pixelsPmm,
                        linestyle='None', marker = m[0], color = c[0])
                    
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
            
            
            #I want to plot the distances from the centreline not as absolute 
            #values. This requieres me to use the boolean 
#            self.
            fig = plt.figure(figsize=(8, 6), dpi=200,); ax = fig.add_subplot(111)      
            tempx=[]; tempy=[]
            for i in range(len(self.volumeFlux)):
                if self.volumeFlux[i] == vf:
                    tempx.append(self.Offset[i]/self.pixelsPmm[i])
                    tempy.append(self.finalDistCentre[i]/self.pixelsPmm[i])
            plt.errorbar(tempx, tempy,
                        #yerr=self.finalDistCentreSTD/self.pixelsPmm, 
                        #xerr=self.OffsetSTD/self.pixelsPmm,
                        linestyle='None', marker = m[0], color = c[0])
                    
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
            
    def _plotAllFinalVsInitialOffset(self):
        '''Plot initial offset from centreline in the channel leading to the 
        semi-sphere versus volume flux Q'''
        m = ana.markerSymbols(); c = ana.coloursForMarker()
        fig = plt.figure(figsize=(8, 6), dpi=200,); ax = fig.add_subplot(111)   
        counter=0
        for vf in self.volumeFluxes:
   
            tempx=[]; tempy=[]
            for i in range(len(self.volumeFlux)):
                if self.volumeFlux[i] == vf:
                    tempx.append(self.Offset[i]/self.pixelsPmm[i])
                    tempy.append(self.finalDistCentre[i]/self.pixelsPmm[i])
            plt.errorbar(np.abs(tempx), tempy,
                        #yerr=self.finalDistCentreSTD/self.pixelsPmm, entries[30]
                        #xerr=self.OffsetSTD/self.pixelsPmm,
                        linestyle='None', marker = m[counter], color = c[counter], 
                        label = "Q = $%.f \ ml/min$" %vf)
            counter +=1
        
        plt.legend(loc='best', fontsize=10, ncol=1)
        plt.title('%s Initial vs Final Offset' %(self.batchName))
        plt.ylabel('Final Offset [$mm$]');   plt.xlabel('Initial Centering [$mm$]')  
        ax = ana._standardPlotSize(ax)        
        
        print("setting fixed y-lims - may need to adjust")
        plt.ylim((7,12))
        ana._saveFig(fig, os.path.join(self.directory, 'RESLT') ,'%s_AllInitiaVsFinalOffset.jpg' %(self.batchName))
            
            
            
    def _plotNormedFinalVsInitialOffset(self):
        '''Plot initial offset from centreline in the channel leading to the 
        semi-sphere versus volume flux Q'''
        m = ana.markerSymbols(); c = ana.coloursForMarker()
        fig = plt.figure(figsize=(8, 6), dpi=200,); ax = fig.add_subplot(111)
        
                
        
        initialCentring=[]; finalDist_normed=[];
        for vft in self.volumeFlux:                
            tempx, tempy =self._returnVolFlux(vft, self.volumeFlux, 
                                         self.Offset/self.pixelsPmm, 
                                         self.finalDistCentre/self.pixelsPmm)
            avefDC, acefDCSTD =self._returnVolFlux(vft, self.volumeFluxes, 
                                         self.aveFinalDistCentre/self.avePPmm, 
                                         self.aveFinalDistCentreSTD/self.avePPmm)
            for kk in range(len(tempx)):
                initialCentring.append(tempx[kk])
                finalDist_normed.append(tempy[kk]/avefDC)
                
#        my, by, ry, errmy, errby   = ana.lsqfity(np.abs(initialCentring),finalDist_normed)
        z = np.polyfit(np.abs(initialCentring),finalDist_normed, 1)
        print('m = %f c= %f' %(z[0], z[1]))
                
        plt.errorbar(initialCentring, 
             finalDist_normed, 
             #yerr=self.finalDistCentreSTD/self.pixelsPmm,
            linestyle='None', marker = m[0], color = c[0],
            markersize=6)
            
        maxx= np.max(np.abs(initialCentring))
#        p = np.poly1d(z)
#        plt.plot([0.0, maxx], [z[0], maxx*z[1]+z[0]], '-', color=c[2])
                    
        plt.title('%s Initial vs Final Offset' %self.batchName)
        plt.ylabel('Final Offset normalized by average offset');   plt.xlabel('Initial Centering [$mm$]')  
        
#        textstring = 'Number of Runs = %d' %(len(self.Offset))
#        plt.text(0.1, 0.1, textstring, fontsize = 12, 
#                     horizontalalignment='left', verticalalignment='center', 
#                    transform = ax.transAxes)
        
        plt.ylim((0.8, 1.4))
                     
        ax = ana._standardPlotSize(ax)        
        ana._saveFig(fig, os.path.join(self.directory, 'RESLT') ,'%s_InitiaVsFinalNormedOffset.jpg' %self.batchName)

            
    def _returnVolFlux(self, tragetVF, vf, x, y):
        '''Returns all x at tragetVF'''
        tempx=[]; tempy=[]
        for i in range(len(vf)):
            if vf[i] == tragetVF:
                tempx.append(x[i])
                tempy.append(y[i])
        return  np.array(tempx), np.array(tempy)
            
    def _plotQVsInitialOffset(self):
        '''Plot final versus initial offset'''
        m = ana.markerSymbols(); c = ana.coloursForMarker()
        fig = plt.figure(figsize=(8, 6), dpi=200,); ax = fig.add_subplot(111)
        
        plt.errorbar(self.volumeFlux, self.Offset/self.pixelsPmm, 
                     #yerr=self.OffsetSTD/self.pixelsPmm,
                    linestyle='None', marker = m[0], color = c[0])
                    
        plt.errorbar(self.volumeFluxes, self.aveOffset/self.avePPmm, 
                     yerr=self.aveOffsetSTD/self.avePPmm, 
                     linestyle='None', marker = m[3], markersize=8, color = c[3])
        plt.title('%s Q vs Centreing' %self.batchName)
        plt.xlabel('Volume Flux [$ml/min$]');   plt.ylabel('Initial Centering [$mm$]')  
        ax = ana._standardPlotSize(ax)        
        
        ana._saveFig(fig, os.path.join(self.directory, 'RESLT') ,'%s_QvsCentering.jpg' %self.batchName)
        
    def _plotQVsAngle(self, fit=False, mode='Linear'):
        '''Plot final versus initial offset'''
        m = ana.markerSymbols(); c = ana.coloursForMarker()
        fig = plt.figure(figsize=(8, 6), dpi=200,); ax = fig.add_subplot(111)
        
        plt.errorbar(self.volumeFlux, self.angle, 
                     #yerr=self.OffsetSTD/self.pixelsPmm,
                    linestyle='None', marker = m[0], color = c[0])
                    
        plt.errorbar(self.volumeFluxes, self.aveAngle, yerr=self.aveAngleSTD , 
                     linestyle='None', marker = m[3], markersize=8, color = c[3])
                     
        if fit:
            x=[]; y=[]
            for ss in range(len(self.finalDistCentre)):
                if self.finalDistCentre[ss] > 0.0:
                    x.append(self.volumeFlux[ss])
                    y.append(self.angle[ss])
            x = np.array(x); y=np.array(y)
            if mode.lower() =='linear':                    
                slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
                plt.plot([0.0, np.max(x)], 
                          [slope*0.0+intercept, slope*np.max(x)+intercept], 
                          '--k', 
                          label='Linear Fit with slope = %f and intercept = %f' %(slope, intercept) )                
                print("In Q vs Angle plot: linear fit results in slope = %f and intercept = %f" %(slope, intercept) )    
                print("std_err: %f" %std_err)
                print("r-squared: %f" %r_value**2)
                
                my, by, ry, smy, sby   = ana.lsqfity(x, y)
                print("2nd way: linear fit results in slope = %f and intercept = %f" %(my, by) )    
                print("ry = %f, smy = %f, sby = %f" %(ry, smy, sby))

            elif mode.lower() == 'sqrt':
                fitParams, fitCovariance = optimize.curve_fit(ana.fitFuncSqrt, x, y)
                error = [] 
                for i in range(len(fitCovariance)):
                    error.append(np.absolute(fitCovariance[i][i])**0.5)
                perr_leastsq = np.array(error)
                #from:
                # http://stackoverflow.com/a/21844726/1826893
                
                #get r^2
                residuals = y - ana.fitFuncSqrt(x, fitParams[0], fitParams[1])
                ss_res = np.sum(residuals**2)
                ss_tot = np.sum((y-np.mean(y))**2)
                r_squared = 1 - (ss_res / ss_tot)

                print("_plotQVsFinalOffset - fitting sqrt: y =  b * np.sqrt(x) + a")
                print("a = %f pm %f \t b = %f pm %f " %(fitParams[0], perr_leastsq[0], fitParams[1], perr_leastsq[1]))
                print("ss_res = %f r^2 = %f" %(ss_res, r_squared))
                xfit = np.arange(0.0, np.max(x))
                yfit = ana.fitFuncSqrt(xfit, fitParams[0], fitParams[1])
                plt.plot(xfit, yfit, '--k')
                plt.plot([0, np.max(x)], [fitParams[0] +fitParams[1],fitParams[0] +fitParams[1]], '-.b')


            elif mode.lower() == 'poly2':
                fitParams, fitCovariance = optimize.curve_fit(ana.fitFuncPoly2, x, y)
                sigma = [fitCovariance[0,0], \
                        fitCovariance[1,1], \
                        fitCovariance[2,2], \
                        ]
                print("_plotQVsAngle - fitting poly2: y = a*x^2 + b*x+c ")
                print("a = %f pm %f \t b = %f pm %f \t c = %f pm %f" %(fitParams[0], sigma[0], fitParams[1], sigma[1], fitParams[2], sigma[2]))
                xfit = np.arange(0.0, np.max(x))
                yfit = ana.fitFuncPoly2(xfit, fitParams[0], fitParams[1] , fitParams[2])
                plt.plot(xfit, yfit, '--k')

            elif mode.lower() == 'exp':
                fitParams, fitCovariance = optimize.curve_fit(ana.fitFuncExp, x, y)
                sigma = [fitCovariance[0,0], \
                        fitCovariance[1,1], \
                        fitCovariance[2,2], \
                        ]

                error = [] 
                for i in range(len(fitCovariance)):
                    error.append(np.absolute(fitCovariance[i][i])**0.5)
                perr_leastsq = np.array(error)
                #from:
                # http://stackoverflow.com/a/21844726/1826893
                
                # find r^2
                residuals = y- ana.fitFuncExp(x, fitParams[0], fitParams[1] , fitParams[2])
                ss_res = np.sum(residuals**2)
                ss_tot = np.sum((y-np.mean(y))**2)
                r_squared = 1 - (ss_res / ss_tot)

                print("_plotQVsAngle - fitting exp: y =  c * np.exp(-b*x) + a")
                print("a = %f pm %f \t b = %f pm %f \t c = %f pm %f" %(fitParams[0], perr_leastsq[0], fitParams[1], perr_leastsq[1], fitParams[2], perr_leastsq[2]))
                print("ss_re = %f r^2 = %f" %(ss_res, r_squared))
                xfit = np.arange(0.0, np.max(x))
                yfit = ana.fitFuncExp(xfit, fitParams[0], fitParams[1] , fitParams[2])
                plt.plot(xfit, yfit, '--k')
                plt.plot([0, np.max(x)], [fitParams[0] +fitParams[2],fitParams[0] +fitParams[2]], '-.b')

            else:
                print("mode specified is not defined")
                    #self.aveAngle, self.aveAngleSTD
        plt.title('%s Q vs Angle' %self.batchName)
        plt.xlabel('Volume Flux [$ml/min$]');   plt.ylabel('Angle [$rad$]')  
        ax = ana._standardPlotSize(ax)        
        
        ana._saveFig(fig, os.path.join(self.directory, 'RESLT') ,'%s_QvsAngle.jpg' %self.batchName)
        
    def _plotQVsMeanHeightEndSS(self):
        ''' consider height of capsule next to obstacle '''        
        m = ana.markerSymbols(); c = ana.coloursForMarker()
        fig = plt.figure(figsize=(8, 6), dpi=200,); ax = fig.add_subplot(111)
        
        plt.errorbar(self.volumeFlux, self.mean_yd_EndSS/self.pixelsPmm, 
                     #yerr=self.OffsetSTD/self.pixelsPmm,
                    linestyle='None', marker = m[0], color = c[0], 
                    label = 'Central Height End Half-Cylinder')
                    
        plt.errorbar(self.volumeFluxes, self.aveMean_yd_EndSS/self.avePPmm, 
                     yerr=self.aveMean_yd_EndSSSTD/self.avePPmm , 
                     linestyle='None', marker = m[1], markersize=8, color = c[1], 
                    label = 'Average Central Height End Half-Cylinder')


                                   
        plt.legend(loc='best', fontsize=6, ncol=1)
        plt.title('%s Q vs Central Width next to Half-Cylinder' %self.batchName)
        plt.xlabel('Volume Flux [$ml/min$]');   plt.ylabel('Height [$mm$]')  
        ax = ana._standardPlotSize(ax)        
        
        ana._saveFig(fig, os.path.join(self.directory, 'RESLT') ,'%s_QvsCentralHeight.jpg' %self.batchName)
        
        print('aveMeanHeightEndSS \taveMinHeightEndSS')
        for gg in range(len(self.aveMeanHeightEndSS)):
            print('%f \t\t\t %f' %(self.aveMeanHeightEndSS[gg], self.aveMinHeightEndSS[gg]))
        
        
    

        
        
        
        
        