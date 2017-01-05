# -*- coding: utf-8 -*-
"""
Collection of functions to analysis result file from tracking capsules in 
top-view images using OpenCV.


Created on 08.02.2016

@author: Edgar Haener
edgar.haner@gmail.com

"""
from __future__ import absolute_import, division, print_function
import matplotlib.pyplot as plt
import numpy as np
import cv2
import shutil
import scipy as sp
#from scipy import interpolate
from scipy      import optimize

import os
#import time

import general as gen
from general import filenameClass

ERROR_CONST=-1
AVERAGINGWINDOW = 6
AVERAGINGREPEATS = 1

FOR_LATEX = False  #sets what save names are



def subtractBackgroundImage(path, background, target, fileType='.png'):
    ''' Takes a folder of images and subtracts an image from all of them
    
    path            path to the folder of images
    background      path to image to be subtracted
    target          path to target folder where images will be saved
    sortfunction    function to read images in folder    
    
    '''
    fileList, leng =sortfunction(path) #, fileType='.jpg')
    print('leng = %d' %leng)
    back = cv2.imread(background,0)
    counterL=-1
    
    if not os.path.exists(target): os.makedirs(target)
    
    print('path = %s' %(path))
    #Iterating over all pictures
    for i in range(leng):
        counterL +=1
        if i%10 == 0: #print info every 10th picture
            print("Current frame " + str(i) + "\t and counter is on " +  str(counterL))# + " and memory used is: ")
        #read image
        readPath=path+fileList[i].fn
    #        print("readfind_max_extendPath: %s" %readPath)
        img =cv2.imread(readPath,0)
        
        
        
        if back.shape != img.shape:
            import warnings
            string = "\n Background image and loaded image do not have the same size. \n" + \
                    "Background image has shape " + str(back.shape) +' and loaded image  ' + str(img.shape) + \
                    "\n Loaded image path: %s \n Skipping Image!" %(readPath)
            warnings.warn(string, RuntimeWarning)
            continue
    
            
        imgC =  cv2.absdiff(back, img)
        img = (255 - imgC)
        
        cv2.imwrite(target+fileList[i].fn, img)
        
def _readImageNames(path, fileType = '.png', prefixleng=10):
    """
    Returns a list of filenameClass, one for each image in directory
    """
    dirs = os.listdir(path)
    numberOfJPG=0

    filenameList=[]
    
    for fname in dirs:
        if (fname[-4:len(fname)] == fileType):
            try:    
                temp=filenameClass(fname, -1, 1 )
                filenameList.append(temp)
                numberOfJPG += 1
            except:
                print('\t Not inluding %s' %fname)    
#    sort by milliseconds
#    filenameList.sort(key=lambda filenameClass: filenameClass.number)
    return filenameList, numberOfJPG
    
    
def plotCentoridOnImage(path1, path2, geo1, geo2, savepath, rotate=-0.3):
    centroid_x1, centroid_y1, _, _, _, _, _, _, _, _, _, _, _, _, _ = readResultsFile(path1, second=False)
    centroid_x2, centroid_y2, _, _, _, _, _, _, _, _, _, _, _, _, _ = readResultsFile(path2, second=False)

    fileList1, lengPic1 =tr.sortPhotosPCO(path1)
    fileList2, lengPic2 =tr.sortPhotosPCO(path2)

    c=2
    for i in range(lengPic1-55):     
        if c%10 ==0:
            print('c = %d' %(c))
        readPath=path1+fileList1[i].fn
        img =cv2.imread(readPath,0)
        
        if rotate != 0:
            img=tr.rotateImage(img, rotate)
        
        img = cv2.cvtColor(img,cv2.COLOR_GRAY2RGB)
        
        for j in range(i):
            if j >0:
                if centroid_x1[j] != -1 and centroid_x1[j-1] != -1:
                    cv2.line(img,(int(round(centroid_x1[j-1])), int(round(centroid_y1[j-1]))),(int(round(centroid_x1[j])), int(round(centroid_y1[j]))), (0,0,255), 3)
        
        img = cv2.resize(img, dsize=(0,0), fx=0.5, fy=0.5) 
        cv2.imwrite(savepath+'Image_'+str(c)+'.jpg', img)
        c=c+1
    
#    gelBeadsXShift=int(round(((geo2[2] - geo1[2])+(geo2[3] - geo1[3]))/2.0))
    gelBeadsXShift=15
    for i in range(25,lengPic2): 
        if c%10 ==0:
            print('c = %d' %(c))
        readPath=path2+fileList2[i].fn
        img =cv2.imread(readPath,0)
        
        img = shiftImage(img, gelBeadsXShift, 0)
        
        if rotate != 0:
            img=tr.rotateImage(img, rotate)
        
        img = cv2.cvtColor(img,cv2.COLOR_GRAY2RGB)
        
        for k in range(lengPic1):
            if k>0:
                xo, xt = int(round(centroid_x1[k])), int(round(centroid_x1[k-1]))
                yo, yt = int(round(centroid_y1[k])), int(round(centroid_y1[k-1]))
                distance = np.sqrt(np.power(xo-xt, 2) + np.power(xo-xt, 2))
                if centroid_x1[k] != -1 and centroid_x1[k-1] != -1 and distance < 25:
                    cv2.line(img,(int(round(centroid_x1[k-1])), int(round(centroid_y1[k-1]))),(int(round(centroid_x1[k])), int(round(centroid_y1[k]))), (0,0,255), 3)
        
        for j in range(i):
            if j>0:
                xo, xt = int(round(centroid_x2[j])), int(round(centroid_x2[j-1]))
                yo, yt = int(round(centroid_y2[j])), int(round(centroid_y2[j-1]))
                distance = np.sqrt(np.power(xo-xt, 2) + np.power(xo-xt, 2))
                if centroid_x2[j] != -1 and centroid_y2[j] != -1 and centroid_x2[j] != 0 and centroid_y2[j] != 0 and distance < 25:
                    cv2.line(img,(int(round(centroid_x2[j-1]))+gelBeadsXShift, int(round(centroid_y2[j-1]))),(int(round(centroid_x2[j]))+gelBeadsXShift, int(round(centroid_y2[j]))), (0,255,0), 3)
        
        img = cv2.resize(img, dsize=(0,0), fx=0.5, fy=0.5) 
        cv2.imwrite(savepath+'Image_'+str(c)+'.jpg', img)
        c=c+1
