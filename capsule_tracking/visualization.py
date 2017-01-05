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
        
def coverImages(ParameterClass, savefolder, color=(255,255,255)):
    """Goes though all the pictures in directory, thresholds them and finds 
    several measurments for the resulting object, which should be a capsule if 
    the threshold has been chossen correctly
    """
    
    if not os.path.exists(savefolder): os.makedirs(savefolder)

    print('Printing Arguments')
    print(' path :  ' + str( ParameterClass.path))
    print('d0 :  ' + str(ParameterClass.d0))
    print('pPmm :  ' + str(ParameterClass.pPmm))
    print('background :  ' + str(ParameterClass.backgroundImage))
    print('rotate :  ' + str(ParameterClass.rotation))
    print('printDebugInfo :  ' + str(ParameterClass.printDebugInfo))

    #for ease of writing
    path = ParameterClass.path; d0=ParameterClass.d0; plot = ParameterClass.plot
    printDebugInfo = ParameterClass.printDebugInfo; pPmm = ParameterClass.pPmm
    
    background= True #bool that states whether we are using background substraction
    
    fileID=gen.find_batchName(path);

    #close all graphs
    plt.close("all")
    
    #'forShow is a boolean variable that governs whether plot is displayed 
    #during program run. Used for debugging only
    forShow=False    
    
    #get sorted list of filenames
    fileList, leng =ParameterClass.listingImageFunc(path) #, fileType='.jpg')
    print('leng = %d' %leng)
    
    #counter of images
    counterL=-1
    
    back = cv2.imread(ParameterClass.backgroundImage,cv2.IMREAD_COLOR)

    #Iterating over all pictures
    for i in range(leng):
        plt.close('all')
        counterL=counterL+1

        if printDebugInfo:
            print('\n\n%d th image (number %d, filename = %s )' %(i, fileList[i].number, fileList[i].fn))
        elif i%10 == 0: #print info every 10th picture
            print("Current frame " + str(i) + "\t and counter is on " +  str(counterL))# + " and memory used is: ")
        #read image
        readPath=path+fileList[i].fn
#        print("readfind_max_extendPath: %s" %readPath)
        img =cv2.imread(readPath,cv2.IMREAD_COLOR)
        t1 = img.copy()
        
        if back.shape != img.shape:
            import warnings
            string = "\n Background image and loaded image do not have the same size. \n" + \
                    "Background image has shape " + str(back.shape) +' and loaded image  ' + str(img.shape) + \
                    "\n Loaded image path: %s \n Skipping Image!" %(readPath)
            warnings.warn(string, RuntimeWarning)
            continue

            
        imgC =  cv2.absdiff(back, img)
        img = (255 - imgC)
        
        if ParameterClass.rotation != 0:
            img=gen.rotateImage(img, ParameterClass.rotation)
            t1=gen.rotateImage(t1, ParameterClass.rotation)
        
        if ParameterClass.flip:
            img=cv2.flip(img,1)
            
        #add a border of black pixels around image in case capsule touches edge 
        bordersize=0
        img=cv2.copyMakeBorder(img,bordersize,bordersize,bordersize,bordersize,
                               cv2.BORDER_CONSTANT,value=color)  

        #Get size of image
        yImg,xImg, _ = img.shape
        
        #Increase Contrast
        img = ParameterClass.coverImgFunc(img, bordersize, 
                                          ParameterClass)
                                          
#        #======================================================================                                          
#        #Option 1: not good but maybe with parameter tweak
#        imghsv = cv2.cvtColor(img, cv2.COLOR_BGR2HSV)
#        imghsv[:,:,2] = [[max(pixel - 25, 0) if pixel < 190 else min(pixel + 25, 255) for pixel in row] for row in imghsv[:,:,2]]
#        img = cv2.cvtColor(imghsv, cv2.COLOR_HSV2BGR)
        
        #======================================================================                                          
        #Option 2
        #-----Converting image to LAB Color model----------------------------------- 
        lab= cv2.cvtColor(img, cv2.COLOR_BGR2LAB)
        cv2.imshow("lab",lab)
        
        #-----Splitting the LAB image to different channels-------------------------
        l, a, b = cv2.split(lab)
        
        #-----Applying CLAHE to L-channel-------------------------------------------
        clahe = cv2.createCLAHE(clipLimit=3.0, tileGridSize=(8,8))
        cl = clahe.apply(l)
        
        #-----Merge the CLAHE enhanced L-channel with the a and b channel-----------
        limg = cv2.merge((cl,a,b))
        
        #-----Converting image from LAB Color model to RGB model--------------------
        img = cv2.cvtColor(limg, cv2.COLOR_LAB2BGR)

        img = ParameterClass.coverImgFunc(img, bordersize, 
                                          ParameterClass, color=color, forPub=True)
        cv2.imwrite(os.path.join(savefolder,'%03d.jpg' %counterL),img)
