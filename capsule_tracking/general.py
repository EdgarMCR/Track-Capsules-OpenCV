"""
Collection of functions to track capsules in top-view images using OpenCV.


Created on 08.02.2016

@author: Edgar Haener
edgar.haner@gmail.com

"""
from __future__ import absolute_import, division, print_function
import matplotlib.pyplot as plt
import numpy as np
import cv2
import scipy as sp
from scipy import interpolate

import os
import time




#------------------------------------------------------------------------------
# CONSTANTS
#------------------------------------------------------------------------------
#maxArea=1.15   #for judging contour quality
#minArea=0.4
#minAreaCH=0.65
#maxCurvatureCutOff = 1.4
#maxPerimeter = 1.7

#when finding capsule contour, set limits for how long expect contour to be
MINLENGTHCONTOUR = 75
MAXLENGTHCONTOUR = 450

UsingOpenCV3 = False
ERROR_CONST=-1

PADDING = 60
BORDERSIZE=5 #border that is added to image 

useDenoising = False

ENHANCE_CONTRAST = False

def setUsingOpenCV3(boolean):
    global UsingOpenCV3
    UsingOpenCV3 = boolean
    
def setUseDenoising(boolean):
    global useDenoising
    useDenoising = boolean
    
def setENHANCE_CONTRAST(boolean):
    global ENHANCE_CONTRAST
    ENHANCE_CONTRAST = boolean
    
    
    

#Parameter Class: A class that contains all the nessecary parameters for the run
#script

class Parameters:
    """Collect all information need from user in one class path, Size, 
    pixels/mm,  rotation
    
    To inizialise provide the following arguments or full path to parameter
    file

    self.baseDirectory      Folder where folders with images are stored
    self.folder             Specific folder in baseDirectory
    self.d0                 Diameter of Capsule
    self.pPmm               Pixels/mm on the images
    self.FPS                Frames per Second at which was recorded
    self.rotation           Rotation of image need to align channel with x-axis
    backgroundImage         Path from base folder to background image that can 
                            be subtracted 
                            
    Runtime arguments
    
    """    
    
    def __init__(self, baseDirectory, folder, d0, pPmm, rotation, FPS, 
                 backgroundImage, plot=False, printDebugInfo=False, 
                 denoising=False, usingOpenCV3=False):
                     
        assert(baseDirectory != None)
        self.baseDirectory = baseDirectory
        self.folder = folder
        self.d0 = d0
        self.pPmm = pPmm
        self.FPS = FPS
        self.rotation = rotation
        self.backgroundImage = backgroundImage
        
        self.updatePath()
        
        self.plot = plot
        self.printDebugInfo=printDebugInfo
        self.denoising=denoising
        
        self.usingOpenCV3=usingOpenCV3
        self.flip = False
        
        #Functions that need to be provided be the specific geometry
        self.listingImageFunc=None #takes path and returns an ordered list and leng 
        self.coverImgFunc=None # Takes image and cuts areas out that are not needed in evaluation
        
    def updatePath(self):
        self.path=os.path.join(os.sep, self.baseDirectory, self.folder, '')
    
    def setFunc(self, a,b):
        '''Set Functions'''
        self.listingImageFunc=a #takes path and returns an ordered list and leng 
        self.coverImgFunc=b # Takes image and cuts areas out that are not needed in evaluation

    
class filenameClass:
    """Stores filename and millisecond time"""
    fn='-1'      #filename
    ms=-1        #millisecondtime   
    number=-1    #number in filename
    timestamp1=-1
    timestamp2=-1
    fps=-1
    mainQ=-1
    sideQ=-1
    
    def __repr__(self):
            return repr((self.fn, self.ms, self.number, self.timestamp1, self.timestamp2, self.fps))
            
    def __init__(self, filename, millisecondTime, number):
        self.fn = filename 
        self.ms = millisecondTime
        self.number=number
        
    def printOut(self):
        print('Filename is \t %s and at \t %.d milliseconds and counter is %d.' %(self.fn, self.ms, self.number))


        
    
class Imagename:
    """Stores filename, counter and millisecond time"""
    fn='-1'      #filename
    ms=-1        #millisecondtime   
    number=-1    #number in filename
    timestamp1=-1
    timestamp2=-1
    fps=-1
    
    def __repr__(self):
            return repr((self.fn, self.ms, self.number, self.timestamp1, self.timestamp2, self.fps))
            
    def __init__(self, filename, millisecondTime, number):
        self.fn = filename
        self.ms = millisecondTime
        self.number=number
        
    def printOut(self):
        print('Filename is \t %s and at \t %.d milliseconds and counter is %d.' %(self.fn, self.ms, self.number))

def _writeLog(runTime, path, background):
    '''Record run times'''
    pathFile = 'capsule_tracking-RunTime.txt' 
    
    b=True
    if background == None:
        b=False
    fileParameters = open(pathFile, 'a')
    fileParameters.write('%f \t%s \t%s \t%r \n' %(runTime, time.strftime("%Y-%m-%d %H:%M:%S") , path, b))        
    fileParameters.close()
    
    
def wholeRun(ParameterClass):
    """ Run the script for finding measurments on all folders.
    
    Makes the assumption that the background image is called
        Background%dFPS.png' %fps
    and is stored in the base directory.
    """
    
    import traceback
    listFolders=os.listdir(ParameterClass.baseDirectory)
    indexDelet=[]
    for i in range(len(listFolders)):
        if os.path.isfile(ParameterClass.baseDirectory+listFolders[i]):
            indexDelet.append(i)

    for i in range(len(indexDelet)-1, -1, -1):
        del listFolders[indexDelet[i]]

    #find FPS info
    fpslist=[]
    for f in listFolders:
        sp1=f.find('FPS')
        if sp1 == -1:
            sp1=f.find('fps')
        sp2=f[:sp1].rfind('_')
        if sp2 == -1 or (len(f[:sp1]) - sp2 ) > 5:
            sp2=f[:sp1].rfind('-')
        
        try:
            FPSString=int(f[sp2+1:sp1])  
    #                print('got fps string: %s' %(FPSString))
        except Exception, err:
            print('Didnt work for \t %s because %s' %(f, err))
            print('f[sp2+1:sp1] = %s \t d1, d2 = %d, %d' %(f[sp2+1:sp1], sp1, sp2))
            print(traceback.format_exc())
            continue
        
        if FPSString not in fpslist:
            fpslist.append(FPSString)
    fpslist.sort()
    for fps in fpslist:
        ParameterClass.FPS = fps
        ParameterClass.backgroundImage = ParameterClass.baseDirectory +'Background%dFPS.png' %fps
        runOneFPS(ParameterClass)
            

def runOneFPS(ParameterClass):
    """ Run the script for finding measurments on all folders with a given
    FPS (and hence volume flux. The convention is to label the folder name with
    both the volum flux 'XXXmlPmin' and frame rate 'XXFPS'. """
    FPS = ParameterClass.FPS
    import traceback
    listFolders=os.listdir(ParameterClass.baseDirectory)
    indexDelet=[]
    for i in range(len(listFolders)):
        if os.path.isfile(ParameterClass.baseDirectory+listFolders[i]):
            indexDelet.append(i)

    for i in range(len(indexDelet)-1, -1, -1):
        del listFolders[indexDelet[i]]
        
    foldersThatWorked=[]
#    print(listFolders)        
    for f in listFolders:
#        print('fname = %s' %(f))
        sp1=f.upper().find('FPS')
        if sp1 == -1:
            sp1=f.find('fps')
        sp2=f[:sp1].rfind('_')
        if sp2 == -1 or (len(f[:sp1]) - sp2 ) > 5:
            sp2=f[:sp1].rfind('-')
#            print('entriesLine[0][sp2+1:sp1] = %s ' %entriesLine[0][sp2+1:sp1])
#            print('reached try statment')
        try:
            FPSString=int(f[sp2+1:sp1])  
    #                print('got fps string: %s' %(FPSString))
        except Exception, err:
            print('Didnt work for \t %s because %s' %(f, err))
            print('f[sp2+1:sp1] = %s \t d1, d2 = %d, %d' %(f[sp2+1:sp1], sp1, sp2))
            print(traceback.format_exc())
            continue
            
        if FPSString != FPS:
                continue
        print('\n Starting \t %s \n' %f)
        ParameterClass.folder = f
        ParameterClass.updatePath()
        track_capsules(ParameterClass=ParameterClass)
        foldersThatWorked.append(f)
#            find_max_extend(directory+f)

    
    print('\nWorked in :')
    for f in foldersThatWorked:
        print(f)
    
def track_capsules(ParameterClass):
    """Goes though all the pictures in directory, thresholds them and finds 
    several measurments for the resulting object, which should be a capsule if 
    the threshold has been chossen correctly
    """
    import timeit
    start_time = start_time = timeit.default_timer()

    print('Printing Arguments')
    print(' path :  ' + str( ParameterClass.path))
    print('d0 :  ' + str(ParameterClass.d0))
    print('pPmm :  ' + str(ParameterClass.pPmm))
    print('background :  ' + str(ParameterClass.backgroundImage))
    print('rotate :  ' + str(ParameterClass.rotation))
    print('printDebugInfo :  ' + str(ParameterClass.printDebugInfo))
    print('denoising   ' + str(ParameterClass.denoising))

    setUseDenoising(ParameterClass.denoising)
    
    
    #for ease of writing
    path = ParameterClass.path; d0=ParameterClass.d0; plot = ParameterClass.plot
    printDebugInfo = ParameterClass.printDebugInfo; pPmm = ParameterClass.pPmm
    
    background= True #bool that states whether we are using background substraction
    
    fileID=find_batchName(path);
    
    #initial 2D projected area
    A0=np.pi*np.power(d0/2,2)


    #close all graphs
    plt.close("all")
    
    #'forShow is a boolean variable that governs whether plot is displayed 
    #during program run. Used for debugging only
    forShow=False    
    
    #get sorted list of filenames
    fileList, leng =ParameterClass.listingImageFunc(path) #, fileType='.jpg')
    print('leng = %d' %leng)
    
    
    #Initializing variable to hold the area, perimeter and vertical distance
    #(y-direction) for each image
    area = np.zeros((leng), dtype=float)
    perimeter = np.zeros((leng), dtype=float)
    verticalDistance = np.zeros((leng), dtype=float)
    horizontalDistance = np.zeros((leng,1), dtype=float)
    width = np.zeros((leng), dtype=float)
    heigth = np.zeros((leng), dtype=float)
    centroid_x = np.zeros((leng), dtype=float)
    centroid_y = np.zeros((leng), dtype=float)
    D = np.zeros((leng), dtype=float)
    Drect = np.zeros((leng), dtype=float) #Taylor deformation parameter found by fitting a rotated rectangle
    angle = np.zeros((leng), dtype=float)
    velx = np.zeros((leng), dtype=float)
    vely = np.zeros((leng), dtype=float)
    xCentralDistance= np.zeros((leng), dtype=float) 
    yCentralDistance= np.zeros((leng), dtype=float) 
    meanC1 = np.zeros((leng), dtype=float) 
    meanC2 = np.zeros((leng), dtype=float) 
    meanC3 = np.zeros((leng), dtype=float) 
    meanC4 = np.zeros((leng), dtype=float)
    
    #counter of images
    counterL=-1
    
    fullMontyFailed=0 #how many images failed
    failedIndexes=[] #index of failed images
    
    back = cv2.imread(ParameterClass.backgroundImage,0)

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
        img =cv2.imread(readPath,0)
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
            img=rotateImage(img, ParameterClass.rotation)
            t1=rotateImage(t1, ParameterClass.rotation)
        
        if ParameterClass.flip:
            img=cv2.flip(img,1)
#        if plot:
#            plt.figure(num='Initial', facecolor='w', edgecolor='k', figsize=(24, 8), dpi=600)
#            a1=plt.subplot(131)
#            plt.imshow(img, cmap='Greys')
            
        #add a border of black pixels around image in case capsule touches edge 
        img=cv2.copyMakeBorder(img,BORDERSIZE,BORDERSIZE,BORDERSIZE,BORDERSIZE,
                               cv2.BORDER_CONSTANT,value=(255, 255, 255))  
                               
#        #crop image to remove black border
#        if plot:
#            savepath1=os.path.join(path,'Check')
#            if not os.path.exists(savepath1): os.makedirs(savepath1)
#            cv2.imwrite(os.path.join(savepath1, 'beforCrop_%d.jpg' %i), img)
#            yImg,xImg = img.shape
#            crop=20
#            img = img[crop:yImg-crop, crop:xImg]   #We first supply the startY and endY coordinates, followed by the startX and endX coordinates to the slice                    
#            cv2.imwrite(os.path.join(savepath1, 'afterCrop_%d.jpg' %i), img)
            
        #Get size of image
        yImg,xImg = img.shape
        
        img = ParameterClass.coverImgFunc(img, BORDERSIZE, ParameterClass)

        contourFound, cntCanny , xCanny,yCanny,wCanny,hCanny, thresh, xBT, yBT, fullMontyFailed, failedIndexes = _findContourBySubtraction(img, fullMontyFailed, failedIndexes, counter = i,threshold=30, plot=plot,  path=path, expectedArea=A0*np.power(pPmm,2), printDebugInfo=printDebugInfo, padding = PADDING, twoCapsules=False)
        

        area[counterL], perimeter[counterL], width[counterL], heigth[counterL], verticalDistance[counterL], \
            horizontalDistance[counterL], angle[counterL], D[counterL], centroid_x[counterL], centroid_y[counterL], ma, MA, box, Drect[counterL], x,y \
            = _assigneValues(pPmm, BORDERSIZE, contourFound, cntCanny , xCanny,yCanny,wCanny,hCanny, thresh, xBT, yBT, background, printDebugInfo)
            
        xCentralDistance[counterL], yCentralDistance[counterL], meanC1[counterL], meanC2[counterL], meanC3[counterL], meanC4[counterL] = _findCentreLengthsandCurvature(cntCanny, centroid_x[counterL], centroid_y[counterL], numPoints=500, plot=plot, path=path, counter=counterL, batchID=fileID)

        #Save outline to file to enable later anlaysis of shape. 
        if np.any(cntCanny) != None:
            savepath1=path+'Outlines'
            if not os.path.exists(savepath1): os.makedirs(savepath1)
            pathFile2=os.path.join(os.sep,  path, 'Outlines', 'OutlineBinary_%s_%d' %(fileID, i))
            np.save(pathFile2, cntCanny) 
            
            if plot:
                t = img.copy()
                t = cv2.cvtColor(t,cv2.COLOR_GRAY2RGB)
                cv2.drawContours(t, [cntCanny], -1, (255,0,0), 1)
                cv2.circle(t, (int(centroid_x[counterL])+BORDERSIZE, int(centroid_y[counterL])+BORDERSIZE), 2, (0,255,0))
                cv2.line(t, (int(centroid_x[counterL])+BORDERSIZE, int(centroid_y[counterL] - yCentralDistance[counterL]/2.0)+BORDERSIZE), (int(centroid_x[counterL])+BORDERSIZE, int(centroid_y[counterL]+yCentralDistance[counterL]/2.0)+BORDERSIZE), (0,0,255))
                cv2.line(t, (int(centroid_x[counterL]- xCentralDistance[counterL]/2.0)+BORDERSIZE, int(centroid_y[counterL] )+BORDERSIZE), (int(centroid_x[counterL]+xCentralDistance[counterL]/2.0)+BORDERSIZE, int(centroid_y[counterL])+BORDERSIZE), (0,255,255))
                savepath=path+'Check'
                if not os.path.exists(savepath): os.makedirs(savepath)
                cv2.imwrite(savepath+'\\AAA_ImageContours_%d.jpg' %i, t)
                
                t1 = cv2.cvtColor(t1,cv2.COLOR_GRAY2RGB)
                
                cntC = cntCanny.copy()
                cntC[:, 0,0] = cntC[:, 0,0] - BORDERSIZE
                cntC[:, 0,1] = cntC[:, 0,1] - BORDERSIZE
                cv2.drawContours(t1, [cntC], -1, (255,0,0), 1)
                cv2.circle(t1, (int(centroid_x[counterL]), int(centroid_y[counterL])), 2, (0,255,0))
                cv2.line(t1, (int(centroid_x[counterL]), int(centroid_y[counterL] - yCentralDistance[counterL]/2.0)), (int(centroid_x[counterL]), int(centroid_y[counterL]+yCentralDistance[counterL]/2.0)), (0,0,255))
                cv2.line(t1, (int(centroid_x[counterL]- xCentralDistance[counterL]/2.0), int(centroid_y[counterL] )), (int(centroid_x[counterL]+xCentralDistance[counterL]/2.0), int(centroid_y[counterL])), (0,255,255))
                savepath=path+'Check'
                if not os.path.exists(savepath): os.makedirs(savepath)
                cv2.imwrite(savepath+'\\BBB_ImageContours_%d.jpg' %i, t1)
                
                
        if plot:    
            savepath=path+'Check' + os.sep
            if not os.path.exists(savepath): os.makedirs(savepath)
                
            if contourFound:                
                #Draw Ellipse
                imC=img.copy()
                filesavepath=savepath+'imgEllipse-'+str(i)+'.jpg'
                cv2.ellipse(imC, (int(y), int(y)), (int(round(ma/2.0)), int(round(MA/2.0))), int(angle[0]), 0, 360, (0,0,0), 2)
                cv2.drawContours(imC,[box],0,(0,0,255),2)
                strW= 'D_{12} ellipse = %f rect = %f' %(D[counterL], Drect[counterL])
                font = cv2.FONT_HERSHEY_SIMPLEX
                cv2.putText(imC,strW,(10,500), font, 0.75,(0,0,0),2)
                cv2.imwrite(filesavepath, imC)
    
            filesavepath=savepath+'Junction_'+str(i)+'.jpg'
#            plt.savefig(filesavepath, dpi=300)
            if(forShow==False):
                plt.clf()
                plt.close("all")
                
                
        del img
        #plt.close("all")
    #finished running over all the images
        
    r = _plotAndWriteToFile(path, d0, pPmm, fileList, area, perimeter, verticalDistance, horizontalDistance, width, heigth, centroid_x, centroid_y, D, angle, velx, vely, Drect, xCentralDistance, yCentralDistance, meanC1, meanC2, meanC3, meanC4, True)
    
    fo = open(path+fileID+'_Failed.txt', "w")
    fo.write('\bFailed to follow proper protocoll %d times. The indexes on which this occured are: \n' %fullMontyFailed)
    for s in failedIndexes :
        fo.write('\t %d' %s)
    fo.close()
    
    print('\b Failed to follow proper protocoll %d times out of %d images.' %(fullMontyFailed,leng))
        
    end_time =  timeit.default_timer()
    _writeLog(end_time - start_time , path, True)
    ParameterClass._writeToFile()
    return r 

def _findContourBySubtraction(img, fullMontyFailed, failedIndexes, counter,threshold=-1, plot=False,  path=None, expectedArea=None, printDebugInfo=None,  padding=30, twoCapsules=False):
    '''
    Find capsule in img by doing an adaptive thresholding and the canny filter 
    the area found
    '''       
    if ENHANCE_CONTRAST:
        img = increaseContrast(img, phi=0.9, theta = 0.95)
        
    if plot:
        savepath=path+'Check\\'
        if not os.path.exists(savepath): os.makedirs(savepath)
        filesavepath=savepath+'imgSubbed-'+str(counter)+'.jpg'
        cv2.imwrite(filesavepath, img)
        
#        cnt, thresh, xBT,yBT,w,h = findContoursByThresholding(img, threshold, plot=plot, counter=i, path=path, expectedArea=A0*np.power(pPmm,2), printDebugInfo=printDebugInfo)
    cnt, thresh, xBT,yBT,w,h, otsu = _findContoursByThresholding(img, adaptive=False,  plot=plot, counter=counter, path=path, expectedArea=expectedArea, printDebugInfo=printDebugInfo, subbed = True, twoCapsules=twoCapsules)
    
    if not twoCapsules:
        if cnt is not None:
            contourFound=True
        else:
            contourFound=False
            fullMontyFailed +=1
            failedIndexes.append(counter)
            if printDebugInfo:
                print('Thresholding failed')
        return contourFound, cnt , xBT,yBT,w,h, thresh, xBT, yBT, fullMontyFailed, failedIndexes
    else:
        contourFound = [np.any(cnt[0]) != None, np.any(cnt[1]) != None]
        if np.any(contourFound) == False:
            if contourFound[0] ==False:
                fullMontyFailed +=1
                failedIndexes.append(counter)
                if printDebugInfo:
                    print('Thresholding failed')
            if contourFound[1] ==False:
                fullMontyFailed +=1
                failedIndexes.append(counter)
                if printDebugInfo:
                    print('Thresholding failed')
        return contourFound, cnt , xBT,yBT,w,h, thresh, xBT, yBT, fullMontyFailed, failedIndexes
        
def _findContoursByThresholding(img,  adaptive, plot=False, counter=0, path=None, expectedArea=None, printDebugInfo=False, subbed=False, twoCapsules=False):
    '''
    Take image, threshold it to a black and white image and find the longest
    contour. Fit a bounding box to this.
    Inputs:
    img         The image to be thresholded
    threshold   The threshold, between 0 and 255 for a fixed threshold, -1 for
                adaptive thresholding
    
    Outputs:
    cnt         the longest contour
    x           corner of bounding box
    y           corner of bounding box
    w           width of bounding box
    h           height of bounding box
    '''
    if subbed:
        color=0
    else:
        color=255
#    print('Thresholding - Type img = ' +str(type(img)))
    
    #convert the image to B&W with the given threshold. 'thresh' is the 
    # the B&W image
    if adaptive:
        blockSize=25
        c=5
        thresh = cv2.adaptiveThreshold(img,255,cv2.ADAPTIVE_THRESH_GAUSSIAN_C,\
                cv2.THRESH_BINARY,blockSize,c)
        threshold = -1
    else:    
        # correct threhsold by removing all pixels equal to 255 and
        # and getting threshold, then apply to image
        tempThresImg = img[img !=  255]

#        thresholdOrg, _ = cv2.threshold(img,0,255,cv2.THRESH_BINARY+cv2.THRESH_OTSU) #using otsu threshold
#        print('Otsu Threshold on whole image = %d' %(thresholdOrg))
        
        thresholdOrg, _ = cv2.threshold(tempThresImg,0,255,cv2.THRESH_BINARY+cv2.THRESH_OTSU) #using otsu threshold
#        print('Otsu Threshold when 255 removed = %d' %(thresholdOrg))
        dif = 255 - thresholdOrg
        mlt=0.0
        mlt=0.75
        threshold = thresholdOrg + int(dif*mlt)
#        threshold = 242
#        print("Threshold is %d" %(threshold))
        ret,thresh = cv2.threshold(img,threshold,255,cv2.THRESH_BINARY)
#        ret,thresh = cv2.threshold(img,245,255,cv2.THRESH_BINARY)
#        threshold2, _ = cv2.threshold(tempThresImg,0,255,cv2.THRESH_BINARY+cv2.THRESH_OTSU)
        if printDebugInfo:
            print('Threshold = %d (befor addition = %d\t' %(threshold, thresholdOrg), end='')
        multiplyer=[0.9, 1.1, 0.8, 1.2]
    
    th=thresh.copy()
    #find the contours in the image
    # Details under: http://docs.opencv.org/modules/imgproc/doc/structural_analysis_and_shape_descriptors.html#findcontours
    #Good tutorial: http://opencvpython.blogspot.co.uk/2012/06/hi-this-article-is-tutorial-which-try.html
    if UsingOpenCV3:
        imgR, contours, hierarchy = cv2.findContours(thresh, cv2.RETR_LIST,cv2.CHAIN_APPROX_SIMPLE) #cv2.RETR_TREE #cv2.cv.CV_RETR_LIST
    else:
        contours, hierarchy = cv2.findContours(thresh, cv2.cv.CV_RETR_LIST,cv2.CHAIN_APPROX_SIMPLE) #cv2.RETR_TREE #cv2.cv.CV_RETR_LIST

    if plot:
        t = th.copy()
        t = cv2.cvtColor(t,cv2.COLOR_GRAY2RGB)
        cv2.drawContours(t, contours, -1, (255,0,0), 3)
        savepath=path+'Check'
        if not os.path.exists(savepath): os.makedirs(savepath)
        cv2.imwrite(savepath+'\\Contours_%d.jpg' %counter, t)
    
    if twoCapsules:
        notFoundOne=True
        notFoundTwo=True
        count2=0
        cnt2 = None; x2,y2,w2,h2 =ERROR_CONST, ERROR_CONST, ERROR_CONST, ERROR_CONST;
    else:
        notFoundOne=True
        notFoundTwo=False
    cnt = None; x,y,w,h =ERROR_CONST, ERROR_CONST, ERROR_CONST, ERROR_CONST;
    count=0
    while notFoundOne  and count < 4:
#        print('Iteration %d of Loop 1 with notFoundOne = %s and notFoundTeo = %s' %(count, notFoundOne, notFoundTwo))
        while notFoundTwo and count2 < 4:
#            print('Before len(contours) = %d \t' %(len(contours)), end = "")
#            print('Iteration %d of Loop 2 with notFoundTwo = %s' %(count2, notFoundTwo))
            if plot:
                t = cv2.cvtColor(th,cv2.COLOR_GRAY2RGB)
                cv2.drawContours(t, contours, -1, (255,0,0), 1)
                countContoursCorrectLength=0
                for ind in range(len(contours)):
                    if len(contours[ind]) > MINLENGTHCONTOUR and len(contours[ind]) < MAXLENGTHCONTOUR:
                        countContoursCorrectLength += 1
                        cv2.drawContours(t, contours, ind, (0,0,255), 2)
                cv2.putText(t,"Number of correct length contours = %d" %countContoursCorrectLength, (50,50), cv2.FONT_HERSHEY_SIMPLEX, 0.5, 255)
                savepath=path+'Check\\Contours'
                if not os.path.exists(savepath): os.makedirs(savepath)
                cv2.imwrite(savepath+'\\Contours_%d_%d.jpg' %(counter, count2), t)
                
            cnt2, x2,y2,w2,h2, para, _ = _checkBWContours(contours, expectedArea, plot, printDebugInfo)                
#            print('After len(contours) = %d \t' %(len(contours)))
            if printDebugInfo:
                savepath=path+'Check'
                if not os.path.exists(savepath): os.makedirs(savepath)
                    
                pathFile=os.path.join(os.sep,  path, 'Check', 'Parameters_%d.txt'%counter)
                if not os.path.isfile(pathFile):
                    fileData = open(pathFile, 'w')
                else:   
                    fileData = open(pathFile, 'a')
                fileData.write("\n\n Capsule 2,Try %d \n" %(count2+1)) 
                for item in para:
                    fileData.write("%s\n" %item)   
                fileData.close()
            
            if x2 == ERROR_CONST:
                if count2 == 0 and useDenoising:
                    img = cv2.fastNlMeansDenoising(img,None,10,7,21)          
#                    cv2.fastNlMeansDenoising(img,None,10,7,21)  
                else:
                    img = cv2.medianBlur(img, ksize=5)
                    
                if adaptive:
                    thresh = cv2.adaptiveThreshold(img,255,cv2.ADAPTIVE_THRESH_GAUSSIAN_C,\
                        cv2.THRESH_BINARY,blockSize,c)
                else:
                    ret, thresh = cv2.threshold(img,multiplyer[count]*threshold,255,cv2.THRESH_BINARY)
                    
                if UsingOpenCV3:
                    imgR, contours, hierarchy = cv2.findContours(thresh, cv2.RETR_LIST,cv2.CHAIN_APPROX_NONE) #cv2.RETR_TREE #cv2.cv.CV_RETR_LIST
                else:
                    contours, hierarchy = cv2.findContours(thresh, 
                           cv2.cv.CV_RETR_LIST,cv2.CHAIN_APPROX_NONE)
#                print('mehhh')
                count2 += 1
            else:
                notFoundTwo=False
                count2 += 1
                break
                
#        #TODO: find better method of seperating the two contours
        if twoCapsules:
            posExcludedContour = [x2, y2, w2, h2]
            cnt, x,y,w,h, para, _ = _checkBWContours(contours, expectedArea, plot, printDebugInfo, posExcludedContour=posExcludedContour)
        else:
            cnt, x,y,w,h, para, _ = _checkBWContours(contours, expectedArea, plot, printDebugInfo)

                
        if printDebugInfo:
                savepath=path+'Check'
                if not os.path.exists(savepath): os.makedirs(savepath)
                    
                pathFile=os.path.join(os.sep,  path, 'Check', 'Parameters_%d.txt'%counter)
                if not os.path.isfile(pathFile):
                    fileData = open(pathFile, 'w')
                else:   
                    fileData = open(pathFile, 'a')
                fileData.write("\n\n Capsule 1,Try %d \n" %(count+1)) 
                for item in para:
                    fileData.write("%s\n" %item)   
                fileData.close()

        if x == ERROR_CONST:
            if count == 0 and useDenoising:
                img = cv2.fastNlMeansDenoising(img,None,10,7,21)                
            else:
                img = cv2.medianBlur(img, ksize=5)
                
            if adaptive:
                thresh = cv2.adaptiveThreshold(img,255,cv2.ADAPTIVE_THRESH_GAUSSIAN_C,\
                    cv2.THRESH_BINARY,blockSize,c)
            else:
                ret, thresh = cv2.threshold(img,multiplyer[count]*threshold,255,cv2.THRESH_BINARY)
                
            if UsingOpenCV3:
                imgR, contours, hierarchy = cv2.findContours(thresh, cv2.RETR_LIST,cv2.CHAIN_APPROX_NONE) #cv2.RETR_TREE #cv2.cv.CV_RETR_LIST
            else:
                contours, hierarchy = cv2.findContours(thresh, 
                       cv2.cv.CV_RETR_LIST,cv2.CHAIN_APPROX_NONE)
#            print('mehhh')
            count += 1
        else:
            notFoundOne=False
            count += 1
            break
            

    if plot:
        imgForPlot=img.copy()
        cv2.rectangle(imgForPlot,(x,y),(x+w,y+h),(color,color,color),2)  
        if twoCapsules:
            cv2.rectangle(imgForPlot,(x2,y2),(x2+w2,y2+h2),(color,color,color),2) 
        
        savepath=path+'Check'+os.sep 
        if not os.path.exists(savepath): os.makedirs(savepath)
        filesavepath=savepath+'ImageWContourFromAdpativeThresholding_'+str(counter)+'.png'
        cv2.imwrite(filesavepath, imgForPlot)
        th = cv2.cvtColor(th,cv2.COLOR_GRAY2RGB)
        cv2.drawContours(th, contours, -1, (255,0,0), 1)
        cv2.rectangle(th,(x,y),(x+w,y+h),(0, 0,255),2)  
        if twoCapsules:
            print('excuting 2nd rectangle')
            cv2.rectangle(th,(x2,y2),(x2+w2,y2+h2),(0, 0,255),2)
        filesavepath=savepath+'Thres_'+str(counter)+'.png'
        cv2.imwrite(filesavepath, th)
        del imgForPlot
    
    if printDebugInfo:
        print('x,y,w,h = %d, %d, %d, %d and number of retries needed = %d in findContoursByThresholding(...)' %(x,y,w,h, count))
        
    if twoCapsules:
        return [cnt2, cnt], th, [x2, x], [y2, y] , [w2, w], [h2, h], threshold
    else:
        return cnt, th, x,y,w,h, threshold
        
def _checkBWContours(contours, expectedArea, plot, printDebugInfo, posExcludedContour=None):
    notFoundMatchingContour=True
    counterWhile=-1
    para=[]
    while notFoundMatchingContour and counterWhile < len(contours):  
        counterWhile += 1
#        print('loop iteration %d and notFoundMatchingContour =  %s' %(counterWhile, notFoundMatchingContour))
        
         #Select longest contour as this should be the capsule
        lengthC=0
        ID=-1
        idCounter=-1
        for xc in contours:
            idCounter=idCounter+1 
            if len(xc) > lengthC:
                lengthC=len(xc)
                ID=idCounter
        
        x,y,w,h = -1, -1, -1, -1
        #if longest contour was found, then ID is the index of it
        if ID != -1:
            if len(contours[ID]) > MINLENGTHCONTOUR and len(contours[ID]) < MAXLENGTHCONTOUR:
                cnt = contours[ID] 
                if expectedArea != None:
                    M = cv2.moments(cnt)
                    if M['m00'] > 0.2*expectedArea and M['m00'] < 1.8*expectedArea:
                        if printDebugInfo:
                            print('Adaptive Threshold: Area = %.f with length=%d, expected Are = %.f' %(M['m00'], len(contours[ID]), expectedArea))
                        
                        #find bounding rectangle of countour
                        x,y,w,h = cv2.boundingRect(cnt)
                        
                        #check that contour has reasonable aspect ratio, max is 6:1
                        if w >= h:
                            aspectRatio= w / h 
                        else:
                            aspectRatio= h/ w
                        
                        if aspectRatio < 7:
                            if np.any(posExcludedContour) == None:
                                notFoundMatchingContour=False
                            else:
                                if abs(x-posExcludedContour[0]) < (w+posExcludedContour[2])/2 and abs(y-posExcludedContour[1]) < (h+posExcludedContour[3])/2:
                                    notFoundMatchingContour=True
                                    para.append('Too close to other contour')
                                    del contours[ID]
                                    cnt=None
                                else:
                                    notFoundMatchingContour=False
#                            print('Match found! (using expected area), notFoundMatchingContour= %s' %notFoundMatchingContour)
                        else:
                            para.append('Aspect Ratio > 7 : %f' %aspectRatio)
                            if printDebugInfo:
                                print(para[-1])
                            del contours[ID]
                            cnt=None
                            
                    else:
                        para.append('Area wrong : expected = %f actual = %f' %(expectedArea, M['m00']))
                        if printDebugInfo:
                                print(para[-1])
                        del contours[ID]
                        cnt=None
                        
                else:
                    #find bounding rectangle of countour
                    x,y,w,h = cv2.boundingRect(cnt)
                    
                    #check that contour has reasonable aspect ratio, max is 6:1
                    if w >= h:
                        aspectRatio= w / h 
                    else:
                        aspectRatio= h/ w
                    
                    if aspectRatio < 7:
                        if np.any(posExcludedContour) == None:
                                notFoundMatchingContour=False
                        else:
                            if abs(x-posExcludedContour[0]) < (w+posExcludedContour[2])/2 and abs(y-posExcludedContour[1]) < (h+posExcludedContour[3])/2:
                                notFoundMatchingContour=True
                                para.append('Too close to other contour')
                                del contours[ID]
                                cnt=None
                            else:
                                notFoundMatchingContour=False
#                        print('Match found! (NOT using expected area)')

                    else:
                        para.append('Aspect Ratio > 7 : %f' %aspectRatio)
                        if printDebugInfo:
                            print(para[-1])
                        del contours[ID]
                        cnt=None
                        
                    
            else:
                para.append('Contour to long/short %d, should be %d < %d ' %(len(contours[ID]), MINLENGTHCONTOUR, MAXLENGTHCONTOUR))
                if printDebugInfo:
                    print(para[-1])
                del contours[ID]
                cnt=None
        else: # ID == -1 i.e. no contours
            para.append('No contour')
            if printDebugInfo:
                print(para[-1])
            cnt=None
            
    para.append('No error to report in checkBWContours')
#    print('para[-1] = %s ' %para[-1])
    return cnt, x,y,w,h, para, ID
    
    
def _assigneValues(pPmm, borderSize, contourFound, cntCanny , xCanny,yCanny,wCanny,hCanny, thresh, xBT, yBT, background, printDebugInfo):
    '''
    Take contour and extract relevant parameters
    '''
    
    if contourFound:  
            #find are of contour and perimeter
            area = cv2.contourArea(cntCanny)
            perimeter = cv2.arcLength(cntCanny,True)
            
            width=wCanny
            heigth=hCanny
            verticalDistance = (hCanny+0.0)/(pPmm+0.0)
            horizontalDistance = (wCanny+0.0)/(pPmm+0.0)      
            
            (x,y),(ma, MA),angle = cv2.fitEllipse(cntCanny)
        
            #Taylor parameter
            D=(MA-ma)/(MA+ma)
            
            rect = cv2.minAreaRect(cntCanny)
            (xr, yr), (w, h), anlgeRect = rect
            if UsingOpenCV3:
                box = cv2.boxPoints(rect)
            else:
                box = cv2.cv.BoxPoints(rect)
            box = np.int0(box)
            D2 = (abs(w-h+0.0))/(w+h+0.0)

            print('(x,y,w,h) = (%d, %d,%d, %d)\tD12 = %.2f, rect D12 = %.2f , (ma, MA) = (%.2f, %.2f)' %(xCanny,yCanny,wCanny,hCanny, D, D2, ma, MA))
        
            #find centroid
            M = cv2.moments(cntCanny)
            try:
                if background:
                    centroid_x = float(M['m10']/M['m00']) - borderSize
                    centroid_y = float(M['m01']/M['m00']) - borderSize
                else:    
                    centroid_x = xBT - borderSize - PADDING + float(M['m10']/M['m00']) #xCanny + int(M['m10']/M['m00'])
                    centroid_y = yBT - borderSize - PADDING + float(M['m01']/M['m00']) #yCanny + int(M['m01']/M['m00'])
                if printDebugInfo:
                    print('x = %d , borderSize = %d, padding = %d, centroid_x = %d ' %(xBT, borderSize, PADDING, centroid_x))
            except: 
                centroid_x = ERROR_CONST
                centroid_y = ERROR_CONST
                    
    else:
        area = ERROR_CONST
        perimeter = ERROR_CONST
        
        width=ERROR_CONST
        heigth=ERROR_CONST
        verticalDistance = ERROR_CONST
        horizontalDistance = ERROR_CONST  
        
        D=ERROR_CONST
        centroid_x = ERROR_CONST
        centroid_y = ERROR_CONST
        angle = ERROR_CONST
        ma = ERROR_CONST
        MA = ERROR_CONST
        box = ERROR_CONST
        D2  = ERROR_CONST
        x, y = ERROR_CONST, ERROR_CONST
        print("No contour!")
        
    return area, perimeter, width, heigth, verticalDistance, horizontalDistance, angle, D, centroid_x, centroid_y, ma, MA, box, D2, x- borderSize,y- borderSize, 

def _findCentreLengthsandCurvature(cnt, centroid_x, centroid_y, numPoints=500, plot=False, path='', counter=0, batchID=''):
    '''Find the measurments of capsule length along centre line
    and the cruvature of each quadrant. 
    
    The approch is to get a parametic interpolation of x and y and use this
    to have a large number of points on the contour. 

    With this, the curvature and lengths are evaluated.     
    
    -Input-
    cnt         contour of capsule to be evaluated
    numPoints   number of points on the interpolated membrane
    
    -Output-
    LfrX        distance alonge centreline, x-direction   
    LfrY        distance alonge centreline, y-direction  
    C1          Curvature first qudrant
    C2          Curvature second qudrant
    C3          Curvature third quadrant
    C4          Curvature fourth qudrant
    '''
    tol = 5.0 # tolerance: within 2*tol of central line(in pixels)
    offsetFromC= 20.0
    xDistance=ERROR_CONST; yDistance = ERROR_CONST
    meanC1,  meanC2,  meanC3,  meanC4 = ERROR_CONST, ERROR_CONST, ERROR_CONST, ERROR_CONST

    if cnt is not None:        
        x, y = np.hsplit(cnt.squeeze(), 2)
        x = x.squeeze() ; y=y.squeeze();
        
        cD, xS, yS = _splineInterpolationOfOutline(x,y,numPoints) #Interpolate data

        B = np.where(np.abs(xS - centroid_x - BORDERSIZE) < tol) #  and xS < centroid_x[ii] + tol)
        if len(B[0]) > 0:
            ymax = np.max(yS[B]) +0.0
            ymin = np.min(yS[B]) +0.0
            yDistance=ymax-ymin
            
        B = np.where(np.abs(yS - centroid_y - BORDERSIZE) < tol)#  and yS < centroid_y[ii] + tol)
        if len(B[0]) > 0:
            xmax = np.max(xS[B]) +0.0 
            xmin = np.min(xS[B]) +0.0
            xDistance=xmax-xmin
            

#        if plot:
#            f=plt.figure(figsize=(6, 6), dpi=200); ax = f.add_subplot(111)
#            plt.plot(x, y, 'ro', label='$L_{fr}$ Daugther')
#            plt.plot(xS, yS, 'kd--', label='$L_{fr}$ Daugther')
#            yi, yx = plt.ylim(); xi, xx = plt.xlim();
#            plt.plot( [xmin, xmax], [centroid_y+BORDERSIZE, centroid_y+BORDERSIZE], '-y')
#            plt.plot([centroid_x+BORDERSIZE, centroid_x+BORDERSIZE],[ymin, ymax], '-g')
#            s = 'Distance y = %f \nDistance x = %f' %(ymax-ymin, xmax-xmin)
#            ax.text(0.1, 0.7, s, fontsize=10, horizontalalignment='left', verticalalignment='center', transform = ax.transAxes)
#            plt.title('Central Distance %s' %batchID)
#            plt.xlabel('X-Position [pixels]'); plt.ylabel('Y-Position [pixels]')
#            savename = path + 'Check' + os.sep + batchID+'_CentralDistance_%d.jpg' %counter
#            plt.savefig(savename, dpi=300,bbox_inches='tight')
#            plt.close(f)
            
            
        #Measure curvature in the 4 quadrants
        #1. Split into quadrants
        xS=xS-centroid_x; yS=yS-centroid_y
        ind1 = ((xS > 0.0)& (yS > 0.0)).nonzero()[0]
        x1=xS[ind1];y1=yS[ind1]
        
        ind2 = ((xS < 0.0)& (yS > 0.0)).nonzero()[0]
        x2=xS[ind2];y2=yS[ind2]

        ind3 = ((xS < 0.0)& (yS < 0.0)).nonzero()[0]
        x3=xS[ind3];y3=yS[ind3]

        ind4 = ((xS > 0.0)& (yS < 0.0)).nonzero()[0]
        x4=xS[ind4];y4=yS[ind4]
        
        meanC1=np.mean(_findMaxCurvatureInContour(x1,y1, n=12))
        meanC2=np.mean(_findMaxCurvatureInContour(x2,y2, n=12))
        meanC3=np.mean(_findMaxCurvatureInContour(x3,y3, n=12))
        meanC4=np.mean(_findMaxCurvatureInContour(x4,y4, n=12))
                
    return xDistance, yDistance,  meanC1,  meanC2,  meanC3,  meanC4
    
def _splineInterpolationOfOutline(x,y,num):
    '''Returns interpolated function
    
    Takes the length along the outline as a parameter and finds spline interpolation
    
    -Input-
    x       x-position
    y       y-position
    num     number of points in spline interpolation
    
    
    -Output-
    cD      cumulative distance where x,y are being evaluated
    xS      interpolated values for x at num points
    yS      interpolated values for y at num points
    '''
    
    #Smooth outline 
    x_smoothed = smooth(x,window_len=11,window='hanning');  y_smoothed = smooth(y,window_len=11,window='hanning');
    x_smoothed = smooth(x_smoothed,window_len=11,window='hanning');  y_smoothed = smooth(y_smoothed,window_len=11,window='hanning');
    x_smoothed = smooth(x_smoothed,window_len=11,window='hanning');  y_smoothed = smooth(y_smoothed,window_len=11,window='hanning');

    #parametize by the lenght along the contour
    N=len(x); n=1;
    distance=np.zeros((N), dtype=float)
    cumulativeDistance = np.zeros((N), dtype=float)
    count=-1
    for jj in range(N):
        count +=1
        
        x2, y2 = x[(jj)%N], y[(jj)%N]
        x1, y1 = x[(jj - n) % N], y[(jj - n) % N] 
        distance[jj] = np.sqrt((x2 - x1)**2 + (y2- y1)**2)
    
    #move the starting point along as not to caputre hole in contour that often appears
    #at the start and end of contour
    offset=20
    xNew = np.zeros((N), dtype=float)
    yNew = np.zeros((N), dtype=float)
    for kk in range(N):
        xNew[kk] = x_smoothed[(kk+offset) %N]
        yNew[kk] = y_smoothed[(kk+offset) %N]
        
        if kk !=0:
            cumulativeDistance[kk] = cumulativeDistance[kk-1] + distance[(kk+offset) %N]
    #interpolate
    xspline = interpolate.splrep(cumulativeDistance, xNew) #scipy sub-module
    yspline = interpolate.splrep(cumulativeDistance, yNew) #scipy sub-module
    
    cD = np.linspace(0.0, cumulativeDistance[-1], num = num)
    xS = sp.interpolate.splev(cD, xspline)
    yS = sp.interpolate.splev(cD, yspline)
    
    return cD, xS, yS
    
def _findMaxCurvatureInContour(x,y, size=None, mmtopixels=None, plot=False, path=None, title='', n=6):
    '''
    This function takes a contour and returns the local curvature between all 
    points, n points appart. 
    
    Consider the angle between the ith and the ith+n point on contour
    
    Taken from http://stackoverflow.com/questions/22029548/is-it-possible-in-opencv-to-plot-local-curvature-as-a-heat-map-representing-an-o
    '''
    assert len(x) == len(y)
    
    N = len(x)
    if N <2: 
        import warnings
        warnings.warn("Length of vector to short. (N=%d, n=%d). Returning Error Code"%(N, n))
        return -1, -1
        
    if (N -n) < 1:
        import warnings
        warnings.warn("Length of vector to short for given window. (N=%d, n=%d). Returning Error Code." %(N, n))
        return -1, -1
    dx =np.gradient(x); dy =np.gradient(y)

    if plot:
        heatmap = np.zeros(size, dtype=np.float)
    
    
#    measure = np.zeros((N), dtype=float)
#    for i in xrange(N):
#        gx1, gy1 = dx[i], dy[i]
#        gx2, gy2 = dx[(i + n) % N], dy[(i + n) % N] 
#    print('N-n = %d' %(N-n))
    measure = np.zeros((N-n), dtype=float)
    for i in xrange(N-n):
        gx1, gy1 = dx[i], dy[i]
        gx2, gy2 = dx[(i + n)], dy[(i + n)] 
    
        # Angle between (gx1, gy1) and (gx2, gy2), normalized by the distance between the points
        cos_angle = gx1 * gx2 + gy1 * gy2
        cos_angle /= (np.linalg.norm((gx1, gy1)) * np.linalg.norm((gx2, gy2)))
        angle = np.arccos(cos_angle) #/ (distance)
        
        if cos_angle < 0.0:
            angle = np.pi - angle

        x1, y1 = x[((2*i + n) // 2) % N], y[((2*i + n) // 2) % N]  # Get the middle point between i and (i + n)
        
        if plot:
            px1 = int(round((mmtopixels*x1)))
            py1 = int(round((mmtopixels*y1)))
            heatmap[py1, px1] = angle  # Use angle between points as score
            
        measure[i] = angle
    
#    if plot:
##        heatmap = cv2.GaussianBlur(heatmap, (3, 3), heatmap.max())
#        f2=plt.figure(figsize=(6, 6), dpi=200); ax = f2.add_subplot(111)
#        plt.imshow(heatmap, cmap=plt.cm.jet)
#        plt.colorbar()
#        plt.title('Curvature '+title)
#        
#        plt.show()
#        savepath=path + 'Check' + os.sep
#        if not os.path.exists(savepath): os.makedirs(savepath)
#        filesavepath=savepath+'HeatmapCurvature_'+title+'.jpg'
#        plt.savefig(filesavepath, dpi=300,bbox_inches='tight')
#        plt.close(f2)
        
        
    return measure

def _plotAndWriteToFile(path, d0, pPmm, fileList, area, perimeter, verticalDistance, horizontalDistance, width, heigth, centroid_x, centroid_y, D, angle, velx, vely,  Drect, xCentralDistance, yCentralDistance, meanC1, meanC2, meanC3, meanC4, constantFrequency=True,  string=''):
    leng = len(area)
    fileID=find_batchName(path);
    
    disp_x  = np.zeros((leng-1), dtype=float)
    disp_y = np.zeros((leng-1), dtype=float)
    vel_x = np.zeros((leng-1), dtype=float)
    vel_y= np.zeros((leng-1), dtype=float)
    
    speed= np.zeros((leng-1), dtype=float)
    speed_ave= np.zeros((leng-1), dtype=float)
    time= np.zeros((leng), dtype=float)
    time[0]=0.0
    
    for i in range(1,leng-1):
        disp_x[i]= centroid_x[i] - centroid_x[i-1] 
        disp_y[i]= centroid_y[i] - centroid_y[i-1] 
        
        dispCutOff=100.0
        if abs(disp_x[i]) > dispCutOff:
            print('disp_x[i] > %f: disp_x[%d]=%e' %(dispCutOff, i,disp_x[i]))
            disp_x[i]=0
        
        if abs(disp_y[i]) > dispCutOff:
            print('disp_y[i] > %f: disp_y[%d]=%e' %(dispCutOff, i,disp_y[i]))
            disp_y[i]=0
        
        if constantFrequency:
            fps=fileList[i].fps+0.0
            deltaT=1.0/fps
            time[i]=time[i-1]+deltaT
            vel_x[i] = (disp_x[i]/pPmm)/deltaT
            vel_y[i] = (disp_y[i]/pPmm)/deltaT            
        else:
            deltaT=((fileList[i].timestamp1-fileList[i-1].timestamp1)/1000.0)
            time[i]=time[i-1]+deltaT
            vel_x[i] = (disp_x[i]/pPmm)/deltaT
            vel_y[i] = (disp_y[i]/pPmm)/deltaT
        
        velocityCutOff=100.0
        if abs(vel_x[i]) > velocityCutOff:
            print('vel_x[i] > %f: vel_x[%d]=%e' %(velocityCutOff, i,vel_x[i]))
            vel_x[i]=0
        
        if abs(vel_y[i]) > velocityCutOff:
            print('vel_y[i] > %f: vel_y[%d]=%e' %(velocityCutOff, i,vel_y[i]))
            vel_y[i]=0
    
    #Running average
    aveWindow=4
    speed = np.sqrt( np.power( vel_x ,2) + np.power(vel_y ,2) )
    speed_ave = np.sqrt( np.power( runningMean(vel_x, aveWindow) ,2) + np.power(runningMean(vel_y, aveWindow) ,2) )
    speed_ave =  runningMean(speed_ave, aveWindow)  
    speed_ave =  runningMean(speed_ave, aveWindow)   
    
    #identify different regions
    dx = gradient(speed_ave, h=1/60)
    zx=np.arange(len(dx))
    
    plt.figure(figsize=(6, 6), dpi=200)
    plt.plot(zx, dx, 'oc',label='Gradient Speed', markersize=5)
    plt.title("dx"+" " + fileID)
    plt.ylabel("dx")
    plt.legend(loc='best')
    plt.savefig(path+fileID+"_GradientSpeed_Graph"+string+".jpg", dpi=300)
    
    
    #Write results to file    
    fo = open(path+fileID+'_Results'+string+'.txt', "w")
    fo.write("# \t Ver Dist \t Hor Dist \t Perimeter \t Area \t centroid_x \t centroid_y \t width[i] \t heigth[i] \t D \t angle \t time (s) \tvel_x (mm/s)\t vel_y (mm/s)\t D_12 from Rectnagle \tCentral x distance \t central y distance \t curvature 1st qudrant \t curvature 2nd qudrant \t curvature 3rd qudrant \t curvature 4th qudrant \t D_Rect\n\n\n")
    for i in range(leng) :
        if i==0:
            fo.write("%d \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t%.4f \t%.4f \t%.4f \t%.4f \t%.4f \t%.4f \n" %(i, verticalDistance[i], horizontalDistance[i], perimeter[i], area[i], centroid_x[i], centroid_y[i], width[i], heigth[i], D[i], angle[i],time[i], 0.0, 0.0, Drect[i], xCentralDistance[i], yCentralDistance[i], meanC1[i], meanC2[i], meanC3[i], meanC4[i]))
        else:
            fo.write("%d \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t%.4f \t%.4f \t%.4f \t%.4f \t%.4f \t%.4f \n" %(i, verticalDistance[i], horizontalDistance[i], perimeter[i], area[i], centroid_x[i], centroid_y[i], width[i], heigth[i], D[i], angle[i],time[i], vel_x[i-1], vel_y[i-1], Drect[i], xCentralDistance[i], yCentralDistance[i], meanC1[i], meanC2[i], meanC3[i], meanC4[i]))
    fo.close()

    #plot all variables
                
    #Remove crazy entries         
    x=np.arange(leng)
    
    #normalize by initial length
    verticalDistance=verticalDistance/(d0)
    horizontalDistance=horizontalDistance/(d0)

    #plost results
    plt.figure(figsize=(6, 6), dpi=200)
    plt.plot(x, verticalDistance, 'oc',label='Vertical', markersize=5)
    plt.plot(x, horizontalDistance, 'sb',label='Horizontal', markersize=5)
    plt.title("Extend of Capsule"+" " + fileID)
    plt.xlabel("Picture # ")
    plt.ylabel("Extend [d0]")
    plt.legend(loc='best')
    plt.savefig(path+fileID+"_Extend_Graph"+string+".jpg", dpi=300)
    
    A0=np.pi*np.power(d0/2,2)
    p0=2.0*np.pi*(d0/2)
    
    area=area/(np.power(pPmm,2))
    area=area/A0
    
    plt.figure(figsize=(6, 6), dpi=200)
    plt.plot(x, area, 'or',label='Area', markersize=5)
    plt.title("Area"+" " + fileID)
    plt.xlabel("Picture # ")
    plt.ylabel("Area / Initial Area")
    plt.savefig(path+fileID+"_Area_Graph"+string+".jpg", dpi=300)
    
    perimeter=perimeter/pPmm
    perimeter=perimeter/p0
    
    plt.figure(figsize=(6, 6), dpi=200)
    plt.plot(x, perimeter, 'sb',label='Perimeter', markersize=5)
    plt.title("Perimeter"+string+" "+" " + fileID)
    plt.xlabel("Picture # ")
    plt.ylabel("Perimeter/ Initial Perimeter")
    plt.savefig(path+fileID+"_Perimeter_Graph"+string+".jpg", dpi=300)
    
    print("Max Values for vertical extend, perimeter and area \t"+ str(np.max(verticalDistance))+"\t"+ str(np.max(perimeter))+"\t"+ str(np.max(area)))
    
    r=np.zeros(3)
    r[0]=np.max(verticalDistance)
    r[1]=np.max(perimeter)
    r[2]=np.max(area)    
    
    #plotVelocity graphs
    # evaluate speed
    
    plt.figure(figsize=(6, 6), dpi=200)
    plt.plot(centroid_x,centroid_y, 'sb',label='Centroid Position', markersize=5)
    plt.title("Centroid Position"+string+" "+" " + fileID)
    plt.xlabel("Z Position [pixels]")
    plt.ylabel("Y Position [pixels]")
    plt.legend(loc='best')
    plt.savefig(path+fileID+"_CentroidPosition_Graph"+string+".jpg", dpi=300)        

    
    x1=np.arange(leng-1)
    plt.figure(figsize=(6, 6), dpi=200)
    plt.plot(x1, vel_x, 'sb',label='v_x', markersize=5)
    plt.plot(x1, vel_y, 'or',label='v_y', markersize=5)
    plt.title("Velocity"+string+" "+" " + fileID)
    plt.xlabel("Picture #  ")
    plt.ylabel("Velocity [mm/s]")
    plt.legend()
    plt.savefig(path+fileID+"_Velocity_Graph"+string+".jpg", dpi=300)
    
    plt.figure(figsize=(6, 6), dpi=200)
    plt.plot(x1, disp_x, 'sb',label='disp_x', markersize=5)
    plt.plot(x1, disp_y, 'or',label='disp_y', markersize=5)
    plt.title("Displacment"+string+" "+" " + fileID)
    plt.xlabel("Picture #  ")
    plt.ylabel("Displacment [pixels]")
    plt.legend()
    plt.savefig(path+fileID+"_Displacment_Graph"+string+".jpg", dpi=300)
    
    plt.figure(figsize=(6, 6), dpi=200)
    plt.plot(x, time, 'sb',label='Time (s)', markersize=5)
    plt.title("Time [s] "+string+" "+" " + fileID)
    plt.xlabel("Picture #  ")
    plt.ylabel("Time [s] ")
    plt.legend()
    plt.savefig(path+fileID+"_Time_Graph"+string+".jpg", dpi=300)
    
    plt.figure(figsize=(6, 6), dpi=200)
    plt.plot(x1, speed_ave, 'sb',label='Speed Averaged', markersize=5)
    plt.title("SpeedAveraged"+string+" "+" " + fileID)
    plt.xlabel("Picture #  ")
    plt.ylabel("Speed[mm/s]")
    plt.legend()
    plt.savefig(path+fileID+"_SpeedAve_Graph"+string+".jpg", dpi=300)
    
    plt.figure(figsize=(6, 6), dpi=200)
    plt.plot(x1, runningMean(disp_x, 6), 'sb',label='running mean disp_x', markersize=5)
    plt.plot(x1, runningMean(disp_y, 6), 'or',label='running mean disp_y', markersize=5)
    plt.title("Displacment Runnung Average"+string+" "+" " + fileID)
    plt.xlabel("Picture #  ")
    plt.ylabel("Displacment [pixels]")
    plt.legend()
    plt.savefig(path+fileID+"_DisplacmentRunnungAverage_Graph"+string+".jpg", dpi=300)
    
    plt.figure(figsize=(6, 6), dpi=200)
    plt.plot(x1, speed, 'or',label='Speed', markersize=5)
    plt.title("Speed"+string+" "+" " + fileID)
    plt.xlabel("Picture #  ")
    plt.ylabel("Speed[mm/s]")
    plt.legend()
    plt.savefig(path+fileID+"_Speed_Graph"+string+".jpg")
    
    plt.figure(figsize=(6, 6), dpi=200)
    plt.plot(x,centroid_x, 'sb',label='Centroid x ', markersize=5)
    plt.plot(x,centroid_y, 'or',label='Centroid y ', markersize=5)
    plt.title("Centroid Position"+string+" "+" " + fileID)
    plt.xlabel("Picture #  ")
    plt.ylabel("Centroid Position [pixels]")
    plt.legend()
    plt.savefig(path+fileID+"_CentroidPosition2_Graph"+string+".jpg", dpi=300) 
    
    plt.figure(figsize=(6, 6), dpi=200)
    plt.plot(x,width, 'sb',label='Width ', markersize=5)
    plt.plot(x,heigth, 'or',label='Heigth', markersize=5)
    plt.title("Size of Particle"+string+" "+" " + fileID)
    plt.xlabel("Picture #  ")
    plt.ylabel("Width / Hiegth [pixels]")
#    plt.ylim((70,82))
    plt.legend(loc='best')
    plt.savefig(path+fileID+"_width-heigth_Graph"+string+".jpg", dpi=300) 
    
    plt.figure(figsize=(6, 6), dpi=200)
    plt.plot(x,D, 'sb',label='Taylor Deformation Parameter ', markersize=5)
    plt.title("Taylor Deformation Parameter"+string+" "+" " + fileID)
    plt.ylabel("Taylor Deformation Parameter")
    plt.xlabel("Picture # ")
#    plt.ylim((70,82))
    plt.legend(loc='best')
    plt.savefig(path+fileID+"_Talyor-Deformation_Graph"+string+".jpg", dpi=300)
    
    plt.figure(figsize=(6, 6), dpi=200)
    plt.plot(x,xCentralDistance, 'sb',label='x-Distance along centre ', markersize=5)
    plt.plot(x,yCentralDistance, 'ro',label='y-Distance along centre ', markersize=5)
    plt.title("Central Distance"+string+" "+" " + fileID)
    plt.xlabel("Picture #")
    plt.ylabel("Central distance [pixels]")
    plt.legend(loc='best')
    plt.savefig(path+fileID+"_CentralDistance"+string+".jpg", dpi=300)
    
    plt.figure(figsize=(6, 6), dpi=200)
    plt.plot(x,xCentralDistance, 'sb',label='x-Distance along centre ', markersize=5)
    plt.plot(x,yCentralDistance, 'ro',label='y-Distance along centre ', markersize=5)
    plt.plot(x, heigth, 'y<',label='Height (y)', markersize=5)
    plt.plot(x, width, 'gh',label='Width (x)', markersize=5)

    plt.title("Size"+string+" "+" " + fileID)
    plt.xlabel("Picture #")
    plt.ylabel("Central distance [pixels]")
    plt.legend(loc='best')
    plt.savefig(path+fileID+"_Size"+string+".jpg", dpi=300)
    
    plt.figure(figsize=(6, 6), dpi=200)
    plt.plot(x,meanC1, 'sb',label='Curvature 1st Quadrant ', markersize=3)
    plt.plot(x,meanC2, 'ro',label='Curvature 2nd Quadrant ', markersize=3)
    plt.plot(x,meanC3, 'g<',label='Curvature 3rd Quadrant ', markersize=3)
    plt.plot(x,meanC4, 'mh',label='Curvature 4th Quadrant ', markersize=3)

    plt.title("Curvature "+string+" "+" " + fileID)
    plt.xlabel("Picture #")
    plt.ylabel("Curvature [arbitrary units]")
    plt.legend(loc='best')
    plt.savefig(path+fileID+"_Size"+string+".jpg", dpi=300)
    return r
        
#==============================================================================
# Miscellaneous Functions
#==============================================================================       
    
def find_batchName(path):
    '''Find the Batch name from the path name. Convention is that it is writen
    in the folder name containing the images, starting with "BatchXXXXX"'''
    dSlash = os.sep
       
    d1=path.rfind(dSlash,1,-3) #find first seperator
    reslt=path[d1+1:-1]
    
    if '\\' in reslt:
        reslt=reslt[:-1]
    
    print('find_batchName: %s ' %reslt)

    return reslt
    
    
def rotateImage(image, angle):
    '''Rotate an image by angle radians in the anti-clockwise direction'''
    try:
        row,col = image.shape
    except:
        row,col, _ = image.shape
    center=tuple(np.array([row,col])/2)
    rot_mat = cv2.getRotationMatrix2D(center,angle,1.0)
    new_image = cv2.warpAffine(image, rot_mat, (col,row))
    return new_image
    
def runningMean(x, N):  
    smoothed=np.zeros((len(x)), dtype=float)   
    smoothed[:-N+1]=np.convolve(x, np.ones((N,))/N, mode='valid')
    smoothed[-N:]=np.average(x)
    return smoothed
    
def gradient(x, h=1.0):
    '''2nd order central difference'''
    #TODO: fix so first point isn't zero and make output lenght == input length
    dx = np.zeros((len(x)-4))
    for i in range(3,len(x)-2):
        dx[i-2] = (1/12*x[i-2]-2/3*x[i-1] + 2/3*x[i+1] - 1/12 * x[i+2])/(h+0.0)
    return dx
    
def _differntiate(x, FPS):
    '''Take an array that contains entries with error and 
    differentiates using numpy gradient (with order 1 at the boundaries).'''
    
    dt = np.zeros((len(x)), dtype=float)
    dt.fill(1.0/FPS)   

    dx = np.gradient(x, dt)
    
    # Go through and set the entries that for which there is insufficient 
    # information to the ERROR_CODE
    for ii in range(len(x)):
        if ii ==0:
            if x[ii] == ERROR_CONST or x[ii+1] == ERROR_CONST:
                dx[ii]=ERROR_CONST
        elif ii == len(x)-1:
            if x[ii] == ERROR_CONST or x[ii-1] == ERROR_CONST:
                dx[ii]=ERROR_CONST
        else:
            if x[ii-1] == ERROR_CONST or x[ii+1] == ERROR_CONST:
                dx[ii]=ERROR_CONST
    return dx
    
    
def smooth(x,window_len=11,window='hanning'):
    """smooth the data using a window with requested size.
    
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
    
    input:
        x: the input signal 
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal
        
    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)
    
    see also: 
    
    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter
 
    TODO: the window parameter could be the window itself if an array instead of a string
    NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
    
    Taken from: http://wiki.scipy.org/Cookbook/SignalSmooth
    """

    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."

    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."


    if window_len<3:
        return x


    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"


    s=np.r_[x[window_len-1:0:-1],x,x[-1:-window_len:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')

    y=np.convolve(w/w.sum(),s,mode='valid')
#    return y
    yreturn = y[(window_len/2-1):-(window_len/2)]
    return yreturn[:-1]

def savitzky_golay(y, window_size, order, deriv=0, rate=1):
    r"""Smooth (and optionally differentiate) data with a Savitzky-Golay filter.
    The Savitzky-Golay filter removes high frequency noise from data.
    It has the advantage of preserving the original shape and
    features of the signal better than other types of filtering
    approaches, such as moving averages techniques.
    Parameters
    ----------
    y : array_like, shape (N,)
        the values of the time history of the signal.
    window_size : int
        the length of the window. Must be an odd integer number.
    order : int
        the order of the polynomial used in the filtering.
        Must be less then `window_size` - 1.
    deriv: int
        the order of the derivative to compute (default = 0 means only smoothing)
    Returns
    -------
    ys : ndarray, shape (N)
        the smoothed signal (or it's n-th derivative).
    Notes
    -----
    The Savitzky-Golay is a type of low-pass filter, particularly
    suited for smoothing noisy data. The main idea behind this
    approach is to make for each point a least-square fit with a
    polynomial of high order over a odd-sized window centered at
    the point.
    Examples
    --------
    t = np.linspace(-4, 4, 500)
    y = np.exp( -t**2 ) + np.random.normal(0, 0.05, t.shape)
    ysg = savitzky_golay(y, window_size=31, order=4)
    import matplotlib.pyplot as plt
    plt.plot(t, y, label='Noisy signal')
    plt.plot(t, np.exp(-t**2), 'k', lw=1.5, label='Original signal')
    plt.plot(t, ysg, 'r', label='Filtered signal')
    plt.legend()
    plt.show()
    References
    ----------
    .. [1] A. Savitzky, M. J. E. Golay, Smoothing and Differentiation of
       Data by Simplified Least Squares Procedures. Analytical
       Chemistry, 1964, 36 (8), pp 1627-1639.
    .. [2] Numerical Recipes 3rd Edition: The Art of Scientific Computing
       W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery
       Cambridge University Press ISBN-13: 9780521880688
    """
    import numpy as np
    from math import factorial
    
    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except ValueError, msg:
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order+1)
    half_window = (window_size -1) // 2
    
    # precompute coefficients
    b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
    
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve( m[::-1], y, mode='valid')
    
def shiftImage(img, xshift, yshift):
    '''Shift Image by xshift, yshift'''
    rows,cols = img.shape
    M = np.float32([[1,0,xshift],[0,1,yshift]])
    dst = cv2.warpAffine(img,M,(cols,rows))
    return dst
    
def increaseContrast(img, phi=0.5, theta = 1.0):
    '''Increase contrast by decreasing intensity such that dark pixels 
    become much darker and brighr pixels become sligthly darker'''
    maxIntensity = 255.0 # depends on dtype of image data
    newImg = (maxIntensity/phi)*(img/(maxIntensity/theta))**2
    return np.array(newImg,dtype=np.uint8)
#    return np.array(newImg)
    

