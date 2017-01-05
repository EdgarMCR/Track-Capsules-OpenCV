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
import matplotlib.patches as mpatches
import numpy as np
import cv2
import shutil
import scipy as sp
#from scipy import interpolate
from scipy      import optimize

import os
#import time

import general as gen

ERROR_CONST=-1
AVERAGINGWINDOW = 6
AVERAGINGREPEATS = 1

FOR_LATEX = False  #sets what save names are

class DataRun:
    """The Data Class holds the general measurments for an experimental run, 
    independent of the geometries it is ran in"""
    
    def __init__(self, path1, pPmm, FPS):
        self.path = path1
        self.pixelsPmm=pPmm
        self.FPS = FPS
        #from file name
        self.volumFlux = ''
        self.bathName = ''
        self.name  = ''

        #read from file
        self.centroid_x = ERROR_CONST
        self.centroid_y = ERROR_CONST
        self.d12 = ERROR_CONST
        self.d12rect = ERROR_CONST
        self.width = ERROR_CONST
        self.height = ERROR_CONST
        self.area = ERROR_CONST
        self.perimeter = ERROR_CONST
        self.curvature1 = ERROR_CONST
        self.curvature2 = ERROR_CONST
        self.curvature3 = ERROR_CONST
        self.curvature4 = ERROR_CONST
        
        #calculated
        self.speed = ERROR_CONST
        self.vel_x = ERROR_CONST
        self.vel_y = ERROR_CONST
        self.accelerationMag = ERROR_CONST
        self.acc_x = ERROR_CONST
        self.acc_y = ERROR_CONST
    
        #read in the data
        self._readResultsFile()
        self._findAccelerationAndVelocityFromPosition()
        self. _getVolumnFlux()
        self._setName()
        
    def _readResultsFile(self):        
        self.bathName=gen.find_batchName(self.path)
        path_to_file = os.path.join(self.path, self.bathName+'_Results.txt')
        
        try:
            numOfLines= _file_len(path_to_file)
        except:
#            print("Didn't work with: %s \n\n gen.find_batchName(self.path) : " %path_to_file)
#            print(gen.find_batchName(self.path))
            pa , folder = os.path.split(self.path)
            path_to_file = os.path.join(self.path, folder+'_Results.txt')
            numOfLines= _file_len(path_to_file)
        
        #This opens a handle to your file, in 'r' read mode
        file_handle = open(path_to_file, 'r')
        
        # Read in all the lines of your file into a list of lines
        lines_list = file_handle.readlines()
        
        cen_x=np.zeros((numOfLines), dtype=float)
        cen_y=np.zeros((numOfLines), dtype=float)
        perimeter=np.zeros((numOfLines,1), dtype=float)
        area=np.zeros((numOfLines,1), dtype=float)
        width =np.zeros((numOfLines), dtype=float)
        height =np.zeros((numOfLines,1), dtype=float)
        d12 =np.zeros((numOfLines,1), dtype=float)
        d12rect =np.zeros((numOfLines), dtype=float)
        centralWidth = np.zeros((numOfLines), dtype=float)
        centralHeight = np.zeros((numOfLines), dtype=float)
        curvature1 = np.zeros((numOfLines), dtype=float)
        curvature2 = np.zeros((numOfLines), dtype=float)
        curvature3 = np.zeros((numOfLines), dtype=float)
        curvature4 = np.zeros((numOfLines), dtype=float)
        
        
        indexC=0
        index=0
        linecounter=-1
        
        for line in lines_list:
            linecounter=linecounter+1
            if linecounter==0 :
                continue
            
            entries=line.split()
            
            if len(entries) < 14:
                continue

            perimeter[indexC] = float(entries[3])
            area[indexC]=float(entries[4])
            cen_x[indexC]=float(entries[5])
            cen_y[indexC]=float(entries[6])
            width[indexC]=float(entries[7])
            height[indexC]=float(entries[8])
            d12[indexC]=float(entries[9])
            indexC += 1
            
            if len(entries) >14:
                d12rect[index] = float(entries[14])
                
            if len(entries) > 20:
                centralWidth[index] = float(entries[15])
                centralHeight[index] = float(entries[16])
                curvature1[index] = float(entries[17])
                curvature2[index] = float(entries[18])
                curvature3[index] = float(entries[19])
                curvature4[index] = float(entries[20])
                  
            index=index+1  
            
        self.centroid_x = cen_x
        self.centroid_y =cen_y
        self.d12 = d12
        self.d12rect = d12rect
        self.width = width
        self.height = height
        self.area = area
        self.perimeter = perimeter
        self.curvature1 = curvature1
        self.curvature2 = curvature2
        self.curvature3 = curvature3
        self.curvature4 = curvature4

    def _checkDisplacement(self):
        '''Check whether centroid position makes sudden jumps.'''
        pass
        

    def _findAccelerationAndVelocityFromPosition(self):
        ''' Differentiate position to get velocity and acceleration'''
        dx =gen._differntiate(self.centroid_x, self.FPS)
        dy =gen._differntiate(self.centroid_y, self.FPS)
        
        d2x =gen._differntiate(dx, self.FPS)
        d2y =gen._differntiate(dy, self.FPS)
    
        self.speed = np.sqrt( np.power( dx ,2) + np.power(dy ,2) )
        self.vel_x = dx
        self.vel_y = dy
        self.accelerationMag = np.sqrt( np.power( d2x ,2) + np.power(d2y ,2) )
        self.acc_x = d2x
        self.acc_y = d2y
        
    def _getVolumnFlux(self):
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

        self.volumFlux  = float(volFlux.replace("p", "."))
        
    def _setName(self):
        '''Set the name from path'''
#        print('Entering _setName(self)')
#        print('self.path = %s' %self.path)
        d1 = self.path.rfind(os.sep,1, -2)
        name = self.path[d1+1:]
#        print('name = %s' %name)
        self.name = name.replace(os.sep, '')
        
    def _centroidsWithoutErrors(self):
        '''return centroids without entries that are empty (i.e. ERROR_CONST).
        Keeping them in general as otherwise break the time dependence, i.e. 
        every measurment is 1/FPS from the next'''
        
        #TODO: remove implicit assumption that the values that are missing in x
        # and y are the same
        return self.centroid_x[self.centroid_x != ERROR_CONST], self.centroid_y[self.centroid_y != ERROR_CONST]
        
    
    def _writeToFile(self, headerString, String):
        '''Generic function that writes string to a file
        
        The path to the folder will be used to extrtact batch name and run name
        This assumes that folder name contains batch name, volume flux, frames
        per seconds and run number
                
        headerString    Name of measurements that are provided in string
        
        String          A string of values, tab seperated that are relevant to self.name
                        the specific geometry
        '''
        
        pathData = _pathExperimentRsltFile(self.path)  
        
        print('pathData = %s' %pathData)       
        
#        print('os.path.isfile(pathData) = %r' %os.path.isfile(pathData))  
        #if we write to file for the first time, write the header
        if not os.path.isfile(pathData): #check whether the file for this run has been started
            if headerString[-2:] != '\n':  headerString = headerString + '\n'
            fileData = open(pathData, 'w')
            fileData.write('Vol Flux \t Frames/second \tname \t pixel/mm \t %s' %headerString)
            fileData.close()
        
        
        dataString = '%.1f \t %.1f \t %s \t %.4f \t %s'  %(self.volumFlux, self.FPS, self.name, self.pixelsPmm, String)
        if String[-2:] != '\n':  dataString = dataString + '\n'
        
        if self.name == '':
            self._setName()
            import warnings
            warnings.warn('self.name not set. Trying to remedy issues', UserWarning)
            
        if (self.name+' ') in open(pathData).read(): #check whether the name is in the file
            f = open(pathData,"r")
            lines = f.readlines()
            f.close()

            fileData = open(pathData, 'w')
            for line in lines:
                if (self.name+' ') in line:
                    fileData.write(dataString)
                else:
                    fileData.write(line)      
            fileData.close()
        else: #if this hasn't been written to file before, simple append it to the end
            fileData = open(pathData, 'a')
            fileData.write(dataString)
            fileData.close()

        self._removeDuplicateLines(self.path)
        self._sortLines(self.path)
    
    def _sortLines(self, path):
        ''' Sorts entries according to volume flux '''
        File = _pathExperimentRsltFile(path)  
        directory = os.path.dirname(File)
        backup = os.path.join(directory, 'tempBackup_Sort.txt')
        shutil.copyfile(File, backup) #create backup
        f = open(File,"r"); lines = f.readlines(); f.close()
    #    os.remove(directory + filename) #delete original
#        print(lines)
        listOfNames=[]
        for line in lines:
            entriesLine=line.split('\t')     
            entriesLine[-1] = entriesLine[-1][:-2] #what am I doing here?
            if entriesLine[0] == 'Vol Flux ':
                header=line
                continue
    #        print(entriesLine)
            name = entriesLine[1].strip() # 2nd entry is the name
            
            #number
            #Batch120615-001-#4-10FPS-5mlPmin-3 
            d=name.rfind('-')
            if d ==-1 or d<len(name)-8:
                d=name.rfind('_')
            number = float(name[d+1:])
            
            #find volumn flux
            d1=name.find('mlPmin')
            d2=name[:d1].rfind('-') 
            if (d1-d2) > 6:
                d2=name[:d1].rfind('_')
            d3=name[d2:d1].find('m')
            if d3==-1:
                volFlux=name[d2+1:d1]
            else:
                volFlux=name[d2+1:d2+d3]
            volFlux = volFlux.replace("p", ".")
    #        print('volFlux = %s \t' %volFlux, end="")
            volFlux=float(volFlux)
    #        print('volFlux = %f' %volFlux)
            
            temp=_entryClass()
            temp.line = line
            temp.fn = name
            temp.number= number
            temp.Q=volFlux 
            listOfNames.append(temp)
            
        newlist = sorted(listOfNames, key=lambda _entryClass: (_entryClass.number, _entryClass.Q ))  
        fileD = open(File, 'w')
        fileD.write(header)
        for cls in newlist:
            fileD.write(cls.line)
        fileD.close()            
        os.remove(backup) #delete backup copy
        
    def _removeDuplicateLines(self, path):
        ''' For results files only '''
        filename = _pathExperimentRsltFile(path)  
        directory = os.path.dirname(filename)
        #create backup
        backup = directory + 'tempBackup_RemoveDuplicateLines.txt'
        shutil.copyfile(filename, backup)
        
        f = open(filename,"r")
        lines = f.readlines()
        f.close()
        
        os.remove(filename) 
        
        listOfNames=[]
        duplicates=0
        fileD = open(filename, 'w')
        
        for line in lines:
            entriesLine=line.split('\t')     
            entriesLine[-1] = entriesLine[-1][:-2] #what am I doing here?
            name = entriesLine[2].strip() # 3rd entry is the name
            
            #TODO: write this that XXX-1 is not in XXX-14
            if not (name+' ') in listOfNames:
                listOfNames.append(name+' ')
                fileD.write(line) 
            else:
                duplicates += 1
                
        fileD.close()            
        os.remove(backup) #delete backup copy
#        print(listOfNames)
        print('%d duplicates removed' %duplicates)
        
    def _findCentralWidthAndHeight(self, indexes, plot=False, pos=[]):
        '''Takes indexes and loads the outlines for these indexes from file. 
        Then interplates the outline and finds the central x and y distance.
        Returns average over all indexes
        
        pos is an add-hoc vector added here to pass information from the 
        semi-sphere / half-cylinder setup to this function. This means the function
        should be in the semi-sphere file but I don't want to copy the stuff over now
        TODO: fix this
        pos containts top, bottom 
        '''
        
        xd=[]; yd=[];
        ydf=[];# forward distance
        for jj in indexes:
            x, y =  self._loadOutline(jj)
            cD, xS, yS = gen._splineInterpolationOfOutline(x,y,1000)
            
            ctr_x = self.centroid_x[jj] + gen.BORDERSIZE
            ctr_y = self.centroid_y[jj] + gen.BORDERSIZE
            
            tol=2 # 2 pixels 
            indexesX  = [i for i in range(len(xS)) if xS[i] < ctr_x+tol  and xS[i] > ctr_x-tol ]
            ymin = np.min(yS[indexesX])
            ymax = np.max(yS[indexesX])
            
            if len(pos) > 0:
                print("top = %d bottom = %d centreling = %.1f" %(pos[0], pos[1], (pos[0]+pos[1]+0.0)/2.0))
                centreline = (pos[0]+pos[1]+0.0)/2.0
                if self.centroid_y[jj] < centreline: #capsule is above obstacle
#                    print("Top of cylinder: self.centroid_y[%d] = %.1f" %(jj, self.centroid_y[jj]))
                    topHalfCylinder = centreline - 4.0*self.pixelsPmm + gen.BORDERSIZE
                    if ymax > topHalfCylinder:
                        print("ymin = %d ymax = %d topHalfCylinder = %d gen.BORDERSIZE = %d" %(ymin, ymax, topHalfCylinder, gen.BORDERSIZE))                        
                        ymax = topHalfCylinder
                        
                else: #capsule below half-cylinder
#                    print("Bottom of cylinder: self.centroid_y[%d] = %.1f" %(jj, self.centroid_y[jj]))                    
                    bottomHalfCylinder = centreline + 4.0*self.pixelsPmm + gen.BORDERSIZE
                    if ymin < bottomHalfCylinder:
                        print("ymin = %d bottomHalfCylinder = %d" %(ymin, bottomHalfCylinder))                          
                        ymin=bottomHalfCylinder
            indexesY  = [i for i in range(len(yS)) if yS[i] < ctr_y+tol  and yS[i] > ctr_y-tol ]
            xmin = np.min(xS[indexesY])
            xmax = np.max(xS[indexesY])
            
            xd.append(xmax - xmin)
            if (ymax - ymin)/self.pixelsPmm > 2.0:
                yd.append(ymax - ymin)
            
            if plot and len(yd) > 0:
                # Plot on image to check
                img = self._loadImage(jj)
                cv2.line(img, (int(ctr_x -gen.BORDERSIZE), int(ctr_y - yd[-1]/2.0 - gen.BORDERSIZE)), (int(ctr_x - gen.BORDERSIZE), int(ctr_y + yd[-1]/2.0) - gen.BORDERSIZE), (255, 0, 0))
                cv2.line(img, (int(ctr_x - xd[-1]/2.0 - gen.BORDERSIZE), int(ctr_y - gen.BORDERSIZE)), (int(ctr_x + xd[-1]/2.0 - gen.BORDERSIZE), int(ctr_y-gen.BORDERSIZE)), (0, 255, 255))
                
                for rr in range(len(x)):
                    cv2.circle(img, (x[rr]-gen.BORDERSIZE,y[rr]-gen.BORDERSIZE), 1, (255, 255, 0))
                
                for rr in indexesY:
                    cv2.circle(img, (int(xS[rr]-gen.BORDERSIZE), int(yS[rr]-gen.BORDERSIZE)), 1, (255, 0, 255))
                
                for rr in indexesX:
                    cv2.circle(img, (int(xS[rr]-gen.BORDERSIZE), int(yS[rr]-gen.BORDERSIZE)), 1, (255, 0, 255))
                dire = os.path.join(self.path, 'CentralDistance')
                if not os.path.exists(dire): os.makedirs(dire)
                    
                sn = os.path.join(dire, '%s_CentralHeight_%04d.png' %(self.name, jj))
                cv2.imwrite(sn, img)
                sn2 = os.path.join('C:\\Users\\mbbxkeh2\\CentralDistancePlots', '%s_CentralHeight_%04d.png' %(self.name, jj))
                cv2.imwrite(sn2, img)
                
                
#            #Trying to get vertical distance forward of centre
#            indexesXforward  = [i for i in range(len(xS)) if xS[i] < ctr_x - xd[jj]/3.0+tol  and xS[i] > ctr_x- xd[jj]/3.0-tol ]
#            ymin = np.min(yS[indexesXforward])
#            ymax = np.max(yS[indexesXforward])
#            ydf.append(ymax - ymin)
        return np.mean(xd), np.mean(yd)  #, np.mean(ydf)
        
    def _findHeightAtXpos(self, indexes, xpos):
        '''Takes indexes and loads the outlines for these indexes from file. 
        Then interplates the outline and finds the  x distance at xpos.
        Returns average over all indexes
        '''
        
        yd=[];
        for jj in indexes:
            x, y =  self._loadOutline(jj)
            cD, xS, yS = gen._splineInterpolationOfOutline(x,y,1000)
            
            tol=2 # 2 pixels 
            indexesX  = [i for i in range(len(xS)) if xS[i] < xpos+tol  and xS[i] > xpos-tol ]
            
            if len(indexesX) != 0 :
                ymin = np.min(yS[indexesX])
                ymax = np.max(yS[indexesX])
                
                yd.append(ymax - ymin)
            
        return np.mean(yd)
        
    def _loadOutline(self, index):
        '''Load the indexth outline and return x,y'''
        fn = 'OutlineBinary_%s_%d.npy' %(self.name, index)
        p = os.path.join(self.path, 'Outlines', fn)
        cnt = np.load(p)
            
        x=np.array(cnt[:,0,0])
        y=np.array(cnt[:,0,1])
        return x,y
        
    def _loadImage(self, index):
        '''Loaf indexth image '''
        fn = '%s_%04d.png' %(self.name, index)
        p = os.path.join(self.path, fn)
        return cv2.imread(p,1)
        
        
    def fitCircle(self, indexes):
#        print('Starting "fitCircle" in analysis.py with len(indexes) = %d' %(len(indexes)))
        d=[]; errd=[]
        
        for ii in indexes:
            x, y = self._loadOutline(ii)
            
#            sp = os.path.join(self.PP.path, 'Outlines')
#            td, terrd = self._fitCircleToXY(x,y, savepath = sp, nr=ii)
            td, terrd = self._fitCircleToXY(x,y)
            d.append(td); errd.append(terrd)
#            print('ii=%d \t d = %f =/- %f ' %(ii, td/self.pixelsPmm, terrd/self.pixelsPmm))
        
        return np.mean(d), np.mean(errd)  #TODO: error analyisis
            
    
    def calc_R(self, x, y, xc, yc):
        """ calculate the distance of each 2D points from the center (xc, yc) """
        return np.sqrt(np.power(x-xc, 2) + np.power(y-yc, 2))
    
    def f_2(self, c, x,y):
        """ calculate the algebraic distance between the data points and the mean circle centered at c=(xc, yc) """
        Ri = self.calc_R(x,y,*c)
        return Ri - np.mean(Ri)
        
    def _fitCircleToXY(self, x,y, savepath = None, nr=0):
        '''From: http://scipy-cookbook.readthedocs.org/items/Least_Squares_Circle.html'''
        
        # coordinates of the barycenter
        x_m = np.mean(x);  y_m = np.mean(y)
        
        center_estimate = x_m, y_m
        center_2, ier = optimize.leastsq(self.f_2, center_estimate, args=(x,y,))
        
        xc_2, yc_2 = center_2
        Ri_2       = self.calc_R(x,y, *center_2)
        R_2        = Ri_2.mean()
        residu_2   = np.sum(np.power(Ri_2 - R_2, 2.0))
        std = np.sqrt(residu_2/(len(Ri_2)+0.0))
        
        if savepath != None:
            fig = plt.figure()
            plt.plot(x,y, 'bs', label='Data')
    
            leng=100; thetaStep = 2.0*np.pi/leng; xr=[]; yr=[]        
            for kk in range(leng):
                theta = kk* thetaStep
                xr.append(R_2*np.cos(theta) + xc_2)
                yr.append(R_2*np.sin(theta) + yc_2)
            
            plt.plot(xr, yr, 'ro', label='fit')
            plt.legend(loc='best')
            plt.axes().set_aspect('equal')
            plt.savefig(os.path.join(savepath,'CircleFitOutline_%d.png' %nr))
            plt.close(fig)

        return 2*R_2, std #do error

def _get_foler_list(path):
    ''' List folders '''
    listFolders=os.listdir(path)
    indexDelet=[]
    for i in range(len(listFolders)):
        if os.path.isfile(path+listFolders[i]):
            indexDelet.append(i)

    for i in range(len(indexDelet)-1, -1, -1):
        del listFolders[indexDelet[i]]
        
    return listFolders
    
def wholeRun(func, path):
    """ Runs func with argument path+fodler. """
    import traceback
    listFolders = _get_foler_list(path)

    #find FPS info
    fpslist=[]
    for f in listFolders:
        sp1=f.upper().find('FPS')
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
        
    for fps in fpslist:
        runOneFPS(func, path, fps)
            

def runOneFPS(func, path, fps):
    """ Runs func for all folder with 'fpsFPS' in there name."""
    FPS = fps
    import traceback

    listFolders = _get_foler_list(path)
        
    foldersThatWorked=[]
#    print(listFolders)        
    for f in listFolders:
#        print('fname = %s' %(f.upper()))
        sp1=f.upper().find('FPS')
        
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
        func(path+f)
#            find_max_extend(directory+f)

    
    print('\nWorked in :')
    for f in foldersThatWorked:
        print(f)

def _pathExperimentRsltFile(path, isDirectory=False):
    '''Defines name of experiment rslt file'''
    
#    print('path: %s' %path)
    if not isDirectory:
        ss1=path.rfind(os.sep, 1,-1)
        directory = path[:ss1+1]
        #TODO: this directory is the actual folder with the images instead the folder
        #containing all the folders with images. Not sure how this possibly could have worked
        #The following line breaks backward compatibility
#        ss=directory.rfind(os.sep, 1,-1)
#        directory = directory[:ss+1]
        
    else:
        directory = path

#    print('directory: %s' %directory)
    ss2=directory.rfind(os.sep, 1,-1)
    ss3=directory.rfind(os.sep, 1,ss2)
#    print('ss2 = %d, ss3 = %d' %(ss2, ss3))
#    print('directory[ss3+1:ss2] = %s' %(directory[ss3+1:ss2]))

    return os.path.join(directory,  directory[ss3+1:ss2]+'-'+directory[ss2+1:].replace(os.sep, '')+'_Results.txt')
        
class ResultsClass:
    """A class to hold and read results from experiments"""
    
    def __init__(self, directory):
        self.directory = directory
        self.volumeFlux =[]
#        self.volFlux2=[]
        self.fps=[]
        self.name =[]
        self.batchName =''
        self.pixelsPmm = []
        
        #read in the data
        self._readData()
        self._setBatchName()
        
        self.directoryName = self._setDirName()

    def _setBatchName(self):
        '''find batch name'''
        s = self.name[0]
        d1=s.find('-')

        #following the first dash, there is either one or three numbers
        d2 = s.find('_', d1+1, d1+5)
        if d2 == -1:
            d2 = s.find('-', d1+1, d1+5)
#            print('d2 = ' + str(d2))
        self.batchName = s[:d2]
        
    def _setDirName(self):
        ''' Find name '''
        pa , folder = os.path.split(self.directory)
        pa , folder = os.path.split(pa)
        return folder
        

    def _readData(self):
        """ Read the generic entries in the experiment file:
        Vol Flux 
        Frames/second 
        name 
        """
        self.filePath =_pathExperimentRsltFile(self.directory, isDirectory=True)
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

            if len(entries)>=4:
                self.volumeFlux.append(float(entries[0]))
                self.fps.append(float(entries[1]))
                self.name.append(entries[2])
                self.pixelsPmm.append(float(entries[3]))
                ii += 1
            
        self.volumeFlux=np.array(self.volumeFlux)
        self.fps=np.array(self.fps)
        self.pixelsPmm=np.array(self.pixelsPmm)
        

    def averageByQ(self, xInput, grZero=False):
        '''     Average the values for each volumn flux and return average and std     '''
        if not hasattr(self, "volumeFluxes"):
            #get all volumn fluxes
            self.volumeFluxes=[]
            for VF in self.volumeFlux:
                if VF not in self.volumeFluxes:
                    self.volumeFluxes.append(VF)
            self.volumeFluxes = np.array(self.volumeFluxes)
        
        lenvfl=len(self.volumeFluxes)
        x=np.zeros((lenvfl))
        xSTD=np.zeros((lenvfl))
        
        counter=0
        for vfl in self.volumeFluxes:
            tempx=[]
#            print('Counter =%d' %counter)
            for i in range(len(self.volumeFlux)):
                if self.volumeFlux[i] == vfl:
                    tempx.append(xInput[i])
            
            if grZero:
                tempx = [kk for kk in tempx if kk > 0.0]
                
            x[counter]=np.average(tempx)
            xSTD[counter]=np.std(tempx, ddof=1)
            
            counter +=1
            
        return x, xSTD
        
    def averageByQWithoutErrors(self, xInput, grZero=False):
        '''     Average the values for each volumn flux and return average and std     '''
        if not hasattr(self, "volumeFluxes"):
            #get all volumn fluxes
            self.volumeFluxes=[]
            for VF in self.volumeFlux:
                if VF not in self.volumeFluxes:
                    self.volumeFluxes.append(VF)
            self.volumeFluxes = np.array(self.volumeFluxes)
        
        lenvfl=len(self.volumeFluxes)
        x=np.zeros((lenvfl))
        xSTD=np.zeros((lenvfl))
        
        counter=0
        for vfl in self.volumeFluxes:
            tempx=[]
#            print('Counter =%d' %counter)
            for i in range(len(self.volumeFlux)):
                if self.volumeFlux[i] == vfl:
                    if xInput[i] != ERROR_CONST:
                        tempx.append(xInput[i])
            
            if grZero:
                tempx = [kk for kk in tempx if kk > 0.0]
                
            x[counter]=np.average(tempx)
            xSTD[counter]=np.std(tempx, ddof=1)
            
            counter +=1
            
        return x, xSTD
        
    def averageBy(self,Q, xInput, grZero=False):
        '''     Average the values for each volumn flux and return average and std     '''

        #get all volumn fluxes
        volumeFluxes=[]
        for VF in Q:
            if VF not in  volumeFluxes:
                volumeFluxes.append(VF)
        volumeFluxes = np.array(volumeFluxes)
        
        lenvfl=len(volumeFluxes)
        x=np.zeros((lenvfl))
        xSTD=np.zeros((lenvfl))
        
        counter=0
        for vfl in volumeFluxes:
            tempx=[]
#            print('Counter =%d' %counter)
            for i in range(len(Q)):
                if Q[i] == vfl:
                    tempx.append(xInput[i])
            
            if grZero:
                tempx = [kk for kk in tempx if kk > 0.0]
                
            x[counter]=np.average(tempx)
            xSTD[counter]=np.std(tempx, ddof=1)
            
            counter +=1
            
        return volumeFluxes, x, xSTD
                
      
    def fitStraightLineScipy(self, x,y, yerr=None):
        m, RChi2, perr, p, rsquared = fitStraightLineScipy(x,y, yerr)
        return m, RChi2, perr, p, rsquared
        
    def _linearFit(self, x,y):
        '''simple linear fit'''
        p1, p2, fres = linearFit( x,y)
        return p1, p2, fres
        

#==============================================================================
# Miscellaneous Functions
#==============================================================================

def linearFit( x,y):
    '''simple linear fit'''
    from scipy.optimize import curve_fit
    
    popt, pcov = curve_fit(func,x, y,p0=(1.0,0.2))
    p1 = popt[0]
    p2 = popt[1]
    residuals = y - func(x,p1,p2)
    fres = sum(residuals**2)
    
    return p1, p2, fres
        
def func(x, p1,p2):
    '''Linear fit function'''
    return p1*x + p2


def fitStraightLineScipy(x,y, yerr=None):
    '''Fit straight line to data with zero intercept'''
    
    x = np.array(x)
    y = np.array(y)
    fitfunc = lambda p, x: np.array(p*x) # Target function
    errfunc = lambda p, x, y: fitfunc(p, x) -  np.array(y) # Distance to the target function
    
    m, pcov = sp.optimize.curve_fit(fitfunc, x, y,  p0=0.75)

    ss_err=(errfunc(m, x, y)**2).sum()
    ss_tot=((y-y.mean())**2).sum()
    rsquared=1-(ss_err/ss_tot)
    
    if yerr: 
        yerr = np.array(yerr)
        chi2 = (errfunc(m,x, y)**2).sum()/(yerr**2).sum()
        p = sp.stats.chisqprob(chi2, len(x)-1)
    else:
        chi2, p = sp.stats.chisquare(f_obs=x, f_exp=fitfunc(m, x), ddof=len(x)-1)

    RChi2 = chi2/(len(x)-1)  #the reduced chi-squared of the fit        
    perr = np.sqrt(np.diag(pcov))

    return m, RChi2, perr, p, rsquared

def _1d_remove_outliers(x,y, window=10, cutoff=1.8, check=False):
    assert len(x) == len(y)
    x_new=[]; y_new=[]    
        
    if check: adist=[]
    for ii in range(1, len(x) - 1):
        #find local average distance        
        dist = np.zeros(2*window)
        counter = 0
        
        start = -window
        end = window
        if ii < window and ii > len(x) - window:
            print("Window to large for this short array")
            start = -ii +2
            dif = window - (len(x) -ii)
            end = window - dif -1            
        elif ii < window:
            start = -ii +2
            end = window +(window-ii) -1
        elif ii > len(x) - window:
            dif = window - (len(x) -ii)
            start = -window -dif +1
            end = window - dif -1            
            
        dist = np.zeros(end-start)
#        print("ii = %d, len(x) = %d, start = %d end= %d" %(ii, len(x), start, end))        
        for hh in range(start, end):
            nr = ii +hh
            dist[counter] = np.abs(x[nr] - x[nr-1])
            counter += 1
        
        ave_dist = np.mean(dist)
        if check: adist.append(ave_dist)
        #check wheter distance is longer than cutoff 
        if(np.abs(x[ii] - x[ii-1]) < cutoff * ave_dist and
            np.abs(x[ii] - x[ii+1]) < cutoff * ave_dist):
            x_new.append(x[ii])
            y_new.append(y[ii])
            
    if check:
        plt.figure()
        plt.plot(np.arange(len(x)), x, 'ro')
        plt.plot(np.arange(1, len(x_new)+1), x_new, 'hc')
        ax = plt.gca()
        ax2 = ax.twinx()
        ax2.plot(np.arange(1,len(adist)+1), adist, 'bs')
        plt.show()

    return np.array(x_new), np.array(y_new)

    

def _remove_outliers(x,y, window=10, cutoff=1.8, check=False, sp=''):
    ''' Reject outliers from a timeseries of x,y coordinates '''
    
    assert len(x) == len(y)
    x_new=[]; y_new=[]    
    
#    distance=np.zeros(len(x)-1)    
#    for ii in range(1, len(x)):
#        distance[ii-1] = _2d_distance(x[ii],y[ii], x[ii-1], y[ii-1])    
#    ave_d = np.mean(distance)
        
    for ii in range(window+1, len(x) - window):
        #find local average distance     
        counter = 0
                
        start = -window
        end = window
        if ii < window and ii > len(x) - window:
            print("Window to large for this short array")
            start = -ii +2
            dif = window - (len(x) -ii)
            end = window - dif -1            
        elif ii < window:
            start = -ii +2
            end = window +(window-ii) -1
        elif ii > len(x) - window:
            dif = window - (len(x) -ii)
            start = -window -dif +1
            end = window - dif -1            
            
        dist = np.zeros(end-start)
#        print("ii = %d, len(x) = %d, start = %d end= %d" %(ii, len(x), start, end))        
        for hh in range(start, end):
            nr = ii +hh
            dist[counter] = _2d_distance(x[nr],y[nr], x[nr-1], y[nr-1])    
            counter += 1
        
        ave_dist = np.mean(dist)
        #check wheter distance is longer than cutoff 
        if(_2d_distance(x[ii],y[ii], x[ii-1], y[ii-1]) < cutoff * ave_dist and
            _2d_distance(x[ii],y[ii], x[ii+1], y[ii+1]) < cutoff * ave_dist):
            x_new.append(x[ii])
            y_new.append(y[ii])
    
    if check and sp != '':
        fig = plt.figure()
        plt.plot(x,y, 'ro', label='Original Data')
        plt.plot(x_new,y_new, 'bs', label='Filter Data')
        plt.legend(loc='best'    )
        plt.title("Window = %d, cutoff = %f" %(window, cutoff))
        plt.savefig(os.path.join('M:\\EdgarHaener\\check\\', sp), dpi=300)
        plt.close(fig)
    return np.array(x_new), np.array(y_new)
    
def _2d_distance(x,y, x2, y2):
    return np.sqrt(np.power(y2 -y, 2.0) + np.power(x2-x, 2.0))
    
def markerSymbols():
    '''A set of distingushable markers'''
    return ['o', 's', 'h', '>', '*', 'D', '^', 'p', 'd', '8',
            'o', 's', 'h', '>', '*', 'D', '^', 'p', 'd', '8',
            'o', 's', 'h', '>', '*', 'D', '^', 'p', 'd', '8',
            'o', 's', 'h', '>', '*', 'D', '^', 'p', 'd', '8',
            'o', 's', 'h', '>', '*', 'D', '^', 'p', 'd', '8',
            'o', 's', 'h', '>', '*', 'D', '^', 'p', 'd', '8',
            'o', 's', 'h', '>', '*', 'D', '^', 'p', 'd', '8',
            'o', 's', 'h', '>', '*', 'D', '^', 'p', 'd', '8',
            'o', 's', 'h', '>', '*', 'D', '^', 'p', 'd', '8',
            'o', 's', 'h', '>', '*', 'D', '^', 'p', 'd', '8',
            'o', 's', 'h', '>', '*', 'D', '^', 'p', 'd', '8',
            'o', 's', 'h', '>', '*', 'D', '^', 'p', 'd', '8',
            'o', 's', 'h', '>', '*', 'D', '^', 'p', 'd', '8',
            'o', 's', 'h', '>', '*', 'D', '^', 'p', 'd', '8',
            'o', 's', 'h', '>', '*', 'D', '^', 'p', 'd', '8',
            'o', 's', 'h', '>', '*', 'D', '^', 'p', 'd', '8']

def coloursForMarker( n=4):
    '''
    Returns n colours that will be dinstingusiable in greyscale
    From: http://colorbrewer2.org/
    '''
    assert(n<11)
    
    if n <=4:

        c1= '#d7191c'
        c2= '#fdae61'
        c3= '#abd9e9'
        c4= '#2c7bb6'
        return [c1, c2, c3, c4]

    
    elif n > 4 and n <= 6:

        c1= '#d73027'
        c2= '#fc8d59'
        c3= '#fee090'
        c4= '#e0f3f8'
        c5= '#91bfdb'
        c6= '#4575b4'
        return [c1, c2, c3, c4, c5, c6]
    
    elif n > 6 and n <=8:
        c1= '#b2182b'
        c2= '#d6604d'
        c3= '#f4a582'
        c4= '#fddbc7'
        c5= '#d1e5f0'
        c6= '#92c5de'
        c7= '#4393c3'
        c8= '#2166ac'
        return [c1, c2, c3, c4, c5, c6, c7, c8]
    
    elif n>8 and n<=10:
        c1= '#a50026'
        c2= '#d73027'
        c3= '#f46d43'
        c4= '#fdae61'
        c5= '#fee090'
        c6= '#e0f3f8'
        c7= '#abd9e9'
        c8= '#74add1'
        c9= '#4575b4'
        c10= '#313695'
        return [c1, c2, c3, c4, c5, c6, c7, c8, c9, c10]
        
def _saveFig(fig, directory, name):
    fig.tight_layout()
    sn = makeSaveNameLatexSafe(name)
    sp = os.path.join(directory, sn)
    if not os.path.exists(directory): os.makedirs(directory)
    plt.savefig(sp, dpi=300)
        
def _standardPlotSize(ax):
    '''Sets standard sizes for axis etc'''
    #fontsize
    ax.tick_params(labelsize=12)
    ax.yaxis.label.set_size(12)
    ax.xaxis.label.set_size(12)
    
    #gridlines
    ax.xaxis.grid(True, color='#D0D0D0')
    ax.yaxis.grid(True, color='#D0D0D0')
    
    xmin, xmax = ax.get_xlim()
    ax.set_xlim([0.95*xmin, xmax*1.05])

    
    return ax
   
def _shiftImage(img, xshift, yshift):
    '''Shift Image by xshift, yshift'''
    rows,cols = img.shape
    M = np.float32([[1,0,xshift],[0,1,yshift]])
    dst = cv2.warpAffine(img,M,(cols,rows))
    return dst
    
    
def _file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1
    
class _entryClass:
    """Store Q and number"""
    line=''     #whole line
    fn='-1'      #filename
    Q=-1        #millisecondtime   
    number=-1    #number in filename
    fps=-1
    
    def __repr__(self):
            return repr((self.fn, self.Q, self.number, self.fps))

def makeSaveNameLatexSafe(sn):
    if FOR_LATEX:
        return _save_for_latex(sn)
    else:
        return sn
    
def _save_for_latex(sn):
     sn=sn.replace('_', '-')
     sn=sn.replace('#', 'nr')
     return sn
        
def _removeERRORCONST(x):
    '''remove entries which are == ERROR_CONST'''
    return x[x != ERROR_CONST]
    
def is_outlier(points, thresh=3.5):
    """
    Returns a boolean array with True if points are outliers and False 
    otherwise.

    Parameters:
    -----------
        points : An numobservations by numdimensions array of observations
        thresh : The modified z-score to use as a threshold. Observations with
            a modified z-score (based on the median absolute deviation) greater
            than this value will be classified as outliers.

    Returns:
    --------
        mask : A numobservations-length boolean array.

    References:
    ----------
        Boris Iglewicz and David Hoaglin (1993), "Volume 16: How to Detect and
        Handle Outliers", The ASQC Basic References in Quality Control:
        Statistical Techniques, Edward F. Mykytka, Ph.D., Editor. 
        
    From:
    http://stackoverflow.com/questions/22354094/pythonic-way-of-detecting-outliers-in-one-dimensional-observation-data
    """
    if len(points.shape) == 1:
        points = points[:,None]
    median = np.median(points, axis=0)
    diff = np.sum((points - median)**2, axis=-1)
    diff = np.sqrt(diff)
    med_abs_deviation = np.median(diff)

    modified_z_score = 0.6745 * diff / med_abs_deviation

    return modified_z_score > thresh
    
def lsqfity(X, Y):
    """
    Calculate a "MODEL-1" least squares fit.

    The line is fit by MINIMIZING the residuals in Y only.

    The equation of the line is:     Y = my * X + by.

    Equations are from Bevington & Robinson (1992)
    Data Reduction and Error Analysis for the Physical Sciences, 2nd Ed."
    pp: 104, 108-109, 199.

    Data are input and output as follows:

    my, by, ry, smy, sby = lsqfity(X,Y)
    X     =    x data (vector)
    Y     =    y data (vector)
    my    =    slope
    by    =    y-intercept
    ry    =    correlation coefficient
    smy   =    standard deviation of the slope
    sby   =    standard deviation of the y-intercept
    
    from: http://stackoverflow.com/questions/27764220/retrieve-the-standard-deviation-of-the-y-intercept

    """

    X, Y = map(np.asanyarray, (X, Y))

    # Determine the size of the vector.
    n = len(X)

    # Calculate the sums.

    Sx = np.sum(X)
    Sy = np.sum(Y)
    Sx2 = np.sum(X ** 2)
    Sxy = np.sum(X * Y)
    Sy2 = np.sum(Y ** 2)

    # Calculate re-used expressions.
    num = n * Sxy - Sx * Sy
    den = n * Sx2 - Sx ** 2

    # Calculate my, by, ry, s2, smy and sby.
    my = num / den
    by = (Sx2 * Sy - Sx * Sxy) / den
    ry = num / (np.sqrt(den) * np.sqrt(n * Sy2 - Sy ** 2))

    diff = Y - by - my * X

    s2 = np.sum(diff * diff) / (n - 2)
    smy = np.sqrt(n * s2 / den)
    sby = np.sqrt(Sx2 * s2 / den)

    return my, by, ry, smy, sby  
    
def fitFuncSqrt(y, a, b):
    return a * np.sqrt(y) + b
    
def fitFuncPoly2(x, a, b, c):
    return a * np.power(x, 2) + b * np.power(x, 1) + c
    
def fitFuncExp(x, a, b, c):
    return c * np.exp(-b*x) + a
    
    
def redchisqg(ydata,ymod,deg=2,sd=None):  
      """  
      From: http://astropython.blogspot.co.uk/2012/02/computing-chi-squared-and-reduced-chi.html
     Returns the reduced chi-square error statistic for an arbitrary model,   
     chisq/nu, where nu is the number of degrees of freedom. If individual   
     standard deviations (array sd) are supplied, then the chi-square error   
     statistic is computed as the sum of squared errors divided by the standard   
     deviations. See http://en.wikipedia.org/wiki/Goodness_of_fit for reference.  
       
     ydata,ymod,sd assumed to be Numpy arrays. deg integer.  
       
     Usage:  
     chisq=redchisqg(ydata,ymod,n,sd)  
     where  
      ydata : data  
      ymod : model evaluated at the same x points as ydata  
      n : number of free parameters in the model  
      sd : uncertainties in ydata  
       
     Rodrigo Nemmen  
     http://goo.gl/8S1Oo  
       """  
      # Chi-square statistic  
      if sd==None:  
           chisq=np.sum((ydata-ymod)**2)  
      else:  
           chisq=np.sum( ((ydata-ymod)/sd)**2 )  
             
      # Number of degrees of freedom assuming 2 free parameters  
      nu=ydata.size-1-deg  
        
      return chisq/nu  
    