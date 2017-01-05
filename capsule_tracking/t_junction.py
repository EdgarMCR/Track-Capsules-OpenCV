# -*- coding: utf-8 -*-
"""
Created on Mon Nov 07 13:08:43 2016

Functions that try and read old T-Junction data

@author: mbbxkeh2
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

#Import folder where t-junction code is stored

import sys
sys.path.append('C:\\Users\\Edgar\\Dropbox\\PhD\\Python\\OpenCV\\T-Junction') #Machine Schuster G.21
sys.path.append('C:\\Users\\mbbxkeh2\\Dropbox\\PhD\\Python\\OpenCV\\T-Junction') #Aland Turing 2.105
sys.path.append('/home/magda/Dropbox/PhD/Python/OpenCV/T-Junction') #Aland Turing 2.105

import ReadOutputFile as TJ_ROF
import track_capsule_TJ as TR

ERROR_CONST=-1

def _find_Q_from_path(path):
    '''Find volumn Flux '''
#    print(path)
    name=path[:-2]
    d1=name.find('mlPmin')
    d2=name[:d1].rfind('_')
#    print('d1=%d d2 = %d' %(d1, d2))
    if d2 ==-1 or (d1-d2)>20:
        d2 = name[:d1].rfind('-')
    d3=name[d2:d1].find('m')
    if d3==-1:
        volFlux=name[d2+1:d1]
    else:
#        print('d3= %d' %d3)
#        print(name[d2+1:d2+d3])
        volFlux=name[d2+1:d2+d3]
#    print('Name: %s \nVolumn Flux = %s' %(name, volFlux))
    volFlux = volFlux.replace("p", ".")
#    print('volFlux = %s' %volFlux)
    volFlux=float(volFlux)
#    print('volFlux = %f' %volFlux)
    return volFlux


def tryout(directory):
    listFolders = ana._get_foler_list(directory)
    for f in listFolders:
        print("directory = %s" %directory)
        print("f = %s" %f)
        path = (os.path.join(directory, f))
        print("path = %s" %path)
        cent_x, cent_y, _, _, _, _, _, _, _, _, _, _, _, _, _ = TJ_ROF.readResultsFile(path, second=False)
        plt.plot(cent_x, cent_y)

def averageTrajectories(directory, savePath, pPmm, batchName, forPub=False):
    ''' Average all trajectories for each volumn flux '''    
    m = ana.markerSymbols()
    #step 1
    # Interpolate them y  is a function of x to get equal sampeling points
    from scipy.interpolate import interp1d
    
    listFolders = ana._get_foler_list(directory)
    minX=0; maxX=10000;
    inter_funcs=[]; Q=[]; xs=[]; ys=[]
    count=0
    for f in listFolders:
        try:
            path = (os.path.join(directory, f))
            cent_x, cent_y, _, _, _, _, _, _, _, _, _, _, _, _, _ = TJ_ROF.readResultsFile(path, second=False)

        except:
            print("\nDidn't work for \n%s\n\n" %f)
            continue
        print(f)
        #remove entries with no values from centroids
        cent_x[cent_x != ERROR_CONST]; cent_y[cent_y != ERROR_CONST]
        
        plt.plot(cent_x, cent_y, marker=m[count],linestyle="None")

#        #check if capsule went right
#        if np.mean(cent_x[-30:-10]) < np.mean(cent_x[10:30]):
#            print("reflecting data")
#            axis_of_reflextion = np.mean(cent_x[10:50])
#            
#            new_cent_x =[]
#            for x in cent_x:
#                dist = x - axis_of_reflextion 
#                new_cent_x.append(axis_of_reflextion +np.abs(dist)) 
#            cent_x = np.array(new_cent_x)


#        x, y = ana._remove_outliers(cent_x,cent_y, window=15, cutoff=1.5)
        x, y = ana._remove_outliers(cent_x,cent_y, window=15, cutoff=1.5, check=False, sp='%s_filtered.jpg' %f)
#        x2, y2 = ana._1d_remove_outliers(cent_x,cent_y, window=30, cutoff=1.5)
#        y2,x2 = ana._1d_remove_outliers(y2,x2, window=30, cutoff=1.5)
        

        try:
            inter_funcs.append(interp1d(x, y))
        except:
            print("couldn't interpolate for %s" %f)
            continue
        
        Q.append(_find_Q_from_path(path))
        xs.append(cent_x)
        ys.append(cent_y)
        if np.min(x) > minX: minX = np.min(x)
        if np.max(x) < maxX: maxX = np.max(x)
        
        count +=1       
            
    Q, inter_funcs = zip(*sorted(zip(Q, inter_funcs)))
      
    print(Q)
    
    qs = np.unique(Q)
        
    x = np.arange(minX, maxX, 1)

    for q in qs:
        fig = plt.figure(figsize=(8, 6), dpi=300); ax = fig.add_subplot(111)
        for ii in range(len(inter_funcs)):    
            if q == Q[ii]:
                plt.plot(xs[ii]/pPmm, ys[ii]/pPmm, marker=m[ii],linestyle="None", label='Q = %.1f' %(Q[ii]))
        plt.title("Trajectories-  " + batchName + " Q = %.1f" %q, fontsize=12)
        plt.xlabel('$X$ [mm]')
        plt.ylabel('$Y$ [mm]')
        plt.ylim((2, 4.5))
        plt.xlim((28, 29.5))
        plt.savefig(("M:\\EdgarHaener\\check\\"+"Trajectory_Q=%.1f" %q +"_" + batchName+".jpg"), dpi=300)
        
    minX = 635
    maxX = 1100
    
    ave_y=[]        
    for q in qs:
        ys=[]  
        for ii in range(len(Q)):
            if q == Q[ii]:
                ys.append(inter_funcs[ii](x))
        y=[]
        for kk in range(len(ys[0])):
            y.append(0.0)
            for jj in range(len(ys)):
                y[kk] += ys[jj][kk]                                
            y[kk] = y[kk] / len(ys)+0.0
        ave_y.append(np.array(y))

    c = ana.coloursForMarker(len(ave_y))
    fig = plt.figure(figsize=(8, 6), dpi=200); ax = fig.add_subplot(111)
    for ii in range(len(ave_y)):    
        plt.plot(x/pPmm, ave_y[ii]/pPmm, marker=m[ii], color=c[ii], label='Q = %.1f' %(qs[ii]))
    
#    if forPub:
#        minylim= 150; maxylim = 450
#        plt.xlim( (0, 900))
#    else:
#        maxylim = DRP.PP.ChannelBottom + 4.0*DRP.PP.pPmm
#        minylim = DRP.PP.ChannelTop - 8.0*.pPmm
#
#                    
#    plt.ylim( ( minylim/DRP.PP.pPmm, maxylim/DRP.PP.pPmm))
    
#    from matplotlib.patches import Rectangle
#    currentAxis = plt.gca()
#    print('EndChannel, ChannelBottom = %d, %d' %(DRP.PP.EndChannel/DRP.PP.pPmm, DRP.PP.ChannelBottom/DRP.PP.pPmm))
#    currentAxis.add_patch(Rectangle((DRP.PP.RightCornerInflow/DRP.PP.pPmm, DRP.PP.ChannelBottom/DRP.PP.pPmm), 
#                                    1.0*DRP.PP.pPmm/DRP.PP.pPmm, (maxylim - DRP.PP.ChannelBottom)/DRP.PP.pPmm, facecolor="grey"))
#
#    currentAxis.add_patch(Rectangle((0, DRP.PP.ChannelBottom/DRP.PP.pPmm), 
#                                    DRP.PP.LeftCornerInflow/DRP.PP.pPmm, (maxylim - DRP.PP.ChannelBottom)/DRP.PP.pPmm, facecolor="grey"))
#    
#    dist = DRP.PP.ChannelTop - minylim
#    currentAxis.add_patch(Rectangle((0, DRP.PP.ChannelTop/DRP.PP.pPmm), 
#                                    DRP.PP.RightCornerInflow/DRP.PP.pPmm, -dist/DRP.PP.pPmm, facecolor="grey"))

    plt.title("Average Trajectories-  " + batchName, fontsize=12)
    plt.xlabel('$X$ [mm]')
    plt.ylabel('$Y$ [mm]')
    plt.legend(loc='best', fontsize=6, ncol=1)
    fig.tight_layout()

    
#    ax_inset=fig.add_axes([0.4,0.6,0.3,0.3])
#    for ii in range(len(ave_y)):    
#        ax_inset.plot(x/pPmm, ave_y[ii]/pPmm, marker=m[ii], color=c[ii], label='Q = %.1f' %(qs[ii]))

#    currentAxis = plt.gca()
#    currentAxis.add_patch(Rectangle((DRP.PP.RightCornerInflow/DRP.PP.pPmm, DRP.PP.ChannelBottom/DRP.PP.pPmm), 
#                                    1.0*DRP.PP.pPmm/DRP.PP.pPmm, (maxylim - DRP.PP.ChannelBottom)/DRP.PP.pPmm, facecolor="grey"))
#
#    currentAxis.add_patch(Rectangle((0, DRP.PP.ChannelBottom/DRP.PP.pPmm), 
#                                    DRP.PP.LeftCornerInflow/DRP.PP.pPmm, (maxylim - DRP.PP.ChannelBottom)/DRP.PP.pPmm, facecolor="grey"))
#    
#    dist = DRP.PP.ChannelTop - minylim
#    currentAxis.add_patch(Rectangle((0, DRP.PP.ChannelTop/DRP.PP.pPmm), 
#                                    DRP.PP.RightCornerInflow/DRP.PP.pPmm, -dist/DRP.PP.pPmm, facecolor="grey"))

#    ax_inset.set_ylim((20,23))
#    ax_inset.set_xlim((8,10))


def plotTrajecotry(path, pPmm, geometryTJ, savepath = '', forPub=False):
    ''' Read results file from image analyisis and plit Trajectory'''
    name = TR.find_batchName(path); name = name.replace('\\', '')
    #centroid_x, centroid_y, v_x, v_y, area, width, d12, height, d12rect, centralWidth, centralHeight, curvature1, curvature2, curvature3, curvature4= readResultsFile(path)
    centroid_x, centroid_y, _, _, _, _, _, _, _, _, _, _, _, _, _ = TJ_ROF.readResultsFile(path)
    
    centreline = (geometryTJ[2] + geometryTJ[3])/2.0
    
    capD = 4.0*pPmm
#    itemindex = np.where(centroid_y > geometryTJ[1] + 2*capD and centroid_y < np.max(centroid_y) + capD)
    itemindex = [i for i,x in enumerate(centroid_y) if (x > geometryTJ[1] + 2*capD and x < np.max(centroid_y) + capD)]
    initial_x = centroid_x[itemindex]
    initial_y = centroid_y[itemindex]
    
    centring = np.mean(initial_x)
    if forPub:
        fig = plt.figure(figsize=(6, 6), dpi=200); ax = fig.add_subplot(111); FS=20; MS=8

    else:
        fig = plt.figure(figsize=(8, 6), dpi=200); ax = fig.add_subplot(111); FS=12; MS=6

    plt.plot(centroid_x/pPmm, centroid_y/pPmm, 'bs', markeredgecolor='b', markersize=MS, label='Centroid position')
    plt.plot(initial_x/pPmm, initial_y/pPmm, 'r<', markeredgecolor='r', markersize=MS, label='Centroid position used for initial offset')

    maxY = np.max(centroid_y)/pPmm                    
    plt.plot([centreline/pPmm, centreline/pPmm], [geometryTJ[1]/pPmm, maxY], 'g--', linewidth =3, label='Centreline')
    plt.plot([centring/pPmm, centring/pPmm], [geometryTJ[1]/pPmm, maxY], 'r-', linewidth =1, label='Initial Centring')   
    plt.xlabel('X [mm]', fontsize=FS); plt.ylabel('Y [mm]', fontsize=FS)
    
    
#    if not forPub:
#        plt.text(0.01, 0.8, 'Final Offset from Centreline = %.3f mm' %(self.distanceFromCentreline/self.PP.pPmm), fontsize = 16, horizontalalignment='left', verticalalignment='center', transform = ax.transAxes)

    from matplotlib.patches import Rectangle
    currentAxis = plt.gca()
    
    ymin, ymax = ax.get_ylim(); xmin, xmax = ax.get_xlim()
    
    XY=[]
    #top rectnalge
    x11 = 0 ; y11 = 0
    x12 = xmax*3 ; y12 = geometryTJ[0]/pPmm
    XY += [x11, y11, x12, y12]
    
    x21 = 0.0 ; y21 = geometryTJ[1]/pPmm
    x22 =  geometryTJ[2]/pPmm; y22 = ymax
    XY += [x21, y21, x22, y22]
    
    x31 = geometryTJ[3]/pPmm; y31 = geometryTJ[1]/pPmm
    x32 = xmax ; y32 = ymax
    XY += [x31, y31, x32, y32] ; XY = np.array(XY)
    

    currentAxis.add_patch(Rectangle((XY[0], XY[1]), XY[2], XY[3], facecolor="grey"))
    currentAxis.add_patch(Rectangle((XY[4], XY[5]), XY[6], XY[7], facecolor="grey"))
    currentAxis.add_patch(Rectangle((XY[8], XY[9]), XY[10], XY[11], facecolor="grey"))


    
    if xmax < geometryTJ[3]/pPmm: xmax = geometryTJ[3]/pPmm + 1.0
    plt.ylim(0.5, 30) ;#plt.ylim(0.5, ymax - 4)
    plt.xlim(6, 35.5); #plt.xlim(6, xmax)
    
    if not forPub:
        plt.title("Centroid position -  " +name, fontsize=12)
        plt.legend(loc='best', fontsize=FS/2, ncol=2)
#            else:
#                plt.legend(loc='best', fontsize=FS/2)
        
#    a = plt.axes([0.175, 0.3, .4, .4], axisbg='w')
#    plt.plot(centroid_x/pPmm, centroid_y/pPmm, 'bs', markeredgecolor='b');
#    currentAxis = plt.gca()
#    currentAxis.add_patch(Rectangle((XY[0], XY[1]), XY[2], XY[3], facecolor="grey"))
#    currentAxis.add_patch(Rectangle((XY[4], XY[5]), XY[6], XY[7], facecolor="grey"))
#    currentAxis.add_patch(Rectangle((XY[8], XY[9]), XY[10], XY[11], facecolor="grey"))
#    plt.xlim(24, 34)
#    plt.ylim(1, 6)
    fig.tight_layout()
    
    if savepath != '':
        if not forPub:
            plt.savefig((savepath+"CentroidPosition_" +name+".jpg"), dpi=300)
        else:
#            plt.gca().invert_yaxis()
            fig.set_size_inches(6.5, 6)
            fig.tight_layout()
            sn = ana._save_for_latex("CentroidPosition_" +name+"_forpub.jpg")
            plt.savefig(savepath+sn, dpi=400)

    
def plotMigrationVelocity(path, pPmm, FPS, geometryTJ, savepath = '', forPub=False, lim=[], anPos=[]):
    ''' Read results file from image analyisis and plit Trajectory'''
    name = TR.find_batchName(path); name = name.replace('\\', '')
    #centroid_x, centroid_y, v_x, v_y, area, width, d12, height, d12rect, centralWidth, centralHeight, curvature1, curvature2, curvature3, curvature4= readResultsFile(path)
    centroid_x, centroid_y,  v_x, v_y, _, width, _, _, _, _, _, _, _, _, _ = TJ_ROF.readResultsFile(path)
    vy = np.gradient(centroid_y, 1/FPS)
    vys = gen.runningMean(vy, 5)
    
    t = np.arange(0, len(centroid_x)); t = t/FPS
#    print(len(v_y)); print(len(vy))
    fig = plt.figure()
    plt.plot(t[:-1], -v_y, 'bs'); plt.plot(t, vys/pPmm, 'ro')
    print(len(t)); print(len(width))
    fig = plt.figure()
    plt.plot(t,width, 'ro')
    capD = 4.0*pPmm

    if forPub:
        fig = plt.figure(figsize=(6, 6), dpi=200); ax = fig.add_subplot(111); FS=20; MS=8
    else:
        fig = plt.figure(figsize=(8, 6), dpi=200); ax = fig.add_subplot(111); FS=12; MS=6

#    plt.plot(t, vy/pPmm, 'bs', markeredgecolor='k', markersize=MS, label='Centroid position')
    plt.plot(t, vys/pPmm, 'bo--', markeredgecolor='k', markersize=MS, label='Centroid position Smoothed')
    if len(anPos)==8:
        ax.annotate('Maximal Deformation', xy=(anPos[0], anPos[1]), xytext=(anPos[2], anPos[3]),
            arrowprops=dict(facecolor='black', shrink=0.05), fontsize=16)
        ax.annotate('Steady State reached', xy=(anPos[4], anPos[5]), xytext=(anPos[6], anPos[7]),
            arrowprops=dict(facecolor='black', shrink=0.05), fontsize=16)
    plt.xlabel('Time [s]', fontsize=FS); plt.ylabel('Migration Velocity [mm/s]', fontsize=FS)
    
    ymin, ymax = ax.get_ylim(); xmin, xmax = ax.get_xlim()    
    
    if len(lim) == 4:
        plt.xlim(lim[0], lim[1]); plt.ylim(lim[2], lim[3])
#    plt.xlim(6, 35.5); #plt.xlim(6, xmax)
    
    if not forPub:
        plt.title("Migration Velocity -  " +name, fontsize=12)
        plt.legend(loc='best', fontsize=FS/2, ncol=2)
    fig.tight_layout()
    
    if savepath != '':
        if not forPub:
            plt.savefig((savepath+"MigrationVelocity_%s.jpg" %(name)), dpi=300)
        else:
#            plt.gca().invert_yaxis()
            fig.set_size_inches(6.5, 6)
            fig.tight_layout()
            sn = ana._save_for_latex("MigrationVelocity_%s_forpub.jpg" %(name))
            plt.savefig(savepath+sn, dpi=400)
            
def plotVelocity(path, pPmm, FPS,  savepath = '', forPub=False, TJlim=[], lim=[], cutoff=0):
    ''' Read results file from image analyisis and plit Trajectory'''
    name = TR.find_batchName(path); name = name.replace('\\', '')
    #centroid_x, centroid_y, v_x, v_y, area, width, d12, height, d12rect, centralWidth, centralHeight, curvature1, curvature2, curvature3, curvature4= readResultsFile(path)
    centroid_x, centroid_y,  v_x, v_y, _, width, _, _, _, _, _, _, _, _, _ = TJ_ROF.readResultsFile(path)
    vy = np.gradient(centroid_y, 1.0/FPS)
    vx = np.gradient(centroid_x, 1.0/FPS)
    vmag = np.sqrt(vy*vy + vx*vx)

    isOutlier=  ana.is_outlier(vmag, thresh=10)
    for jj in range(20, len(vmag)):
        if isOutlier[jj]:
            dv = vmag[jj-1] -vmag[jj-2]
            vmag[jj] = vmag[jj-1] +dv
    vmag = gen.runningMean(vmag, 5)
#    vmag = np.sqrt(v_y*v_y*v_x*v_x)
    
    t = np.arange(0, len(centroid_x)); t = t/FPS

#    print(len(v_y)); print(len(vy))
    fig = plt.figure()
    plt.plot(t[:-1], -v_y, 'bs'); plt.plot(t, vy/pPmm, 'ro')
    plt.plot(t[:-1], -v_x, 'ys'); plt.plot(t, vx/pPmm, 'go')
    plt.ylim(-1, 1)
#    print(len(t)); print(len(width))
#    fig = plt.figure()
#    plt.plot(t,width, 'ro')
    if cutoff!= 0:    vmag = vmag[:-cutoff]; t = t[:-cutoff]
    if forPub:
        fig = plt.figure(figsize=(6, 6), dpi=200); ax = fig.add_subplot(111); FS=20; MS=8
    else:
        fig = plt.figure(figsize=(8, 6), dpi=200); ax = fig.add_subplot(111); FS=12; MS=6

#    plt.plot(t, vy/pPmm, 'bs', markeredgecolor='k', markersize=MS, label='Centroid position')
    plt.plot(t, vmag/pPmm, 'bo--', markeredgecolor='k', markersize=MS, label='Centroid position Smoothed')
    if len(TJlim) == 2:
        plt.axvspan(TJlim[0], TJlim[1], color='y', alpha=0.5)  
#        ymin, ymax = ax.get_ylim(); #xmin, xmax = ax.get_xlim()   
        plt.text(TJlim[0] + 0.1 * (TJlim[1]-TJlim[0]),0.95*(lim[3]- lim[2]), "T-Junction", fontsize=15)
        plt.text(lim[0]+0.1*(TJlim[0]-lim[0]) ,0.3*(lim[3]- lim[2]), "Main \nChannel", fontsize=15)
        plt.text(1.05*TJlim[1] ,0.3*(lim[3]- lim[2]), "Daughter \nChannel", fontsize=15)
#        plt.text(1,0.2, "T-Junction", fontsize=15)
        
    plt.xlabel('Time [s]', fontsize=FS); plt.ylabel('Velocity [mm/s]', fontsize=FS)
    
     
    
    if len(lim) == 4:
        plt.xlim(lim[0], lim[1]); plt.ylim(lim[2], lim[3])
#    plt.xlim(6, 35.5); #plt.xlim(6, xmax)
    
    if not forPub:
        plt.title("Velocity -  " +name, fontsize=12)
        plt.legend(loc='best', fontsize=FS/2, ncol=2)
    fig.tight_layout()
    
    if savepath != '':
        if not forPub:
            plt.savefig((savepath+"Velocity_%s.jpg" %(name)), dpi=300)
        else:
#            plt.gca().invert_yaxis()
            fig.set_size_inches(6.5, 6)
            fig.tight_layout()
            sn = ana._save_for_latex("Velocity_%s_forpub.jpg" %(name))
            plt.savefig(savepath+sn, dpi=400)
            
def smoothAndInterpolate(x):
    isOutlier=  ana.is_outlier(x, thresh=10)
    for jj in range(20, len(x)):
        if isOutlier[jj]:
            dv = x[jj-1] -x[jj-2]
            x[jj] = x[jj-1] +dv    
    return gen.runningMean(x, 5)
    
def plotAcceleration(path, pPmm, FPS,  savepath = '', forPub=False, TJlim=[], lim=[], cutoff=0):
    ''' Read results file from image analyisis and plit Trajectory'''
    name = TR.find_batchName(path); name = name.replace('\\', '')
    #centroid_x, centroid_y, v_x, v_y, area, width, d12, height, d12rect, centralWidth, centralHeight, curvature1, curvature2, curvature3, curvature4= readResultsFile(path)
    centroid_x, centroid_y,  v_x, v_y, _, width, _, _, _, _, _, _, _, _, _ = TJ_ROF.readResultsFile(path)
    vy = np.gradient(centroid_y, 1.0/FPS)
    vx = np.gradient(centroid_x, 1.0/FPS)

    vy = smoothAndInterpolate(vy); vx = smoothAndInterpolate(vx);
    ay = np.gradient(vy, 1.0/FPS)
    ax = np.gradient(vx, 1.0/FPS)
    amag = ay + ax

    isOutlier=  ana.is_outlier(amag, thresh=10)
    for jj in range(20, len(amag)):
        if isOutlier[jj]:
            dv = amag[jj-1] -amag[jj-2]
            amag[jj] = amag[jj-1] +dv
    amag = gen.runningMean(amag, 5)
#    vmag = np.sqrt(v_y*v_y*v_x*v_x)
    
    t = np.arange(0, len(centroid_x)); t = t/FPS

    if cutoff!= 0:    amag = amag[:-cutoff]; t = t[:-cutoff]
    if forPub:
        fig = plt.figure(figsize=(6, 6), dpi=200); ax = fig.add_subplot(111); FS=20; MS=8
    else:
        fig = plt.figure(figsize=(8, 6), dpi=200); ax = fig.add_subplot(111); FS=12; MS=6

#    plt.plot(t, vy/pPmm, 'bs', markeredgecolor='k', markersize=MS, label='Centroid position')
    plt.plot(t, amag/pPmm, 'ro--', markeredgecolor='k', markersize=MS, label='Acceleration Smoothed')
    if len(TJlim) == 2:
        plt.axvspan(TJlim[0], TJlim[1], color='y', alpha=0.5)  
#        ymin, ymax = ax.get_ylim(); #xmin, xmax = ax.get_xlim()   
        plt.text(TJlim[0] + 0.1 * (TJlim[1]-TJlim[0]),0.95*(lim[3]- lim[2]), "T-Junction", fontsize=15)
        plt.text(lim[0]+0.1*(TJlim[0]-lim[0]) ,0.3*(lim[3]- lim[2]), "Main \nChannel", fontsize=15)
        plt.text(1.05*TJlim[1] ,0.3*(lim[3]- lim[2]), "Daughter \nChannel", fontsize=15)
#        plt.text(1,0.2, "T-Junction", fontsize=15)
        
    plt.xlabel('Time [s]', fontsize=FS); plt.ylabel('Acceleration [mm/s${}^2$]', fontsize=FS)

    if len(lim) == 4:
        plt.xlim(lim[0], lim[1]); plt.ylim(lim[2], lim[3])
#    plt.xlim(6, 35.5); #plt.xlim(6, xmax)
    
    if not forPub:
        plt.title("Velocity -  " +name, fontsize=12)
        plt.legend(loc='best', fontsize=FS/2, ncol=2)
    fig.tight_layout()
    
    if savepath != '':
        if not forPub:
            plt.savefig((savepath+"Velocity_%s.jpg" %(name)), dpi=300)
        else:
#            plt.gca().invert_yaxis()
            fig.set_size_inches(6.5, 6)
            fig.tight_layout()
            sn = ana._save_for_latex("Velocity_%s_forpub.jpg" %(name))
            plt.savefig(savepath+sn, dpi=400)