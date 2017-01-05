# -*- coding: utf-8 -*-
"""
Created on Mon Jun 08 17:27:02 2015

@author: mbbxkeh2

Compare different Runs, i.e. deal with several results classes
"""

from __future__ import absolute_import, division, print_function
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import os


import analysis as ana
import semi_sphere as ss


#=========================
#Constant
SIGNATURE= 'Edgar Haener, ALT 2.105'
HFONT = {'fontname':'Helvetica'}
HFONT = {'fontname':'Comic Sans MS'}
ERRFORCE=1.0

PUBMS=10     #marker size for pulbications
PUBFS=18    #fontsize axis for publication

def plotQVsFinalOffset(instances, savePath, forPub=False, names=[], savename='_'):
        '''Plot final distance from centreline versus volume flux Q'''
        leng = len(instances)
        assert leng <= 10
        c =  ana.coloursForMarker(n=leng)    
        m =ana.markerSymbols()
        if forPub:
            fig = plt.figure(figsize=(6, 6), dpi=200,); FS=PUBFS; MS=PUBMS
        else:
            fig = plt.figure(figsize=(8, 6), dpi=200,); FS=10; MS=6
            
        ax = fig.add_subplot(111)
        
        for ii in range(leng):
            if len(names) != leng:
                nam = instances[ii].directoryName
            else:
                nam = names[ii]
                
            plt.errorbar(instances[ii].volumeFlux, 
                         instances[ii].finalDistCentre/instances[ii].pixelsPmm, 
                         #yerr=self.finalDistCentreSTD/self.pixelsPmm,
                        linestyle='None', marker = m[ii], color = c[ii],
                        markersize=4, label = nam)
                    
            plt.errorbar(instances[ii].volumeFluxes, 
                         instances[ii].aveFinalDistCentre/instances[ii].avePPmm, 
                         yerr=instances[ii].aveFinalDistCentreSTD/instances[ii].avePPmm, 
                         linestyle='None', marker = m[ii], markersize=8, 
                         color = c[ii], label = nam + ' Mean')
                         
        plt.legend(loc='best', fontsize=6, ncol=1)
        plt.title('Q vs Offset in Semi Cylinder Setup');    plt.xlabel('Volume Flux [$ml/min$]');   plt.ylabel('Final Offset [$mm$]')  
        ax = ana._standardPlotSize(ax)        
        x1,x2,y1,y2 = plt.axis()
        plt.axis((x1,x2,5,15))
        
        _addSchematic(ax)
        
        ana._saveFig(fig, os.path.join(savePath, 'RESLT') ,'QvsOffset%s.jpg' %savename)
        
def plotQVsFinalOffsetAve(instances, savePath, forPub=False, names=[], savename='_'):
        '''Plot final distance from centreline versus volume flux Q'''
        leng = len(instances)
        assert leng <= 10
        c =  ana.coloursForMarker(n=leng)    
        m =ana.markerSymbols()
        if forPub:
            FS=PUBFS; MS=PUBMS
        else:
            FS=10; MS=8
            
        fig = plt.figure(figsize=(8, 6), dpi=200,); ax = fig.add_subplot(111)
        
        for ii in range(leng):
            if len(names) != leng:
                nam = instances[ii].directoryName
            else:
                nam = names[ii]
                
#            plt.errorbar(instances[ii].volumeFlux, 
#                         instances[ii].finalDistCentre/instances[ii].pixelsPmm, 
#                         #yerr=self.finalDistCentreSTD/self.pixelsPmm,
#                        linestyle='None', marker = m[ii], color = c[ii],
#                        markersize=4, label = nam)
                    
            plt.errorbar(instances[ii].volumeFluxes, 
                         instances[ii].aveFinalDistCentre/instances[ii].avePPmm, 
                         yerr=instances[ii].aveFinalDistCentreSTD/instances[ii].avePPmm, 
                         linestyle='None', marker = m[ii], markersize=MS, ecolor='k',
                         color = c[ii], label = nam + ' Mean')
                         
        if not forPub:
            plt.legend(loc ='best', fontsize=6, ncol=1)
            plt.title('Q vs Offset in Semi Cylinder Setup');    
            _addSchematic(ax)

        ax = ana._standardPlotSize(ax)              
        plt.xlabel('Volume Flux [ml/min]', fontsize=FS);   plt.ylabel('Final Offset $L_{fo}$ [mm]', fontsize=FS)        
        x1,x2,y1,y2 = plt.axis()


        plt.axis((x1,x2,5,15))
        plt.show()             
        
        ana._saveFig(fig, os.path.join(savePath, 'RESLT') ,'AveQvsOffset%s.jpg' %savename)

def plotVelocitiesVsFinalOffsetAve(instances, savePath, forPub=False, names=[], savename='_'):
        '''Plot final distance from centreline versus volume flux Q'''
        leng = len(instances)
        assert leng <= 10
        c =  ana.coloursForMarker(n=leng)    
        m =ana.markerSymbols()
        if forPub:
            fig = plt.figure(figsize=(6, 6), dpi=200,); FS=PUBFS; MS=PUBMS
        else:
            fig = plt.figure(figsize=(8, 6), dpi=200,); FS=10; MS=6
            
        ax = fig.add_subplot(111)
        
        for ii in range(leng):
            if len(names) != leng:
                nam = instances[ii].directoryName
            else:
                nam = names[ii]
                
#            plt.errorbar(instances[ii].volumeFlux, 
#                         instances[ii].finalDistCentre/instances[ii].pixelsPmm, 
#                         #yerr=self.finalDistCentreSTD/self.pixelsPmm,
#                        linestyle='None', marker = m[ii], color = c[ii],
#                        markersize=4, label = nam)
                    
            plt.errorbar(instances[ii].aveInitialVelocities, 
                         instances[ii].aveFinalDistCentre/instances[ii].avePPmm, 
                         xerr = instances[ii].aveInitialVelocitiesSTD,
                         yerr = instances[ii].aveFinalDistCentreSTD/instances[ii].avePPmm, 
                         linestyle='None', marker = m[ii], markersize=8, ecolor='k',
                         color = c[ii], label = nam + ' Mean')
                         
        plt.legend(loc='best', fontsize=6, ncol=1)
        plt.title('Velocity vs Offset in Semi Cylinder Setup');    plt.xlabel('Velocity in initial channel [$mm/s$]');   plt.ylabel('Final Offset [$mm$]')  
        ax = ana._standardPlotSize(ax)        
        x1,x2,y1,y2 = plt.axis()

        if forPub:
            plt.axis((0,50,5,15))
        else:
            plt.axis((x1,x2,5,15))
            
        _addSchematic(ax)
        
        ana._saveFig(fig, os.path.join(savePath, 'RESLT') ,'AveVelocityvsOffset%s.jpg' %savename)
        
def plotQVsFinalOffset_specificQ(instances, savePath, VolFlux=20.0, forPub=False, names=[], savename='_'):
        '''Plot final distance from centreline versus volume flux Q'''
        leng = len(instances)
        assert leng <= 10
        c =  ana.coloursForMarker(n=leng)    
        m =ana.markerSymbols()
        if forPub:
            FS=PUBFS; MS=PUBMS
        else:
             FS=10; MS=6
        fig = plt.figure(figsize=(8, 6), dpi=200,); ax = fig.add_subplot(111)
        
        for ii in range(leng):
            if len(names) != leng:
                nam = instances[ii].directoryName
            else:
                nam = names[ii]
            
            q=[]; pos=[]                
            for kk in range(len(instances[ii].volumeFlux)):
                if instances[ii].volumeFlux[kk] == VolFlux:
                    q.append(instances[ii].volumeFlux[kk])
                    pos.append(instances[ii].finalDistCentre[kk]/instances[ii].pixelsPmm[kk])
                    
            plt.errorbar(q, 
                         pos,
                         linestyle='None', marker = m[ii], markersize=MS, 
                         color = c[ii], label = nam)
                         
        if not forPub:
            plt.legend(loc='best', fontsize=6, ncol=1)
            plt.title('Q = %.0f vs Offset in Half Cylinder Setup' %VolFlux);    
            _addSchematic(ax)
        ax = ana._standardPlotSize(ax)        
        plt.xlabel('Volume Flux = %.1f [$ml/min$]' %VolFlux);   plt.ylabel('Final Offset [$mm$]')  
        
        x1,x2,y1,y2 = plt.axis()
        plt.axis((x1,x2,5,15))
        
        ana._saveFig(fig, os.path.join(savePath, 'RESLT') ,'Q%.0f_0vsOffset%s.jpg' %(VolFlux,savename))
        
        fig = plt.figure(figsize=(8, 6), dpi=200,); ax = fig.add_subplot(111)
            
        ax = fig.add_subplot(111)
        x = range(leng)    
        plt.xticks(x, names, rotation=0, fontsize=FS)
        
        for ii in range(leng):
            if len(names) != leng:
                nam = instances[ii].directoryName
            else:
                nam = names[ii]
            
            pos=[]; err_pos=[]                
            for kk in range(len(instances[ii].volumeFlux)):
                if instances[ii].volumeFlux[kk] == VolFlux:
                    pos.append(instances[ii].finalDistCentre[kk]/instances[ii].pixelsPmm[kk])
                    err_pos.append(instances[ii].finalDistCentreSTD[kk]/instances[ii].pixelsPmm[kk])
                    
            xtemp = np.arange(len(pos))
#            print('len(xtemp) = %d' %len(xtemp))            
            xtemp.fill(x[ii])
#            print('len(pos) = %d' %len(pos))
#            print('len(xtemp) = %d' %len(xtemp))
            ax.errorbar(xtemp, 
                         pos, yerr = err_pos,
                         linestyle='None', marker = m[ii], markersize=MS, 
                         ecolor = 'k', color = c[ii], label = nam)
        plt.xlim((min(x)-0.5, max(x)+0.5))
#        plt.legend(loc='best', fontsize=6, ncol=1)
        if not forPub:
            plt.title('Q=%.0f in Half Cylinder Setup' %VolFlux);  
        
        plt.ylabel('Final Offset $L_{fo}$ [mm]', fontsize = FS)
        plt.xlabel('Capsule / Gel Bead', fontsize = FS)
#        ax.tick_params(labelsize=FS)
        
        ana._saveFig(fig, os.path.join(savePath, 'RESLT') ,'FinalOffset_fixedQ_%s.jpg' %savename)
        

        
def _plotQVsAngle(instances, savePath, forPub=False, names=[], savename='_'):
        '''Plot angle versus volume flux'''
        leng = len(instances)
        assert leng <= 10
        c =  ana.coloursForMarker(n=leng)    
        m =ana.markerSymbols()
        if forPub:
            fig = plt.figure(figsize=(6, 6), dpi=200,); FS=PUBFS; MS=PUBMS
        else:
            fig = plt.figure(figsize=(8, 6), dpi=200,); FS=10; MS=6
            
        ax = fig.add_subplot(111)
        
        for ii in range(leng):
            if len(names) != leng:
                nam = instances[ii].directoryName
            else:
                nam = names[ii]
                
            plt.errorbar(instances[ii].volumeFlux, 
                         instances[ii].angle, 
                         #yerr=self.finalDistCentreSTD/self.pixelsPmm,
                        linestyle='None', marker = m[ii], color = c[ii],
                        markersize=4, label = nam)
                    
            plt.errorbar(instances[ii].volumeFluxes, 
                         instances[ii].aveAngle, 
                         yerr=instances[ii].aveAngleSTD, 
                         linestyle='None', marker = m[ii], markersize=8, 
                         color = c[ii], label = nam + ' Mean')
                         
        plt.legend(loc='best', fontsize=6, ncol=1)
        plt.title('Q vs Angle in Semi Cylinder Setup');    plt.xlabel('Volume Flux [$ml/min$]');   plt.ylabel('Angle in Diffuser [$rad$]')  
        ax = ana._standardPlotSize(ax)  
        
        _addSchematic(ax)
#        x1,x2,y1,y2 = plt.axis()
#        plt.axis((x1,x2,0,y2))
        
        ana._saveFig(fig, os.path.join(savePath, 'RESLT') ,'QvsAngle%s.jpg' %savename)
        
def plotQVsAngleAve(instances, savePath, forPub=False, names=[], savename='_'):
        '''Plot angle versus volume flux'''
        leng = len(instances)
        assert leng <= 10
        c =  ana.coloursForMarker(n=leng)    
        m =ana.markerSymbols()
        if forPub:
            FS=PUBFS; MS=PUBMS
        else:
            FS=10; MS=8
            
        fig = plt.figure(figsize=(8, 6), dpi=200,); ax = fig.add_subplot(111)
        
        for ii in range(leng):
            if len(names) != leng:
                nam = instances[ii].directoryName
            else:
                nam = names[ii]
                    
            plt.errorbar(instances[ii].volumeFluxes, 
                         instances[ii].aveAngle* 180 / np.pi, 
                         yerr=instances[ii].aveAngleSTD* 180 / np.pi, 
                         linestyle='None', marker = m[ii], markersize=MS, ecolor='k',
                         color = c[ii], label = nam + ' Mean')
                         
        if not forPub:
            plt.legend(loc='best', fontsize=6, ncol=1)
            plt.title('Q vs Angle in Semi Cylinder Setup');    
            _addSchematic(ax)
        ax = ana._standardPlotSize(ax)             
        plt.xlabel('Volume Flux [ml/min]', fontsize=FS);   plt.ylabel(r'Angle $\theta$ [degrees]', fontsize=FS)  
        
        ana._saveFig(fig, os.path.join(savePath, 'RESLT') ,'AveQvsAngle%s.jpg' %savename)
        
def plotVelocityVsAngleAve(instances, savePath, forPub=False, names=[], savename='_'):
        '''Plot angle versus volume flux'''
        leng = len(instances)
        assert leng <= 10
        c =  ana.coloursForMarker(n=leng)    
        m =ana.markerSymbols()
        if forPub:
            fig = plt.figure(figsize=(6, 6), dpi=200,); FS=PUBFS; MS=PUBMS
        else:
            fig = plt.figure(figsize=(8, 6), dpi=200,); FS=10; MS=6
            
        ax = fig.add_subplot(111)
        
        for ii in range(leng):
            if len(names) != leng:
                nam = instances[ii].directoryName
            else:
                nam = names[ii]
                
#            plt.errorbar(instances[ii].volumeFlux, 
#                         instances[ii].angle, 
#                         #yerr=self.finalDistCentreSTD/self.pixelsPmm,
#                        linestyle='None', marker = m[ii], color = c[ii],
#                        markersize=4, label = nam)

            plt.errorbar(instances[ii].aveInitialVelocities, 
                         instances[ii].aveAngle* 180 / np.pi, 
                         xerr = instances[ii].aveInitialVelocitiesSTD,
                         yerr=instances[ii].aveAngleSTD* 180 / np.pi, 
                         linestyle='None', marker = m[ii], markersize=8, ecolor='k',
                         color = c[ii], label = nam + ' Mean')
                         
        plt.legend(loc='best', fontsize=6, ncol=1)
        plt.title('Velocity vs Angle in Semi Cylinder Setup');    plt.xlabel('Velocity [$mm/s$]');   plt.ylabel('Angle in Diffuser [degrees]')  
        ax = ana._standardPlotSize(ax)        
        
        if forPub:
            x1,x2,y1,y2 = plt.axis()
            plt.axis((0,50,0.14,0.36))
            
        _addSchematic(ax)
        
        ana._saveFig(fig, os.path.join(savePath, 'RESLT') ,'AveVelocityvsAngle%s.jpg' %savename)


def plotAveAngleVsInitialOffset(instances, savePath, forPub=False, names=[], savename=''):
    '''Plot initial offset from centreline in the channel leading to the 
    semi-cylinder versus angle in the diffuser'''     
    leng = len(instances)
    assert leng <= 10
    c =  ana.coloursForMarker(n=leng)    
    m =ana.markerSymbols()
    if forPub:
        fig = plt.figure(figsize=(6, 6), dpi=200,); FS=PUBFS; MS=PUBMS
    else:
        fig = plt.figure(figsize=(8, 6), dpi=200,); FS=10; MS=6
        
    ax = fig.add_subplot(111)
    
    for ii in range(leng):
        if len(names) != leng:
            nam = instances[ii].directoryName
        else:
            nam = names[ii]
            
        plt.errorbar(instances[ii].Offset/instances[ii].pixelsPmm, 
                     instances[ii].angle*180 / np.pi, 
#                         yerr=instances[ii].aveAngleSTD, 
                    linestyle='None', marker = m[ii], color = c[ii],
                    markersize=4, label = nam)
                
#        plt.errorbar(instances[ii].aveOffset/instances[ii].avePPmm, 
#                     instances[ii].aveFinalDistCentre/instances[ii].avePPmm, 
#                     yerr=instances[ii].aveOffsetSTD/instances[ii].avePPmm, 
#                     linestyle='None', marker = m[ii], markersize=8, 
#                     color = c[ii], label = nam + ' Mean')
                     
    plt.legend(loc='best', fontsize=6, ncol=1)
    plt.title('Initial centring vs Angle in Semi Cylinder Setup');    plt.xlabel('Initial Centring [$mm$]');   plt.ylabel('(Absolute) Angle in Diffuser [degrees]')  
    ax = ana._standardPlotSize(ax)        
    x1,x2,y1,y2 = plt.axis()
    plt.axis((x1,x2,5,y2))
    
    _addSchematic(ax)
    
    ana._saveFig(fig, os.path.join(savePath, 'RESLT') ,'InitialVsFinalPos%s.jpg' %savename)

        
def plotFinalVsInitialOffset(instances, savePath, forPub=False, names=[], savename='_'):
    '''Plot initial offset from centreline in the channel leading to the 
    semi-cylinder versus volume flux Q'''     
    leng = len(instances)
    assert leng <= 10
    c =  ana.coloursForMarker(n=leng)    
    m =ana.markerSymbols()
    if forPub:
        fig = plt.figure(figsize=(6, 6), dpi=200,); FS=PUBFS; MS=PUBMS
    else:
        fig = plt.figure(figsize=(8, 6), dpi=200,); FS=10; MS=6
        
    ax = fig.add_subplot(111)
    
    for ii in range(leng):
        if len(names) != leng:
            nam = instances[ii].directoryName
        else:
            nam = names[ii]
            
        plt.errorbar(instances[ii].Offset/instances[ii].pixelsPmm, 
                     instances[ii].finalDistCentre/instances[ii].pixelsPmm, 
                     #yerr=self.finalDistCentreSTD/self.pixelsPmm,
                    linestyle='None', marker = m[ii], color = c[ii],
                    markersize=4, label = nam)
                
        plt.errorbar(instances[ii].aveOffset/instances[ii].avePPmm, 
                     instances[ii].aveFinalDistCentre/instances[ii].avePPmm, 
                     yerr=instances[ii].aveOffsetSTD/instances[ii].avePPmm, 
                     linestyle='None', marker = m[ii], markersize=8, 
                     color = c[ii], label = nam + ' Mean')
                     
    plt.legend(loc='best', fontsize=6, ncol=1)
    plt.title('Initial centring vs Final Offset in Semi Cylinder Setup');    plt.xlabel('Initial Centring [$mm$]');   plt.ylabel('(Absolute) Final Offset [$mm$]')  
    ax = ana._standardPlotSize(ax)        
    x1,x2,y1,y2 = plt.axis()
    plt.axis((x1,x2,5,y2))
    
    _addSchematic(ax)
    
    ana._saveFig(fig, os.path.join(savePath, 'RESLT') ,'InitialVsFinalPos%s.jpg' %savename)
    
def plotFinalNormedVsInitialOffset(instances, force, savePath, forPub=False, names=[], savename='_'):
    '''Plot initial offset from centreline in the channel leading to the 
    semi-cylinder versus volume flux Q'''     
    leng = len(instances)
    assert leng <= 10
    c =  ana.coloursForMarker(n=leng)    
    m =ana.markerSymbols()
    if forPub:
        FS=PUBFS; MS=PUBMS
    else:
        FS=10; MS=6
        
    fig = plt.figure(figsize=(8, 6), dpi=200,);    ax = fig.add_subplot(111)
    
    std_initialCentring=[]; std_finalDistance=[];
    for ii in range(leng):
        if len(names) != leng:
            nam = instances[ii].directoryName
        else:
            nam = names[ii]
        
        initialCentring=[]; finalDist_normed=[];
        err_initialCentring=[]; err_finalDist_normed=[];
        for vft in instances[ii].volumeFlux:                
            tempx, tempy =_returnVolFlux(vft, instances[ii].volumeFlux, 
                                         instances[ii].Offset/instances[ii].pixelsPmm, 
                                         instances[ii].finalDistCentre/instances[ii].pixelsPmm)
            err_tempx, err_tempy =_returnVolFlux(vft, instances[ii].volumeFlux, 
                                         instances[ii].OffsetSTD/instances[ii].pixelsPmm, 
                                         instances[ii].finalDistCentreSTD/instances[ii].pixelsPmm)
            avefDC, acefDCSTD =_returnVolFlux(vft, instances[ii].volumeFluxes, 
                                         instances[ii].aveFinalDistCentre/instances[ii].avePPmm, 
                                         instances[ii].aveFinalDistCentreSTD/instances[ii].avePPmm)
            for kk in range(len(tempx)):
                initialCentring.append(tempx[kk])
                finalDist_normed.append(tempy[kk]/avefDC)

                err_initialCentring.append(err_tempx[kk])
                err_finalDist_normed.append(err_tempy[kk]/avefDC)
                                   
        
        aveinitialCentring = np.mean(initialCentring)
        aveinitialCentringSTD = np.std(initialCentring)
        avefinalDist_normed = np.mean(finalDist_normed)
        avefinalDist_normedSTD = np.std(finalDist_normed)
        
        std_initialCentring.append(aveinitialCentringSTD)
        std_finalDistance.append(avefinalDist_normedSTD)
        
        if ii!=0 and ii!=1:            
            plt.errorbar(initialCentring, 
                         finalDist_normed, 
#                         xerr= err_initialCentring,
#                         yerr = err_finalDist_normed,
                         linestyle='None', marker = m[ii], color = c[ii],
                        markersize=MS-2, label = nam)
            if not forPub:
                plt.errorbar(aveinitialCentring, 
                         avefinalDist_normed, 
                         yerr=avefinalDist_normedSTD, 
                         xerr=aveinitialCentringSTD,
                         linestyle='None', marker = m[ii], markersize=12, 
                         color = c[ii], label = nam + ' Mean')
    if not forPub:
        plt.legend(loc='best', fontsize=6, ncol=1)
        plt.title('Initial centring vs Final Offset in Semi Cylinder Setup');    
        _addSchematic(ax)        
    ax = ana._standardPlotSize(ax)            
    plt.xlabel('Initial Centring [mm]', fontsize=FS);   plt.ylabel('Normalized Final Offset', fontsize=FS)  

#    x1,x2,y1,y2 = plt.axis()
    plt.axis((-1.0,1.0,0.75,1.25))

    ana._saveFig(fig, os.path.join(savePath, 'RESLT') ,'InitialVsFinalNormedPos%s.jpg' %savename)
    
    # Plot Standard Deviation of Capsules
    if forPub:
        fig = plt.figure(figsize=(6, 6), dpi=200,); FS=PUBFS; MS=PUBMS
    else:
        fig = plt.figure(figsize=(8, 6), dpi=200,); FS=10; MS=6
        
    ax = fig.add_subplot(111)
    x = range(leng)    
    plt.xticks(x, names, rotation=45)
    ax.plot(x, std_initialCentring,linestyle='None', marker = 's', color='b', label='Standard Deviation of Initial Centring', markersize=8)
    ax.plot(x, std_finalDistance, linestyle='None',marker = 'o', color='r', label='Standard Deviation of Final Distance from Centreline', markersize=8)
    
    plt.legend(loc='best', fontsize=6, ncol=1)
    plt.title('Variation in Semi Cylinder Setup');  plt.ylabel('Standard Deviation')  
#    ax = ana._standardPlotSize(ax)            
#    _addSchematic(ax)
    
    ana._saveFig(fig, os.path.join(savePath, 'RESLT') ,'VariationSTD_%s.jpg' %savename)

def _returnVolFlux(tragetVF, vf, x, y):
    '''Returns all x at tragetVF'''
    tempx=[]; tempy=[]
    for i in range(len(vf)):
        if vf[i] == tragetVF:
            tempx.append(x[i])
            tempy.append(y[i])
    return  np.array(tempx), np.array(tempy)
    
def plotAbsFinalVsInitialOffset(instances, savePath, forPub=False):
    '''Plot initial offset from centreline in the channel leading to the 
    semi-cylinder versus volume flux Q'''      
    
    leng = len(instances)
    assert leng <= 10
    c =  ana.coloursForMarker(n=leng)    
    m =ana.markerSymbols()
    if forPub:
        FS=PUBFS; MS=PUBMS
    else:
        FS=10; MS=6
        
    fig = plt.figure(figsize=(8, 6), dpi=200,); ax = fig.add_subplot(111)
    
    for ii in range(leng):
        p1, p2, fres = instances[ii]._linearFit(np.abs(instances[ii].Offset/instances[ii].pixelsPmm),instances[ii].finalDistCentre/instances[ii].pixelsPmm)    

    
        plt.errorbar(np.abs(instances[ii].Offset/instances[ii].pixelsPmm), 
                     instances[ii].finalDistCentre/instances[ii].pixelsPmm,
                    #yerr=self.finalDistCentreSTD/self.pixelsPmm, 
                    #xerr=self.OffsetSTD/self.pixelsPmm,
                    markersize=MS,
                    linestyle='None', marker = m[ii], color = c[ii], 
                    label=instances[ii].directoryName)
                
        plt.plot([0.0, np.max(np.abs(instances[ii].Offset/instances[ii].pixelsPmm))], 
                         [instances[ii].func(0.0, p1,p2), instances[ii].func(np.max(np.abs(instances[ii].Offset/instances[ii].pixelsPmm)), p1,p2)], 
                          color = c[ii], 
                          label='%s: Linear Best Fit with m= %.4f and c = %.4f (fres = %.4f)' %(instances[ii].directoryName,p1, p2, fres))
    if not forPub:
        plt.legend(loc='best', fontsize=6, ncol=1)
        plt.title('Initial vs Final Offset')
        ax = ana._standardPlotSize(ax)        
    plt.ylabel('(Absolute) Final Offset [mm]', fontsize=FS);   plt.xlabel('(Absolute) Initial Centering [mm]', fontsize=FS)  
    
    ana._saveFig(fig, os.path.join(savePath, 'RESLT') ,'AbsInitiaVsFinalOffset.jpg')
    
    vfs=[]
    for ii in range(leng):
        vfs.extend(instances[ii].volumeFlux)       
    vfs = list(set(vfs)) # remove duplicate entries
    
   
    counter=0
    for vft in vfs:
        numbRuns=[]
        fig = plt.figure(figsize=(8, 6), dpi=200,); ax = fig.add_subplot(111)  
        for ii in range(leng):
            if vft in instances[ii].volumeFlux:
                
                tempx, tempy =_returnVolFlux(vft, instances[ii].volumeFlux, 
                                             instances[ii].Offset/instances[ii].pixelsPmm, 
                                             instances[ii].finalDistCentre/instances[ii].pixelsPmm)
                p1, p2, fres = instances[ii]._linearFit(np.abs(tempx),tempy)
                
                plt.errorbar(np.abs(tempx), tempy,
                            #yerr=self.finalDistCentreSTD/self.pixelsPmm, 
                            #xerr=self.OffsetSTD/self.pixelsPmm,
                            linestyle='None', marker = m[ii], color = c[ii], 
                            label=instances[ii].directoryName)
                            
                labelString = '%s: Linear Best Fit with m= %.4f and c = %.4f (RSS = %.4f)' %(instances[ii].directoryName, p1, p2, fres)
        
                plt.plot([0.0, np.max(np.abs(tempx))], [instances[ii].func(0.0, p1,p2), 
                          instances[ii].func(np.max(np.abs(tempx)), p1,p2)], 
                          color=c[ii], label=labelString)
                numbRuns.append(len(tempx))
        if not forPub:
            plt.legend(loc='best', fontsize=6, ncol=1)

            textstring = 'Number of Runs = '
            for d in numbRuns:
                textstring += '%d, ' %d
        
            plt.text(0.1, 0.6, textstring, fontsize = 16, 
                 horizontalalignment='left', verticalalignment='center', 
                 transform = ax.transAxes)
            plt.title('Initial vs Final Offset, Q = %.f' %(vft))
        ax = ana._standardPlotSize(ax)        
        plt.ylabel('(Absolute) Final Offset [mm]', fontsize=FS);   plt.xlabel('(Absolute) Initial Centering [mm]', fontsize=FS)  
                 
        
        ana._saveFig(fig, os.path.join(savePath, 'RESLT') ,'Q=%.f_AbsInitiaVsFinalOffset.jpg' %(vft))
        counter +=1

def plotInitialOffsetVsQ(instances, savePath, forPub=False, names=[], savename='_'):
    '''Plot initial offset from centreline in the channel leading to the 
    semi-cylinder versus volume flux Q'''      
    
    leng = len(instances)
    assert leng <= 10
    c =  ana.coloursForMarker(n=leng)    
    m =ana.markerSymbols()
    if forPub:
        fig = plt.figure(figsize=(6, 6), dpi=200,); FS=PUBFS; MS=PUBMS
    else:
        fig = plt.figure(figsize=(8, 6), dpi=200,); FS=10; MS=6
        
    ax = fig.add_subplot(111)
    
    for ii in range(leng):
        if len(names) != leng:
            nam = instances[ii].directoryName
        else:
            nam = names[ii]
#        p1, p2, fres = instances[ii]._linearFit(np.abs(instances[ii].Offset/instances[ii].pixelsPmm),instances[ii].finalDistCentre/instances[ii].pixelsPmm)        
        plt.errorbar(instances[ii].volumeFlux, 
                     instances[ii].Offset/instances[ii].pixelsPmm,
                    #yerr=self.finalDistCentreSTD/self.pixelsPmm, 
                    #xerr=self.OffsetSTD/self.pixelsPmm,
                    linestyle='None', marker = m[ii], color = c[ii], 
                    label=nam)
                
#        plt.plot([0.0, np.max(np.abs(instances[ii].Offset/instances[ii].pixelsPmm))], 
#                         [instances[ii].func(0.0, p1,p2), instances[ii].func(np.max(np.abs(instances[ii].Offset/instances[ii].pixelsPmm)), p1,p2)], 
#                          color = c[ii], 
#                          label='%s: Linear Best Fit with m= %.4f and c = %.4f (fres = %.4f)' %(instances[ii].directoryName,p1, p2, fres))
    
    plt.legend(loc='best', fontsize=6, ncol=1)
    plt.title('Initial vs Final Offset')
    plt.ylabel('Initial Offset [$mm$]');   plt.xlabel('Volume Flux [$ml/min$]')  
    
#    textstring = 'Number of Runs = %d' %(len(self.Offset))
#    plt.text(0.1, 0.1, textstring, fontsize = 12, 
#                 horizontalalignment='left', verticalalignment='center', 
#                 transform = ax.transAxes)
    ax = ana._standardPlotSize(ax)  

    _addSchematic(ax)
    
    ana._saveFig(fig, os.path.join(savePath, 'RESLT') ,'InitalCentringQ%s.jpg' %savename)
    
def plotAveInitialOffsetVsQ(instances, savePath, forPub=False, names=[], savename='_'):
    '''Plot initial offset from centreline in the channel leading to the 
    semi-cylinder versus volume flux Q'''      
    
    leng = len(instances)
    assert leng <= 10
    c =  ana.coloursForMarker(n=leng)    
    m =ana.markerSymbols()
    if forPub:
        fig = plt.figure(figsize=(6, 6), dpi=200,); FS=PUBFS; MS=PUBMS
    else:
        fig = plt.figure(figsize=(8, 6), dpi=200,); FS=10; MS=6
        
    ax = fig.add_subplot(111)
    
    for ii in range(leng):
        if len(names) != leng:
            nam = instances[ii].directoryName
        else:
            nam = names[ii]
#        p1, p2, fres = instances[ii]._linearFit(np.abs(instances[ii].Offset/instances[ii].pixelsPmm),instances[ii].finalDistCentre/instances[ii].pixelsPmm)        
        plt.errorbar(instances[ii].volumeFluxes, 
                     instances[ii].aveOffset/instances[ii].avePPmm,
                    yerr=instances[ii].aveOffsetSTD/instances[ii].avePPmm, 
                    #xerr=self.OffsetSTD/self.pixelsPmm,
                    linestyle='None', marker = m[ii], color = c[ii], 
                    label=nam)
                
#        plt.plot([0.0, np.max(np.abs(instances[ii].Offset/instances[ii].pixelsPmm))], 
#                         [instances[ii].func(0.0, p1,p2), instances[ii].func(np.max(np.abs(instances[ii].Offset/instances[ii].pixelsPmm)), p1,p2)], 
#                          color = c[ii], 
#                          label='%s: Linear Best Fit with m= %.4f and c = %.4f (fres = %.4f)' %(instances[ii].directoryName,p1, p2, fres))
    
    plt.legend(loc='best', fontsize=6, ncol=1)
    plt.title('Initial Centring')
    plt.ylabel('Initial Offset [$mm$]');   plt.xlabel('Volume Flux [$ml/min$]')  
    
#    textstring = 'Number of Runs = %d' %(len(self.Offset))
#    plt.text(0.1, 0.1, textstring, fontsize = 12, 
#                 horizontalalignment='left', verticalalignment='center', 
#                 transform = ax.transAxes)
    ax = ana._standardPlotSize(ax)  

    _addSchematic(ax)
    
    ana._saveFig(fig, os.path.join(savePath, 'RESLT') ,'AveInitalCentringQ%s.jpg' %savename)
    
def QtoCa(Q, F):
    '''Convert volume flux to pseud capillary number for the semi-cylinder setup
    
    Q   -   Array of volumn fluxes
    F   -   Measure of force
    '''
    #5'000 cSt Silicon Oil
    nu = 0.005 #m^2 / s
    rho = 997.0 #kg/m^3
    mu = nu/rho    
    
    #for the area, what cross - section to use? Best the channel, as that is where the deformation happens
    area = 4e-3 * 16e-3
    
    distance = 4e-3 #diameter or approximate capsules diamter
    
    Ca = mu * Q * distance / (area*F)
    
    return Ca
    
def plotCaVsAngleAve(instances, force, savePath, forPub=False, names=[], savename='_'):
    '''Plot angle versus volume flux'''
    leng = len(instances)
    assert leng <= 10
    c =  ana.coloursForMarker(n=leng)    
    m =ana.markerSymbols()
    if forPub:
        FS=PUBFS; MS=PUBMS
    else:
        FS=10; MS=8
        
    fig = plt.figure(figsize=(8, 6), dpi=200,); ax = fig.add_subplot(111)
    
    for ii in range(leng):
        if len(names) != leng:
            nam = instances[ii].directoryName
        else:
            nam = names[ii]
            
#            plt.errorbar(instances[ii].volumeFlux, 
#                         instances[ii].angle, 
#                         #yerr=self.finalDistCentreSTD/self.pixelsPmm,
#                        linestyle='None', marker = m[ii], color = c[ii],
#                        markersize=4, label = nam)
        ca = QtoCa(instances[ii].volumeFluxes, force[ii])
        plt.errorbar(ca, 
                     instances[ii].aveAngle * 180 / np.pi, 
                     yerr=instances[ii].aveAngleSTD* 180 / np.pi, 
                     linestyle='None', marker = m[ii], markersize=MS, ecolor="k", 
                     color = c[ii], label = nam + ' Mean')
    if not forPub:
        plt.legend(loc='best', fontsize=6, ncol=1)
        plt.title('Ca vs Angle in Semi Cylinder Setup');    
        _addSchematic(ax)    
    ax = ana._standardPlotSize(ax)   
    plt.xlabel('Pseudo-Capillary Number', fontsize=FS);   plt.ylabel(r'Angle $\theta$ [degrees]', fontsize=FS)  
     
    
    ana._saveFig(fig, os.path.join(savePath, 'RESLT') ,'AveCavsAngle%s.jpg' %savename)
    
def plotCaVsFinalOffsetAve(instances, force, savePath, forPub=False, names=[], savename='_'):
    '''Plot final distance from centreline versus volume flux Q'''
    leng = len(instances)
    assert leng <= 10
    c =  ana.coloursForMarker(n=leng)    
    m =ana.markerSymbols()
    if forPub:
        FS=PUBFS; MS=PUBMS
    else:
        FS=10; MS=8
        
    fig = plt.figure(figsize=(8, 6), dpi=200,); ax = fig.add_subplot(111)
    
    for ii in range(leng):
        if len(names) != leng:
            nam = instances[ii].directoryName
        else:
            nam = names[ii]
            
#            plt.errorbar(instances[ii].volumeFlux, 
#                         instances[ii].finalDistCentre/instances[ii].pixelsPmm, 
#                         #yerr=self.finalDistCentreSTD/self.pixelsPmm,
#                        linestyle='None', marker = m[ii], color = c[ii],
#                        markersize=4, label = nam)
        ca = QtoCa(instances[ii].volumeFluxes, force[ii])        
        plt.errorbar(ca, 
                     instances[ii].aveFinalDistCentre/instances[ii].avePPmm, 
                     yerr=instances[ii].aveFinalDistCentreSTD/instances[ii].avePPmm, 
                     linestyle='None', marker = m[ii], markersize=MS, ecolor="k", 
                     color = c[ii], label = nam + ' Mean')
    if not forPub:                     
        plt.legend(loc='best', fontsize=6, ncol=1)
        plt.title('Ca vs Offset in Semi Cylinder Setup');    
        _addSchematic(ax)    
#    print('plotCaVsFinalOffsetAve: forPub = %r \t FS = %d' %(forPub, FS))
    ax = ana._standardPlotSize(ax)         
    plt.xlabel('Pseudo-Capillary Number', fontsize=FS);   plt.ylabel(r'Final Offset $L_{fo}$ [mm]', fontsize=FS)  
       
    x1,x2,y1,y2 = plt.axis()
    plt.axis((x1,x2,5,15)) 
    
    ana._saveFig(fig, os.path.join(savePath, 'RESLT') ,'AveCavsOffset%s.jpg' %savename)
    
def signature(ax):
    ax.text(1.01, -0.05, SIGNATURE, fontsize = 3, horizontalalignment='left', verticalalignment='center', transform = ax.transAxes, **HFONT)
    
        
def _addSchematic(ax):
    '''Add schematic in top right'''

    from matplotlib.offsetbox import OffsetImage, AnnotationBbox
    
    
    fn= 'M:\\EdgarHaener\\Images\\Semi-Sphere.png'
    
    from matplotlib._png import read_png
    arr_lena = read_png(fn)

    imagebox = OffsetImage(arr_lena, zoom=0.2)

    xy = (0.3, 0.85)
    
    ab = AnnotationBbox(imagebox, xy,
                        xybox=(0, 0),
                        xycoords=ax.transAxes,
                        boxcoords="offset points")
#                        bbox=dict(facecolor='none', edgecolor='red'))
#                        pad=0.5,
#                        arrowprops=dict(arrowstyle="->",
#                                        connectionstyle="angle,angleA=0,angleB=90,rad=3")
#                        )

    ax.add_artist(ab)
    
def plotCaVsFinalOffsetAveSelected(instances, force, savePath, initialCentringLimit=0.3, forPub=False, names=[], savename='_'):
    '''Plot final distance from centreline versus volume flux Q'''
    leng = len(instances)
    assert leng <= 10
    c =  ana.coloursForMarker(n=leng)    
    m =ana.markerSymbols()
    if forPub:
        fig = plt.figure(figsize=(6, 6), dpi=200,); FS=PUBFS; MS=PUBMS
    else:
        fig = plt.figure(figsize=(8, 6), dpi=200,); FS=10; MS=6
        
    ax = fig.add_subplot(111)
    
    nrPloted=0
    for ii in range(leng):
        if len(names) != leng:
            nam = instances[ii].directoryName
        else:
            nam = names[ii]
            
#            plt.errorbar(instances[ii].volumeFlux, 
#                         instances[ii].finalDistCentre/instances[ii].pixelsPmm, 
#                         #yerr=self.finalDistCentreSTD/self.pixelsPmm,
#                        linestyle='None', marker = m[ii], color = c[ii],
#                        markersize=4, label = nam)
        
        vf=[]; fd =[]; fdSTD=[]   ; offs=[]     
        
        #select values that are within initialCentringLimit        
        for jj in range(len(instances[ii].Offset)):
            
            if abs(instances[ii].Offset[jj]/instances[ii].pixelsPmm[jj]) < initialCentringLimit:
                vf.append(instances[ii].volumeFlux[jj])
                fd.append(instances[ii].finalDistCentre[jj]/instances[ii].pixelsPmm[jj])
                fdSTD.append(instances[ii].finalDistCentreSTD[jj]/instances[ii].pixelsPmm[jj])
                offs.append(instances[ii].Offset[jj]/instances[ii].pixelsPmm[jj])
            else:
                offs.append('x')
        nrPloted += len(fd)                
#        print("instances[ii].Offset/instances[ii].pixelsPmm = ", end='')
#        print(instances[ii].Offset/instances[ii].pixelsPmm)
#        print("offs = ", end='')
#        print(offs)
        
        vfs, aveFD, aveFDstd  = instances[ii].averageBy(vf, fd)
        print("aveFD = ", end='')
        print(aveFD)
        ca = QtoCa(vfs, force[ii]) 
#        print(cap)
#        print('%s \t len(cap) = %d' %(nam, len(cap)))
        if len(aveFD) != 0:
            plt.errorbar(ca, aveFD, yerr=aveFDstd, 
                         linestyle='None', marker = m[ii], markersize=8, 
                         color = c[ii], label = nam + ' Selected')

    textstring = 'Number of Runs = %d' %(nrPloted)
    plt.text(0.1, 0.1, textstring, fontsize = 12, 
                 horizontalalignment='left', verticalalignment='center', 
                 transform = ax.transAxes)
    
    plt.legend(loc='best', fontsize=6, ncol=1)
    plt.title('Ca vs Offset in Semi Cylinder Setup - Initial Centring withing %.2f mm' %(initialCentringLimit));    plt.xlabel('Pseudo-Capillary Number');   plt.ylabel('Final Offset [$mm$]')  
    ax = ana._standardPlotSize(ax)        
    x1,x2,y1,y2 = plt.axis()

    if forPub:
        plt.axis((0,50,5,15))
    else:
        plt.axis((x1,x2,5,15))
        
    _addSchematic(ax)
    ana._saveFig(fig, os.path.join(savePath, 'RESLT') ,'AveCavsOffset_Selected_%s.jpg' %savename)



def plotNormedFinalVsInitialOffset(instances, savePath, forPub=False, names=[], savename='_'):
    '''Plot initial offset from centreline in the channel leading to the 
    semi-sphere versus volume flux Q'''
    leng = len(instances)
    assert leng <= 10
    c =  ana.coloursForMarker(n=leng)    
    m =ana.markerSymbols()
    if forPub:
        FS=PUBFS; MS=PUBMS
    else:
        FS=10; MS=6

    fig = plt.figure(figsize=(8, 6), dpi=200,); ax = fig.add_subplot(111)
    
    for ii in range(leng):   
        if len(names) != leng:
            nam = instances[ii].directoryName
        else:
            nam = names[ii]
                                         
        initialCentring=[]; finalDist_normed=[];
        for vft in instances[ii].volumeFlux:                
            tempx, tempy =instances[ii]._returnVolFlux(vft, instances[ii].volumeFlux, 
                                         instances[ii].Offset/instances[ii].pixelsPmm, 
                                         instances[ii].finalDistCentre/instances[ii].pixelsPmm)
            avefDC, acefDCSTD =instances[ii]._returnVolFlux(vft, instances[ii].volumeFluxes, 
                                         instances[ii].aveFinalDistCentre/instances[ii].avePPmm, 
                                         instances[ii].aveFinalDistCentreSTD/instances[ii].avePPmm)
            for kk in range(len(tempx)):
                initialCentring.append(tempx[kk])
                finalDist_normed.append(tempy[kk]/avefDC)
                
        plt.errorbar(initialCentring, 
             finalDist_normed, 
             #yerr=self.finalDistCentreSTD/self.pixelsPmm,
            linestyle='None', marker = m[ii], color = c[ii],
            markersize=MS, label=nam)
    if not forPub:            
        plt.legend(loc='best', fontsize=6, ncol=1)
        plt.title('Initial vs Normalized Final Offset')
        _addSchematic(ax)        
    ax = ana._standardPlotSize(ax)        
    plt.ylabel('Normalized Final Offset', fontsize=FS);   plt.xlabel('Initial Centring [mm]', fontsize=FS)  
       
    plt.ylim((0.8, 1.4))
    plt.show()
#    x1,x2,y1,y2 = plt.axis()

    ana._saveFig(fig, os.path.join(savePath, 'RESLT') ,'InitialVsNormedFinalOffset_%s.jpg' %savename)
    
def plotQvsMeanHeightAtSS(instances, savePath, forPub=False, names=[], savename='_'):
    '''Plot initial offset from centreline in the channel leading to the 
    semi-sphere versus volume flux Q'''
    leng = len(instances)
    assert leng <= 10
    c =  ana.coloursForMarker(n=leng)    
    m =ana.markerSymbols()
    if forPub:
        fig = plt.figure(figsize=(8, 6), dpi=200,); FS=PUBFS; MS=PUBMS
    else:
        fig = plt.figure(figsize=(8, 6), dpi=200,); FS=10; MS=6
        
    ax = fig.add_subplot(111)
    
    for ii in range(leng):   
        if len(names) != leng:
            nam = instances[ii].directoryName
        else:
            nam = names[ii]
                                                        
#        plt.errorbar(instances[ii].volumeFluxes, 
#             instances[ii].aveMeanHeightEndSS/instances[ii].avePPmm, 
#             yerr=instances[ii].aveMeanHeightEndSSSTD/instances[ii].avePPmm,
#            linestyle='None', marker = m[ii], color = c[ii],
#            markersize=7, label=names[ii] + ' Contour')

        plt.errorbar(instances[ii].volumeFluxes, 
             instances[ii].aveMean_yd_EndSS/instances[ii].avePPmm, 
             yerr=instances[ii].aveMean_yd_EndSSSTD/instances[ii].avePPmm,
            linestyle='None', marker = m[ii], color = c[ii],
            markersize=MS, label=names[ii])
    if not forPub:
        plt.legend(loc='best', fontsize=6, ncol=2)
        plt.title('Height vs Q')
    plt.ylabel('Deformed Height $L_h$ [mm]', fontsize=PUBFS);   plt.xlabel('Volume Flux [ml/min]', fontsize=PUBFS)  
   
    plt.ylim((2, 4))

    
    if not forPub:
        _addSchematic(ax)
        ax = ana._standardPlotSize(ax)        
    ana._saveFig(fig, os.path.join(savePath, 'RESLT') ,'QVsHeight_%s.jpg' %savename)


    
        