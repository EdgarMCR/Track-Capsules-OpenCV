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

def plotQVsFinalOffset(instances, savePath, forPub=False, names=[], savename='', TJ=False):
        '''Plot final distance from centreline versus volume flux Q'''
        leng = len(instances)
        assert leng <= 10
        c =  ana.coloursForMarker(n=leng*2)    
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
                    
            plt.errorbar(instances[ii].volumeFlux2, 
                         instances[ii].finalDistCentre/instances[ii].pixelsPmm, 
                         #yerr=self.finalDistCentreSTD/self.pixelsPmm,
                        linestyle='None', marker = m[ii], color = c[ii],
                        markersize=4, label = nam)
                    
            plt.errorbar(instances[ii].volumeFluxes2, 
                         instances[ii].aveFinalDistCentre/instances[ii].avePPmm, 
                         yerr=instances[ii].aveFinalDistCentreSTD/instances[ii].avePPmm, 
                         linestyle='None', marker = m[leng+ii], markersize=8, 
                         color = c[leng+ii], label = nam + ' Mean')
                     
        plt.legend(loc='best', fontsize=6, ncol=1)
        plt.title('Q vs Offset in Pinching %s' %(savename.replace('_', ' ')));    plt.xlabel('Pinching Flow Volume Flux [$ml/min$]');   plt.ylabel('Final Offset [$mm$]')  
        ax = ana._standardPlotSize(ax)        
        x1,x2,y1,y2 = plt.axis()
        plt.axis((x1,x2,0,y2))
        
        signature(ax)
        _addSchematic(ax, TJ=TJ)
        
        ana._saveFig(fig, os.path.join(savePath, 'RESLT') ,'QvsOffset_Pinching_%s.jpg' %(savename))

def plotQVsGapAndD12(instances, savePath, forPub=False, names=[], savename='', TJ=False):
        '''Plot final distance from centreline versus volume flux Q'''
        leng = len(instances)
        assert leng <= 10
        c =  ana.coloursForMarker(n=leng*2)    
        m =ana.markerSymbols()
        fig = plt.figure(figsize=(8, 6), dpi=200,);
        if forPub:
            FS=PUBFS; MS=PUBMS
        else:
             FS=12; MS=8
            
        ax = fig.add_subplot(111)
        
        for ii in range(leng):
            if len(names) != leng:
                nam = instances[ii].directoryName
            else:
                nam = names[ii]
                    
            plt.errorbar(instances[ii].volumeFluxes2, 
                         instances[ii].aveGap/instances[ii].avePPmm, 
                         yerr=instances[ii].aveGapSTD/instances[ii].avePPmm, 
                         linestyle='None', marker = m[ii], markersize=MS, 
                         color = c[ii], label = nam + ' Gap')
        
        plt.xlabel('Pinching Volume Flux $Q$ [ml/min]');   plt.ylabel('Gap Capsule - Wall [mm]')  
        if not forPub:
            ax = ana._standardPlotSize(ax)        
            plt.legend(loc='best', fontsize=6, ncol=1)
            plt.title('Gap $D_{12}$ Vs $L_{fo}$ %s' %(savename.replace('_', ' ')));            
            signature(ax)
            _addSchematic(ax, TJ=TJ)      
        
        ax2 = ax.twinx()
        for ii in range(leng):                   
            if len(names) != leng:
                nam = instances[ii].directoryName
            else:
                nam = names[ii]
            ax2.errorbar(instances[ii].volumeFluxes2, 
                         instances[ii].aveGap_d12, 
                         yerr=instances[ii].aveGap_d12STD, 
                         linestyle='None', marker = m[leng+ii], markersize=MS, 
                         color = c[leng+ii], label = nam + ' $D_{12}$')
        plt.legend(loc='center right', fontsize=6, ncol=1)
        plt.ylabel('Gap Deformation $D_{12}$');
        ax2.set_ylim((0.0, 0.55))
        ax.set_ylim((-0.2, 1.2))
#        x1,x2,y1,y2 = plt.axis()
#        plt.axis((x1,x2,0,y2))


        if not forPub:
            ana._saveFig(fig, os.path.join(savePath, 'RESLT') ,'GapD12vsQ_Pinching_%s.jpg' %(savename))
        else:
            ana._saveFig(fig, os.path.join(savePath, 'RESLT') ,'GapD12vsQ_Pinching_forPub_%s.jpg' %(savename))  
            
def plotCaVsGapAndQvsD12(instances, force, savePath, forPub=False, names=[], savename=''):
        '''Plot final distance from centreline versus volume flux Q'''
        leng = len(instances)
        assert leng <= 10
        c =  ana.coloursForMarker(n=leng)    
        m =ana.markerSymbols()
        
        
        
        
        fig = plt.figure(figsize=(4, 3.6), dpi=200,);        ax = fig.add_subplot(111)
        if forPub:
            FS=PUBFS; MS=PUBMS
        else:
             FS=12; MS=8
            
     
        for ii in range(leng):
            if len(names) != leng:
                nam = instances[ii].directoryName
            else:
                nam = names[ii]
                
            ca = QtoCa(instances[ii].volumeFluxes2, force[ii])
            plt.errorbar(ca, 
                         instances[ii].aveGap/instances[ii].avePPmm, 
                         yerr=instances[ii].aveGapSTD/instances[ii].avePPmm, 
                         linestyle='None', marker = m[ii], markersize=MS, 
                         color = c[ii], label = nam + ' $L_{gap}$')
        
        plt.xlabel('Pseudo-$Ca$');   plt.ylabel('Gap Width $L_{gap}$ [mm]')  
        if not forPub:
            ax = ana._standardPlotSize(ax)        
            plt.legend(loc='best', fontsize=6, ncol=1)
            plt.title(' $L_{gap}$ Vs $Ca$ - %s' %(savename.replace('_', ' ')));            
            signature(ax)
            _addSchematic(ax)      
        
#        x1,x2,y1,y2 = plt.axis()
#        plt.axis((x1,x2,0,y2))

        if not forPub:
            ana._saveFig(fig, os.path.join(savePath, 'RESLT') ,'D12vsCa_Pinching_%s.jpg' %(savename))
        else:
            ana._saveFig(fig, os.path.join(savePath, 'RESLT') ,'D12vsCa_Pinching_forPub_%s.jpg' %(savename)) 

        fig = plt.figure(figsize=(4, 3.5), dpi=200,);        ax = fig.add_subplot(111)
        if forPub:
            FS=PUBFS; MS=PUBMS
        else:
             FS=12; MS=8
        
        for ii in range(leng):                   
            if len(names) != leng:
                nam = instances[ii].directoryName
            else:
                nam = names[ii]
            ca = QtoCa(instances[ii].volumeFluxes2, force[ii])
            plt.errorbar(ca,
                         instances[ii].aveGap_d12, 
                         yerr=instances[ii].aveGap_d12STD, 
                         linestyle='None', marker = m[ii], markersize=MS, 
                         color = c[ii], label = nam + ' $D_{12}$')
  
        plt.xlabel('Pseudo-$Ca$');   plt.ylabel('Taylor Deformation Parameter $D_{12}$')  
        if not forPub:
            ax = ana._standardPlotSize(ax)        
            plt.legend(loc='best', fontsize=6, ncol=1)
            plt.title(' $D_{12}$ Vs $Ca$ - %s' %(savename.replace('_', ' ')));            
            signature(ax)
            _addSchematic(ax)      
        
#        x1,x2,y1,y2 = plt.axis()
#        plt.axis((x1,x2,0,y2))

        if not forPub:
            ana._saveFig(fig, os.path.join(savePath, 'RESLT') ,'GapvsCa_Pinching_%s.jpg' %(savename))
        else:
            ana._saveFig(fig, os.path.join(savePath, 'RESLT') ,'GapvsCa_Pinching_forPub_%s.jpg' %(savename))             



def plotQVsGapAndQvsD12(instances, savePath, forPub=False, names=[], savename='', TJ=False, top=False, force=[]):
        '''Plot final distance from centreline versus volume flux Q'''
        leng = len(instances)
        assert leng <= 10
        c =  ana.coloursForMarker(n=leng)    
        m =ana.markerSymbols()
        
        fig = plt.figure(figsize=(4, 3.6), dpi=200,);        ax = fig.add_subplot(111)
        if forPub:
            FS=PUBFS; MS=PUBMS
        else:
             FS=12; MS=8
            
     
        for ii in range(leng):
            if len(names) != leng:
                nam = instances[ii].directoryName
            else:
                nam = names[ii]
                
            if force:
                x= QtoCa(instances[ii].volumeFluxes2, force[ii])
            else:
                x= instances[ii].volumeFluxes2
            if top:
                plt.errorbar(x, 
                             instances[ii].aveTopGap/instances[ii].avePPmm, 
                             yerr=instances[ii].aveTopGapSTD/instances[ii].avePPmm, 
                             linestyle='None', marker = m[ii], markersize=MS, 
                             color = c[ii], label = nam + ' $L_{gap}$')
            else:
                plt.errorbar(x, 
                             instances[ii].aveGap/instances[ii].avePPmm, 
                             yerr=instances[ii].aveGapSTD/instances[ii].avePPmm, 
                             linestyle='None', marker = m[ii], markersize=MS, 
                             color = c[ii], label = nam + ' $L_{gap}$')
        
        if force: plt.xlabel('Pinching Pseudo $Ca$ [ml/min]'); qs = 'Ca'
        else: plt.xlabel('Pinching Volume Flux $Q_{pin}$ [ml/min]')  ; qs = 'Q'
        plt.ylabel('Gap Width $L_{gap}$ [mm]')  
        if not forPub:
            ax = ana._standardPlotSize(ax)        
            plt.legend(loc='best', fontsize=6, ncol=1)
            plt.title(' $L_{gap}$ Vs Q - %s' %(savename.replace('_', ' ')));            
            signature(ax)
            _addSchematic(ax, TJ=TJ)      
        
#        x1,x2,y1,y2 = plt.axis()
#        plt.axis((x1,x2,0,y2))

        if not forPub:
            ana._saveFig(fig, os.path.join(savePath, 'RESLT') ,'D12vs%s_Pinching_%s.jpg' %(qs,savename))
        else:
            ana._saveFig(fig, os.path.join(savePath, 'RESLT') ,'D12vs%s_Pinching_forPub_%s.jpg' %(qs, savename)) 

        fig = plt.figure(figsize=(4, 3.5), dpi=200,);        ax = fig.add_subplot(111)
        if forPub:
            FS=PUBFS; MS=PUBMS
        else:
             FS=12; MS=8
        
        for ii in range(leng):                   
            if len(names) != leng:
                nam = instances[ii].directoryName
            else:
                nam = names[ii]
            if force:
                x= QtoCa(instances[ii].volumeFluxes2, force[ii])
            else:
                x= instances[ii].volumeFluxes2
            plt.errorbar(x, 
                         instances[ii].aveGap_d12, 
                         yerr=instances[ii].aveGap_d12STD, 
                         linestyle='None', marker = m[ii], markersize=MS, 
                         color = c[ii], label = nam + ' $D_{12}$')
  
        if force: plt.xlabel('Pinching Pseudo $Ca$ [ml/min]')
        else: plt.xlabel('Pinching Volume Flux $Q_{pin}$ [ml/min]')
            
        plt.ylabel('Taylor Deformation Parameter $D_{12}$')  
        if not forPub:
            ax = ana._standardPlotSize(ax)        
            plt.legend(loc='best', fontsize=6, ncol=1)
            plt.title(' $D_{12}$ Vs Q - %s' %(savename.replace('_', ' ')));            
            signature(ax)
            _addSchematic(ax, TJ=TJ)      
        
        x1,x2,y1,y2 = plt.axis()
        plt.axis((x1,x2,0,y2))

        if not forPub:
            ana._saveFig(fig, os.path.join(savePath, 'RESLT') ,'Gapvs%s_Pinching_%s.jpg' %(qs,savename))
        else:
            ana._saveFig(fig, os.path.join(savePath, 'RESLT') ,'Gapvs%s_Pinching_forPub_%s.jpg' %(qs, savename)) 
  
def plotGapD12VsFinalOffsetAve(instances, savePath, forPub=False, names=[], savename='', pnb=0):
        '''Plot final distance from centreline versus volume flux Q'''
        leng = len(instances)
        assert leng <= 10
        c =  ana.coloursForMarker(n=leng+pnb)    
        m =ana.markerSymbols()
        fig = plt.figure(figsize=(6, 5), dpi=200,);
        if forPub:
            FS=PUBFS; MS=PUBMS
        else:
             FS=12; MS=8
            
        ax = fig.add_subplot(111)
        
        for ii in range(leng):
            if len(names) != leng:
                nam = instances[ii].directoryName
            else:
                nam = names[ii]
                    
            plt.errorbar(instances[ii].aveGap_d12, 
                         instances[ii].aveFinalDistCentre/instances[ii].avePPmm, 
                         xerr= instances[ii].aveGap_d12STD,
                         yerr=instances[ii].aveFinalDistCentreSTD/instances[ii].avePPmm, 
                         linestyle='None', marker = m[ii+pnb], markersize=MS, 
                         color = c[ii+pnb], label = nam + ' Mean')
        
        plt.xlabel('Gap Deformation $D_{12}$');   plt.ylabel('Final Distance from Centreline $L_{fo}$ [$mm$]')  
#        x1,x2,y1,y2 = plt.axis()
#        plt.axis((x1,x2,0,y2))
        if not forPub:
            ax = ana._standardPlotSize(ax)        
            plt.legend(loc='best', fontsize=6, ncol=1)
            plt.title('Gap $D_{12}$ Vs $L_{fo}$ %s' %(savename.replace('_', ' ')));            
            signature(ax)
            _addSchematic(ax)      
            ana._saveFig(fig, os.path.join(savePath, 'RESLT') ,'GapD12vsOffset_Pinching_%s.jpg' %(savename))
        else:
            ana._saveFig(fig, os.path.join(savePath, 'RESLT') ,'GapD12vsOffset_Pinching_forPub_%s.jpg' %(savename))       
            
def plotTopGapVsGapD12Ave(instances, savePath, forPub=False, names=[], savename='', pnb=0):
        '''Plot final distance from centreline versus volume flux Q'''
        leng = len(instances)
        assert leng <= 10
        c =  ana.coloursForMarker(n=leng+pnb)    
        m =ana.markerSymbols()
        fig = plt.figure(figsize=(6, 5), dpi=200,);
        if forPub:
            FS=PUBFS; MS=PUBMS
        else:
             FS=12; MS=8
            
        ax = fig.add_subplot(111)
        
        for ii in range(leng):
            if len(names) != leng:
                nam = instances[ii].directoryName
            else:
                nam = names[ii]
                    
            plt.errorbar(instances[ii].aveGap_d12, 
                         instances[ii].aveTopGap/instances[ii].avePPmm, 
                         xerr= instances[ii].aveGap_d12STD,
                         yerr=instances[ii].aveTopGapSTD/instances[ii].avePPmm, 
                         linestyle='None', marker = m[ii+pnb], markersize=MS, 
                         color = c[ii+pnb], label = nam + ' Mean')
        
        plt.xlabel('Gap Deformation $D_{12}$');   plt.ylabel('Gap $L_{gap}$ [mm]')  
#        x1,x2,y1,y2 = plt.axis()
#        plt.axis((x1,x2,0,y2))
        if not forPub:
            ax = ana._standardPlotSize(ax)        
            plt.legend(loc='best', fontsize=6, ncol=1)
            plt.title('Gap $D_{12}$ Vs $L_{gap}$ %s' %(savename.replace('_', ' ')));            
            signature(ax)
            _addSchematic(ax)      
            ana._saveFig(fig, os.path.join(savePath, 'RESLT') ,'GapVsD12_Pinching_%s.jpg' %(savename))
        else:
            ana._saveFig(fig, os.path.join(savePath, 'RESLT') ,'GapVsD12_Pinching_forPub_%s.jpg' %(savename))   

        
def plotQVsFinalOffsetAve(instances, savePath, forPub=False, names=[], savename='', TJ=False):
        '''Plot final distance from centreline versus volume flux Q'''
        leng = len(instances)
        assert leng <= 10
        c =  ana.coloursForMarker(n=leng)    
        m =ana.markerSymbols()
        fig = plt.figure(figsize=(8, 6), dpi=200,);
        if forPub:
            FS=PUBFS; MS=PUBMS
        else:
             FS=12; MS=8
            
        ax = fig.add_subplot(111)
        
        
        
        for ii in range(leng):
            if len(names) != leng:
                nam = instances[ii].directoryName
            else:
                nam = names[ii]
                    
#            plt.errorbar(instances[ii].volumeFlux2, 
#                         instances[ii].finalDistCentre/instances[ii].pixelsPmm, 
#                         #yerr=self.finalDistCentreSTD/self.pixelsPmm,
#                        linestyle='None', marker = m[ii], color = c[ii],
#                        markersize=4, label = nam)
                    
            plt.errorbar(instances[ii].volumeFluxes2, 
                         instances[ii].aveFinalDistCentre/instances[ii].avePPmm, 
                         yerr=instances[ii].aveFinalDistCentreSTD/instances[ii].avePPmm, 
                         linestyle='None', marker = m[ii], markersize=MS, 
                         color = c[ii], label = nam + ' Mean')
                     
        plt.xlabel('Pinching Flow Volume Flux $Q_{pin}$ [ml/min]', fontsize=FS);   plt.ylabel('Final Offset $L_{fo}$ [mm]', fontsize=FS)        
        x1,x2,y1,y2 = plt.axis()
        plt.axis((0,x2,0,y2))
        
        if not forPub:
            ax = ana._standardPlotSize(ax)        
            plt.legend(loc='best', fontsize=6, ncol=1)
            plt.title('Q vs Mean Offset in Pinching %s' %(savename.replace('_', ' ')));            
            signature(ax)
            _addSchematic(ax, TJ=TJ)      
            ana._saveFig(fig, os.path.join(savePath, 'RESLT') ,'QvsMeanOffset_Pinching_%s.jpg' %(savename))
        else:
            ana._saveFig(fig, os.path.join(savePath, 'RESLT') ,'QvsMeanOffset_Pinching_forPub_%s.jpg' %(savename))        
        

def _returnVolFlux(tragetVF, vf, x, y):
    '''Returns all x at tragetVF'''
    tempx=[]; tempy=[]
    for i in range(len(vf)):
        if vf[i] == tragetVF:
            tempx.append(x[i])
            tempy.append(y[i])
    return  np.array(tempx), np.array(tempy)
    
def plotAbsFinalVsInitialOffset(instances, savePath, forPub=False, TJ=False):
    '''Plot initial offset from centreline in the channel leading to the 
    semi-sphere versus volume flux Q'''     
    
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
#        p1, p2, fres = instances[ii]._linearFit(np.abs(instances[ii].Offset/instances[ii].pixelsPmm),instances[ii].finalDistCentre/instances[ii].pixelsPmm)    

    
        plt.errorbar(np.abs(instances[ii].Offset/instances[ii].pixelsPmm), 
                     instances[ii].finalDistCentre/instances[ii].pixelsPmm,
                    #yerr=self.finalDistCentreSTD/self.pixelsPmm, 
                    #xerr=self.OffsetSTD/self.pixelsPmm,
                    linestyle='None', marker = m[ii], color = c[ii], 
                    label=instances[ii].directoryName)
                
#        plt.plot([0.0, np.max(np.abs(instances[ii].Offset/instances[ii].pixelsPmm))], 
#                         [instances[ii].func(0.0, p1,p2), instances[ii].func(np.max(np.abs(instances[ii].Offset/instances[ii].pixelsPmm)), p1,p2)], 
#                          color = c[ii], 
#                          label='%s: Linear Best Fit with m= %.4f and c = %.4f (fres = %.4f)' %(instances[ii].directoryName,p1, p2, fres))
    
    plt.legend(loc='best', fontsize=6, ncol=1)
    plt.title('Initial vs Final Offset')
    plt.ylabel('(Absolute) Final Offset [$mm$]');   plt.xlabel('(Absolute) Initial Centering [$mm$]')  
    
#    textstring = 'Number of Runs = %d' %(len(self.Offset))
#    plt.text(0.1, 0.1, textstring, fontsize = 12, 
#                 horizontalalignment='left', verticalalignment='center', 
#                 transform = ax.transAxes)
    ax = ana._standardPlotSize(ax)        
    
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
#                p1, p2, fres = instances[ii]._linearFit(np.abs(tempx),tempy)
                
                plt.errorbar(np.abs(tempx), tempy,
                            #yerr=self.finalDistCentreSTD/self.pixelsPmm, 
                            #xerr=self.OffsetSTD/self.pixelsPmm,
                            linestyle='None', marker = m[ii], color = c[ii], 
                            label=instances[ii].directoryName)
                            
#                labelString = '%s: Linear Best Fit with m= %.4f and c = %.4f (RSS = %.4f)' %(instances[ii].directoryName, p1, p2, fres)
        
#                plt.plot([0.0, np.max(np.abs(tempx))], [instances[ii].func(0.0, p1,p2), 
#                          instances[ii].func(np.max(np.abs(tempx)), p1,p2)], 
#                          color=c[ii], label=labelString)
                numbRuns.append(len(tempx))
        
        plt.legend(loc='best', fontsize=6, ncol=1)

        textstring = 'Number of Runs = '
        for d in numbRuns:
            textstring += '%d, ' %d
        
        plt.text(0.1, 0.6, textstring, fontsize = 16, 
                 horizontalalignment='left', verticalalignment='center', 
                 transform = ax.transAxes)
        plt.title('Initial vs Final Offset, Q = %.f' %(vft))
        plt.ylabel('(Absolute) Final Offset [$mm$]');   plt.xlabel('(Absolute) Initial Centering [$mm$]')  
        ax = ana._standardPlotSize(ax)      
        
        signature(ax)
        _addSchematic(ax, TJ=TJ)
        
        ana._saveFig(fig, os.path.join(savePath, 'RESLT') ,'Q=%.f_AbsInitiaVsFinalOffset.jpg' %(vft))
        counter +=1

def plotInitialOffsetVsQ(instances, savePath, forPub=False, TJ=False):
    '''Plot initial offset from centreline in the channel leading to the 
    semi-sphere versus volume flux Q'''      
    
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
#        p1, p2, fres = instances[ii]._linearFit(np.abs(instances[ii].Offset/instances[ii].pixelsPmm),instances[ii].finalDistCentre/instances[ii].pixelsPmm)        
        plt.errorbar(np.abs(instances[ii].volumeFlux), 
                     instances[ii].Offset/instances[ii].pixelsPmm,
                    #yerr=self.finalDistCentreSTD/self.pixelsPmm, 
                    #xerr=self.OffsetSTD/self.pixelsPmm,
                    linestyle='None', marker = m[ii], color = c[ii], 
                    label=instances[ii].directoryName)
                
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

    signature(ax)
    _addSchematic(ax, TJ=TJ)    
    
    ana._saveFig(fig, os.path.join(savePath, 'RESLT') ,'InitiaVsQ.jpg')
    
    
def plotQVsGapSpaceing(instances, savePath, forPub=False, names=[], TJ=False):
        '''Plot final distance from centreline versus volume flux Q'''
        leng = len(instances)
        assert leng <= 10
        c =  ana.coloursForMarker(n=2.0*leng)    
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
                    
#            plt.errorbar(instances[ii].volumeFlux2, instances[ii].gapSpace/instances[ii].pixelsPmm, 
#                         #yerr=self.finalDistCentreSTD/self.pixelsPmm,
#                        linestyle='None', marker = m[ii], color = c[ii], 
#                        label = nam)
                        
            plt.errorbar(instances[ii].volumeFluxes2,  instances[ii].aveGap/instances[ii].avePPmm, 
                         yerr=instances[ii].aveGapSTD/instances[ii].avePPmm, 
                         linestyle='None', marker = m[ii], markersize=8, color = c[ii],
                         label = nam + ' Bottom Gap')
                         
            plt.errorbar(instances[ii].volumeFluxes2,  instances[ii].aveTopGap/instances[ii].avePPmm, 
                         yerr=instances[ii].aveTopGapSTD/instances[ii].avePPmm, 
                         linestyle='None', marker = m[ii+leng], markersize=8, color = c[ii],
                         label = nam + ' Top Gap')

        plt.title('Q vs Gap Spacing'); plt.xlabel('Pinching Flow Volume Flux [$ml/min$]');   plt.ylabel('Distance betweeen Capsule and Wall [$mm$]') 
                     
        plt.legend(loc='best', fontsize=6, ncol=1)
        
        ax = ana._standardPlotSize(ax)        
        x1,x2,y1,y2 = plt.axis()
        plt.axis((x1,x2,0,y2))
        
        signature(ax)
        _addSchematic(ax, TJ=TJ)
        
        ana._saveFig(fig, os.path.join(savePath, 'RESLT') ,'QvsWallGap_Pinching.jpg')

  
        
def QtoCa(Q, F):
    '''Convert volume flux to pseud capillary number
    
    Q   -   Array of volumn fluxes
    F   -   Measure of force
    '''
    nu = 0.005 #m^2 / s
    rho = 997.0 #kg/m^3
    mu = nu/rho    
    
    #for the area, what cross - section to use? Best the channel, as that is where the deformation happens
    area = 4e-3 * 4e-3
    
    distance = 4e-3 #diameter or approximate capsules diamter
    
    Ca = mu * Q * distance / (area*F)
    
    return Ca
    
    
def plotCaVsFinalOffsetAve(instances, force, savePath, forPub=False, names=[], savename='', TJ=False):
        '''Plot final distance from centreline versus volume flux pseud Ca'''
        leng = len(instances)
        assert leng <= 10
        c =  ana.coloursForMarker(n=leng)    
        m =ana.markerSymbols()
        fig = plt.figure(figsize=(8, 6), dpi=200,);
        if forPub:
            FS=PUBFS; MS=PUBMS
        else:
             FS=12; MS=8
            
        ax = fig.add_subplot(111)
        
        
        
        for ii in range(leng):
            if len(names) != leng:
                nam = instances[ii].directoryName
            else:
                nam = names[ii]
                    
#            plt.errorbar(instances[ii].volumeFlux2, 
#                         instances[ii].finalDistCentre/instances[ii].pixelsPmm, 
#                         #yerr=self.finalDistCentreSTD/self.pixelsPmm,
#                        linestyle='None', marker = m[ii], color = c[ii],
#                        markersize=4, label = nam)
            ca = QtoCa(instances[ii].volumeFluxes2, force[ii])
#            ca = instances[ii].volumeFluxes2/force[ii]         
            plt.errorbar(ca, 
                         instances[ii].aveFinalDistCentre/instances[ii].avePPmm, 
                         yerr=instances[ii].aveFinalDistCentreSTD/instances[ii].avePPmm, 
                         linestyle='None', marker = m[ii], markersize=MS, 
                         color = c[ii], label = nam + ' Mean')
                     
        plt.xlabel('Pinching Flow Pseudo-Cappilary Number', fontsize=FS);   plt.ylabel('Final Offset $L_{fo}$ [mm]', fontsize=FS)  

        x1,x2,y1,y2 = plt.axis()
        plt.axis((x1,x2,0,y2))
        if not forPub:
            ax = ana._standardPlotSize(ax)                
            signature(ax)
            _addSchematic(ax, TJ=TJ)
            plt.legend(loc='best', fontsize=6, ncol=1)
            plt.title('Pseud-Ca vs Mean Offset in Pinching %s' %(savename.replace('_', ' ')));            
            ana._saveFig(fig, os.path.join(savePath, 'RESLT') ,'CavsMeanOffset_Pinching_%s.jpg' %(savename))
        else:
            ana._saveFig(fig, os.path.join(savePath, 'RESLT') ,'CavsMeanOffset_Pinching_forPub_%s.jpg' %(savename))
            
            
            
def fit_exponentiol_LfoVsCa(instances, force, savePath, forPub=False, names=[], savename='', TJ=False):
        '''Plot final distance from centreline versus volume flux pseud Ca'''
        leng = len(instances)
        assert leng <= 10
        c =  ana.coloursForMarker(n=leng)    
        m =ana.markerSymbols()
        fig = plt.figure(figsize=(8, 6), dpi=200,);
        if forPub:
            FS=PUBFS; MS=PUBMS
        else:
             FS=12; MS=8
            
        ax = fig.add_subplot(111)
        
        
        xCa=[]; yLfo=[]; eyLfo=[]
        for ii in range(leng):
            if len(names) != leng:
                nam = instances[ii].directoryName
            else:
                nam = names[ii]
                    
            ca = QtoCa(instances[ii].volumeFluxes2, force[ii])
#            ca = instances[ii].volumeFluxes2/force[ii]         
            plt.errorbar(ca, 
                         instances[ii].aveFinalDistCentre/instances[ii].avePPmm, 
                         yerr=instances[ii].aveFinalDistCentreSTD/instances[ii].avePPmm, 
                         linestyle='None', marker = m[ii], markersize=MS, 
                         color = c[ii], label = nam + ' Mean')
#            print(len(ca))
#            xCa = xCa + ca.squeeze()
            xCa.extend(ca)
            yLfo.extend(instances[ii].aveFinalDistCentre/instances[ii].avePPmm)
            eyLfo.extend(instances[ii].aveFinalDistCentreSTD/instances[ii].avePPmm)
            
#        print("xCa:")
#        print(xCa)
#        print("yLfo:")
#        print(yLfo)
        
        for ii in range(len(xCa)-1, -1, -1):
            if np.isnan(xCa[ii]) or np.isnan(yLfo[ii]):
                del xCa[ii]
                del yLfo[ii]
        
        
        print(len(yLfo))
        print(len(xCa))
        from scipy.optimize import curve_fit
        popt, pcov = curve_fit(ana.fitFuncExp, xCa, yLfo,p0=(8.1, 300.0, 1.6))
        a = popt[0]; b=popt[1]; c=popt[2]
        print("a/b/c = %f, %f, %f \tlen(xCa) = %d" %(a,b,c, len(xCa)))

#        fity=[]        
#        for ca in xCa:
#            fity.append(ana.fitFuncExp(ca,a, b, c))
#        print(len(fity))
        cp = 8.1; bp=300; ap = 1.6  # I don't know hwy the fit suddenly changed! Here the old values that five a better fit
        residuals = np.array(yLfo) -  ana.fitFuncExp(np.array(xCa),a, b, c)
        fres = np.sum(np.power(residuals,2))
        ss_tot = np.sum(np.power((yLfo-np.mean(yLfo)), 2))
        r_squared = 1.0 - (fres / ss_tot)
        
        error = [] 
        for i in range(len(popt)):
            try:
                error.append(np.absolute(pcov[i][i])**0.5)
            except:
                error.append( 0.00 )
        perr_curvefit = np.array(error)
        print("Best match for a, b, c = %f +/- %f, %f +/- %f, %f +/- %f with R^2 = %f. \n(Using y = c * np.exp(-b*x) + a)" %(a, perr_curvefit[0], b,perr_curvefit[1],  c, perr_curvefit[2], r_squared))
        
        x = np.arange(0.0, np.max(xCa), 1e-5)
        y = ana.fitFuncExp(x,ap, bp, cp)
        plt.plot(x, y, 'k--', label='Exponential fit y = %f * np.exp(-%f * x) + %f' %(a,b,c))
        plt.xlabel('Pinching Flow Pseudo-Cappilary Number', fontsize=FS);   plt.ylabel('Final Offset $L_{fo}$ [mm]', fontsize=FS)  

        x1,x2,y1,y2 = plt.axis()
        plt.axis((x1,x2,0,y2))
        if not forPub:
            ax = ana._standardPlotSize(ax)                
            signature(ax)
            _addSchematic(ax, TJ=TJ)
            plt.legend(loc='best', fontsize=6, ncol=1)
            plt.title('Pseud-Ca vs Mean Offset in Pinching %s' %(savename.replace('_', ' ')));            
            ana._saveFig(fig, os.path.join(savePath, 'RESLT') ,'CavsMeanOffset_Fitted_Pinching_%s.jpg' %(savename))
        else:
            ana._saveFig(fig, os.path.join(savePath, 'RESLT') ,'CavsMeanOffset_Fitted_Pinching_forPub_%s.jpg' %(savename))
            
    
def signature(ax):
    ax.text(1.01, -0.05, SIGNATURE, fontsize = 3, horizontalalignment='left', verticalalignment='center', transform = ax.transAxes, **HFONT)


def _addSchematic(ax, TJ=False):
    '''Add schematic in top right'''

    from matplotlib.offsetbox import OffsetImage, AnnotationBbox
    
    if TJ:
        fn = 'M:\\EdgarHaener\\Images\\T-JunctionMode.png'
    else:
        fn= 'M:\\EdgarHaener\\Images\\PinchingMode.png'
    
    from matplotlib._png import read_png
    arr_lena = read_png(fn)

    imagebox = OffsetImage(arr_lena, zoom=0.175)

    xy = (0.8, 0.2)
    
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
        
        