# -*- coding: utf-8 -*-
"""
Created on Sun Dec 15 00:48:38 2019

@author: wpiper
"""

import numpy as np
import matplotlib.pyplot as plt
#import matplotlib.patches as mpatches
from scipy.signal import butter, filtfilt, medfilt
from scipy.stats import linregress
from scipy.optimize import curve_fit
from matplotlib.backends.backend_pdf import PdfPages

def PlotSaveAll(filelist=[], pdffileoutput='', denoise=True, denoiseHz = 2, baselinemethod='expdecay', motioncontrol=True):
    objectList = []
    headerstring = 'denoise='+str(denoise)+' , denoiseHz='+str(denoiseHz)+' , baselinemethod='+str(baselinemethod)+' , motioncontrol='+str(motioncontrol)
    if type(filelist)==str: filelist=[filelist]
    with PdfPages(pdffileoutput) as pp:
        firstPage = plt.figure(figsize=(11.69,8.27))
        firstPage.clf()
        txt = 'Fiber Photometry data - WPiper'
        firstPage.text(0.5,0.5,txt, transform=firstPage.transFigure, size=24, ha="center")
        firstPage.text(0.15,0.35,headerstring, transform=firstPage.transFigure, size=10)        
        pp.savefig()
        plt.close()
        for file in filelist:
            X = FullTrace(csvfile=file, trim='end')
            X.preprocess(figsave=True,denoise = denoise, denoiseHz = denoiseHz, baselinemethod=baselinemethod, motioncontrol=motioncontrol)
            figX = X.plot(headerstring=headerstring)
            pp.savefig(figX)
            if 'PTC_walter_285' in file:
                for i in [1,3,5,7,9]: 
                    CohenFig = CohensD_prepoststim(X,i,duration=30)
                    pp.savefig(CohenFig)
            if 'LTMUUS_walter_285' in file:
                for i in [1,2,3,4,5,6]: 
                    CohenFig = CohensD_prepoststim(X,i,duration=30)
                    pp.savefig(CohenFig)
            objectList.append(X)
    return(objectList)
    
        

class FullTrace:
    def __init__(self,csvfile,trim='None',ch1=1,ch2=3):
        try: Protpos = csvfile.index('LTMUUS')
        except: 
            try: Protpos = csvfile.index('PTC')
            except:
                try: Protpos = csvfile.index('U-US')
                except: Protpos=-28
        self.ProtType = csvfile[Protpos:Protpos+4] #Infer protocol type from first 4 characters of file name
        IDpos = csvfile.index('_28')
        self.AnimalID = csvfile[IDpos+1:IDpos+7] #Pull animal ID from file name.
        self.Time, self.GreenCh, self.RedCh, self.EventsCh = np.genfromtxt(csvfile, delimiter=',', skip_header=2, usecols=(0,ch1,ch2,5), unpack=True)
#change columns 1 and 3 back to 2 and 4 when done with this test
        TimeSlope = self.Time[1:] - self.Time[:-1]
        self.samplingrate = 1/np.mean(TimeSlope)
        self.samplingrateplusminus = [(1/np.min(TimeSlope))-self.samplingrate, 1/np.max(TimeSlope)-self.samplingrate]   
        print('Sampling Rate: '+str(self.samplingrate)+' Hz +/-'+str(self.samplingrateplusminus)+' Hz')
        EventSlope = self.EventsCh[1:] - self.EventsCh[:-1]
        EventIndexArray = np.argwhere(EventSlope)+1
        #print(EventIndexArray)
        self.EventList = []
        if EventSlope.max() > 0:
            self.EventList.append(['Whole Protocol', [int(EventIndexArray[0]),int(EventIndexArray[-1])]])
        #print(self.EventList)
        
        if self.ProtType == 'PTC_': # For December2019 animals
            counter = 0
            countbool = False
            for t in range(1,EventIndexArray.shape[0],1):
                timegap = self.Time[EventIndexArray[t]] - self.Time[EventIndexArray[t-1]]
                #print(str(timegap)+' timegap')
                if 27.5 < timegap < 28.5:
                    counter += 1
                    str_label = 'CS'+str(counter)
                    self.EventList.append([str_label, [int(EventIndexArray[t-1]),int(EventIndexArray[t+2])]])
                elif 0.95 < timegap < 1.05 and countbool==True:
                    str_label = 'US'+str(counter)
                    self.EventList.append([str_label, [int(EventIndexArray[t]),int(EventIndexArray[t+1])]])
                    countbool=False
                elif 0.95 < timegap < 1.05 and countbool==False:     
                    print(self.ProtType)
                    countbool=True
                elif 2.8 < timegap < 3.2 and countbool==True:
                    self.EventList[0][1][1] = int(EventIndexArray[t])
                elif 2.8 < timegap < 3.2 and countbool==False:
                    self.EventList[0][1][0] = int(EventIndexArray[t-1])
                    self.EventList[0][0] = 'Protocol whole'
                    countbool=True
        elif self.ProtType == 'LTMU': # For December2019 animals
            counter = 0 
            countbool=False
            for t in range(1,EventIndexArray.shape[0],1):
                timegap = self.Time[EventIndexArray[t]] - self.Time[EventIndexArray[t-1]]
                #print(str(timegap)+' timegap')
                if 29.5 < timegap < 30.5:
                    counter += 1
                    str_label = 'CS'+str(counter)
                    self.EventList.append([str_label, [int(EventIndexArray[t-1]),int(EventIndexArray[t])]])
                elif 4.5 < timegap < 5.5:
                    str_label = 'U-US'
                    self.EventList.append([str_label, [int(EventIndexArray[t-1]),int(EventIndexArray[t])]])
                elif 3.8 < timegap < 4.2 and countbool==True:
                    self.EventList[0][1][1] = int(EventIndexArray[t])
                elif 3.8 < timegap < 4.2 and countbool==False:
                    self.EventList[0][1][0] = int(EventIndexArray[t-1])
                    self.EventList[0][0] = 'Protocol whole'
                    countbool=True
        elif self.ProtType == 'U-US': # For October2019 animals
            counter = False
            for t in range(1,EventIndexArray.shape[0],1):
                timegap = self.Time[EventIndexArray[t]] - self.Time[EventIndexArray[t-1]]
                if 6.93 < timegap < 7.07 and counter == False:
                    self.EventList[0][0] = 'Protocol whole'
                    self.EventList[0][1][0] = int(EventIndexArray[t-1])
                    counter = True
                elif 4.95 < timegap < 5.05:
                    self.EventList.append(['U-US',[int(EventIndexArray[t-1]),int(EventIndexArray[t])]])
                elif 6.93 < timegap < 7.07 and counter == True:
                    self.EventList[0][1][1] = int(EventIndexArray[t])

       
        self.EventIndexArray = EventIndexArray                
        print(self.EventList)
        print(self.AnimalID)
        print(self.ProtType)
        
        if trim == 'end':
            self.Time = self.Time[:self.EventList[0][1][1]+2500]
            self.RedCh = self.RedCh[:self.EventList[0][1][1]+2500]
            self.GreenCh = self.GreenCh[:self.EventList[0][1][1]+2500]
            self.EventsCh = self.EventsCh[:self.EventList[0][1][1]+2500]
        elif trim == 'LTMUUS_LTM':
            self.Time = self.Time[:self.EventList[-1][1][0]]
            self.RedCh = self.RedCh[:self.EventList[-1][1][0]]
            self.GreenCh = self.GreenCh[:self.EventList[-1][1][0]]
            self.EventsCh = self.EventsCh[:self.EventList[-1][1][0]]
        elif trim == 'LTMUUS_UUS':
            self.Time = self.Time[self.EventList[-2][1][1]:]
            self.RedCh = self.RedCh[self.EventList[-2][1][1]:]
            self.GreenCh = self.GreenCh[self.EventList[-2][1][1]:]
            self.EventsCh = self.EventsCh[self.EventList[-2][1][1]:]
    
    def preprocess(self, baselinemethod='expdecay', denoise=True, denoiseHz=2, motioncontrol=True, figsave=False):
        if self.RedCh.min() == np.nan:
            self.RedCh = np.nan_to_num(self.RedCh,False)
            print('Missing values were in red channel!')
        if self.GreenCh.min() == np.nan:
            self.GreenCh = np.nan_to_num(self.GreenCh,False)
            print('Missing values were in green channel!')
        plt.rcParams['figure.figsize'] = [8, 8]
        #Preprocessing steps:
        plt.close()
        #prefiglist=[]
        plt.plot(self.Time, self.GreenCh, 'g', linewidth=0.5)
        plt.plot(self.Time, self.RedCh, 'r', linewidth=0.5)
        plt.title('Raw data for '+str(self.AnimalID)+str(self.ProtType))
        if figsave==True:
            plt.savefig(str(self.AnimalID)+str(self.ProtType)+'_raw_'+'.png')
        plt.show()
        plt.close()
        if denoise == True:
            # Median filtering to remove electrical artifact.
            kernelsize=int(round(self.samplingrate/100))
            if (kernelsize%2)==0:
                kernelsize=kernelsize+1
            print('kernelsize for median filter: '+str(kernelsize))
            Green_denoised = medfilt(self.GreenCh, kernel_size=kernelsize)
            Red_denoised = medfilt(self.RedCh, kernel_size=kernelsize)
            # Lowpass filter - zero phase filtering (with filtfilt) is used to avoid distorting the signal.
            b,a = butter(2, denoiseHz, btype='low', fs=self.samplingrate) #Lowpass butterworth filter with 10Hz
            Green_denoised = filtfilt(b,a, Green_denoised)
            Red_denoised = filtfilt(b,a, Red_denoised)
            plt.plot(self.Time, Green_denoised, 'g', linewidth=0.5)
            plt.plot(self.Time, Red_denoised, 'r', linewidth=0.5)
            plt.title('Denoised data for '+str(self.AnimalID)+str(self.ProtType))
            if figsave==True:
                plt.savefig(str(self.AnimalID)+str(self.ProtType)+'_denoised_'+'.png')
            plt.show()
            plt.close()
        elif denoise == False:
            Green_denoised = self.GreenCh
            Red_denoised = self.RedCh
            
        if baselinemethod == 'mean':
            self.Green_lfr = Green_denoised
            self.Red_lfr = Red_denoised
        elif baselinemethod == 'Polynomial':
            #Photobleaching correction with low-order polynomial
            # Fit 4th order polynomial to GCaMP signal.
            coefs_Green = np.polyfit(self.Time, Green_denoised, deg=4)
            Green_polyfit = np.polyval(coefs_Green, self.Time)
            # Fit 4th order polynomial to TdTomato signal.
            coefs_Red = np.polyfit(self.Time, Red_denoised, deg=4)
            Red_polyfit = np.polyval(coefs_Red, self.Time)
            # Plot fits
            plt.plot(self.Time, Green_denoised, 'g', label='NE_GRAB')
            plt.plot(self.Time, Green_polyfit,'k', linewidth=1.5) 
            plt.plot(self.Time, Red_denoised, 'r', label='red channel (for motion correction)')
            plt.plot(self.Time, Red_polyfit,'k', linewidth=1.5) 
            plt.title('Polynomial fit as baseline. Animal: '+str(self.AnimalID)+' '+str(self.ProtType))
            plt.xlabel('Time (seconds)')
            plt.ylabel('Volts (at photodetector)')
            plt.legend()
            if figsave==True:
                plt.savefig(str(self.AnimalID)+str(self.ProtType)+'_polyfit_'+'.png')
            plt.show()
            plt.clf()
            self.Green_lfr = Green_denoised - Green_polyfit + np.mean(Green_denoised)
            self.Red_lfr = Red_denoised - Red_polyfit + np.mean(Red_denoised)
        elif baselinemethod == 'expdecay':
            
            # The exponential curve we are going to fit.
            def exp_func(x, a, b, c):
               return a*np.exp(-b*x) + c
            
            # Fit curve to GFP-like signal.
            GFP_parms, parm_cov = curve_fit(exp_func, self.Time, Green_denoised, p0=[1,1e-3,1],bounds=([0,0,0],[4,0.1,4]), maxfev=1000)
            Green_expfit = exp_func(self.Time, *GFP_parms)
            
            # Fit curve to red signal. Not sure I should do this given there was no red fluorophore...
            Red_parms, parm_cov = curve_fit(exp_func, self.Time, Red_denoised, p0=[1,1e-3,1],bounds=([0,0,0],[4,0.1,4]), maxfev=1000)
            Red_expfit = exp_func(self.Time, *Red_parms)
            
            plt.plot(self.Time, Green_denoised, 'g', label='NE_GRAB')
            plt.plot(self.Time, Green_expfit,'k', linewidth=1.5) 
            plt.plot(self.Time, Red_denoised, 'r', label='mCherry')
            plt.plot(self.Time, Red_expfit,'k', linewidth=1.5) 
            plt.title('Exponential fit to bleaching for '+str(self.AnimalID)+str(self.ProtType))
            plt.xlabel('Time (seconds)');
            if figsave==True:
                plt.savefig(str(self.AnimalID)+str(self.ProtType)+'_expfit_'+'.png')
            plt.show()
            plt.clf()
            
            #substract exponential decay fit from data...
            Green_lfr = Green_denoised - Green_expfit
            Red_lfr = Red_denoised - Red_expfit
            plt.plot(self.Time, Green_lfr    , 'g', label='NE_GRAB')
            plt.plot(self.Time, Red_lfr-0.1, 'r', label='mCherry')
            #plt.title('Bleaching correction by subtraction of exponential decay fit')
            plt.title('No bleaching correction for '+str(self.AnimalID)+str(self.ProtType))
            plt.xlabel('Time (seconds)');
            if figsave==True:
                plt.savefig(str(self.AnimalID)+str(self.ProtType)+'_expcorrected_'+'.png')
            plt.show()
            plt.clf()
            self.Green_BL = Green_expfit
            self.Red_BL = Red_expfit
            self.Green_lfr = Green_lfr
            self.Red_lfr = Red_lfr
            
        if motioncontrol == True:
            slope, intercept, r_value, p_value, std_err = linregress(x=self.Red_lfr, y=self.Green_lfr)
            print('r_value: '+str(r_value))
            est_motion = intercept + slope * self.Red_lfr #
            plt.plot(self.Time, self.Green_lfr, 'b', linewidth=0.5, label='NE_GRAB pre-motion corrected', alpha=0.5)
            self.Green_lfr = self.Green_lfr - est_motion #+ np.mean(Green_denoised)
            plt.plot(self.Time, self.Green_lfr, 'g', linewidth=0.5, label='NE_GRAB motion corrected', alpha=0.8)
            plt.plot(self.Time, est_motion-0.01, 'r', linewidth=0.5, label='estimated motion', alpha=0.5)
            plt.title('Processing steps for #'+self.AnimalID+' '+self.ProtType, fontsize=14)
            plt.xlabel('Time (seconds)')
            plt.legend()
            if figsave==True:
                plt.savefig(str(self.AnimalID)+str(self.ProtType)+'_motion_'+'.png')
            plt.show()
            plt.clf()
        
        if baselinemethod == 'mean':
            self.Green_BL = np.mean(self.Green_lfr)*np.ones(self.Green_lfr.shape)
        elif baselinemethod == 'Polynomial':
            #Photobleaching correction with low-order polynomial
            # Fit 4th order polynomial to GCaMP signal.
            coefs_Green = np.polyfit(self.Time, self.Green_lfr, deg=4)
            Green_polyfit = np.polyval(coefs_Green, self.Time)
            self.Green_BL = Green_polyfit
            plt.plot(self.Time, self.Green_BL, 'k')
            if figsave==True:
                plt.savefig(str(self.AnimalID)+str(self.ProtType)+'_polyagain_'+'.png')
            plt.show()
            plt.clf()
            print(Green_polyfit.shape)
        
        self.GFP_dF_F = (np.divide((self.Green_lfr+self.Green_BL),self.Green_BL))-1
        print(self.GFP_dF_F.shape)
        plt.plot(self.Time, self.GFP_dF_F, color='g', linewidth=0.5, label='GFP_dF_F')
        plt.plot(self.Time, self.Green_lfr, 'y', linewidth=0.5, label='lfr')
        plt.plot(self.Time, self.Green_BL, 'k', linewidth=0.5, label='baseline')
        plt.xlabel('Time(sec)')
        plt.ylabel('dF/F')
        plt.title('GFP_dF_F for '+str(self.AnimalID)+str(self.ProtType))
        plt.legend()
        if figsave==True:    
            plt.savefig(str(self.AnimalID)+str(self.ProtType)+'_dF_'+'.png')
        plt.show()
        plt.clf()
        self.GFP_dF_F = self.GFP_dF_F *100
        
    def plot(self,buffer_left=20,buffer_right=60,headerstring=''):
        nplots = len(self.EventList)

        # Index, Time, Time_adj, Index_adj
        PlotList = np.zeros([nplots, 4])
        self.NameList = self.EventList[0][:]
        print(self.NameList)
        RectList = np.zeros([nplots-1,2])
        for i in range(nplots):
            s = self.EventList[i][1][0]
            e = self.EventList[i][1][1]
            ts = self.Time[s]
            te = self.Time[e]
            PlotList[i,0] = s
            PlotList[i,1] = e
            PlotList[i,2] = np.argmin(abs(self.Time-(ts-buffer_left)))
            PlotList[i,3] = np.argmin(abs(self.Time-(te+buffer_right)))   
            if i > 0: RectList[i-1,:] = ts,te
        self.PlotList = PlotList
        print(PlotList)     
        print(RectList)
        
        if self.ProtType == 'PTC_':
            plt.clf()
            plt.rcParams['figure.figsize'] = [12, 15]
            fig, axes = plt.subplots(5,1,constrained_layout=True)
            for i in range(len(axes)):
                
                z = i*2 
                try:
                    s = PlotList[z+1,2]
                    e = PlotList[z+1,3]
                    s=int(s)
                    e=int(e)
                    axes[i].plot(self.Time[s:e], self.GFP_dF_F[s:e], linewidth=0.5, color='g')
                    axes[i].plot(self.Time[s:e], self.EventsCh[s:e], color='y',alpha=0.8)
                    #axes[i,j].plot(self.Time[s:e], self.EventsCh[s:e]+np.min(self.GFP_dF_F[s:e]), color='y')
                    axes[i].set_ylabel('delta F (% from baseline)')
                    axes[i].set_xlabel('Time(sec)')
                    axes[i].set_title(self.AnimalID+' '+self.ProtType+' CSUS'+str(i+1))
                    if z == -1:
                        for nprow in range(RectList.shape[0]):
                            axes[i].axvspan(xmin=RectList[nprow,0],xmax=RectList[nprow,1], color='b', alpha=0.1)
                    else:
                        axes[i].axvspan(xmin=RectList[z+1,0],xmax=RectList[z+1,1], color='r', alpha=0.1)
                        axes[i].axvspan(xmin=RectList[z,0],xmax=RectList[z,1], color='b', alpha=0.1)
                except:
                    pass
            plt.show()
            return(fig)
        elif self.ProtType == 'HABI':
            plt.clf()
            plt.rcParams['figure.figsize'] = [12, 12]
            s = PlotList[0,0]
            e = PlotList[0,1]
            s=int(s)
            e=int(e)
            fig, axes = plt.subplots(1,1)
            axes.plot(self.Time[s:e], self.GFP_dF_F[s:e], linewidth=0.5, color='g')
            #axes.plot(self.Time[s:e], self.EventsCh[s:e]+np.min(self.GFP_dF_F[s:e]), color='y', alpha=0.5, linewidth=0.5)
            axes.set_ylabel('delta F (% from baseline)')
            axes.set_xlabel('Time(sec)')
            axes.set_title(self.AnimalID+' '+self.ProtType)
            plt.show()
            return(fig)
        else:
            plt.clf()
            plt.rcParams['figure.figsize'] = [12, 12] # Make default figure size larger.
            fig, axes = plt.subplots(nplots, 1, constrained_layout=True)
            #plt.text(0.01,0.99,headerstring,transform=fig.transFigure, size=10)
 #filelist, pdffileoutput, denoise=True, denoiseHz = 2, baselinemethod='mean', motioncontrol=False
            for z in range(1,nplots,1):
                s = PlotList[z,2]
                e = PlotList[z,3]
                s=int(s)
                e=int(e)
                if z==0:
                    linewidth=0.5
                else:
                    linewidth=1
                axes[z].plot(self.Time[s:e], self.GFP_dF_F[s:e], linewidth=linewidth, color='g')
                axes[z].plot(self.Time[s:e], self.EventsCh[s:e], color='y',alpha=0.8)
                #axes[z].plot(self.Time[s:e], self.EventsCh[s:e]+np.min(self.GFP_dF_F[s:e]), color='y')
                if z == 0:
                    for nprow in range(RectList.shape[0]):
                        if self.ProtType=='U-US':
                            axes[z].axvspan(xmin=RectList[nprow,0],xmax=RectList[nprow,1], color='r', alpha=0.1)
                        if self.ProtType=='PTC_':
                            axes[z].axvspan(xmin=RectList[nprow,1]-1,xmax=RectList[nprow,1], color='r', alpha=0.1)
                            axes[z].axvspan(xmin=RectList[nprow,0],xmax=RectList[nprow,1]-1, color='b', alpha=0.1)
                else:
                    if self.ProtType=='U-US' or self.ProtType=='LTMU':
                        if z == 6:
                            axes[z].axvspan(xmin=RectList[z-1,0],xmax=RectList[z-1,1], color='r', alpha=0.1)
                        else:
                            axes[z].axvspan(xmin=RectList[z-1,0],xmax=RectList[z-1,1], color='b', alpha=0.1)
                    if self.ProtType=='PTC_':
                        axes[z].axvspan(xmin=RectList[z-1,1]-1,xmax=RectList[z-1,1], color='r', alpha=0.1)
                        axes[z].axvspan(xmin=RectList[z-1,0],xmax=RectList[z-1,1]-1, color='b', alpha=0.1)
                #axes[z].axvspan(xmin=self.Time[int(PlotList[z,0])], xmax=self.Time[int(PlotList[z,1])],color='y',alpha=0.3)
                #axes[z].vlines(x=self.Time[s:e][self.EventsCh[s:e]==True], ymin=np.min(self.GFP_dF_F), ymax=np.max(self.GFP_dF_F), color='y')
                axes[z].set_ylabel('dF/F (%)')
                axes[z].set_xlabel('Time(sec)')
                axes[z].set_title(self.AnimalID+' '+self.ProtType+' '+str(z))
            plt.show()
            return(fig)

def AUC(TraceObj, EventListNum, duration=30):
    stimtime = TraceObj.Time[TraceObj.EventList[EventListNum][1][0]]
    stimidx = TraceObj.EventList[EventListNum][1][0]
    
    stimofftime = TraceObj.Time[TraceObj.EventList[EventListNum][1][1]]
    stimoffidx = TraceObj.EventList[EventListNum][1][1]
    
    if duration=='matchstim':
        duration = stimofftime - stimtime 
    
    pretime = TraceObj.Time[TraceObj.EventList[EventListNum][1][0]] - duration
    preidx = np.argmin(abs(TraceObj.Time-(pretime)))
    
    posttime = TraceObj.Time[TraceObj.EventList[EventListNum][1][1]] + duration
    postidx = np.argmin(abs(TraceObj.Time-(posttime)))
    
    #print(pretime),print(preidx)
    #print(stimtime),print(stimidx)
    #print(stimofftime),print(stimoffidx)
    #print(posttime),print(postidx)
    #print(stimofftime - stimtime)
    
    preGauc=np.trapz(TraceObj.GreenCh[preidx:stimidx], x=TraceObj.Time[preidx:stimidx])
    stimGauc=np.trapz(TraceObj.GreenCh[stimidx:stimoffidx], x=TraceObj.Time[stimidx:stimoffidx])
    postGauc=np.trapz(TraceObj.GreenCh[stimoffidx:postidx], x=TraceObj.Time[stimoffidx:postidx])
    
   # npre = len(preG)
    #nstim = len(stimG)
    #npost = len(postG)
    
    preRauc=np.trapz(TraceObj.RedCh[preidx:stimidx], x=TraceObj.Time[preidx:stimidx])
    stimRauc=np.trapz(TraceObj.RedCh[stimidx:stimoffidx], x=TraceObj.Time[stimidx:stimoffidx])
    postRauc=np.trapz(TraceObj.RedCh[stimoffidx:postidx], x=TraceObj.Time[stimoffidx:postidx])
    """
    plt.rcParams['figure.figsize']=[14,4]
    
    fig1, (axpre,axstim,axpost) = plt.subplots(nrows=1,ncols=3,constrained_layout=True)
    axpre.plot(TraceObj.Time[preidx:stimidx], TraceObj.GreenCh[preidx:stimidx], color='g', linewidth=0.5)
    axpre.plot(TraceObj.Time[preidx:stimidx], 0.056+0.004*TraceObj.EventsCh[preidx:stimidx], color='y')
    axpre.set_ylim(ymin=np.min([TraceObj.GreenCh[preidx:postidx].min(),0.056+0.004*TraceObj.EventsCh[preidx:postidx].min()]), ymax=np.max([TraceObj.GreenCh[preidx:postidx].max(),0.056+0.004*TraceObj.EventsCh[preidx:postidx].max()]))
    axpre.set_title('Pre-stim AUC='+str(preGauc/preGauc))
    axpre.set_ylabel('deltaF/F (%)')
    axpre.set_xlabel('Time(sec)')
    axstim.plot(TraceObj.Time[stimidx-1:stimoffidx+1], TraceObj.GreenCh[stimidx-1:stimoffidx+1], color='g', linewidth=0.5)
    axstim.plot(TraceObj.Time[stimidx-1:stimoffidx+1], 0.056+0.004*TraceObj.EventsCh[stimidx-1:stimoffidx+1], color='y')
    axstim.set_ylim(ymin=np.min([TraceObj.GreenCh[preidx:postidx].min(),0.056+0.004*TraceObj.EventsCh[preidx:postidx].min()]), ymax=np.max([TraceObj.GreenCh[preidx:postidx].max(),0.056+0.004*TraceObj.EventsCh[preidx:postidx].max()]))
    axstim.set_title('Stim AUC='+str(stimGauc/preGauc))
    axstim.set_xlabel('Time(sec)')
    axpost.plot(TraceObj.Time[stimoffidx:postidx], TraceObj.GreenCh[stimoffidx:postidx], color='g', linewidth=0.5)
    axpost.plot(TraceObj.Time[stimoffidx:postidx], 0.056+0.004*TraceObj.EventsCh[stimoffidx:postidx], color='y')
    axpost.set_ylim(ymin=np.min([TraceObj.GreenCh[preidx:postidx].min(),0.056+0.004*TraceObj.EventsCh[preidx:postidx].min()]), ymax=np.max([TraceObj.GreenCh[preidx:postidx].max(),0.056+0.004*TraceObj.EventsCh[preidx:postidx].max()]))
    axpost.set_title('Post-stim AUC='+str(postGauc/preGauc))
    axpost.set_xlabel('Time(sec)')
    fig1.suptitle('Green Channel for Event#'+str(EventListNum))
    plt.show()
    plt.close()
    
    fig2, (ax2pre,ax2stim,ax2post) = plt.subplots(nrows=1,ncols=3,constrained_layout=True)
    ax2pre.plot(TraceObj.Time[preidx:stimidx], TraceObj.RedCh[preidx:stimidx], color='r', linewidth=0.5)
    ax2pre.plot(TraceObj.Time[preidx:stimidx], 0.056+0.004*TraceObj.EventsCh[preidx:stimidx], color='y')
    ax2pre.set_ylim(ymin=np.min([TraceObj.RedCh[preidx:postidx].min(),0.056+0.004*TraceObj.EventsCh[preidx:postidx].min()]), ymax=np.max([TraceObj.RedCh[preidx:postidx].max(),0.056+0.004*TraceObj.EventsCh[preidx:postidx].max()]))
    ax2pre.set_title('Pre-stim AUC='+str(preRauc/preRauc))
    ax2pre.set_ylabel('deltaF/F (%)')
    ax2pre.set_xlabel('Time(sec)')
    ax2stim.plot(TraceObj.Time[stimidx-1:stimoffidx+1], TraceObj.RedCh[stimidx-1:stimoffidx+1], color='r', linewidth=0.5)
    ax2stim.plot(TraceObj.Time[stimidx-1:stimoffidx+1], 0.056+0.004*TraceObj.EventsCh[stimidx-1:stimoffidx+1], color='y')
    ax2stim.set_ylim(ymin=np.min([TraceObj.RedCh[preidx:postidx].min(),0.056+0.004*TraceObj.EventsCh[preidx:postidx].min()]), ymax=np.max([TraceObj.RedCh[preidx:postidx].max(),0.056+0.004*TraceObj.EventsCh[preidx:postidx].max()]))
    ax2stim.set_title('Stim AUC='+str(stimRauc/preRauc))
    ax2stim.set_xlabel('Time(sec)')
    ax2post.plot(TraceObj.Time[stimoffidx:postidx], TraceObj.RedCh[stimoffidx:postidx], color='r', linewidth=0.5)
    ax2post.plot(TraceObj.Time[stimoffidx:postidx], 0.056+0.004*TraceObj.EventsCh[stimoffidx:postidx], color='y')
    ax2post.set_ylim(ymin=np.min([TraceObj.RedCh[preidx:postidx].min(),0.056+0.004*TraceObj.EventsCh[preidx:postidx].min()]), ymax=np.max([TraceObj.RedCh[preidx:postidx].max(),0.056+0.004*TraceObj.EventsCh[preidx:postidx].max()]))
    ax2post.set_title('Post-stim AUC='+str(postRauc/preRauc))
    ax2post.set_xlabel('Time(sec)')
    fig2.suptitle('Red Channel for Event#'+str(EventListNum))
    plt.show()
    plt.close()
    """
    AUClist=np.array([preGauc,stimGauc,postGauc,preRauc,stimRauc,postRauc])
    if duration=='matchstim':
        normAUCGlist=np.divide(AUClist[1:3],preGauc)
        normAUCRlist=np.divide(AUClist[4:],preRauc)
    else:
        normAUCGlist=np.divide(AUClist[2],preGauc)
        normAUCRlist=np.divide(AUClist[5],preRauc)
    return(AUClist,normAUCGlist,normAUCRlist)
    
    


def LTMUUSAUC(TraceObj):
    InterestList=[]
    for i in [1,2,3,4,5,6]: 
        (AUClist,normAUCGlist,normAUCRlist)=AUC(TraceObj,i,duration=30)
        #print('Green normed AUC for event #'+str(i)+': '+ str(normAUCGlist[0])+', '+str(normAUCGlist[1]))
        #print('Red normed AUC for event #'+str(i)+': '+ str(normAUCRlist[0])+', '+str(normAUCRlist[1]))
        InterestList.append((AUClist,normAUCGlist,normAUCRlist))
    for contents in InterestList:
        print('Green then Red (% change from pre-to-stim then pre-to-post)')
        print((contents[1]-1)*100)
        print((contents[2]-1)*100)
    return(InterestList)

def PTCAUC(TraceObj):
    InterestList=[]
    for i in [1,3,5,7,9]: 
        (AUClist,normAUCGlist,normAUCRlist)=AUC(TraceObj,i,duration=30)
        #print('Green normed AUC for event #'+str(i)+': '+ str(normAUCGlist[0])+', '+str(normAUCGlist[1]))
        #print('Red normed AUC for event #'+str(i)+': '+ str(normAUCRlist[0])+', '+str(normAUCRlist[1]))
        InterestList.append((AUClist,normAUCGlist,normAUCRlist))
    for contents in InterestList:
        print('Green then Red (% change from pre-to-stim then pre-to-post)')
        print((contents[1]-1)*100)
        print((contents[2]-1)*100)
    return(InterestList)


def CohensD_prepoststim(TraceObj, EventListNum, duration=30):
    # For calculating averages and cohen's d:
    
    
    stimtime = TraceObj.Time[TraceObj.EventList[EventListNum][1][0]]
    stimidx = TraceObj.EventList[EventListNum][1][0]
    
    stimofftime = TraceObj.Time[TraceObj.EventList[EventListNum][1][1]]
    stimoffidx = TraceObj.EventList[EventListNum][1][1]
    
    if duration=='matchstim':
        duration = stimofftime - stimtime 
    
    pretime = TraceObj.Time[TraceObj.EventList[EventListNum][1][0]] - duration
    preidx = np.argmin(abs(TraceObj.Time-(pretime)))
    
    posttime = TraceObj.Time[TraceObj.EventList[EventListNum][1][1]] + duration
    postidx = np.argmin(abs(TraceObj.Time-(posttime)))
    
    print(pretime),print(preidx)
    print(stimtime),print(stimidx)
    print(stimofftime),print(stimoffidx)
    print(posttime),print(postidx)
    print(stimofftime - stimtime)
    
    pre=TraceObj.GFP_dF_F[preidx:stimidx]
    stim=TraceObj.GFP_dF_F[stimidx:stimoffidx]
    post=TraceObj.GFP_dF_F[stimoffidx:postidx]
    
    npre = len(pre)
    nstim = len(stim)
    npost = len(post)
    
    poolstd = ((pre.var()*(npre-1)) + (stim.var()*(nstim-1)))/(npre+nstim-2)
    CohensD12 = (stim.mean()-pre.mean())/poolstd
    #print('Pooled st.dev.(stim-pre): '+str(poolstd))
    print('Cohen\'s D (pre-to-stim): '+str(CohensD12))

    poolstd = ((stim.var()*(nstim-1)) + (post.var()*(npost-1)))/(nstim+npost-2)
    CohensD23 = (post.mean()-stim.mean())/poolstd    
    #print('Pooled st.dev. (post-stim): '+str(poolstd))
    print('Cohen\'s D (stim-to-post): '+str(CohensD23))

    poolstd = ((pre.var()*(npre-1)) + (post.var()*(npost-1)))/(npre+npost-2)
    CohensD13 = (post.mean()-pre.mean())/poolstd    
    #print('Pooled st.dev. (post-pre): '+str(poolstd))
    print('Cohen\'s D (pre-to-post): '+str(CohensD13))
    
    plt.rcParams['figure.figsize']=[14,8]
    
    fig, (axpre,axstim,axpost) = plt.subplots(nrows=1,ncols=3,constrained_layout=True)
    axpre.plot(TraceObj.Time[preidx:stimidx], TraceObj.GFP_dF_F[preidx:stimidx], color='g', linewidth=0.5)
    axpre.plot(TraceObj.Time[preidx:stimidx], TraceObj.EventsCh[preidx:stimidx], color='y')
    axpre.set_ylim(ymin=np.min([TraceObj.GFP_dF_F[preidx:postidx].min(),-0.1]), ymax=np.max([TraceObj.GFP_dF_F[preidx:postidx].max(),1.1]))
    axpre.set_title('Pre-stim mean: %f' % pre.mean())
    axpre.set_ylabel('deltaF/F (%)')
    axpre.set_xlabel('Time(sec)')
    axstim.plot(TraceObj.Time[stimidx-1:stimoffidx+1], TraceObj.GFP_dF_F[stimidx-1:stimoffidx+1], color='g', linewidth=0.5)
    axstim.plot(TraceObj.Time[stimidx-1:stimoffidx+1], TraceObj.EventsCh[stimidx-1:stimoffidx+1], color='y')
    axstim.set_ylim(ymin=np.min([TraceObj.GFP_dF_F[preidx:postidx].min(),-0.1]), ymax=np.max([TraceObj.GFP_dF_F[preidx:postidx].max(),1.1]))
    axstim.set_title('Stim mean: %f' % stim.mean())
    axstim.set_xlabel('Time(sec)')
    axpost.plot(TraceObj.Time[stimoffidx:postidx], TraceObj.GFP_dF_F[stimoffidx:postidx], color='g', linewidth=0.5)
    axpost.plot(TraceObj.Time[stimoffidx:postidx], TraceObj.EventsCh[stimoffidx:postidx], color='y')
    axpost.set_ylim(ymin=np.min([TraceObj.GFP_dF_F[preidx:postidx].min(),-0.1]), ymax=np.max([TraceObj.GFP_dF_F[preidx:postidx].max(),1.1]))
    axpost.set_title('Post-stim mean: %f' % post.mean())
    axpost.set_xlabel('Time(sec)')
    
    fig.suptitle('Cohen\'s D (pre-to-stim): '+str(CohensD12)+'\n'+'Cohen\'s D (stim-to-post): '+str(CohensD23)+'\n'+'Cohen\'s D (pre-to-post): '+str(CohensD13))
    plt.show()
    return(fig)


path='C:\\Users\\waltp\\OneDrive\\Documents\\Dec2019 for transfer\\Fiber Photometry data\\December14 2019\\'
date='_dec14_2019'
IDlist = ['295','296']#,'298','299','300','301','302','303','305','306']
filelist=[]
ProtocolType = 'PTC'
for ID in IDlist: filelist.append(path+ProtocolType+'_walter_285'+str(ID)+date+'.csv')
#ObjectList2=PlotSaveAll(filelist,pdffileoutput='NE2h LTMUUS with CohensD.pdf')
AllList=[]
for file in filelist:
    TraceObj = FullTrace(file)
    InterestList = PTCAUC(TraceObj)
    AllList.append(InterestList)
"""
path='C:\\Users\\waltp\\OneDrive\\Documents\\Dec2019 for transfer\\Fiber Photometry data\\December15 2019\\'
date='_dec15_2019'
IDlist = ['295','296','298','301','302','303']
filelist=[]
ProtocolType = 'LTMUUS'
for ID in IDlist: filelist.append(path+ProtocolType+'_walter_285'+str(ID)+date+'.csv')
#ObjectList2=PlotSaveAll(filelist,pdffileoutput='NE2h LTMUUS with CohensD.pdf')
for file in filelist:
    TraceObj = FullTrace(file)
    InterestList = LTMUUSAUC(TraceObj)
    AllList.append(InterestList)
"""

"""
PTC303 = FullTrace(csvfile='C:\\Users\\wpiper\\Desktop\\Fiber Photometry data\\December14 2019\\PTC_walter_285303_dec14_2019.csv')
PTC303.preprocess(denoiseHz=2,motioncontrol=True,baselinemethod='expdecay')
fig = PTC303.plot()
"""
"""
filelist=[]
IDlist=[295,296,298,301,302,303]
for ID in IDlist: filelist.append('C:\\Users\\wpiper\\Desktop\\Fiber Photometry data\\December15 2019\\LTMUUS_walter_285'+str(ID)+'_dec15_2019.csv')

PlotSaveAll(filelist,pdffileoutput='NE2h_LTMUUS.pdf', motioncontrol=True)
"""
"""
#PTC_walter_285303_dec14_2019.csv
path='C:\\Users\\wpiper\\Desktop\\Fiber Photometry data\\December14 2019\\'
date='_dec14_2019'
IDlist=[295,296,298,299,300,301,302,303,305,306]
"""


"""
objectlist=[]
pdffileoutput='Dec15 UUSonly automated analysis_nomotionctrl.pdf'
#with PdfPages(pdffileoutput) as pp:
for file in filelist:
    #LTML=FullTrace(csvfile=file,trim='LTMUUS_LTM')
    LTMU=FullTrace(csvfile=file,trim='LTMUUS_UUS')
    #LTML.preprocess(denoiseHz=2, baselinemethod='expdecay', motioncontrol=True)
    LTMU.preprocess(denoiseHz=2, baselinemethod='expdecay', motioncontrol=True)
    #objectlist.append([LTMU,LTML,file])
    plt.close()
    plt.plot(LTMU.Time,LTMU.GFP_dF_F,color='g')
    plt.plot(LTMU.Time,LTMU.EventsCh,color='y')
    plt.show()

    buffer_left=20
    buffer_right=60
    s = LTMU.EventList[-1][1][0] - LTMU.EventList[-2][1][1]
    e = LTMU.EventList[-1][1][1] - LTMU.EventList[-2][1][1]
    ts = LTMU.Time[s]
    te = LTMU.Time[e]
    sbuf = np.argmin(abs(LTMU.Time-(ts-buffer_left)))
    ebuf = np.argmin(abs(LTMU.Time-(te+buffer_right)))
    
    plt.close()
    plt.plot(LTMU.Time[sbuf:ebuf],LTMU.GFP_dF_F[sbuf:ebuf],color='g',linewidth=0.5)
    plt.plot(LTMU.Time[sbuf:ebuf],LTMU.EventsCh[sbuf:ebuf],color='y',linewidth=1)
    plt.ylabel('delta F (% from baseline)')
    plt.xlabel('Time(sec)')
    plt.title(LTMU.AnimalID+' '+LTMU.ProtType+' U-US')
    pp.savefig()
    #plt.savefig(str(LTMU.AnimalID)+str(LTMU.ProtType)+'_UUS_'+'.png')
    plt.show()
    
    buffer_left=20
    buffer_right=60
    s=LTMUUS298.EventList[-1][1][0]
    e=LTMUUS298.EventList[-1][1][1]
    ts=LTMUUS298.Time[s]
    te=LTMUUS298.Time[e]
    sbuf = np.argmin(abs(LTMUUS298.Time-(ts-buffer_left)))
    ebuf = np.argmin(abs(LTMUUS298.Time-(te+buffer_right)))
    
    plt.close()
    plt.plot(LTMUUS298.Time[sbuf:ebuf],LTMUUS298.GFP_dF_F[sbuf:ebuf],color='g',linewidth=0.5)
    plt.plot(LTMUUS298.Time[sbuf:ebuf],LTMUUS298.EventsCh[sbuf:ebuf],color='y')
    plt.show()
    plt.close()
    plt.plot(X.Time[sbuf:ebuf],X.GFP_dF_F[sbuf:ebuf],color='g',linewidth=0.5)
    plt.plot(X.Time[sbuf:ebuf],X.EventsCh[sbuf:ebuf],color='y')
    plt.xlabel('Time(sec)')
    plt.ylabel('deltaF/F (%)')
    plt.title(X.AnimalID+' '+X.ProtType+' U-US')
    plt.show()
"""

"""
ptc303 = FullTrace(csvfile=path+'December14 2019\\PTC_walter_285303_dec14_2019.csv')
tst303 = FullTrace(csvfile=path+'December14 2019\\tst_walter_285303_dec14_2019.csv')  
ltm303 = FullTrace(csvfile=path+'December15 2019\\LTMUUS_walter_285303_dec15_2019.csv')
uus303 = FullTrace(csvfile=path+'December15 2019\\LTMUUS_walter_285303_dec15_2019_4.csv')
list303 = [ptc303,tst303,ltm303,uus303]
for all303 in list303:
    print(all303.GreenCh.mean())
    print(all303.GreenCh.min())
    print(all303.GreenCh.max())
    print(all303.GreenCh.std())

plt.close()
fig,((ax1),(ax2),(ax3),(ax4))=plt.subplots(nrows=1,ncols=4)
ax1.plot(ptc303.Time[28800:32000],ptc303.GreenCh[28800:32000])
ax1.set_ylim(ymin=0,ymax=ptc303.GreenCh[28800:32000].max())
ax2.plot(tst303.Time[8800:12000],tst303.GreenCh[8800:12000])
ax2.set_ylim(ymin=0,ymax=tst303.GreenCh[8800:12000].max())
ax3.plot(ltm303.Time[180000:200000],ltm303.GreenCh[180000:200000])
ax3.set_ylim(ymin=0,ymax=ltm303.GreenCh[180000:200000].max())
ax4.plot(uus303.Time[130000:150000],uus303.GreenCh[130000:150000])
ax4.set_ylim(ymin=0,ymax=uus303.GreenCh[130000:150000].max())
plt.show()
plt.close()
"""
"""
ptc298 = FullTrace(csvfile=path+'December14 2019\\PTC_walter_285298_dec14_2019.csv')
ltm298 = FullTrace(csvfile=path+'December15 2019\\LTMUUS_walter_285298_dec15_2019.csv')
plt.close()
fig,((ax1),(ax3))=plt.subplots(nrows=1,ncols=2)
ax1.plot(ptc298.Time[28800:32000],ptc298.GreenCh[28800:32000])
ax1.set_ylim(ymin=0,ymax=ptc298.GreenCh[28800:32000].max())
ax3.plot(ltm298.Time[180000:200000],ltm298.GreenCh[180000:200000])
ax3.set_ylim(ymin=0,ymax=ltm298.GreenCh[180000:200000].max())
plt.show()
plt.close()

ptc296 = FullTrace(csvfile=path+'December14 2019\\PTC_walter_285296_dec14_2019.csv')
ltm296 = FullTrace(csvfile=path+'December15 2019\\LTMUUS_walter_285296_dec15_2019.csv')
plt.close()
fig,((ax1),(ax3))=plt.subplots(nrows=1,ncols=2)
ax1.plot(ptc296.Time[28800:32000],ptc296.GreenCh[28800:32000])
ax1.set_ylim(ymin=0,ymax=ptc296.GreenCh[28800:32000].max())
ax3.plot(ltm296.Time[180000:200000],ltm296.GreenCh[180000:200000])
ax3.set_ylim(ymin=0,ymax=ltm296.GreenCh[180000:200000].max())
plt.show()
plt.close()

ptc301 = FullTrace(csvfile=path+'December14 2019\\PTC_walter_285301_dec14_2019.csv')
ltm301 = FullTrace(csvfile=path+'December15 2019\\LTMUUS_walter_285301_dec15_2019.csv')
plt.close()
fig,((ax1),(ax3))=plt.subplots(nrows=1,ncols=2)
ax1.plot(ptc301.Time[28800:32000],ptc301.GreenCh[28800:32000])
ax1.set_ylim(ymin=0,ymax=ptc301.GreenCh[28800:32000].max())
ax3.plot(ltm301.Time[180000:200000],ltm301.GreenCh[180000:200000])
ax3.set_ylim(ymin=0,ymax=ltm301.GreenCh[180000:200000].max())
plt.show()
plt.close()


fac=[0.1,1]
for fac in fac:
    dist25=int(round(320*fac))
    dist4=int(round(2000*fac))
    s25=28800
    e25=s25+dist25
    s4=180000
    e4=s4+dist4
    # s25:e25
    # s4:e4
    
    plt.close()
    fig,((ax1),(ax2),(ax3))=plt.subplots(nrows=3,ncols=1)
    ax1.plot(ptc303.Time[s25:e25],ptc303.GreenCh[s25:e25],color='g')
    ax1.set_ylim(ymin=ptc303.GreenCh[s25:e25].min(),ymax=ptc303.GreenCh[s25:e25].max())
    ax2.plot(ltm303.Time[s4:e4],ltm303.GreenCh[s4:e4],color='g')
    ax2.set_ylim(ymin=ltm303.GreenCh[s4:e4].min(),ymax=ltm303.GreenCh[s4:e4].max())
    ax3.plot(uus303.Time[s4:e4],uus303.GreenCh[s4:e4],color='g')
    ax3.set_ylim(ymin=uus303.GreenCh[s4:e4].min(),ymax=uus303.GreenCh[s4:e4].max())
    plt.show()
    plt.close()


fac=[0.03,1]
for fac in fac:
    dist25=int(round(320*fac))
    dist4=int(round(2000*fac))
    s25=28800
    e25=s25+dist25
    s4=180000
    e4=s4+dist4
    # s25:e25
    # s4:e4
    
    plt.close()
    fig,((ax1),(ax2),(ax3))=plt.subplots(nrows=3,ncols=1)
    ax1.plot(ptc303.Time[s25:e25],ptc303.RedCh[s25:e25],color='r')
    ax1.set_ylim(ymin=ptc303.RedCh[s25:e25].min(),ymax=ptc303.RedCh[s25:e25].max())
    ax2.plot(ltm303.Time[s4:e4],ltm303.RedCh[s4:e4],color='r')
    ax2.set_ylim(ymin=ltm303.RedCh[s4:e4].min(),ymax=ltm303.RedCh[s4:e4].max())
    ax3.plot(uus303.Time[s4:e4],uus303.RedCh[s4:e4],color='r')
    ax3.set_ylim(ymin=uus303.RedCh[s4:e4].min(),ymax=uus303.RedCh[s4:e4].max())
    plt.show()
    plt.close()

#ptc303.preprocess(baselinemethod='expdecay',denoiseHz=2,motioncontrol=False)
#ltm303.preprocess(baselinemethod='expdecay',denoiseHz=2,motioncontrol=False)
#uus303.preprocess(baselinemethod='expdecay',denoiseHz=2,motioncontrol=False)


fac=[0.1,1]
for fac in fac:
    dist25=int(round(320*fac))
    dist4=int(round(2000*fac))
    s25=28800
    e25=s25+dist25
    s4=180000
    e4=s4+dist4
    # s25:e25
    # s4:e4
    
    plt.close()
    fig,((ax1),(ax2),(ax3))=plt.subplots(nrows=3,ncols=1)
    ax1.plot(ptc303.Time[s25:e25],ptc303.Green_lfr[s25:e25],color='g')
    ax1.set_ylim(ymin=ptc303.Green_lfr[s25:e25].min(),ymax=ptc303.Green_lfr[s25:e25].max())
    ax2.plot(ltm303.Time[s4:e4],ltm303.Green_lfr[s4:e4],color='g')
    ax2.set_ylim(ymin=ltm303.Green_lfr[s4:e4].min(),ymax=ltm303.Green_lfr[s4:e4].max())
    ax3.plot(uus303.Time[s4:e4],uus303.Green_lfr[s4:e4],color='g')
    ax3.set_ylim(ymin=uus303.Green_lfr[s4:e4].min(),ymax=uus303.Green_lfr[s4:e4].max())
    plt.show()
    plt.close()


fac=[0.03,1]
for fac in fac:
    dist25=int(round(320*fac))
    dist4=int(round(2000*fac))
    s25=28800
    e25=s25+dist25
    s4=180000
    e4=s4+dist4
    # s25:e25
    # s4:e4
    
    plt.close()
    fig,((ax1),(ax2),(ax3))=plt.subplots(nrows=3,ncols=1)
    ax1.plot(ptc303.Time[s25:e25],ptc303.Red_lfr[s25:e25],color='r')
    ax1.set_ylim(ymin=ptc303.Red_lfr[s25:e25].min(),ymax=ptc303.Red_lfr[s25:e25].max())
    ax2.plot(ltm303.Time[s4:e4],ltm303.Red_lfr[s4:e4],color='r')
    ax2.set_ylim(ymin=ltm303.Red_lfr[s4:e4].min(),ymax=ltm303.Red_lfr[s4:e4].max())
    ax3.plot(uus303.Time[s4:e4],uus303.Red_lfr[s4:e4],color='r')
    ax3.set_ylim(ymin=uus303.Red_lfr[s4:e4].min(),ymax=uus303.Red_lfr[s4:e4].max())
    plt.show()
    plt.close()

#-----------------------------------------
#-----------------------------------------
#-----------------------------------------
#-----------------------------------------
#-----------------------------------------
#-----------------------------------------
#-----------------------------------------
#-----------------------------------------
#-----------------------------------------
#-----------------------------------------
  
rcParams['figure.figsize']=[14,3]  
plt.close()  

fac=[3,10,100]
for fac in fac:
    dist4=uus303.samplingrate*fac
    s4=int(round(uus303.EventList[1][1][0] - dist4/2))
    e4=int(round(uus303.EventList[1][1][0] + dist4))
    if s4 < 0:
        s4=0
    if e4 > len(uus303.Time):
        e4=len(uus303.Time)
        # s4:e4
    print(dist4), print(s4), print(e4)
    plt.close()
    fig,(ax1)=plt.subplots(nrows=1,ncols=1)

    plt.plot(uus303.Time[s4:e4],uus303.Red_lfr[s4:e4],color='r', linewidth=0.5)
    plt.plot(uus303.Time[s4:e4],uus303.EventsCh[s4:e4]*(uus303.Red_lfr[s4:e4].max()-uus303.Red_lfr[s4:e4].min())*(0.99)+uus303.Red_lfr[s4:e4].min(),color='y')
    plt.ylim(ymin=uus303.Red_lfr[s4:e4].min(),ymax=uus303.Red_lfr[s4:e4].max())
    plt.show()
    plt.close()
    
    
fac=[3,10,100]
for fac in fac:
    dist4=uus303.samplingrate*fac
    s4=int(round(uus303.EventList[1][1][0] - dist4/2))
    e4=int(round(uus303.EventList[1][1][0] + dist4))
    if s4 < 0:
        s4=0
    if e4 > len(uus303.Time):
        e4=len(uus303.Time)
        # s4:e4
    print(dist4), print(s4), print(e4)
    plt.close()
    fig,(ax1)=plt.subplots(nrows=1,ncols=1)

    ax1.plot(uus303.Time[s4:e4],uus303.Green_lfr[s4:e4],color='b', linewidth=0.5)
    ax1.plot(uus303.Time[s4:e4],(uus303.EventsCh[s4:e4]*(uus303.Green_lfr[s4:e4].max()-uus303.Green_lfr[s4:e4].min())*(0.99))+uus303.Green_lfr[s4:e4].min(),color='y')
    ax1.set_ylim(ymin=uus303.Green_lfr[s4:e4].min(),ymax=uus303.Green_lfr[s4:e4].max())
    plt.show()
    plt.close()

fac=[3,10,100]
for fac in fac:
    dist4=uus303.samplingrate*fac
    s4=int(round(uus303.EventList[1][1][0] - dist4/2))
    e4=int(round(uus303.EventList[1][1][0] + dist4))
    if s4 < 0:
        s4=0
    if e4 > len(uus303.Time):
        e4=len(uus303.Time)
        # s4:e4
    print(dist4), print(s4), print(e4)
    plt.close()
    fig,(ax1)=plt.subplots(nrows=1,ncols=1)

    ax1.plot(uus303.Time[s4:e4],uus303.GFP_dF_F[s4:e4],color='g', linewidth=0.5)
    ax1.plot(uus303.Time[s4:e4],(uus303.EventsCh[s4:e4]*(uus303.GFP_dF_F[s4:e4].max()-uus303.GFP_dF_F[s4:e4].min())*(0.99))+uus303.GFP_dF_F[s4:e4].min(),color='y')
    ax1.set_ylim(ymin=uus303.GFP_dF_F[s4:e4].min(),ymax=uus303.GFP_dF_F[s4:e4].max())
    plt.show()
    plt.close()
"""







"""
plt.close()
fig,((ax1),(ax2),(ax3),(ax4))=plt.subplots(nrows=1,ncols=4)
ax1.plot(ptc303.Time[28800:32000],ptc303.GreenCh[28800:32000])
ax1.set_ylim(ymin=0,ymax=ptc303.GreenCh[28800:32000].max())
ax2.plot(tst303.Time[8800:12000],tst303.GreenCh[8800:12000])
ax2.set_ylim(ymin=0,ymax=tst303.GreenCh[8800:12000].max())
ax3.plot(ltm303.Time[180000:200000],ltm303.GreenCh[180000:200000])
ax3.set_ylim(ymin=0,ymax=ltm303.GreenCh[180000:200000].max())
ax4.plot(uus303.Time[130000:150000],uus303.GreenCh[130000:150000])
ax4.set_ylim(ymin=0,ymax=uus303.GreenCh[130000:150000].max())
plt.show()
plt.close()

plt.close()
plt.plot(ltm303.Time[180000:200000],ltm303.GreenCh[180000:200000])
plt.ylim(ymin=0,ymax=ltm303.GreenCh[180000:200000].max())
plt.show()
plt.close()
plt.plot(uus303.Time[180000:200000],uus303.GreenCh[180000:200000])
plt.ylim(ymin=0,ymax=uus303.GreenCh[180000:200000].max())
plt.show()
"""
#filename = ["U-US_walter_282038_oct20_2019.csv", "U-US_walter_282039_oct20_2019.csv" , "U-US_walter_282040_oct20_2019.csv","LTM_walter_282038_oct20_2019.csv","LTM_walter_282039_oct20_2019.csv","LTM_walter_282040_oct20_2019.csv", "PTC_walter_282038_oct19_2019.csv","PTC_walter_282039_oct19_2019.csv","PTC_walter_282040_oct19_2019.csv","HABITUATION_walter_282038_oct18_2019.csv","HABITUATION_walter_282039_oct18_2019.csv","HABITUATION_walter_282040_oct18_2019.csv"] # <--- SET CSV FILES TO BE PROCESSED