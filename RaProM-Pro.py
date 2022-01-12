##SCRIPT FOR READING AND PROCESSING DATA FROM MRR PRO IN A FOLDER
##The script is compatible with spectrum_raw and spectrum_reflectivity

import numpy as np
import calendar
import datetime
import miepython as mp
from math import e
from netCDF4 import Dataset,num2date,date2num
import netCDF4 as nc4
import glob
import os
import sys


def Aliasing(matriu,fNy,he,temps):#analyse to correct the aliasing
    numberDoppler=len(matriu[0])
    speed=np.arange(0,fNy*numberDoppler,fNy)
    speeddeal=np.arange(-1*fNy*numberDoppler,2*fNy*numberDoppler,fNy)
    etaN=np.copy(matriu)

    etaN_da=[]
    
    vec_null=np.copy(matriu[0])*np.nan
    for i in range(len(matriu)):
        v1=matriu[i]
        
        if i==0 or i==len(matriu)-1:
            if i==0:
                vect=np.concatenate([vec_null,matriu[i],matriu[i+1]])
            if i==len(matriu)-1:
                vect=np.concatenate([matriu[i-1],matriu[i],vec_null])
        else:
            vect=np.concatenate([matriu[i-1],matriu[i],matriu[i+1]])

        etaN_da.append(vect)



            
######CREATE THE SPEED PROFILE
    W=[];Imax=[];Imin=[];Wmax=[];Sig=[];LeNnan=[]

    for i in range(len(matriu)):
        vect=matriu[i]


        PT3=np.nansum(vect)
        w3=np.nansum(np.prod([vect,speed],axis=0))/PT3#estimated velocity
        sigma3=np.sqrt(np.nansum(np.prod([vect,np.power(speed-w3,2)],axis=0))/PT3)# spectral witdh
        W.append(w3)
        Sig.append(sigma3)
        vect2=[];vect3=[]

        vect[vect == np.inf] = np.nan#replace inf values
####        DELETE ISOLATED VALUES
        if np.isnan(vect).all():
            co=numberDoppler/2
        else:
            for j in range(len(vect)):#adapt the vector to delete the unique elements 
                if j==0 or j>=len(vect)-2:
                    if j==0:
                        if ~np.isnan(vect[j]) and np.isnan(vect[j+1]):
                            vect2.append(np.nan)
                        else:
                            vect2.append(vect[j])
                    else:
                        if ~np.isnan(vect[j]) and np.isnan(vect[j-1]):
                            vect2.append(np.nan)
                        else:
                            vect2.append(vect[j])
                        
                else:

                    if ~np.isnan(vect[j]) and np.isnan(vect[j-1]) and np.isnan(vect[j+1]):
                        vect2.append(np.nan)
                    else:
                        vect2.append(vect[j])
######            FOUND THE POSITION OF VALUES

            for j in range(len(vect2)):
                if ~np.isnan(vect2[j]):
                    vect3.append(j)
            First=vect3[0]
            Darrer=vect3[-1]
            co=0
            for j in range(Darrer-First):
                if np.isnan(vect2[j+First]):
                    co+=1
            
        LeNnan.append(co)

        if np.isnan(vect).all():
            Imin.append(np.nan)
            Imax.append(np.nan)

        else:
            Imin.append(First)
            Imax.append(Darrer)       

########     CALCULATE THE W i new matrix DEALISED
    ReVect=[]
    ReVect2=[]
    nul1=np.nan*np.arange(0,fNy*numberDoppler/2,fNy)
    nul2=np.nan*np.arange(0,fNy*numberDoppler*3/2,fNy)
    Nul1=np.nan*np.arange(0,fNy*numberDoppler*5/4,fNy)
    Nul2=np.nan*np.arange(0,fNy*numberDoppler*3/4,fNy)
    W_da=[];Sig_da=[]
    for i in range(len(matriu)):
        Diff=np.diff(W)
        vect=etaN_da[i]
        if i<=2 or np.isnan(Imin[i]) or np.isnan(Imax[i]):
            W_da.append(W[i])
            
            deal=np.concatenate([vec_null,vect[int(numberDoppler):int(numberDoppler*2)],vec_null])
            if np.isnan(deal).all():
                VectRe=vect*np.nan
            else:
                
                I=np.nanargmax(np.asarray(deal))
                Indvel=speeddeal*(vect/vect)
                VectRe,inewV=group(vect,I,3,Indvel)
            
        else:
            if Imax[i]-Imin[i]>(numberDoppler/2) and LeNnan[i]>numberDoppler/4:#apply 2 condition, one for if exist aliasing and the second to avoid a large vector
                if W[i-1]>(fNy*numberDoppler/2):
                    vect2=vect[int(numberDoppler*5/4):int(numberDoppler*9/4)]
                    speed2=speeddeal[int(numberDoppler*5/4):int(numberDoppler*9/4)]
                    PT4=np.nansum(vect2)
                    w4=np.nansum(np.prod([vect2,speed2],axis=0))/PT4#estimated velocity
                    sigma4=np.sqrt(np.nansum(np.prod([vect2,np.power(speed2-w4,2)],axis=0))/PT4)# spectral witdh
                    W_da.append(w4)

                    deal=np.concatenate([Nul1,vect2,Nul2])
                    if np.isnan(deal).all():
                        VectRe=vect*np.nan
                    else:
                            
                        I=np.nanargmax(np.asarray(deal))
                        Indvel=speeddeal*(vect/vect)
                        VectRe,inewV=group(vect,I,3,Indvel)
                else:
                    if W[i-1]<(fNy*numberDoppler/3):
                        vect2=vect[int(numberDoppler*3/4):int(numberDoppler*5/4)]
                        speed2=speeddeal[int(numberDoppler*3/4):int(numberDoppler*5/4)]
                        PT4=np.nansum(vect2)
                        w4=np.nansum(np.prod([vect2,speed2],axis=0))/PT4#estimated velocity
                        sigma4=np.sqrt(np.nansum(np.prod([vect2,np.power(speed2-w4,2)],axis=0))/PT4)# spectral witdh

                        
                        W_da.append(w4)
                        deal=np.concatenate([Nul2,vect2,Nul1])
                        if np.isnan(deal).all():
                            VectRe=vect*np.nan
                        else:
                                    
                            I=np.nanargmax(np.asarray(deal))
                            Indvel=speeddeal*(vect/vect)
                            VectRe,inewV=group(vect,I,3,Indvel)
                    else:
                        vect2=vect[int(numberDoppler):int(numberDoppler*2)]
                        speed2=speeddeal[int(numberDoppler):int(numberDoppler*2)]
                        PT4=np.nansum(vect2)
                        w4=np.nansum(np.prod([vect2,speed2],axis=0))/PT4#estimated velocity
                        sigma4=np.sqrt(np.nansum(np.prod([vect2,np.power(speed2-w4,2)],axis=0))/PT4)# spectral witdh
                        W_da.append(w4)
                        deal=np.concatenate([vec_null,vect2,vec_null])
                        if np.isnan(deal).all():
                            VectRe=vect*np.nan
                        else:
                                        
                            I=np.nanargmax(np.asarray(deal))
                            Indvel=speeddeal*(vect/vect)
                            VectRe,inewV=group(vect,I,3,Indvel)
                W[i]=w4
                
            else:
                if W[i]>(fNy*numberDoppler/2):
                    if W[i-1]>(fNy*numberDoppler/2):
                        vect2=vect[int(numberDoppler*5/4):int(numberDoppler*9/4)]
                        speed2=speeddeal[int(numberDoppler*5/4):int(numberDoppler*9/4)]
                        PT4=np.nansum(vect2)
                        w4=np.nansum(np.prod([vect2,speed2],axis=0))/PT4#estimated velocity
                        sigma4=np.sqrt(np.nansum(np.prod([vect2,np.power(speed2-w4,2)],axis=0))/PT4)# spectral witdh

                        W_da.append(w4)
                        deal=np.concatenate([Nul1,vect2,Nul2])
                        if np.isnan(deal).all():
                            VectRe=vect*np.nan
                        else:
                                
                            I=np.nanargmax(np.asarray(deal))
                            Indvel=speeddeal*(vect/vect)
                            VectRe,inewV=group(vect,I,3,Indvel)
                    else:
                        if W[i-1]<(fNy*numberDoppler*1/4):
                            vect2=vect[int(numberDoppler*3/4):int(numberDoppler*5/4)]
                            speed2=speeddeal[int(numberDoppler*3/4):int(numberDoppler*5/4)]
                            PT4=np.nansum(vect2)
                            w4=np.nansum(np.prod([vect2,speed2],axis=0))/PT4#estimated velocity
                            sigma4=np.sqrt(np.nansum(np.prod([vect2,np.power(speed2-w4,2)],axis=0))/PT4)# spectral witdh

                            W_da.append(w4)
                            deal=np.concatenate([Nul2,vect2,Nul1])
                            if np.isnan(deal).all():
                                VectRe=vect*np.nan
                            else:
                                        
                                I=np.nanargmax(np.asarray(deal))
                                Indvel=speeddeal*(vect/vect)
                                VectRe,inewV=group(vect,I,3,Indvel)
                        else:
                            vect2=vect[int(numberDoppler):int(numberDoppler*2)]
                            speed2=speeddeal[int(numberDoppler):int(numberDoppler*2)]
                            PT4=np.nansum(vect2)
                            w4=np.nansum(np.prod([vect2,speed2],axis=0))/PT4#estimated velocity
                            sigma4=np.sqrt(np.nansum(np.prod([vect2,np.power(speed2-w4,2)],axis=0))/PT4)# spectral witdh
                            W_da.append(w4)
                            deal=np.concatenate([vec_null,vect2,vec_null])
                            if np.isnan(deal).all():
                                VectRe=vect*np.nan
                            else:
                                        
                                I=np.nanargmax(np.asarray(deal))
                                Indvel=speeddeal*(vect/vect)
                                VectRe,inewV=group(vect,I,3,Indvel)
                    W[i]=w4
                else:
                    W_da.append(W[i])
                
                    deal=np.concatenate([vec_null,vect[int(numberDoppler):int(numberDoppler*2)],vec_null])
                    if np.isnan(deal).all():
                        VectRe=vect*np.nan
                    else:
                        
                        I=np.nanargmax(np.asarray(deal))
                        Indvel=speeddeal*(vect/vect)
                        VectRe,inewV=group(vect,I,3,Indvel)
                    



        ReVect.append(VectRe)
    

    return ReVect

def anchor(signal, weight):
    buffer = []
    last = signal[0]
    for i in signal:
        smoothed_val = last * weight + (1 - weight) * i
        buffer.append(smoothed_val)
        last = smoothed_val
    return buffer

def Inter1D(vector):
    y=np.asarray(vector)



    indx=np.argwhere(~np.isnan(y))



    if len(indx)>5:
        nou=[];noux=[];Nanx=[]

        for i in range(int(indx[-1]-indx[0])):
            if ~np.isnan(y[indx[0]+i]):
                
            
                nou.append(float(y[indx[0]+i]))
                noux.append(i)
            Nanx.append(i)
              

        y2= np.interp(Nanx, noux, nou)
        Inici=np.ones(int(indx[0]))*np.nan
        Fi=np.ones(len(y)-int(indx[-1]))*np.nan

        y3=np.concatenate((Inici,y2))

        y4=np.concatenate((y3,Fi))
    else:
        y4=np.copy(y)

    return y4

    
    
def PrepType(dm,nw):
##    convert the matrix dm, nw in  linear vector
    Nw=[];Dm=[]
    if np.isnan(dm).all() or np.isnan(nw).all():
        dm_axes=[]
        nw_axes=[]
        Matrix=[]
    else:
        
        for i in range(len(nw)):
            vector=nw[i]
            for j in range(len(vector)):
                if ~np.isnan(nw[i][j]) and ~np.isnan(dm[i][j]): 
                    Nw.append(nw[i][j])
                    Dm.append(dm[i][j])
        
        Dm=np.asarray(Dm,dtype=float)
        Nw=np.asarray(Nw,dtype=float)
        Dm2=[];Nw2=[];Y2=[]
        for i in range(len(Dm)):
            if ~np.isnan(Dm[i])and ~np.isnan(Nw[i]):
                y=6.3-1.6*Dm[i]
                ThuInd=Nw[i]-y#thurai(2016) index
                Dm2.append(round(Dm[i],2))
                Nw2.append(round(Nw[i],2))
                Y2.append(round(ThuInd,2))
                           

        Dm_axes=np.arange(np.min(Dm2),np.max(Dm2),.01)
        Nw_axes=np.arange(np.min(Nw2),np.max(Nw2),.01)
        dm_axes=[];nw_axes=[]
        for i in range(len(Dm_axes)):
            dm_axes.append(round(Dm_axes[i],2))
        for i in range(len(Nw_axes)):
            nw_axes.append(round(Nw_axes[i],2))

        Matrix=np.ones((len(dm_axes),len(nw_axes)))*np.nan
    ##    Create the matrix results Stra, Trans and Convective
        for i in range(len(dm_axes)):
            for j in range(len(Dm2)):

                if dm_axes[i]==Dm2[j]:

                    for k in range(len(nw_axes)):

                        if nw_axes[k]==Nw2[j]:

                            Matrix[i][k]=Y2[j]

    ##    Matrix is the matrix with the values de index Thurai(2016)
    ##now give the values from precypitation type where convective is 1, transition is 0 and stratiform is -1
        for i in range(len(Matrix)):
            for j in range(len(Matrix[i])):
                if ~np.isnan(Matrix[i][j]):
                    if abs(Matrix[i][j])<=0.3:
                        Matrix[i][j]=0.#transtition
                    else:
                        if Matrix[i][j]<-0.3:
                            Matrix[i][j]=-5.#stratiform
                        if Matrix[i][j]>0.3:
                            Matrix[i][j]=5.#convective
        

    return dm_axes,nw_axes,Matrix


def foundNonan(vector):#return the lats value that is not a nan
    a=np.asarray(vector)
    c=[]
    if np.isnan(a).all():
        c.append(2)
    else:
        for i in range(len(a)):
            nou=len(a)-i-1
                 
            if ~np.isnan(a[nou]):
                c.append(nou)
    return c[0]
    
def BB2(V,ZE,he,SK,KUR,last_bot,last_top,last_peak):#the input are fall speed, equivalent reflectivity and height 

    indNan_b=np.argwhere(~np.isnan(last_bot))
    if len(indNan_b)==0:
        last_bb_bot=np.nan
    else:
        last_bb_bot=last_bot[int(indNan_b[-1])]
    indNan_t=np.argwhere(~np.isnan(last_top))
    if len(indNan_t)==0:
        last_bb_top=np.nan
    else:
        last_bb_top=last_top[int(indNan_t[-1])]

    indNan_t=np.argwhere(~np.isnan(last_peak))
    if len(indNan_t)==0:
        last_bb_peak=np.nan
    else:
        last_bb_peak=last_peak[int(indNan_t[-1])]

    Deltah=he[3]-he[2]
    tiBB=[];h_BB_bot=[];h_BB_top=[];h_BB_peak=[]
    vector=ZE;vectorV=V
    firstV=foundNonan(vector)
    indexLimitBB=int(2000/Deltah)
    nouInd=3
    vect_2Sk=[]#it will the vector to find the top and bottom
    for i in range(len(SK)):
        if SK[i]>=0:
            vect_2Sk.append(1)
        if SK[i]<0:
            vect_2Sk.append(-1)
        if np.isnan(SK[i]):
            vect_2Sk.append(np.nan)
    
    
    if firstV>=8:
        FirstHe=firstV-3
    else:
        FirstHe=firstV
    if len(vector[:FirstHe])<2:
        Grad=np.nan*np.ones(len(vector[:FirstHe]))
        GradW=np.nan*np.ones(len(vectorV[:FirstHe]))
        DiffZ=np.nan*np.ones(len(vector[:FirstHe]))
        GradSk=np.nan*np.ones(len(SK[:FirstHe]))

    else:
        
        Grad=np.gradient(vector[:FirstHe])
        GradW=np.gradient(vectorV[:FirstHe])
        DiffZ=np.diff(vector[:FirstHe])
        GradSk=np.gradient(SK[:FirstHe])
    
    
    ######    Detection Cha 2009
    if np.isnan(Grad[nouInd:]).all() or np.nanargmax(Grad[nouInd:])==0 or np.nanargmin((Grad[nouInd:]))==0 or np.isnan(ZE[2])or np.nanmax(vector[nouInd:])<=15:#Ze less 15 dBZ is weak rain
        hBBtop=np.nan
        hBBbottom=np.nan
        hBBPEAK=np.nan
    else:
        
        idexTop=np.nanargmin((Grad[nouInd:]))


        if idexTop==1:
            if len(((Grad[nouInd:])))>3:
                nouInd=3
                idexTop=np.nanargmin((Grad[nouInd:]))
        PeakIndex=nouInd+np.nanargmax(vector[nouInd:])
        Hpeak=he[PeakIndex]
        Hbot=he[nouInd+np.nanargmax(Grad[nouInd:])]#avoid the firsts three heights gates
        Htop=he[nouInd+np.nanargmin((Grad[nouInd:]))]

        if Hbot>Htop and ~np.isnan(Grad[nouInd:idexTop]).all():
            Hbot=he[nouInd+np.nanargmax(Grad[nouInd:idexTop])]
        if Htop-Hbot>2000.:#bad election from idexTop, here establish the limit from BB in 2 km
            idexTop=np.nanargmin((Grad[nouInd:indexLimitBB+np.nanargmax(Grad[nouInd:idexTop])]))

        hBBbottom=np.nan
        hBBtop=np.nan
        hBBPEAK=np.nan
        if Hbot<Hpeak and Hpeak<Htop:#condition if the BB exists
            indexBot=np.nan;indexTop=np.nan
            if vect_2Sk[PeakIndex]>=0:#the BB exists
                for j in range(PeakIndex):
                    if vect_2Sk[PeakIndex-j]<0:
                        hBBbottom=he[PeakIndex-j]
                        indexBot=PeakIndex-j
                        break
                for j in range(len(he)-PeakIndex):
                    if vect_2Sk[PeakIndex+j]<0:
                        hBBtop=he[PeakIndex+j]
                        indexTop=PeakIndex+j
                        break
                if np.isnan(indexBot) or np.isnan(indexTop):
                    hBBPEAK=np.nan
                else:
                    
                    hBBPEAK=he[indexBot+np.nanargmax(SK[indexBot:indexTop])]
                    if hBBPEAK>hBBtop:
                        hBBPEAK=hBBtop
                    if hBBPEAK<hBBbottom:
                        hBBPEAK=hBBottom
                    

    if ~np.isnan(hBBbottom) and ~np.isnan(last_bb_bot):
        if abs(hBBbottom-last_bb_bot)>300:
            hBBbottom=np.nan
    if ~np.isnan(hBBtop) and ~np.isnan(last_bb_top):
        if abs(hBBtop-last_bb_top)>300:
            hBBtop=np.nan
    if ~np.isnan(hBBPEAK) and ~np.isnan(last_bb_peak):
        if abs(hBBPEAK-last_bb_peak)>300:
            hBBPEAK=np.nan
    if np.isnan(hBBPEAK) and ~np.isnan(hBBtop) and np.isnan(hBBbottom):
        hBBtop=np.nan
        
        
            
    h_BB_bot.append(hBBbottom)
    h_BB_top.append(hBBtop)
    h_BB_peak.append(hBBPEAK)

    return hBBbottom,hBBtop,hBBPEAK    



def group(a,indexcentral,Nnan,d):
    
    d=np.asarray(d)
    a=np.asarray(a)
    b=np.where(a>=0)
    acut=np.asarray(a[NbinsM:2*NbinsM])
    bcut=np.where(acut>=0)
    

    c=b[0]
    ccut=bcut[0]

    if indexcentral<=(NbinsM+(NbinsM/4)) or indexcentral>=(NbinsM+(3*NbinsM/4)):

        for i in range(np.size(b)-1):
            if c[i]-indexcentral<=0 and c[i+1]-indexcentral>=0:
                index=c[i+1]
                break
            else:
                index=indexcentral

        cond=1
        cont=0
        incr1=0#starts at 0

        while cond:
            if cont>=Nnan or (index)+incr1>=(len(a)-1):
                
                cond=0
            else:
            
                if np.isnan(a[index+incr1]):
                    cont+=1
                else:
                    cont=0


            incr1+=1

        cont=0;
        incr2=1#starta at 1
        
        cond=1
        while cond:
            if cont>=Nnan or index-incr2<=0:
                cond=0
            

            if np.isnan(a[index-incr2]):
                cont+=1
            
            else:
                cont=0


            incr2+=1
        vf2=np.copy(a);xf2=np.copy(d)
        vf2[0:index-incr2+1]=np.nan
        vf2[index+incr1:]=np.nan
        xf2[0:index-incr2+1]=np.nan
        xf2[index+incr1:]=np.nan
            

    else:
        
        for i in range(np.size(bcut)-1):
            if ccut[i]-indexcentral<=0 and ccut[i+1]-indexcentral>=0:
                index=ccut[i+1]
                break
            else:
                index=indexcentral-NbinsM

        cond=1
        cont=0
        incr1=0#starts at 0

        while cond:
            if cont>=Nnan or index+incr1>=len(acut)-1:
                cond=0
            else:
                
                if np.isnan(acut[index+incr1]):
                    cont+=1
                else:
                    cont=0


            incr1+=1

        cont=0;
        incr2=1#starts at 1
        
        cond=1
        while cond:
            if cont>=Nnan or index-incr2<=0:
                cond=0
            
            if np.isnan(acut[index-incr2]):
                cont+=1
            
            else:
                cont=0


            incr2+=1
        blanckv=np.nan*np.ones(len(acut))
        vf1=np.copy(acut);xf1=np.copy(d[NbinsM:2*NbinsM])
        vf1[0:index-incr2+1]=np.nan
        vf1[index+incr1:]=np.nan
        vff1=np.concatenate((blanckv,vf1))
        vf2=np.concatenate((vff1,blanckv))
                               
        xf1[0:index-incr2+1]=np.nan
        xf1[index+incr1:]=np.nan
        xff1=np.concatenate((blanckv,xf1))
        xf2=np.concatenate((xff1,blanckv))
        


    return vf2,xf2
    
    



def Process(matrix,he,temps,D,cte,neta,deltavel,code,Noi_spe_ref):#This function is the core of the signal processing
    matrix=np.asarray(matrix)
    lenHei=len(he)
    Doppler_bins=len(matrix[0])
    roW=10**6 #water density g/m3
    Cfact=2#value from cover factor, is the number multiplicate to sigma, Initially I considered as 2.
    ##Found the parameters dv in function of the height (mrr physics equation)
    dv=[]
    for i in range(len(he)):

        dv.append(1+3.68*10**-5*he[i]+1.71*10**-9*he[i]**2)
    if code==0:#from spectrum_raw
        etan=np.copy(matrix)
        etaN=np.multiply(etan,cte)
        etaV=etaN/deltavel#convert eta(n) in eta(v)
        
        
        
        Noise=np.multiply(neta,cte)#noise from eta (n) units m-1
        
        Snr=[]
        for j in range(len(Noise)):
            if np.isnan(Noise[j]):
                Snr.append(np.nan)
            else:
                Snr.append(10.*np.log10(np.nansum(etan[j])/neta[j]))
    if code==1:#from spectrum_reflectivity
        etaV=np.copy(matrix)
        Snr=np.copy(Noi_spe_ref)
        Noise=np.copy(Snr)*np.nan
    state=[]
    zewater=[];Ni=[];VT=[];Z=[];Z_da=[];Vhail=[]
    Vec_Deal=Aliasing(etaV,deltavel,he,temps)#function to avoid the aliasing#################

    PIAind=[]
    Z_pol_h=[];Z_all=[];RR_all=[];LWC_all=[];N_all=[];dm_all=[];nw_all=[]
    speeddeal=np.arange(-Nbins*fNy,2*Nbins*fNy,fNy)
    

    for m in range(len(etaV)): 
        vect1=Vec_Deal[m]
        vect2=vect1[Doppler_bins:2*Doppler_bins]

        
        proba=np.where(~np.isnan(vect2))
        leN=len(vect2)
        

        zewater.append(10**18*lamb**4*np.nansum(vect2)*deltavel/K2w)#Rayleight approach
        nde=[];vt=[];velHail=[]
        for n in range(len(etaV[0])):#Calculate the Ze from every gate without PIA                    
            value=6.18*vect2[n]*dv[m]*e**(-1*0.6*D[m][n])#D in mm
            
            value3=(9.65-10.3*e**(-1*0.6*D[m][n]))*dv[m]

            velHail.append(13.96*np.sqrt(10*D[m][n]))#vel from Hail Ulbrich and atlas 1982
            
            sbk=SigmaScatt[m][n]
            
            value2=10**6*(value/sbk)#(10**6)*(value/sbk)#N in m-3 mm-1
            
            nde.append(value2)#units mm-1 m-3 N(D)
            
            vt.append(value3)#terminal speed in function heigh and diameter
            
        Vhail.append(velHail)    
        VT.append(vt)    
        Ni.append(nde)
        ########    CALCULATE THE PIA INDEPENDENT THE STATE
        dif=[]#diference between diameters for N
        dif2=[]#diference between diameters for Z

        for n in range(len(D[m])):
            if n==0 or n==len(D[m])-1:
                if n==0:
                    dif2.append(D[m][n+1]-D[m][n])
                    dif.append(D[m][n+1]-D[m][n])
                if n==len(D[m])-1:
                    dif.append(abs(D[m][n-1]-D[m][n]))
                    dif2.append(abs(D[m][n]-D[m][n-1]))
            else:
                dif2.append(abs((D[m][n+1]-D[m][n])))
                dif.append(abs((D[m][n+1]-D[m][n-1]))/2.)
        Z_pol_h.append(10*np.log10(np.nansum(np.prod([nde,pow(np.asarray(D[m]),6),dif],axis=0))))

        DeltaAlt=he[3]-he[2]

        if m==0:
            PIAind.append(1.)
        else:
            Np=np.multiply(nde,PIAind[-1])
            Pro=[]
            for k in range(len(Np)):
                pro=SigmaExt[m][k]*Np[k]*dif[k]
                            
                Pro.append(pro)
                           
            kp=np.nansum(Pro)*10**-6
                        
            num=2.*kp*DeltaAlt

            N=-1.*np.multiply(Np,np.log(1-num)/num)
            Pro2=[]
            for k in range(len(N)):
                pro2=SigmaExt[m][k]*N[k]*dif[k]
                Pro2.append(pro2)
                            
            Kr=np.nansum(Pro2)*10**-6

            pia=PIAind[-1]*e**(-2.*Kr*DeltaAlt)

            if pia>=10. or num==0.:
                if num==0:
                    pia=np.nan
                else:
                    pia=10.

            PIAind.append(pia)
        vel=np.copy(speeddeal)
        PT=np.nansum(vect1)
        w=np.nansum(np.prod([vect1,vel],axis=0))/PT#estimated velocity

        value=np.nansum(np.prod([np.power(D[m],6),nde,dif2],axis=0))
        value2=np.nansum(np.prod([np.power(D[m],3),nde,dif2],axis=0))
        value3=np.nansum(np.prod([np.power(D[m],3),nde,dif2,w],axis=0))
        value4=np.nansum(np.prod([np.power(D[m],4),nde,dif2],axis=0))
                
        if np.nansum(nde)<=0.:
            N_all.append(np.nan)
        else:
            N_all.append(np.log10(np.nansum(nde)))
        

        if value<=0. or np.isnan(value):
            Z_all.append(np.nan)

        else:
            Z_all.append(10*np.log10(value))

        if value2<=0. or np.isnan(value2):
            LWC_all.append(np.nan)
        else:
            LWC_all.append(roW*value2*(np.pi/6.)*(10**-9))

        if value3<=0. or np.isnan(value3):
            RR_all.append(np.nan)
        else:
            RR_all.append(value3*(np.pi/6.)*(10**-9)*1000.*3600.)

        if value4<=0 or np.isnan(value2):
            dm_all.append(np.nan)
            nw_all.append(np.nan)
        else:
            dm_all.append(value4/value2)
            nw_all.append(np.log10(256.*(roW*value2*(np.pi/6.))/ (np.pi*roW*(value4/value2)**4)))#units m-3 mm-1
        

        

    z,lwc,rr,ze=Parameters(Ni,D,VT,0)
    z_f=np.copy(z)
    lwc_f=np.copy(lwc)
    rr_f=np.copy(rr)

    
    
    #stratiform case (M-P)
    vwaterR=2.65*np.power(zewater,.114)#values a,b from Atlas et al. 1973
    vsnowR=.817*np.power(zewater,.063)#values a,b from Atlas et al. 1973
    vwaterMiestr=2.65*np.power(ze,.114)
    #Thunderstorm Rain (S-S)
    vwaterMieconv=4.13*np.power(ze,.062)
    vwaterMie=np.nanmean([vwaterMieconv,vwaterMiestr],axis=0)
    
    speeddeal=np.arange(-Nbins*fNy,2*Nbins*fNy,fNy)
    NewM=[];state=[];mov=[];VerTur=[]
    W=[];Sig=[];Sk=[];lwc=[];rr=[];Z_da=[];SnowRate=[];N_da=[];Kurt=[];dm=[];nw=[];N_D_T=[]
    if np.sum(np.nanmean([vwaterR,vsnowR],axis=0))!=0:
        

        DealMatrix=[]
        
        for o in range(len(matrix)):
            ReVect=Vec_Deal[o]
            if np.isnan(ReVect).all():
                S=np.nan
                L=np.nan
                w3=np.nan
                sigma3=np.nan
                
            else:
                                        
                PT3=np.nansum(ReVect)
                w3=np.nansum(np.prod([ReVect,speeddeal],axis=0))/PT3#estimated velocity
                sigma3=np.sqrt(np.nansum(np.prod([ReVect,np.power(speeddeal-w3,2)],axis=0))/PT3)# spectral witdh
                    
                S=(w3-(dv[o]*vsnowR[o]))
                L=(w3-(dv[o]*vwaterMie[o]))

            if abs(S)<abs(L):#snow case

                if abs(S)<=(Cfact*abs(sigma3)) and abs(L)>(Cfact*abs(sigma3)):#case not liquid, possible snow

                    state.append(-10)#snow

                    if S<0 :
                        
                        mov.append(-1)#mvt upward
                    else:
                        mov.append(1)#mvt downpward
                else:
                    state.append(0)#mixed
                    mov.append(-1)
                    
            if abs(S)>=abs(L):#rain case
                        
                if abs(S)>(Cfact*abs(sigma3)) and abs(L)<=(Cfact*abs(sigma3)):#case liquid
                    

                    state.append(10)#rain

                    if L<0:
                        
                        mov.append(-1)
                    else:
                        mov.append(1)
                else:
                    state.append(0)#mixed
                    mov.append(1)
                    
            if np.isnan(L) and np.isnan(S):
                state.append(np.nan)
                mov.append(np.nan)



            if np.isnan(S) and ~np.isnan(L):#case liquid, but possible wrong election
                state.append(10)#rain
                if L<0:
                    mov.append(-1)
                else:
                    mov.append(1)


                        

                    
            if ~np.isnan(S) and np.isnan(L):#case not liquid, but possible wrong election
                state.append(-10)#snow
                if S<0:
                    mov.append(-1)
                else:
                    mov.append(1)
                        



                        
                        
            
            NewM.append(ReVect)
        

            

        for m in range(len(state)):#to avoid sporadic values
            if m!=0 and m!=len(state)-1:
                s1=state[m-1]
                s2=state[m]
                s3=state[m+1]
                if s2==0 and s1==-10 and s3==-10:
                    state[m]=-10
                if s2==0 and s1==10 and s3==10:
                    state[m]=10
                if s2==20 and s1==10 and s3==10:
                    state[m]=10
                if s2==20 and s1==-10 and s3==-10:
                    state[m]=-10
                if s2==-15 and s1==-10 and s3==-10:
                    state[m]=-10
        Mwater=[]
        Msnow=[]
        Mmixed=[]
        Mhail=[]
        MDriz=[]
        Munk=[]
        ZE=[]
        Mgrau=[]
        

        vel=np.copy(speeddeal)
        Nde=[]
        PIA=[]
############        create the vector diff Ze to detect drizzle
        Vector=[]
        for m in range(len(NewM)):
            vector=NewM[m]
            vector2=vector[Nbins:int(Nbins*2)]
            ValueZeD=(10**18*lamb**4*deltavel*np.nansum(vector2))/(np.pi**5*K2w)
            Vector.append(ValueZeD)
        Zediff=np.diff(Vector)
                      
                      

        for m in range(len(NewM)):

            vector=NewM[m]
            NullVector=NewM[m]*np.nan
            
            

            PT=np.nansum(vector)
            w=np.nansum(np.prod([vector,vel],axis=0))/PT#estimated velocity
            sigma=np.sqrt(np.nansum(np.prod([vector,np.power(vel-w,2)],axis=0))/PT)# spectral witdh
            sk=np.nansum(np.prod([vector,np.power(vel-w,3)],axis=0))/(PT*pow(sigma,3))# skewnes
            Kur=np.nansum(np.prod([vector,np.power(vel-w,4)],axis=0))/(PT*pow(sigma,4))# Kurtosis
            ValueZe=(10**18*lamb**4*deltavel*np.nansum(NewM[m]))/(np.pi**5*K2w)

            
            if ValueZe<=0 or np.isnan(ValueZe) or np.isinf(ValueZe):
                ZE.append(np.nan)
            else:
                ZE.append(10*np.log10(ValueZe))
            if w==0.:
                w=np.nan

            
            W.append(w)
            
            Sig.append(sigma)
            Sk.append(sk)
            Kurt.append(Kur)

            
            
            
            if state[m]==10:#rain case
                
                Mwater.append(NewM[m])
                SnowRate.append(np.nan)

                dif=[]#diference between diameters for N
                dif2=[]#diference between diameters for Z
                nde=[]
                indexFinded=[]
                for n in range(len(D[m])):
                    if n==0 or n==len(D[m])-1:
                        if n==0:
                            dif2.append(D[m][n+1]-D[m][n])
                            dif.append(D[m][n+1]-D[m][n])
                        if n==len(D[m])-1:
                            dif.append(abs(D[m][n-1]-D[m][n]))
                            dif2.append(abs(D[m][n]-D[m][n-1]))
                    else:
                        dif2.append(abs((D[m][n+1]-D[m][n])))
                        dif.append(abs((D[m][n+1]-D[m][n-1]))/2.)
                        condFH=speed[n]-w
                        

                        
                
                    EtaV=NewM[m][Nbins:int(Nbins*2)]#interval for water choosed
                    value=6.18*EtaV[n]*dv[m]*e**(-1*0.6*D[m][n])
                    s=SigmaScatt[m][n]
                    value2=(10**6)*(value/s)#N in m-3 mm-1
                    nde.append(value2)#units mm-1 m-3
                N_D_T.append(nde)
                #Calculate the diamater from the mean vel found
                diaWork=np.copy(D[m])
                

                if np.isnan(NewM[m][Nbins:int(Nbins*2)]).all():
                    diamHail=3
                else:
                        
                    diamHail=diaWork[np.nanargmax(NewM[m][Nbins:int(Nbins*2)])]#the diamHail is obtained with the max array position in eta over diameters                 
                    
                

                LastN=nde

                value=np.nansum(np.prod([np.power(D[m],6),LastN,dif2],axis=0))
                value2=np.nansum(np.prod([np.power(D[m],3),LastN,dif2],axis=0))
                value3=np.nansum(np.prod([np.power(D[m],3),LastN,dif2,w],axis=0))
                value4=np.nansum(np.prod([np.power(D[m],4),LastN,dif2],axis=0))
                
                if np.nansum(nde)<=0.:
                    N_da.append(np.nan)
                else:
                    N_da.append(np.log10(np.nansum(LastN)))

                if diamHail>=5:#Hail case
                    
                    Mhail.append(NewM[m])
                    state[m]=-20.
                    
                else:
                    Mhail.append(NullVector)
                

                if ~np.isnan(sk) :
                    

                    if m<len(Zediff):
                        if sk<=-0.5 and Zediff[m]>=1.:#New criteria from empric values. More information in--> https://doi.org/10.1175/JTECH-D-18-0158.1 (It's necessary to check again)

                            
                            MDriz.append(NewM[m])
                            state[m]=5.
                        else:
                            MDriz.append(NullVector)
                    else:
                        MDriz.append(NullVector)
                    
                else:
                    MDriz.append(NullVector)
                        
                        
                        
                if value<=0. or np.isnan(value):
                    Z_da.append(np.nan)
                else:
                    Z_da.append(10*np.log10(value))
                if value2==0.:
                    lwc.append(np.nan)
                else:
                    lwc.append(roW*value2*(np.pi/6.)*(10**-9))
                if value3==0.:
                    rr.append(np.nan)
                else:
                    rr.append(value3*(np.pi/6.)*(10**-9)*1000.*3600.)
                if value4==0:
                    dm.append(np.nan)
                    nw.append(np.nan)
                else:
                    dm.append(value4/value2)
                    nw.append(np.log10(256.*(roW*value2*(np.pi/6.))/ (np.pi*roW*(value4/value2)**4)))#units m-3 mm-1
                if mov[m]==-1:#case rain and upward
                      VerTur.append((2.6*np.power(ValueZe,.107))-w)
                if mov[m]==1:#case rain and downward
                      VerTur.append(w-(2.6*np.power(ValueZe,.107)))
                if np.isnan(mov[m]):#case rain and upward
                      VerTur.append(np.nan)
               
            else:
                Mwater.append(NullVector)
                
                
            if state[m]==-10:#Snow case
                
                Msnow.append(NewM[m])


                if ValueZe<=0 or np.isnan(ValueZe):
                    
                    SnowRate.append(np.nan)
                else:
                    if ~np.isnan(sk):

                        if sk>=0. and w>2.:
                            state[m]=-15.#graupel
                    
                    SnowRate.append(np.power(ValueZe/56.,1/1.2))#following Matrosov (2007) constants - https://link.springer.com/article/10.1007/s00703-011-0142-z#CR15
                Z_da.append(np.nan)
                lwc.append(np.nan)
                nw.append(np.nan)
                dm.append(np.nan)
                rr.append(np.nan)
                N_da.append(np.nan)
                N_D_T.append(np.ones(len(D[m]))*np.nan)
                if mov[m]==-1:#case snow and upward
                      VerTur.append((.817*np.power(ValueZe,.063))-w)
                if mov[m]==1:#case snow and downward
                      VerTur.append(w-(.817*np.power(ValueZe,.063)))
                if np.isnan(mov[m]):
                      VerTur.append(np.nan)
                    
                
                
            else:
                Msnow.append(NullVector)
                


            if state[m]==0:#Mixed case
                
                
                
                if Snr[m]<=15.:
                    state[m]=-10.#snow
                    SnowRate.append(np.power(ValueZe/56.,1/1.2))
                                        
                else:
                    if m<len(Zediff):
                        if ~np.isnan(sk):
                            if sk<=-0.5 and Zediff[m]>1.0:
                                state[m]=5.#drizzle
            ####################criteria from doi:10.5194/acp-16-2997-2016
                            if sk>=0. and w>2.:
                                state[m]=-15.#graupel
                    SnowRate.append(np.nan)
                
                        
                    
                Mmixed.append(NewM[m])
                Value=10**18*lamb**4*deltavel*np.nansum(NewM[m])/(np.pi**5*K2w)
                
                Z_da.append(np.nan)
                lwc.append(np.nan)
                dm.append(np.nan)
                nw.append(np.nan)
                rr.append(np.nan)
                N_da.append(np.nan)
                N_D_T.append(np.ones(len(D[m]))*np.nan)
                

                if mov[m]==-1:#case mixed and upward
                      VerTur.append((.817*np.power(ValueZe,.063))-w)
                if mov[m]==1:#case mixed and downward
                      VerTur.append(w-(.817*np.power(ValueZe,.063)))
                if np.isnan(mov[m]):
                      VerTur.append(np.nan)
                            
            else:
                Mmixed.append(NullVector)
                
            if np.isnan(state[m]):
                Mwater.append(NullVector)
                Msnow.append(NullVector)
                Mmixed.append(NullVector)
                Mhail.append(NullVector)
                MDriz.append(NullVector)
                Munk.append(NullVector)
                Mgrau.append(NullVector)
                Z_da.append(np.nan)
                lwc.append(np.nan)
                nw.append(np.nan)
                dm.append(np.nan)
                rr.append(np.nan)
                N_da.append(np.nan)
                N_D_T.append(np.ones(len(D[m]))*np.nan)
                SnowRate.append(np.nan)

                VerTur.append(np.nan)


                
            if state[m]==20:#cas unknown
                Munk.append(NewM[m])
                
                
                Value=10**18*lamb**4*deltavel*np.nansum(NewM[m])/(np.pi**5*K2w)#use the rayleight estimation
                
                Z_da.append(np.nan)
                lwc.append(np.nan)
                nw.append(np.nan)
                dm.append(np.nan)
                rr.append(np.nan)
                N_da.append(np.nan)
                N_D_T.append(np.ones(len(D[m]))*np.nan)


                SnowRate.append(np.nan)
                Nde.append(np.ones(shape=(len(matrix[1])))*np.nan)
                if mov[m]==-1:#case hail and upward, using the same for snow
                      VerTur.append((.817*np.power(ValueZe,.063))-w)
                if mov[m]==1:#case hail and downward, using the same for snow
                      VerTur.append(w-(.817*np.power(ValueZe,.063)))
                if np.isnan(mov[m]):
                      VerTur.append(np.nan)
                            
            else:
                Munk.append(NullVector)
                            
        for m in range(len(state)):
            if m!=0 and m!=len(state)-1:
                s1=state[m-1]
                s2=state[m]
                s3=state[m+1]
                if s2==-20 and s1==10 and s3==10:
                    state[m]=10        
                
                

        

    else:#There is not Signal
        
        blanck=np.nan*np.ones(shape=(len(matrix)))
        state=blanck
        NewM=(np.ones(shape=(len(he),len(speeddeal)))*np.nan)
        Z_da=blanck
        lwc=blanck
        nw=blanck
        dm=blanck
        rr=blanck
        N_da=blanck
        W=blanck
        Sig=blanck
        Sk=blanck
        Kurt=blanck
        SnowRate=blanck

        N_D_T=(matrix*np.nan)
        ZE=blanck
        mov=blanck
        VerTur=blanck
        Snr=blanck
        PIA=blanck
    
    return state,NewM,Z_da,lwc,rr,SnowRate,W,Sig,Sk,Noise,N_da,N_D_T,ZE,mov,VerTur,Snr,Kurt,PIAind,nw,dm,z_f,lwc_f,rr_f,Z_pol_h,Z_all,RR_all,LWC_all,dm_all,nw_all,N_all


       
    
def group_consecutives(vals, step=1):
    """Return list of consecutive lists of numbers from vals (number list)."""
    run = []
    result = [run]
    expect = None
    for v in vals:
        if (v == expect) or (expect is None):
            run.append(v)
        else:
            run = [v]
            result.append(run)
        expect = v + step
    return result

def Parameters(n,d,v,da):#the differences between diameter aren't constant
    Z=[];lwc=[];rr=[];ze=[]
    
    roW=10**6 #water density g/m3
    for i in range(len(n)):
        D=d[i]
        N=n[i]
        
        w=v[i]
        dif=[]
        
        
        if da==1:#in dealiased axes
            for m in range(len(D)):
                if m==0 or m==len(D)-1:
                    if m==0:
                        dif.append(d[i][m+1]-d[i][m])
                    if m==len(D)-1:
                        dif.append(abs(d[i][m-1]-d[i][m]))
                else:
                    if m<len(w)/2.:
                        dif.append((d[i][m+1]-d[i][m]))
                    else:
                        dif.append(abs((d[i][m]-d[i][m+1])))

        else:
            
            for m in range(len(D)):

                if m==0 or m==len(D)-1:
                    if m==0:
                        dif.append(d[i][m+1]-d[i][m])
                    if m==len(D)-1:
                        dif.append(d[i][m]-d[i][m-1])
                else:
                    dif.append((d[i][m+1]-d[i][m]))
        value=np.nansum(np.prod([np.power(D,6),N,dif],axis=0))
        value2=np.nansum(np.prod([np.power(D,3),N,dif],axis=0))
        value3=np.nansum(np.prod([np.power(D,3),N,dif,w],axis=0))
        if value==0. or np.isnan(value):
            Z.append(np.nan)
            ze.append(np.nan)
        else:
            Z.append(10*np.log10(value))
            ze.append(value)
        if value2==0. or np.isnan(value2):
            lwc.append(np.nan)
        else:
            lwc.append(roW*value2*(np.pi/6.)*(10**-9))
        if value3==0. or np.isnan(value3):
            rr.append(np.nan)
        else:
            rr.append(value3*(np.pi/6.)*(10**-9)*1000.*3600.)

    return Z,lwc,rr,ze


def CorrectWithBBMatrix(estat,Z,Lwc,Rr,SnowRate,Hcolum,bb_bot,bb_top,Ze,z,lwc,rr,sk):#return a matrix with all corrections
    for j in range(len(estat)):
        
    
        for i in range(len(Hcolum)):
            
            if ~np.isnan(bb_bot[j]):
                if (estat[j][i]==-10 or estat[j][i]==-15 or estat[j][i]==0) and Hcolum[i]<bb_bot[j]:#snow, graupel or mixed below the bb is not possible, so exist an upward
                    if sk[j][i]<-.5:
                        estat[j][i]=5.
                    else:
                        estat[j][i]=10.# change to rain, must change the lwr and rr, and Z--> write code
                    SnowRate[j][i]=np.nan
                    Z[j][i]=z[j][i]
                    Lwc[j][i]=lwc[j][i]
                    Rr[j][i]=rr[j][i]
                                    
            if ~np.isnan(bb_top[j]):
                if (estat[j][i]==10 or estat[j][i]==5) and Hcolum[i]>bb_top[j]:#rain beyond the bb is not possible, so exists downward
                    if sk[j][i]>0:
                        estat[j][i]==-15
                    else:
                        estat[j][i]==-10# change to rain, must change the lwr and rr, and Z--> write code
                    Z[j][i]=np.nan
                    Lwc[j][i]=np.nan
                    Rr[j][i]=np.nan

                    if np.isnan(Ze[j][i]):
                        SnowRate[j][i]=np.nan
                    else:
                        ValueZe=np.power(10,Ze[j][i]/10.)
                        SnowRate[j][i]=(np.power(ValueZe/56.,1/1.2))
                    

    

    return estat,Z,Lwc,Rr,SnowRate






def ScatExt(diameter,longW):#for 1 height gate
    ag_lam=longW*1000.#mm Convert lamb from m to mm

    ag_mre=6.417
    ag_mim=2.758
    m = ag_mre + 1.0j * ag_mim
    
    scatt=[];extinct=[]
    for i in range(len(diameter)):
        r=diameter[i]/2.

        if np.isnan(r):
            scatt.append(np.nan)
            extinct.append(np.nan)

        else:
            
            x = 2*np.pi*r/ag_lam;#must be non dimensional, and the units from ag_lam and r must be the same
            qext, qsca, qback, g = mp.mie(m,x)
            absorb  = (qext - qsca) * np.pi * r**2
            scatt.append(qsca * np.pi * r**2)
            extinct.append(qext* np.pi * r**2)
    return scatt,extinct
            
def MrrProNoise2(vector,he,DifRange,limitTime):
    v1=np.copy(vector)
    v2=np.copy(vector)

    meanv=np.nanmean(v1)
    varv=np.nanvar(v1)
    quo=meanv**2/varv
    rate=np.nanmax(v1)/np.nanmean(v1)
    ValLim=limitTime

    if np.isnan(vector).all() or (quo>ValLim and rate<1.3):
        soroll=np.nan
        v10=np.ones(len(vector))*np.nan
    else:
        condition=1
        while condition:
            quo2=np.nanmean(v1)**2/np.nanvar(v1)
            rate2=np.nanmax(v1)/np.nanmean(v1)

            if (quo2>ValLim and rate2<1.3) or np.nanvar(v1)==0.:
                if np.nanvar(v1)==0.:
                    soroll=np.nanmin(vector)
                    
                else:
                    number=len(vector)-np.count_nonzero(~np.isnan(v1))
                    soroll=np.nanmean(v1)
                condition=0
            np.put(v1,np.nanargmax(v1),np.nan)
        v8=vector-soroll
        v9=np.ones(len(v8))
        v9[v8<(np.nanmax(v8)*0.1)]=np.nan
        
        v10=np.multiply(v9,v8)
    return v10,soroll



def Continuity(vector,matrix,deltaH):#vector is the noise and matrix is the signal
    if deltaH<100:
        for i in range(len(vector)-4):
            
            elem=np.ones(len(matrix[i]))
            s1=vector[i]
            s2=vector[i+1]
            s3=vector[i+2]
            s4=vector[i+3]
            
            if np.isnan(s1) and ~np.isnan(s2) and ~np.isnan(s3) and np.isnan(s4):
                
                vector[i+1]=np.nan
                matrix[i+1]=elem*np.nan
                vector[i+2]=np.nan
                matrix[i+2]=elem*np.nan
            if np.isnan(s1) and ~np.isnan(s2) and np.isnan(s3) and np.isnan(s4):
                
                vector[i+1]=np.nan
                matrix[i+1]=elem*np.nan
                vector[i+2]=np.nan
                matrix[i+2]=elem*np.nan
            if np.isnan(s1) and np.isnan(s2) and ~np.isnan(s3) and np.isnan(s4):
                
                vector[i+1]=np.nan
                matrix[i+1]=elem*np.nan
                vector[i+2]=np.nan
                matrix[i+2]=elem*np.nan

    else:#delete spirous
        
        for i in range((len(vector)-3)):
            elem=np.ones(len(matrix[i]))
            s1=vector[i]
            s2=vector[i+1]
            s3=vector[i+2]
            if np.isnan(s1) and ~np.isnan(s2) and np.isnan(s3):
                vector[i+1]=np.nan
                matrix[i+1]=elem*np.nan


    return vector,matrix


def date2unix(date):
    return calendar.timegm(date.timetuple())
def unix2date(unix):
    return datetime.datetime.utcfromtimestamp(unix)


np.warnings.filterwarnings('ignore')#to avoid the error messages

########INCLUDE THE OPTIONS IN EXECUTATION
if len(sys.argv)==1:
    option=0
c_opt=0;c1=0;c2=0;c3=0;h0_opt=np.nan
if len(sys.argv)>1:
    for i in sys.argv:
   
        if i=='-spe3D':
            #print('Your chosen option is to save the corrected spectral reflectivity values\n')
            option=1
            c_opt+=1
            c1=1
        
        if i=='-dsd3D':
            #print('Your chosen option is to save the 3D DSD\n')
            option=2
            c_opt+=1
            c2=1
        if i[0:2]=='-h':
            #print('The first height has been changed\n')
            h0_opt=float(i[2:])
            option=0
            c_opt+=1
            c3=1
    if c_opt!=len(sys.argv)-1:

        print('Please check the syntaxis, possible command line arguments available are (more than one is possible, in any order):\n -spe3D -dsd3D -hxxx \nwhere\n-spe3D : saves corrected spectral reflectivity values\n-dsd3D : saves 3D DSD\n -hxxx : forces the antenna height is at xxx meters above sea level.\n')
        sys.exit()
    if c1==1 and c2==1:
        option=3
    if c1==1 and c2==0 and c3==0:
        print('\nYour chosen option is to save the corrected spectral reflectivity values\n')
    if c1==0 and c2==1 and c3==0:
        print('\nYour chosen option is to save the 3D DSD\n')
    if c1==0 and c2==0 and c3==1:
        print('\nThe antenna height has been changed\n')
    if c1==0 and c2==1 and c3==1:
        print('\nThe first height has been changed and the 3D DSD is saved\n')
    if c1==1 and c2==0 and c3==1:
        print('\nThe antenna height has been changed and the corrected spectral reflectivity is saved\n')
    if c1==1 and c2==1 and c3==1:
        print('\nThe antenna height has been changed and the corrected spectral reflectivity and 3D DSD are saved\n')
    if c1==1 and c2==1 and c3==0:
        print('\nThe corrected spectral reflectivity and 3D DSD are saved\n')



    

##########option=1 #the spectra afetr noise ande dealiasing is saved
##########option=2 #the 3D_DSD is saved, but the user should create the matrix of dropsize to each height
##########option=3 #the two option are applied



print('Insert the path of directory containing the netcdf file(s) to be processed. For instance d:\MrrProdata/')
Root=input()  #input from the user 
os.chdir(Root)     

folder=Root
dircf=glob.glob(Root+'*.nc')
dircf=np.sort(dircf)
Ndir=[]
for i in dircf:
    if i[-12:-3]!='processed':
        Ndir.append(i)
dircf=np.copy(Ndir)

##get the main parameters
TimeInteFiles=[];Cons=[];TimeFinal=np.nan
Count_Cons=0;Cons.append(0)
for i in dircf:
    
    f=nc4.Dataset(i,'r')
    for k in f.variables:
        if k=='spectrum_reflectivity':
            Code_spectrum=1
        if k=='spectrum_raw':
            Code_spectrum=0
    if Code_spectrum==0:
        raw=f.variables['spectrum_raw'][:,:]#log_attenuated_power  in dB
        Ntime,NheiM,NbinsM=np.shape(raw)
    if Code_spectrum==1:
        raw2=f.variables['spectrum_reflectivity'][:,:]#log_attenuated_power  in dB in 1m2/m3
        Ntime,NheiM,NbinsM=np.shape(raw2)
    
    Latitude=f.variables['latitude'][:]
    Longitude=f.variables['longitude'][:]
    #print(Latitude,Longitude)

    TF=f.variables['transfer_function'][:]
    CC=f.variables['calibration_constant'][:]
    Range=f.variables['range'][:]
    TimeVector=f.variables['time'][:]
    for k in f.ncattrs():
        
        if k=='title':
            new_title='Processed data from '+str(f.getncattr(k))
        
        if k=='instrument_name':
            new_Instr='Processed data from '+str(f.getncattr(k))
        if k=='site_name':
            Site=str(f.getncattr(k))

        
    
    f.close()
    if ~np.isnan(TimeFinal):
        valor=TimeVector[0]-TimeFinal

        if valor<3600:#the file is consecutive supposing that file are hourly
            Cons.append(1)#1 is consecutive
        else:
            Cons.append(0)#0 is non consecutive
    else:
        Cons.append(0)
    
    if len(TimeVector)>1:
        TimeInteFiles.append(round(np.min(np.diff(TimeVector)),0))
        TimeFinal=TimeVector[-1]
    else:
        TimeFinal=np.nan
FullIntervalTime= round(np.min(TimeInteFiles),0)#minimum time resolution of all files
if FullIntervalTime==1 or FullIntervalTime==5 or FullIntervalTime==10 or FullIntervalTime==60:
    TimeInt=FullIntervalTime
else:
    print('The time resolution found in files is not the usual one. Please insert the time resolution of the files (in seconds), for instance: 10')
    val=input()  #input from the user
    FullIntervalTime=int(val)
    TimeInt=FullIntervalTime

######################NOTE THAT THE RANGE FROM NETCDF STARTS AT THE FIRST HEIGHT GATE
DeltaH=Range[3]-Range[2]
if np.isnan(h0_opt):
    Hcolum=Range
else:
    Hcolum=np.arange(h0_opt,len(Range)*DeltaH+h0_opt,DeltaH)

FTcolum=np.asarray(TF)
C=CC



#features from radar Mrr Pro

velc=299792458.#light speed 
lamb=velc/(24.23*1e9)  #The frequency of radar is 24.23 GHz
ag_lam=lamb
fsampling=500000#Hz valor samplig frequency
fNy=fsampling*lamb/(2*2*NheiM*NbinsM) 
K2w=0.92
Deltaf=fsampling/(2*NheiM*NbinsM) 
CetNtoetaV=2./(Deltaf*lamb)
Deltav=Deltaf*lamb/2.
AveSpectre=TimeInt*fsampling/(2*NheiM*NbinsM)


##Found the parameters dv in function of the height (mrr physics equation)
dv=[]
for i in range(len(Hcolum)):
    
    dv.append(1+3.68*10**-5*Hcolum[i]+1.71*10**-9*Hcolum[i]**2)

speed=np.arange(0,NbinsM*fNy,fNy)
speed21=np.arange(0,NbinsM*fNy/2,fNy)
speed22=np.arange(-1*NbinsM*fNy/2,0,fNy)
speed2=np.concatenate((speed21,speed22),axis=0)
speed3=np.arange(-1*NbinsM*fNy,2*NbinsM*fNy,fNy)

##Found the diameters in function of height and speed
D=[]
for i in range(len(dv)):
    d=[]
    for j in range(len(speed)):
            
        b=speed[j]/dv[i]
    
        if b>=0.002 and b<=9.37:#Condition of diameter is good for 0.109 mm< D< 6 mm
            d.append(np.log((9.65-b)/10.3)*(-1/0.6))
        else:
            d.append(np.nan)
    D.append(d)



    
##Found the scatter and extint sections in function of height and speed
SigmaScatt=[]
SigmaExt=[]
for i in range(len(D)):#entry in height
    sig1,sig2=ScatExt(D[i],lamb)
    SigmaScatt.append(sig1)
    SigmaExt.append(sig2)

##start the processing
print(datetime.datetime.now())
countfile=0
CountNumberofFile=0;Count_Cons=0
BB_last_bot=np.nan;BB_last_peak=np.nan;BB_last_top=np.nan
for i in dircf:
    condTime=0
    if len(dircf)==0:
        print('There are no netcdf files in this folder! Please check your folder')
        break
    else:
        if len(dircf)>=2 and CountNumberofFile==0:
            print(' There are '+str(len(dircf))+' netcdf files in this folder')
            CountNumberofFile+=1
        else:
            if CountNumberofFile==0:
                print('There is only 1 netcdf file in the folder')
        print('Processing file ',i)
    ##UPLOAD ALL RAW VALUES
    Nw_2=[];Dm_2=[];NewMatrix_full=[]
    Raw=[];Time=[]
    EmptyRaw=np.asarray(np.ones([NheiM,NbinsM])*np.nan)


    
    
    
    f=nc4.Dataset(i,'r')
    for k in f.variables:
        if k=='spectrum_reflectivity':
            Code_spectrum=1
        if k=='spectrum_raw':
            Code_spectrum=0
    if Code_spectrum==0:
        raw=f.variables['spectrum_raw'][:,:]#log_attenuated_power  in dB
        Ntime,Nhei,Nbins=np.shape(raw)
        Snr_Refl=[]
        
    if Code_spectrum==1:
        raw2=f.variables['spectrum_reflectivity'][:,:]#log_attenuated_power  in dB in 1m2/m3
        Ntime,Nhei,Nbins=np.shape(raw2)
        Snr_Refl=f.variables['SNR'][:,:]

    t=f.variables['time'][:]
    Range=f.variables['range'][:]
    TF=f.variables['transfer_function'][:]
    CC=f.variables['calibration_constant'][:]
    
    f.close()
    TimeInt=FullIntervalTime
    
    DeltaH=Range[3]-Range[2]
    if np.isnan(h0_opt):
        Hcolum=Range
    else:
        Hcolum=np.arange(h0_opt,len(Range)*DeltaH+h0_opt,DeltaH)
    FTcolum=np.asarray(TF)
    C=CC



    #features from radar Mrr Pro
    

    velc=299792458.#light speed 
    lamb=velc/(24.23*1e9)  #La frequency of radar is 24.23 GHz
    ag_lam=lamb
    fsampling=500000#Hz valor samplig frequency
    fNy=fsampling*lamb/(2*2*Nhei*Nbins) 
    K2w=0.92
    Deltaf=fsampling/(2*Nhei*Nbins) 
    CetNtoetaV=2./(Deltaf*lamb)
    Deltav=Deltaf*lamb/2.
    AveSpectre=TimeInt*fsampling/(2*Nhei*Nbins)


    ##Found the parameters dv in function of the height (mrr physics equation)
    dv=[]
    for j in range(len(Hcolum)):
        
        dv.append(1+3.68*10**-5*Hcolum[j]+1.71*10**-9*Hcolum[j]**2)

    speed=np.arange(0,NbinsM*fNy,fNy)
    speed21=np.arange(0,NbinsM*fNy/2,fNy)
    speed22=np.arange(-1*NbinsM*fNy/2,0,fNy)
    speed2=np.concatenate((speed21,speed22),axis=0)
    speed3=np.arange(-1*NbinsM*fNy,2*NbinsM*fNy,fNy)

    ##Found the diameters in function of height and speed
    D=[]
    for k in range(len(dv)):
        d=[]
        for j in range(len(speed)):
                
            b=speed[j]/dv[k]
        
            if b>=0.002 and b<=9.37:#Condition of diameter is good for 0.109 mm< D< 6 mm
                d.append(np.log((9.65-b)/10.3)*(-1/0.6))
            else:
                d.append(np.nan)
        D.append(d)



        
    ##Found the scatter and extint sections in function of height and speed
    SigmaScatt=[]
    SigmaExt=[]
    for j in range(len(D)):#entry in height
        sig1,sig2=ScatExt(D[j],lamb)
        SigmaScatt.append(sig1)
        SigmaExt.append(sig2)

    

    if (3600/TimeInt)-len(t)!=0:
        print('Attention, there are ',(3600/TimeInt)-len(t),' time steps missing in this file')
    if (3600/TimeInt)-len(t)>=353:
        print('The file has not enought values, less of 1 minute')
        
    #Ntime,Nhei,Nbins=np.shape(raw)
    if NheiM!=Nhei or NbinsM!=Nbins:
        print('Attention, there are netcdf files with different configuration!')

    count=0
    Min_strats=unix2date(int(round(t[0]))).minute
    Hour_strats=unix2date(int(round(t[0]))).hour
    Sc_strats=unix2date(int(round(t[0]))).second
    DD_strats=unix2date(int(round(t[0]))).day
    MM_strats=unix2date(int(round(t[0]))).month
    YY_strats=unix2date(int(round(t[0]))).year

    if Min_strats!=0 or Sc_strats!=0:
        print('Attention, this file started late ')
        datIni=datetime.datetime(year = YY_strats, month = MM_strats, day = DD_strats, hour = Hour_strats, minute = 0, second = 0)
        DatIni=date2unix(datIni)
    else:
        DatIni=int(round(t[0]))#date2unix(datIni)
    print('started at',unix2date(DatIni))
    


    TimeCorr=np.arange(DatIni,DatIni+(3600),TimeInt)#The files are of one hour
    t=np.asarray(t,dtype='int')
##    LOOP TO CORRECT THE GAPS EQUIPEMNT
    ind=0#index for temps
    for j in range(len(TimeCorr)):
        
        if TimeCorr[j]<=t[ind] and ind<len(t)-1:    
            if TimeCorr[j]==t[ind]:
                Time.append(t[ind])
                if Code_spectrum==0:
                    Raw.append(raw[ind])
                if Code_spectrum==1:
                    Raw.append(raw2[ind])
                ind+=1

            else:
                
                Time.append(TimeCorr[j])
                Raw.append(EmptyRaw)
        else:
            Time.append(TimeCorr[j])
            Raw.append(EmptyRaw)
        
        
        
    dataset=Dataset(i[:-3]+'-processed'+'.nc','w',format='NETCDF4')
    dataset.description=new_title
    dataset.author='Albert Garcia Benad'+u'\xed'
    dataset.orcid='0000-0002-5560-4392 '
    dataset.instrument_details=new_Instr
    if not Site:
        dataset.site='Undefined'
    else:
        dataset.site=Site
    if not Latitude:
        dataset.latitude='Undefined'

    else:
        dataset.latitude=Latitude
    if not Longitude:
        dataset.longitude='Undefined'
    else:
        dataset.longitude=Longitude

    dataset.createDimension('DropSize',len(D[0]))

    dataset.createDimension('Height',len(Hcolum))

    dataset.createDimension('Speed',len(speed3))

    dataset.createDimension('PIA_Height',len(Hcolum)+1)

    dataset.createDimension('BB_Height',1)

    dataset.createDimension('time',None)
    dataset.createDimension('time_utc',None)
    
    nc_times=dataset.createVariable('time','float64',('time',))
    nc_Format_times=dataset.createVariable('time_utc', 'float64', ('time_utc',))
    
    
    nc_ranges_H=dataset.createVariable('Height','f',('Height',))
    nc_ranges_V=dataset.createVariable('Speed','f',('Speed',))
    nc_ranges_DropSize=dataset.createVariable('DropSize','f',('DropSize',))
    nc_ranges_H_PIA=dataset.createVariable('PIA_Height','f',('PIA_Height',))
    nc_ranges_H_BB=dataset.createVariable('BB_Height','f',('BB_Height',))

    nc_times.units = 'UNIX Time Stamp, SECONDS SINCE 1970-01-01'
    nc_times.description='Time in unix format'


    nc_Format_times.units='seconds since 1970-01-01'
    nc_Format_times.calendar='standard'
    nc_Format_times.decription='time UTC'

    nc_ranges_H.units = 'm'
    nc_ranges_H.description = 'Heights in meters'

    nc_ranges_V.units = 'm/s'
    nc_ranges_V.description = 'speed without aliasing'

    nc_ranges_DropSize.units = 'mm'
    nc_ranges_DropSize.description = 'Size of the water drops'

    nc_ranges_H_PIA.units = 'm'
    nc_ranges_H_PIA.description = 'Heights in meters a.s.l.'

    nc_ranges_H_BB.units = 'm'
    nc_ranges_H_BB.description = 'Heights in meters a.s.l.'
                    
      
    nc_ranges_H[:]=np.array(Hcolum,dtype='f4')
    
    nc_ranges_V[:]=np.array(speed3,dtype='f4')
        
    #constant value to conevrt S/TF to eta(n)
    Cte=DeltaH*C/10**20

    PotCorrSum=[]
    Timecount=0
    countwork=0

    bb_bot_full=[];bb_top_full=[];bb_peak_full=[]
    for i in range(len(Time)):
        tst1=(datetime.datetime.now()).timestamp()
        
        Pot=[];NewNoise=[]

        
        if Code_spectrum==0:
            for k in range(len(Hcolum)):
                COL=np.asarray(Raw[i][k])#len of 64 doppler bins
                COL=np.power(10,COL/10.)#convert dB power to units unknow
                COL2,Noise=MrrProNoise2(COL,k,DeltaH,TimeInt)
                
                NewNoise.append(Noise*(k)**2/FTcolum[k])
                Pot.append((COL2*((k)**2))/FTcolum[k])#col is the value s(n,i) spectrum power in dB so change the dB power to units
            Snr_Refl_2=[]
        if Code_spectrum==1:
            for k in range(len(Hcolum)):
                COL=np.asarray(Raw[i][k])#len of 64 doppler bins
                Pot.append(np.power(10,COL/10.))#convert dB power to units m2/m3 and the noise has been substracted
            if i<=len(Snr_Refl)-1:
                Snr_Refl_2=Snr_Refl[i]
            else:
                Snr_Refl_2=Snr_Refl_2*np.nan


        if len(Pot)==0:
            proeta=np.ones([NheiM,NbinsM])*np.nan
        else:
            
            NewNoise,Pot=Continuity(NewNoise,Pot,DeltaH)#delete puntual values, considered bad signals
            proeta=Pot

        nc_times[Timecount:Timecount+1]=Time[i]#
        nc_Format_times[Timecount:Timecount+1]=date2num(unix2date(Time[i]),units=nc_Format_times.units,calendar=nc_Format_times.calendar)
        
        nc_ranges_DropSize[:]=np.array(np.ma.masked_invalid(D[0]),dtype='f')
        ncShape2D = ('time','Height',)
        ncShape2D_utc = ('time_utc','Height',)
        ncShape2D_BB = ('time_utc','BB_Height',)
        ncShape3D = ('time_utc','Height','DropSize',)
        ncShape3D2 = ('time_utc','Height','Speed',)
        if Timecount==0:

            if option==1:
                nc_etaV=dataset.createVariable('spe_3D','f',ncShape3D2)
                nc_etaV.description='Spectral reflectivity Distribution in function of time and height after noise and dealisaing'
                nc_etaV.units=' mm-1'
            if option==2:
                nc_3D=dataset.createVariable('dsd_3D','f',ncShape3D)
                nc_3D.description='Drop Size Distribution in function of time, height and diameter'
                nc_3D.units='log10(m-3 mm-1)'
            if option==3:
                nc_etaV=dataset.createVariable('spe_3D','f',ncShape3D2)
                nc_etaV.description='Spectral reflectivity Distribution in function of time and heightafter noise and dealisaing'
                nc_etaV.units=' mm-1'

                nc_3D=dataset.createVariable('dsd_3D','f',ncShape3D)
                nc_3D.description='Drop Size Distribution in function of time, height and diameter'
                nc_3D.units='log10(m-3 mm-1)'
                
                
            
                    ##################create the netcdf############

            nc_w=dataset.createVariable('W','f',ncShape2D_utc)
            nc_w.description='fall speed with aliasing correction'
            nc_w.units='m s-1'
                    

            nc_sig=dataset.createVariable('spectral width','f',ncShape2D_utc)
            nc_sig.description='spectral of the spectral reflectivity width with dealiasing'
            nc_sig.units='m s-1'
                    

            nc_sk=dataset.createVariable('Skewness','f',ncShape2D_utc)
            nc_sk.description='skewness of the spectral reflectivity with dealiasing'
            nc_sk.units='none'

            nc_kur=dataset.createVariable('Kurtosis','f',ncShape2D_utc)
            nc_kur.description='Kurtosis of the spectral reflectivity with dealiasing'
            nc_kur.units='none'

            nc_PIA=dataset.createVariable('DBPIA','f',ncShape2D_utc)
            nc_PIA.description='Path Integrated Attenuation, expressed as 10*log(PIA) without type hydrometeor consideration'
            nc_PIA.units='dB'

            nc_state=dataset.createVariable('Type','f',ncShape2D_utc)
            nc_state.description='Indicate the type from hydrometeor as unknown (20), rain (10), drizzle (5), mixed (0), snow (-10), graupel (-15) and hail (-20)'
            nc_state.units=''

            nc_LWC_all=dataset.createVariable('LWC_all','f',ncShape2D_utc)
            nc_LWC_all.description='Liquid water content suposing that all hidrometeors are in liquid phase'
            nc_LWC_all.units='g m-3'
                        

            nc_RR_all=dataset.createVariable('RR_all','f',ncShape2D_utc)
            nc_RR_all.description='Rain Rate suposing that all hidrometeors are in liquid phase'
            nc_RR_all.units='mm hr-1'
                    

            nc_LWC=dataset.createVariable('LWC','f',ncShape2D_utc)
            nc_LWC.description='Liquid water content'
            nc_LWC.units='g m-3'
                    

            nc_RR=dataset.createVariable('RR','f',ncShape2D_utc)
            nc_RR.description='Rain Rate'
            nc_RR.units='mm hr-1'

            nc_SnowR=dataset.createVariable('SR','f',ncShape2D_utc)
            nc_SnowR.description='Snow Rate'
            nc_SnowR.units='mm hr-1'

            nc_Z_DA=dataset.createVariable('Za','f',ncShape2D_utc)
            nc_Z_DA.description='Attenuated Reflectivity considering only liquid drops'
            nc_Z_DA.units='dBZ'

            nc_Z_all=dataset.createVariable('Z_all','f',ncShape2D_utc)
            nc_Z_all.description='Reflectivity suposing that all hidrometeors are in liquid phase'
            nc_Z_all.units='dBZ'
                    

            nc_Z_da=dataset.createVariable('Z','f',ncShape2D_utc)
            nc_Z_da.description='Reflectivity considering only liquid drops'
            nc_Z_da.units='dBZ'

            nc_Z_e=dataset.createVariable('Ze','f',ncShape2D_utc)
            nc_Z_e.description='Attenuated Equivalent Reflectivity'
            nc_Z_e.units='dBZ'

            nc_Z_ea=dataset.createVariable('Zea','f',ncShape2D_utc)
            nc_Z_ea.description='Equivalent Reflectivity attenuated'
            nc_Z_ea.units='dBZ'

            #nc_VerMov=dataset.createVariable('Vmov','f',ncShape2D_utc)
            #nc_VerMov.description='Verical movement +1 downward -1 upward'
            #nc_VerMov.units='None'
                    
            nc_N_da=dataset.createVariable('N(D)','f',ncShape2D_utc)
            nc_N_da.description='Drop Size Distribution'
            nc_N_da.units='log10(m-3 mm-1)'

            nc_N_all=dataset.createVariable('N(D)_all','f',ncShape2D_utc)
            nc_N_all.description='Drop Size Distribution suposing that all hidrometeors are in liquid phase'
            nc_N_all.units='log10(m-3 mm-1)'

                    
            nc_SNR=dataset.createVariable('SNR','f',ncShape2D_utc)
            nc_SNR.description='Signal noise relation from signal without dealiasing'
            nc_SNR.units='dB'

            nc_Noi=dataset.createVariable('Noise','f',ncShape2D_utc)
            nc_Noi.description='Noise'
            nc_Noi.units='m-1'

            nc_nw=dataset.createVariable('Nw','f',ncShape2D_utc)
            nc_nw.description='Intercept parameter of the gamma distribution normalized to the liquid water content'
            nc_nw.units='log10(mm-1 m-3)'

            nc_dm=dataset.createVariable('Dm','f',ncShape2D_utc)
            nc_dm.description='mean mas-weighted raindrop diameter'
            nc_dm.units='mm'

            nc_NW=dataset.createVariable('Nw_all','f',ncShape2D_utc)
            nc_NW.description='Intercept parameter of the gamma distribution normalized to the liquid water content suposing that all hidrometeors are liquid pahse'
            nc_NW.units='log10(mm-1 m-3)'

            nc_DM=dataset.createVariable('Dm_all','f',ncShape2D_utc)
            nc_DM.description='mean mas-weighted raindrop diameter  content suposing that all hidrometeors are liquid phase'
            nc_DM.units='mm'


            #nc_VelTur=dataset.createVariable('Fall speed variability','f',ncShape2D_utc)
            #nc_VelTur.description='Estimate the fall speed variability'
            #nc_VelTur.units='m/s'

            nc_bb_bot=dataset.createVariable('BB_bottom','f',ncShape2D_BB)
            nc_bb_bot.description='height from BB bottom in meters a.g.l.'
            nc_bb_bot.units='m'

            nc_bb_top=dataset.createVariable('BB_top','f',ncShape2D_BB)
            nc_bb_top.description='height from BB Top in meters a.g.l.'
            nc_bb_top.units='m'

            nc_bb_peak=dataset.createVariable('BB_peak','f',ncShape2D_BB)
            nc_bb_peak.description='height from BB Peak in meters a.g.l.'
            nc_bb_peak.units='m'


        estat,NewMatrix,z_da,Lwc,Rr,SnowRate,w,sig,sk,Noi,DSD,NdE,Ze,Mov,velTur,snr,kur,PiA,NW,DM,z_P,lwc_P,rr_P,Z_h,Z_all,RR_all,LWC_all,dm_all,nw_all,N_all=Process(proeta,Hcolum,Time[i],D,Cte,NewNoise,Deltav,Code_spectrum,Snr_Refl_2)
        if Timecount==0 or Timecount==1:
            if Timecount==0:
                bb_bot,bb_top,bb_peak=BB2(w,Ze,Hcolum,sk,kur,np.ones(2)*np.nan,np.ones(2)*np.nan,np.ones(2)*np.nan)#the bb is necessary calcultae after determinate w
                if Cons[Count_Cons]==1:
                    
                    if ~np.isnan(BB_last_bot) and ~np.isnan(bb_bot):
                        bb_bot=BB_last_bot
                        bb_peak=BB_last_peak
                        bb_top=BB_last_top
                        BB_last_bot=np.nan
                        BB_last_peak=np.nan
                        BB_last_top=np.nan
            if Timecount==1:
                bb_bot,bb_top,bb_peak=BB2(w,Ze,Hcolum,sk,kur,np.ones(2)*bb_bot_full,np.ones(2)*bb_top_full,np.ones(2)*bb_peak_full)#the bb is necessary calcultae after determinate w
        else:
            bb_bot,bb_top,bb_peak=BB2(w,Ze,Hcolum,sk,kur,bb_bot_full,bb_top_full,bb_peak_full)#the bb is necessary calcultae after determinate w

        pIA=10*np.log10(PiA)
        ZeCorrec=[]
        ZaCorrec=[]
        ZaCorrec_all=[]
            
        
        bb_bot_full.append(bb_bot);bb_top_full.append(bb_top);bb_peak_full.append(bb_peak)
        for j in range(len(Ze)):
            ZaCorrec_all.append(Z_all[j]-pIA[j])
            if estat[j]==10 or estat[j]==5:#PIA IS APPLIED ony for drizzle and rain
                ZeCorrec.append(Ze[j]-pIA[j])
                ZaCorrec.append(z_da[j]-pIA[j])
            else:
                ZeCorrec.append(Ze[j])
                ZaCorrec.append(np.nan)
        
        if option==1:
            nc_etaV[Timecount,:,:]=np.array(np.ma.masked_invalid(NewMatrix),dtype='f')
            
                
        if option==2:
            nc_3D[Timecount,:,:]=np.array(np.ma.masked_invalid(np.log10(NdE)),dtype='f')
            
                
        if option==3:
            nc_etaV[Timecount,:,:]=np.array(np.ma.masked_invalid(NewMatrix),dtype='f')
            nc_3D[Timecount,:,:]=np.array(np.ma.masked_invalid(np.log10(NdE)),dtype='f')
            
        if Timecount==0:
            estat_full=estat;sk_full=sk;kur_full=kur;PIA_full=pIA;
            w_full=w;LWC_full=Lwc;RR_full=Rr;Z_da_full=z_da;Z_ea_full=Ze;Z_e_full=ZeCorrec;Z_a_full=ZaCorrec
            sig_full=sig;VerMov_full=Mov;SnowR_full=SnowRate;Noi_full=Noi
            SNR_full=snr;N_da_full=DSD;VelTur_full=velTur;nw_full=NW;dm_full=DM

            z_all_full=ZaCorrec_all;lwc_all=LWC_all;rr_all=RR_all;DM_all=dm_all;NW_all=nw_all;n_all=N_all

            

                
            Z_P=z_P;LWC_P=lwc_P;RR_P=rr_P
        else:

            estat_full=np.vstack((estat_full,estat));sk_full=np.vstack((sk_full,sk));kur_full=np.vstack((kur_full,kur));PIA_full=np.vstack((PIA_full,pIA));Z_a_full=np.vstack((Z_a_full,ZaCorrec))
            w_full=np.vstack((w_full,w));LWC_full=np.vstack((LWC_full,Lwc));RR_full=np.vstack((RR_full,Rr));Z_da_full=np.vstack((Z_da_full,z_da));Z_ea_full=np.vstack((Z_ea_full,Ze));Z_e_full=np.vstack((Z_e_full,ZeCorrec))
            sig_full=np.vstack((sig_full,sig));VerMov_full=np.vstack((VerMov_full,Mov));SnowR_full=np.vstack((SnowR_full,SnowRate));Noi_full=np.vstack((Noi_full,Noi))
            SNR_full=np.vstack((SNR_full,snr));N_da_full=np.vstack((N_da_full,DSD));VelTur_full=np.vstack((VelTur_full,velTur));nw_full=np.vstack((nw_full,NW));dm_full=np.vstack((dm_full,DM))

                
            z_all_full=np.vstack((z_all_full,ZaCorrec_all));lwc_all=np.vstack((lwc_all,LWC_all));rr_all=np.vstack((rr_all,RR_all))
            NW_all=np.vstack((NW_all,nw_all));DM_all=np.vstack((DM_all,dm_all));n_all=np.vstack((n_all,N_all))
                
            Z_P=np.vstack((Z_P,z_P));LWC_P=np.vstack((LWC_P,lwc_P));RR_P=np.vstack((RR_P,rr_P))

            

            

        

        if ~np.isnan(DM).all():
            Nw_2.append(NW)
            Dm_2.append(DM)



        Timecount+=1
        tst2=(datetime.datetime.now()).timestamp()
        DeltaTimeWork=round((3600./float(TimeInt))*float(tst2-tst1),0)
        
        if condTime==0:
            print('The estimated time to finish this file is ',DeltaTimeWork,' seconds')
            condTime+=1
        if countwork==3:
            sys.stdout.write('   \b—\r')
            
        if countwork==0:
            sys.stdout.write('   \b|\r')
        if countwork==1:
            sys.stdout.write('   \b/\r')
        if countwork==2:
            sys.stdout.write('   \b—\r')
        if countwork==5:
            sys.stdout.write('   \b\\\r')
            countwork=0
        
        countwork+=1
        

        
########    generate a smooth signal from BB    
    bb_bot_full3=Inter1D(bb_bot_full)
    bb_top_full3=Inter1D(bb_top_full)
    bb_peak_full3=Inter1D(bb_peak_full)

    bb_bot_full2=anchor(bb_bot_full3,.95)
    bb_top_full2=anchor(bb_top_full3,.95)
    bb_peak_full2=anchor(bb_peak_full3,.95)
######        CHECK THE LEVEL AND DIMENSION
    for j in range(len(bb_bot_full2)):
        if bb_peak_full2[j]>bb_top_full2[j]:
            bb_peak_full2[j]=bb_top_full2[j]-DeltaH
        if bb_peak_full2[j]<bb_bot_full2[j]:
            bb_peak_full2[j]=bb_bot_full2[j]+DeltaH

       
        
        if np.isnan(bb_peak_full2[j]) and ~np.isnan(bb_bot_full2[j]) and np.isnan(bb_top_full2[j]):
            bb_bot_full2[j]=np.nan
        if np.isnan(bb_peak_full2[j]) and np.isnan(bb_bot_full2[j]) and ~np.isnan(bb_top_full2[j]):
            bb_top_full2[j]=np.nan
        if ~np.isnan(bb_peak_full2[j]) and np.isnan(bb_bot_full2[j]) and ~np.isnan(bb_top_full2[j]):
            bb_top_full2[j]=np.nan
            bb_peak_full2[j]=np.nan
        if ~np.isnan(bb_peak_full2[j]) and ~np.isnan(bb_bot_full2[j]) and np.isnan(bb_top_full2[j]):
            bb_bot_full2[j]=np.nan
            bb_peak_full2[j]=np.nan
        if ~np.isnan(bb_peak_full2[j]) and np.isnan(bb_bot_full2[j]) and np.isnan(bb_top_full2[j]):
            bb_bot_full2[j]=np.nan
            bb_top_full2[j]=np.nan
        if np.isnan(bb_peak_full2[j]) and ~np.isnan(bb_bot_full2[j]) and ~np.isnan(bb_top_full2[j]):
            bb_peak_full2[j]=bb_bot_full2[j]+((bb_top_full2[j]-bb_bot_full2[j])/2)

########    Correct the differnt values in function of the BB
    estat_full,Z_da_full,LWC_full,RR_full,SnowR_full=CorrectWithBBMatrix(estat_full,Z_da_full,LWC_full,RR_full,SnowR_full,Hcolum,bb_bot_full2,bb_top_full2,Z_ea_full,Z_P,LWC_P,RR_P,sk_full)
    

##    if option==1:
##        nc_etaV[:,:,:]=np.array(np.ma.masked_invalid(NewMatrix_full),dtype='f')
##    if option==2:
##        nc_3D[:,:,:]=np.array(np.ma.masked_invalid(np.log10(NdE_full)),dtype='f')
##    if option==3:
##        nc_etaV[:,:,:]=np.array(np.ma.masked_invalid(NewMatrix_full),dtype='f')
##        nc_3D[:,:,:]=np.array(np.ma.masked_invalid(np.log10(NdE_full)),dtype='f')
        
    nc_state[:,:]=np.array(np.ma.masked_invalid(estat_full),dtype='f')
    nc_w[:,:]=np.array(np.ma.masked_invalid(w_full),dtype='f')
    nc_sig[:,:]=np.array(np.ma.masked_invalid(sig_full),dtype='f')
    nc_sk[:,:]=np.array(np.ma.masked_invalid(sk_full),dtype='f')
    nc_kur[:,:]=np.array(np.ma.masked_invalid(kur_full),dtype='f')

    nc_PIA[:,:]=np.array(np.ma.masked_invalid(PIA_full),dtype='f')
                
    nc_LWC[:,:]=np.array(np.ma.masked_invalid(LWC_full),dtype='f')
    nc_RR[:,:]=np.array(np.ma.masked_invalid(RR_full),dtype='f')
                
    nc_Z_da[:,:]=np.array(np.ma.masked_invalid(Z_da_full),dtype='f')
    nc_Z_DA[:,:]=np.array(np.ma.masked_invalid(Z_a_full),dtype='f')
    nc_Z_ea[:,:]=np.array(np.ma.masked_invalid(Z_ea_full),dtype='f')
    nc_Z_e[:,:]=np.array(np.ma.masked_invalid(Z_e_full),dtype='f')
    nc_Z_all[:,:]=np.array(np.ma.masked_invalid(z_all_full),dtype='f')
    nc_LWC_all[:,:]=np.array(np.ma.masked_invalid(lwc_all),dtype='f')
    nc_RR_all[:,:]=np.array(np.ma.masked_invalid(rr_all),dtype='f')
    #nc_VerMov[:,:]=np.array(np.ma.masked_invalid(Mov),dtype='f')
                
    nc_SnowR[:,:]=np.array(np.ma.masked_invalid(SnowR_full),dtype='f')
    nc_Noi[:,:]=np.array(np.ma.masked_invalid(Noi_full),dtype='f')
    nc_SNR[:,:]=np.array(np.ma.masked_invalid(SNR_full),dtype='f')
    nc_N_da[:,:]=np.array(np.ma.masked_invalid(N_da_full),dtype='f')
    nc_N_all[:,:]=np.array(np.ma.masked_invalid(n_all),dtype='f')


    #nc_VelTur[:,:]=np.array(np.ma.masked_invalid(VelTur_full),dtype='f')
    nc_nw[:,:]=np.array(np.ma.masked_invalid(nw_full),dtype='f')
    nc_dm[:,:]=np.array(np.ma.masked_invalid(dm_full),dtype='f')

    nc_NW[:,:]=np.array(np.ma.masked_invalid(NW_all),dtype='f')
    nc_DM[:,:]=np.array(np.ma.masked_invalid(DM_all),dtype='f')
        

    nc_bb_bot[:,:]=np.array(np.ma.masked_invalid(bb_bot_full2),dtype='f')
    nc_bb_top[:,:]=np.array(np.ma.masked_invalid(bb_top_full2),dtype='f')
    nc_bb_peak[:,:]=np.array(np.ma.masked_invalid(bb_peak_full2),dtype='f')
    
    ##########        parameters of control continuity BB
    if ~np.isnan(bb_bot_full2).all():
            
        for m in range(len(bb_bot_full2)):
            if ~np.isnan(bb_bot_full2[-(m+1)]):
                indLastNoNan=m
                break
            if m==len(bb_bot_full2)-2:
                indLastNoNan=10
                break

        if indLastNoNan<=5:#only checked the last 5 elements
            BB_last_bot=bb_bot_full2[-(indLastNoNan+1)]
            BB_last_peak=bb_peak_full2[-(indLastNoNan+1)]
            BB_last_top=bb_top_full2[-(indLastNoNan+1)]

    if len(Dm_2)>2:
        dm_ax,nw_ax,PrepTypeC=PrepType(Dm_2,Nw_2)
        dataset.createDimension('Dm_ax',len(dm_ax))
        dataset.createDimension('Nw_ax',len(nw_ax))

        nc_ranges_Dm=dataset.createVariable('Dm_ax','f',('Dm_ax',))
        nc_ranges_Nw=dataset.createVariable('Nw_ax','f',('Nw_ax',))
        
        nc_ranges_Dm.description = 'mean diameter axes to Rainfall type'
        nc_ranges_Dm.units = '(mm)'
        
        nc_ranges_Nw.description = 'Intecept parameter axes to Rainfall type'
        nc_ranges_Nw.units = 'log(m-3 mm-1)'

        nc_ranges_Dm[:]=np.array(dm_ax,dtype='f4')
        nc_ranges_Nw[:]=np.array(nw_ax,dtype='f4')



        nc_TypePrecipitation=dataset.createVariable('TyPrecipi','f',('Dm_ax','Nw_ax',))
        nc_TypePrecipitation.description='Rainfall type where the value 5 is convective, 0 is transition and -5 is stratiform'
        nc_TypePrecipitation.units='none'

        nc_TypePrecipitation[:,:]=np.array(np.ma.masked_invalid(PrepTypeC),dtype='f')

        DM_ax,NW_ax,PrepTypeC_all=PrepType(DM_all,NW_all)
        dataset.createDimension('DM_ax',len(DM_ax)) 
        dataset.createDimension('NW_ax',len(NW_ax))

        nc_ranges_DM=dataset.createVariable('DM_ax','f',('DM_ax',))
        nc_ranges_NW=dataset.createVariable('NW_ax','f',('NW_ax',))
            
        nc_ranges_DM.description = 'mean diameter axes to Rainfall type suposing all hydrometeors are kiquid phase'
        nc_ranges_DM.units = '(mm)'
            
        nc_ranges_NW.description = 'Intecept parameter axes to Rainfall type  suposing all hydrometeors are kiquid phase'
        nc_ranges_NW.units = 'log(m-3 mm-1)'

        nc_ranges_DM[:]=np.array(DM_ax,dtype='f4')
        nc_ranges_NW[:]=np.array(NW_ax,dtype='f4')
        nc_TypePrecipitation_all=dataset.createVariable('TyPrecipi_all','f',('DM_ax','NW_ax',))
            
        nc_TypePrecipitation_all.description='Rainfall type where the value 5 is convective, 0 is transition and -5 is stratiform,  suposing all hydrometeors are kiquid phase'
        nc_TypePrecipitation_all.units='none'

        nc_TypePrecipitation_all[:,:]=np.array(np.ma.masked_invalid(PrepTypeC_all))
            
    
    print(datetime.datetime.now())
    print('\n')


    dataset.close()
    Count_Cons=Count_Cons+1
            


