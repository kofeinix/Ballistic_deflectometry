# -*- coding: utf-8 -*-
"""
Created on Sat Sep 11 23:03:36 2021

@author: Юрий
"""

from PIL import Image, ImageDraw
import matplotlib.pyplot as plt
import numpy as np
import math
from math import sqrt, pi, sin, cos
import random
import time
from numba import jit, njit, prange
import os
import matplotlib.cm as cm 
import scipy
import My_Colormaps
os.environ["NUMBA_WARNINGS"] = "0"

def field_energy(B,E):
    wB=0
    for i in range(0,len(B)-1,1):
        wB=(((B[i,3])**2+(B[i,4])**2+(B[i,5])**2)*(step)**3)/(2.51*10**(-6)) + wB
    wE=0    
    for i in range(0,len(E)-1,1):
        wE=(((E[i,3])**2+(E[i,4])**2+(E[i,5])**2)*(step)**3)*(4.425*10**(-12)) + wE

save_path='K:\\Work\\Python\\GitHub\\Ballistic\\output_files\\'


print("loading magnetic field...")
B1=np.genfromtxt('K:\\Work\\Python\\GitHub\\Ballistic\\input_files\\60-70-0-70.txt')
B2=np.loadtxt('K:\\Work\\Python\\GitHub\\Ballistic\\input_files\\xyz05.txt')                           #contains x y z coordinates 

print("filtering magnetic field...")
#Radia sometimes gives some inf and NANs - should get rid of it by filtering.
for i in range(0,len(B1)-1,1):
    if math.isfinite(B1[i,0])==False:  
        B1[i,0]=(B1[i-1,0]+(B1[i-1,0]-B1[i-2,0]))    
    if math.isfinite(B1[i,1])==False:
        B1[i,1]=(B1[i-1,1]+(B1[i-1,1]-B1[i-2,1]))  
    if math.isfinite(B1[i,2])==False :
        B1[i,2]=(B1[i-1,2]+(B1[i-1,2]-B1[i-2,2]))  

print("making one mag field array...")

#Making one big file with columns: x y z Bx By Bz
B=np.zeros((len(B1),6))
B[:,0]=B2[:,0]
B[:,1]=B2[:,1]
B[:,2]=B2[:,2]
B[:,3]=B1[:,0]
B[:,4]=B1[:,1]
B[:,5]=B1[:,2]

#Convert it from default Radia mm into meters
for i in range(0,len(B),1):
    for j in range(0,3,1):
        B[i,j]=round((B[i,j]/1000),6)

#Reading electric field files.
#Used version of Comsol provides data in format: x y z Ei
print("loading electric field...")
E1=np.loadtxt('K:\\Work\\Python\\GitHub\\Ballistic\\input_files\\ExCCW.txt')
E2=np.loadtxt('K:\\Work\\Python\\GitHub\\Ballistic\\input_files\\EyCCW.txt')
E3=np.loadtxt('K:\\Work\\Python\\GitHub\\Ballistic\\input_files\\EzCCW.txt')


print("making one el field array...")
#Making one big file with columns: x y z Ex Ey Ez
E=np.zeros((len(E1),6))
E[:,0]=np.round(E1[:,0],6)
E[:,1]=np.round(E1[:,1],6)
E[:,2]=np.round(E1[:,2],6)
E[:,3]=E1[:,3]
E[:,4]=E2[:,3]
E[:,5]=E3[:,3]

#Checking for artifacts
for i in range(0,len(E),1):
    if math.isfinite(E[i,3])==False or math.isfinite(E[i,4])==False or math.isfinite(E[i,5])==False:    
        E[i,3]=E[i,4]=E[i,5]=1.0
        
print('deleting temp')                 
#Delete temporary files
del(E1,E2,E3,B1,B2)

#Define parameters that were used to create field files. If try to use other coordinates - it won't work. Why?
#At the initial version the current coordinates of the proton were used for SEARCHing for the corresponding field value in arrays.
#However that was too slow using any method availiable
#Thus now it's not a search, but a calculation of the correct row in the array using known coordinate limits, known sorting of x y z, and the step value.

xmax=600e-6
xmin=-600e-6
step=10e-6
ymax=500e-6
ymin=-500e-6
zmin=-500e-6
zmax=1000e-6

print('defining coordinate binding')
# Search for a row, where some given coordinate will change to a larger coordinate
def find_max(E,coord,pos,coordmin,coordmax,minim,maxim):
    j=maxim
    i=minim
    m = int((j+i)/2)
    mmax=0
    while i < j:
        if coord>coordmax:
            mmax=-1
            break
        if coord<coordmin:
            mmax=-1
            break;
        if coord>E[m,pos]:
            i = m+1
            m = int(round((i+j)/2))
            #print('1',i,j,m)
    
        if coord<E[m,pos]:
            j=m-1
            m = int(round((i+j)/2))
            #print('2',i,j,m)
        if coord==E[m,pos]:
            if maxim-m<10:
                mmax=maxim
                break
            if E[m-1,pos]==E[m,pos]:
                i=m
                m = int(round((i+j)/2))+1
                #print('4',i,j,m)
            if E[m-1,pos]<E[m,pos]:
                mmax=m-1
                #print('3')
                break

    return(mmax) 
    
# Search for a row, where some given coordinate will change to a smaller coordinate
def find_min(E,coord,pos,coordmin,coordmax,minim,maxim):
    j=maxim
    i=minim
    m = int((j+i)/2)
    mmin=0
    while i < j:
        if coord>coordmax:
            mmin=-1
            break
        if coord<coordmin:
            mmin=-1
            break;
        if coord>E[m,pos]:
            i = m+1
            m = int((i+j)/2)
            #print('1',i,j,m)
        if coord<E[m,pos]:
            j=m-1
            m = int((i+j)/2)
            #print('2',i,j,m)
        if coord==E[m,pos]:
            if E[m+1,0]==E[m,0]:
                j=m
                m = int((i+j)/2)-1
               # print('4',i,j,m)
            if E[m+1,pos]>E[m,pos]:
                mmin=m+1 
                #print('3',i,j,m)
                break
    return(mmin)

imaxx=round((xmax-xmin)/step)+1             #How much steps in x direction
imaxy=round((ymax-ymin)/step)+1             #How much steps in y direction
imaxz=round((zmax-zmin)/step)+1             #How much steps in z direction
M=np.zeros((imaxx,2))                       #Here "row, when x changes" is stored
N=np.zeros((imaxy,2))                       #Here "row, when y changes for one fixed x" is stored
for i in range(0,imaxx,1):
    x=xmin+step*i 
    M[i,1]=find_max(B,round(x,6),0,xmin,xmax,0,(len(B)-1))
    M[i,0]=find_min(B,round(x,6),0,xmin,xmax,0,(len(B)-1))
for i in range(0,imaxy,1):
    y=ymin+i*step
    N[i,0]=find_min(B,round(y,6),1,ymin,ymax,0,M[0,1])
    N[i,1]=find_max(B,round(y,6),1,ymin,ymax,0,M[0,1])

#Function to change the electric field (linear dependence so it's ok)
def coeff(E,ke):
    E1=np.zeros_like(E)
    E1[:,0]=E[:,0]
    E1[:,1]=E[:,1]
    E1[:,2]=E[:,2]
    E1[:,3]=E[:,3]*ke
    E1[:,4]=E[:,4]*ke
    E1[:,5]=E[:,5]*ke
    return(E1)

#Function to plot the field    
#arguments - 1. field array. 2. which plane. 3. coordinate along axis perpendicular to chosen cross-section. 4. which field component to plot. 5-6. upper and lower limitsof plot (colors), 7-8. coord and numbers for naming

def PlotField(Field, AxisA, coord, Component,lowlimit,uplimit,position,kk):
    if AxisA=='XY':  
        emptest=np.zeros((int((xmax-xmin)/step),int((ymax-ymin)/step)))
        x=xmin
        z=coord
        y=ymin
        for i in range(0,(int((xmax-xmin)/step)),1):
            y=ymin
            x=x+step
        
            for j in range(0,(int((ymax-ymin)/step)),1):
                y=y+step
                ix=int(round((x-xmin)/step))
                iy=int(round((y-ymin)/step))
                iz=int(round((z-zmin)/step))
                mminx=M[int(ix),0]
                mminy=N[int(iy),0]
                mminz=mminy+mminx
                rez=int(mminz+iz)
                Ex=Field[rez,3]
                Ey=Field[rez,4]
                Ez=Field[rez,5]
                Emod=Ex**2+Ey**2+Ez**2
                if Component=='Abs':
                    emptest[i,j]=sqrt(Emod)
                elif Component=='x':
                    emptest[i,j]=Ex
                elif Component=='y':
                    emptest[i,j]=Ey
                elif Component=='z':
                    emptest[i,j]=Ez
                else:
                    print('BAAAKA')
                    break
        plt.imshow(np.flipud(emptest),cmap=My_Colormaps.my_dir,vmin=lowlimit, vmax=uplimit)
        cbar=plt.colorbar()
        plt.title('FieldXY, '+(Component)+'-Component,z='+str(round(position,6)))
        plt.savefig(save_path+'FieldXY'+str(kk)+'.png',bbox_inches='tight',dpi=200)
        plt.close()

    if AxisA=='XZ':  
        emptest=np.zeros((int((xmax-xmin)/step),int((zmax-zmin)/step)))
        x=xmin
        y=coord
        z=zmin
        for i in range(0,(int((xmax-xmin)/step)),1):
            z=zmin
            x=x+step
        
            for j in range(0,(int((zmax-zmin)/step)),1):
                z=z+step
                ix=int(round((x-xmin)/step))
                iy=int(round((y-ymin)/step))
                iz=int(round((z-zmin)/step))
                mminx=M[int(ix),0]
                mminy=N[int(iy),0]
                mminz=mminy+mminx
                rez=int(mminz+iz)
                Ex=Field[rez,3]
                Ey=Field[rez,4]
                Ez=Field[rez,5]
                Emod=Ex**2+Ey**2+Ez**2
                if Component=='Abs':
                    emptest[i,j]=sqrt(Emod)
                elif Component=='x':
                    emptest[i,j]=Ex
                elif Component=='y':
                    emptest[i,j]=Ey
                elif Component=='z':
                    emptest[i,j]=Ez
                else:
                    print('BAAAKA')
                    break     
        plt.imshow(np.flipud(emptest),cmap=My_Colormaps.my_dir,vmin=lowlimit, vmax=uplimit)
        cbar=plt.colorbar()
        plt.title('FieldXz, '+(Component)+'-Component,y='+str(round(position,6)))
        plt.savefig(save_path+'FieldXZ'+str(kk)+'.png',bbox_inches='tight',dpi=200)
        plt.close()



    
#Some constants and particle parameters. dt - time step. T - particle energy in MeV. Vmaxall - convertion of this energy into speed (relativistic)
#nup - number of particles that didn't pass through the field area. U is initial potential of target, used in COMSOL!
print('initializing parameters')

q=1.6e-19 #proton charge
m=1.67e-27 #proton mass
c=3e8   #speed of light
dt=1e-12 #time step
E0=m*c**2  #rest energy
T=3.6e6  #initial particle energy
Vmaxall=c*sqrt(1-(E0**2)/((T*1.6e-19+E0)**2))  #particle speed
nup=0
U=1000 #potential from COMSOL

#numba is used to accelerate calculations
@njit(parallel=True)
def calculation(B,E,Npart,Edelete1, interpolation):
    nup=0
    Picture=np.zeros((5000,5000)) 
    for i in prange (0,Npart,1):  # CAN BE SWAPPED TO SIMPLE RANGE
        if (i % 100000==0):
            print('Progress report..',((i/Npart)*100)) # messed up for prange, but useful for range
        teta=random.gauss(0,10)                             #particle launch angle
        phi=random.randrange(0,3600,1)                      #uniform distribution over phi
        Vy=-Vmaxall*cos(teta*pi/180)                        #speed projections
        Vx=Vmaxall*sin(teta*pi/180)*cos(phi*pi/1800)
        Vz=Vmaxall*sin(teta*pi/180)*sin(phi*pi/1800)
        x=0#random.uniform(-2e-6,2e-6)                      #can be not a point-like source, but more physical
        y=3000e-6                                           #Initial coordinate of a particle source (TNSA foil). They fly towards x=y=z=0, where the field is centered. 
        z=0#random.uniform(-2e-6,2e-6)
        
        ttomesh=abs((y-1310e-6)/Vy)                     #Time from TNSA foil to the mesh.
        xmesh=x+Vx*ttomesh                              #Coordinates of the mesh to stop a particle 
        zmesh=z+Vz*ttomesh
        for n in range(-10000,10000,84):
            if (n-18)*(2e-7)<=xmesh<=(n+18)*(2e-7) or (n-18)*(2e-7)<=zmesh<=(n+18)*(2e-7):
                x=1e6
                z=1e6
                break
        

        ttofield=abs((y-ymax)/Vy)                       #Simple fly of a particle to the start of the field area 
        y=ymax
        x=x+Vx*ttofield
        z=z+Vz*ttofield
        #Then use equation of motion when inside the field area
        ax=ay=az=0
        while y>ymin:
            if interpolation==True:
                ix1=int(math.floor((x-xmin)/step))
                iy1=int(math.floor((y-ymin)/step))
                iz1=int(math.floor((z-zmin)/step))
                
                ix2=ix1+1
                iy2=iy1+1
                iz2=iz1+1
                # this is Trilinear interpolation of the field at the point where particle is. 
                # The results does not differ much though, thus this can be modified and deleted. 
    
                if ix1>=0 and ix1<imaxx and iy1>=0 and iy1<imaxy and iz1>=0 and iz1<imaxz and ix2>=0 and ix2<imaxx and iy2>=0 and iy2<imaxy and iz2>=0 and iz2<imaxz :
                    
                    mminx1=M[int(ix1),0]
                    mminy1=N[int(iy1),0]
                    mminx2=M[int(ix2),0]
                    mminy2=N[int(iy2),0]
    
                    mminz11=mminy1+mminx1
                    rez=int(mminz11+iz1)
                    B111=B[rez,3:6]
                    E111=E[rez,3:6]
                    x1,y1,z1=B[rez,0:3]
    
                    mminz12=mminy1+mminx2
                    rez=int(mminz12+iz1)
                    B211=B[rez,3:6]
                    E211=E[rez,3:6]
                   
                    rez=int(mminz12+iz2) 
                    B212=B[rez,3:6]
                    E212=E[rez,3:6]
                    
                    mminz22=mminy2+mminx2
                    rez=int(mminz22+iz1) 
                    B221=B[rez,3:6]
                    E221=E[rez,3:6]
                    
                    rez=int(mminz22+iz2) 
                    B222=B[rez,3:6]
                    E222=E[rez,3:6]
                    x2,y2,z2=B[rez,0:3]
    
                    rez=int(mminz11+iz2) 
                    B112=B[rez,3:6]
                    E112=E[rez,3:6]
    
                    mminz21=mminy2+mminx1
                    rez=int(mminz21+iz1) 
                    B121=B[rez,3:6]
                    E121=E[rez,3:6]
                    
                    rez=int(mminz21+iz2) 
                    B122=B[rez,3:6]
                    E122=E[rez,3:6]
                    
                    if y1!=y2 and x1!=x2 and z1!=z2:                
                        B111a=B111*(x2-x)*(y2-y)*(z2-z) 
                        B211a=B211*(x-x1)*(y2-y)*(z2-z)  
                        B212a=B212*(x-x1)*(y2-y)*(z-z1)  
                        B221a=B221*(x-x1)*(y-y1)*(z2-z)  
                        B222a=B222*(x-x1)*(y-y1)*(z-z1)  
                        B112a=B112*(x2-x)*(y2-y)*(z-z1)  
                        B121a=B121*(x2-x)*(y-y1)*(z2-z) 
                        B122a=B122*(x2-x)*(y-y1)*(z-z1)
                        Bx,By,Bz=(B111a+B211a+B212a+B221a+B222a+B112a+B121a+B122a)/((x2-x1)*(y2-y1)*(z2-z1))
                        
                        E111a=E111*(x2-x)*(y2-y)*(z2-z) 
                        E211a=E211*(x-x1)*(y2-y)*(z2-z)  
                        E212a=E212*(x-x1)*(y2-y)*(z-z1)  
                        E221a=E221*(x-x1)*(y-y1)*(z2-z)  
                        E222a=E222*(x-x1)*(y-y1)*(z-z1)  
                        E112a=E112*(x2-x)*(y2-y)*(z-z1)  
                        E121a=E121*(x2-x)*(y-y1)*(z2-z) 
                        E122a=E122*(x2-x)*(y-y1)*(z-z1)
                        Ex,Ey,Ez=(E111a+E211a+E212a+E221a+E222a+E112a+E121a+E122a)/((x2-x1)*(y2-y1)*(z2-z1))
    
                    else:
                        Bx,By,Bz=B222
                        Ex,Ey,Ez=E222
                        #(E[rez,3], E[rez,4], E[rez,5])
                        #print('interpolation is bad')
                    #Edelete1 is a specialy made value, that allow to delete particles that try to pass through the target bulk.
                    if (E111[0]==Edelete1 or E111[1]==Edelete1 or E111[2]==Edelete1 or
                        E211[0]==Edelete1 or E211[1]==Edelete1 or E211[2]==Edelete1 or
                        E212[0]==Edelete1 or E212[1]==Edelete1 or E212[2]==Edelete1 or
                        E221[0]==Edelete1 or E221[1]==Edelete1 or E221[2]==Edelete1 or 
                        E222[0]==Edelete1 or E222[1]==Edelete1 or E222[2]==Edelete1 or
                        E112[0]==Edelete1 or E112[1]==Edelete1 or E112[2]==Edelete1 or
                        E121[0]==Edelete1 or E121[1]==Edelete1 or E121[2]==Edelete1 or
                        E122[0]==Edelete1 or E122[1]==Edelete1 or E122[2]==Edelete1):
                        x=1e6
                        z=1e6
                        nup=nup+1
                        break
    
    
                else:
                    Bx,By,Bz=(0,0,0)
                    Ex,Ey,Ez=(0,0,0)

            if interpolation==False:
                ix=int(round((x-xmin)/step))
                iy=int(round((y-ymin)/step))
                iz=int(round((z-zmin)/step))
                if ix>=0 and ix<imaxx and iy>=0 and iy<imaxy and iz>=0 and iz<imaxz:
                    mminx=M[int(ix),0]
                    mminy=N[int(iy),0]
                    mminz=mminy+mminx
                    rez=int(mminz+iz)
                    Bx,By,Bz=(B[rez,3], B[rez,4], B[rez,5])
                    Ex,Ey,Ez=(E[rez,3], E[rez,4], E[rez,5])
                else:
                    Bx,By,Bz=(0,0,0)
                    Ex,Ey,Ez=(0,0,0)
                if Ex==Edelete1 and Ey==Edelete1 and Ez==Edelete1:
                    x=1e6
                    z=1e6
                    nup=nup+1
                    #print(nup)
                    break
                
            Vx=Vx+ax*dt
            Vy=Vy+ay*dt
            Vz=Vz+az*dt        
            ax=(q*(Vy*Bz-Vz*By)/m) + q*Ex/m
            ay=(-q*(Vx*Bz-Vz*Bx)/m) + q*Ey/m
            az=(q*(Vx*By-Vy*Bx)/m) + q*Ez/m
            x=x+Vx*dt+(ax*dt**2)/2
            y=y+Vy*dt+(ay*dt**2)/2
            z=z+Vz*dt+(az*dt**2)/2
            #if speed is derected at the other side, the particle is reflected.
            if Vy>0:
                x=1e6
                z=1e6
                nup=nup+1
                break
        #simple fly to the detector. 111e-3 is a distance from the target at (0,0,0) to the detector (RCF)
        ttorcf=abs((111e-3-y)/Vy)
        xrcf=x+Vx*ttorcf
        zrcf=z+Vz*ttorcf
        if abs(xrcf)>2499e-5 or abs(zrcf)>2499e-5:
            xrcf=2499e-5
            zrcf=2499e-5
        Vfin=sqrt(Vx**2+Vy**2+Vz**2)
        Efin=E0*(1/sqrt(1-(Vfin**2)/(c**2))-1)
        Tfin=(Efin/(1.6e-19))*100/T -100
        #Check the difference in initial and final energies. RCF can detect only a small range thus only particles with small difference from initial are recorded.
        if (abs(Tfin)<5):
            uu=Picture[int(xrcf*1e5)+2500,int(zrcf*1e5)+2500]
            if (uu<8):
                Picture[int(xrcf*1e5)+2500,int(zrcf*1e5)+2500]=uu+1
    return(Picture)
   
        

ix=int(round((-xmin)/step))
iy=int(round((-ymin)/step))
iz=int(round((-zmin)/step))
mminx=M[int(ix),0]
mminy=N[int(iy),0]
mminz=mminy+mminx
rez=int(mminz+iz)

# very important - this is a field in point (0,0,0)! and this is z component!
centr=B[rez,5]
#value to delete particles that try to fly the target bulk
Edelete=1.0
Picture=np.zeros((5000,5000))
Bznach=[300, 600]
Eznach=[25000, 100000]

print('starting calculations...')
for tt in range(0,len(Bznach),1):
    B1=coeff(B,Bznach[tt]/centr)
    for ttt in range (0, len(Eznach),1):
        E1=coeff(E,Eznach[ttt]/1000)
        U1=U*(Eznach[ttt]/1000)
        Edelete1=Edelete*(Eznach[ttt]/1000)
        Picture=calculation((B1),E1,100000000,Edelete1, interpolation=True)
        print(U1,B1[rez,5] )
        plt.imsave(save_path+'B='+str(Bznach[tt])+'U='+str(U1)+'.png', np.array(Picture).reshape(5000,5000),cmap=cm.Greys)

print('plotting field')
kk=0
for i in np.arange(ymin,ymax,step):
    #arguments - 1. field array. 2. which plane. 3. coordinate along axis perpendicular to chosen cross-section. 4. which field component to plot. 5-6. upper and lower limitsof plot (colors), 7-8. coord and numbers for naming
    PlotField(B,'XZ', i,'x',-1000,1000,i,kk)
    kk=kk+1

