#! /usr/bin/python3

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
from scipy import interpolate
from scipy.interpolate import RegularGridInterpolator
import sys


fname=sys.argv[1] #result file for plotting
J=int(sys.argv[2]) #order of statistics

title=r'$\bar{\xi}_'+str(J)+'$, example run - extended plot'
contours_vals=[3e-3,3e-2,1e-1,3e-1,1,3]  #contours
ticks=['-50','-20','-10','-5','2','5','10','20','50']


scales=np.unique(np.genfromtxt(fname, usecols=(0),unpack=True)) #all scales present here
[naxa,naxb]=[len(scales),len(scales)]  #no. of r_par and r_perp grid points (still script for square only)
grid_1D=len(scales) #grid points in 1D

dax=np.log10(scales[-1]/scales[0])/(grid_1D-1) #step of power (starting from 0, thus -1), helpful for grid <-> vaue conversion



#2D interpolation
def Interp_2D(data, xvals, yvals, x0,y0):
    interp=RegularGridInterpolator((xvals,yvals),data,method='cubic')
    return interp([x0,y0])




#position on the heatmap for given [rpar,rperp]=[axa,axb] (accounting for mirrored versions)
def Position(axa,axb):

    #right-top corner of plot
    posx=np.log10(axa/scales[0])/dax +naxa-1
    posy=naxa-1- np.log10(axb/scales[0])/dax
    
    return [posx,posy,
    2*naxa-2-posx,posy,
    posx,naxa-1+np.log10(axb/scales[0])/dax,
    2*naxa-2-posx, naxa-1+np.log10(axb/scales[0])/dax]



#inverse: obtaining [axa,axb] from position on the grid (may be not integers)
#only from right-upper mirrored part
def Vals(posx,posy):
    axa=10**(np.log10(scales[0]) +(posx-naxa+1)*dax)
    axb=10**(np.log10(scales[0])+ np.fabs(naxa-1-posy)*dax )
    
    if np.fabs(scales[0]-axa)<1e-5: axa=scales[0]
    if np.fabs(scales[0]-axb)<1e-5: axb=scales[0]
    if np.fabs(scales[-1]-axa)<1e-5: axa=scales[-1]
    if np.fabs(scales[-1]-axb)<1e-5: axb=scales[-1]
    
    return [axa,axb]




#lines of constant volume for given effective radius (from definition 4/3 pi R_eff^3 = 4/3 pi * r_par * r_perp^2)
#nn - number of points that form the line
def Constvol(R_eff,nn):
    xx,yy=np.zeros((4,nn)),np.zeros((4,nn)) #line for each mirrored part
    #a,b axes used (need to not exceed scales range)
    b=np.logspace(np.log10(max((R_eff**3/max(scales))**0.5,scales[0])),np.log10(min((R_eff**3/min(scales))**0.5,scales[-1])),nn)
    a=R_eff**3 /b**2 #condition for constant volume
    
    xx[0,:],yy[0,:],xx[1,:],yy[1,:],xx[2,:],yy[2,:],xx[3,:],yy[3,:]=Position(a,b) #converting [a,b]=[rpar,rperp] into positions on the heatmap
    
    return xx,yy
    


#plotting lines of constant volumes
def Plot_constvol(ax):
    nnReff=5 #desired number of constant volume lines
    R_eff=np.logspace(np.log10(2.*min(scales)),np.log10(0.9*max(scales)),nnReff) #effective radii for that lines
    for i in range(0,nnReff): #each constant volume line:
        xx,yy=Constvol(R_eff[i],500) #constant volume line for this effective radius, made of 500 points
        for j in range(0,4): #each mirrored part
            ax.plot(xx[j,:],yy[j,:],color='darkblue',linewidth=1.5)

    #lines of rpar=rperp
    arg=np.arange(0,2*len(scales)-1)
    ax.plot(arg,arg,linewidth=1,linestyle='--',color='black',alpha=0.5)
    ax.plot(arg,2*len(scales)-2-arg,linewidth=1,linestyle='--',color='black',alpha=0.5)
    
    return




#lines of maximum value for given volume
def Plot_maxconstvol(data,ax):
    nnReff=50 #number of effective radii considered
    R_eff=np.logspace(np.log10(1.3*min(scales)),np.log10(0.9*max(scales)),nnReff) #effective radii
    
    maxx,maxy=np.empty(0),np.empty(0) #storing max localisations
    for i in range(0,nnReff): #for each effective radii (constant volume)
        print('Searching volumes: ',i,'/',nnReff)
        xx,yy=Constvol(R_eff[i],100) #positions on the heatmap of points with effective radii of R_eff[i]
        maxval=-1e6 #for searching of maximal value
        maxvalloc=[-1,-1] #for location of maximal value
        
        for j in range(0,len(xx[0])): #searching for max: (points along current i-th volume)
            if xx[0,j]<naxa or xx[0,j]>2*naxa-2 or yy[0,j]<0 or yy[0,j]>2*naxb-1: continue
            axa,axb=Vals(xx[0,j],yy[0,j]) #[rpar,rperp] values for this position

            val=Interp_2D(data,scales,scales,axa,axb) #Xi_J value here
            if val>maxval: #condition if its larger than maximum yet found
                maxval=val
                maxvalloc=[xx[0,j],yy[0,j]] #location of maximum on the grid
        
        axa_max,axb_max=Vals(maxvalloc[0],maxvalloc[1]) #[rpar,rperp] for maximized signal at i-th constant volume
        if axa_max==min(scales) or axa_max==max(scales) or axb_max==min(scales) or axb_max==max(scales):
            continue #on the border - max should be outside of range
        maxx=np.append(maxx,maxvalloc[0]) #adding this maximized point to the array
        maxy=np.append(maxy,maxvalloc[1]) #adding this maximized point to the array
    
    #plotting all mirrored maximized sets:
    ax.plot(maxx,maxy,linewidth=1.5,color='orangered')
    ax.plot(2.*naxa-2-maxx,maxy,linewidth=1.5,color='orangered')
    ax.plot(maxx,2.*naxb-2-maxy,linewidth=1.5,color='orangered')
    ax.plot(2.*naxa-2-maxx,2.*naxb-2-maxy,linewidth=1.5,color='orangered')
    return







#collecting data:-------------------------------------------------------------------------------------------------------------------------------
val_onequarter=np.zeros((naxa,naxb)) #data - one quarter of the plot (independent)
val=np.zeros((2*naxa-1,2*naxb-1)) #data - full, containing mirrored val_onequarter data

val_onequarter_1D=np.genfromtxt(fname, usecols=(2*J-2),unpack=True) #1D array of results
val_onequarter=val_onequarter_1D.reshape([naxb,naxa]) #reshaping into 2D grid

#mirrored versions:
val[:naxa-1,naxb-1:]=val_onequarter[::-1,:][:naxa-1,:] #left top
val[naxa-1:,naxb-1:]=val_onequarter[:,:] #right top
val[:naxa-1,:naxb-1]=val_onequarter[::-1,::-1][:naxa-1,:naxb-1] #left bottom
val[naxa-1:,:naxb-1]=val_onequarter[:,::-1][:,:naxb-1] #right bottom

#plotting ------------------------------------------------------------------------------------------------------------------------------

#axes ticks:
tickvals=np.zeros(len(ticks))
for i in range(0,len(ticks)): tickvals[i]=float(ticks[i])

tlocs=np.zeros(len(tickvals))
for i in range(0,len(tickvals)): #location of ticks
    tlocs[i]=len(scales) -1 + np.sign(tickvals[i])* np.log10(np.abs(tickvals[i])/scales[0])/dax


fig, ax = plt.subplots(figsize=(6,5))
ax.set_box_aspect(1) #square size
ax.grid(linestyle='--',color='grey',alpha=.7) #grid

ax.set_xticks(tlocs,ticks,fontsize=13)
ax.set_yticks(tlocs,ticks[::-1],fontsize=13)

ax.set_xlabel('$r_{\parallel}$' +' [Mpc]',fontsize=14)
ax.set_ylabel('$r_{\perp}$' +' [Mpc]',fontsize=14)
ax.set_title(title,fontsize=15,alpha=.7)


#making heatmap (logarithmic colorscale)
cf=plt.imshow(np.transpose(val[:,::-1]), cmap='cubehelix',interpolation = 'bicubic',aspect=1,norm=mpl.colors.LogNorm(vmin=np.amin(val),vmax=np.amax(val)))


#adding contours
contours = plt.contour(np.arange(0,2*grid_1D-1,1),np.arange(0,2*grid_1D-1,1), np.transpose(val[:,::-1]),contours_vals, linestyles='--',colors='cyan')
plt.clabel(contours, inline=True, fontsize=11)
        

#colorbar
cb=fig.colorbar(cf,ax=ax)
cb.set_label(r'$\bar{\xi}_'+str(J)+'$', fontsize=13) #label
cb.ax.tick_params(labelsize=13) #ticks label size


Plot_constvol(ax) #pllotting constant volumes lines
Plot_maxconstvol(val_onequarter,ax) #plotting lines of maximized Xi_J at constant volume

plt.savefig('Heatmap_Xi_'+str(J)+'_extended.png',bbox_inches='tight')