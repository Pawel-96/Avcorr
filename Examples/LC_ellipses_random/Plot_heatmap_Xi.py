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

title=r'$\bar{\xi}_'+str(J)+'$, example run'
contours_vals=[3e-3,3e-2,1e-1,3e-1,1,3]  #contours
ticks=['-50','-20','-10','-5','2','5','10','20','50']


scales=np.unique(np.genfromtxt(fname, usecols=(0),unpack=True)) #all scales present here
[naxa,naxb]=[len(scales),len(scales)]  #no. of r_par and r_perp grid points (still script for square only)
grid_1D=len(scales) #grid points in 1D

dax=np.log10(scales[-1]/scales[0])/(grid_1D-1) #step of power (starting from 0, thus -1), helpful for grid <-> vaue conversion


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

plt.savefig('Heatmap_Xi_'+str(J)+'.png',bbox_inches='tight')