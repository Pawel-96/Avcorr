#! /usr/bin/python3

import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import sys

fname=sys.argv[1] #result file for plotting
J=int(sys.argv[2]) #order of statistics

color='blue' #line color
title=r'$\bar{\xi}_'+str(J)+'$, example run'
ALPHA=0.2 #alpha for error filling
LW=2 #line width



def Logerr(data,error): #to have nice symmetric errors on logscale
    if data<=0: return ['nan',0]
    lgerr=math.fabs(error/(1.0*data*math.log(10))) #logarithmic error
    err_ranges=[data-10**(math.log10(data)-lgerr),10**(math.log10(data)+lgerr)-data]
    return err_ranges

#plotting--------------------------------------------------------------------------------------------

fig,ax =plt.subplots() #figsize=(16,8))

ax.grid(color='grey',linestyle='--',alpha=.7)
ax.set_xlabel('$R[Mpc]$',fontsize=13)
ax.set_ylabel(r'$\bar{\xi}_'+str(J)+'(R)$',fontsize=13)
ax.set_yscale('log')
ax.set_xscale('log')
ax.set_title(title, fontsize=18,alpha=0.5)
ax.tick_params(axis='both', which='major', labelsize=13)


#reading data
scale,value,error=np.genfromtxt(fname,usecols=(0, 2*J-3, 2*J-2),unpack=True) #format is: scale, Xi2, u_Xi2, Xi3, uXi3,...

nnR=len(scale) #number of scales considered
errr=np.zeros((2,nnR)) #for nicely looking symmetric logscale error
for m in range(0,nnR): errr[:,m]=Logerr(value[m],error[m]) #setting logscale errors

ax.fill_between(scale,value-errr[0,:],value+errr[1,:],color=color,alpha=ALPHA)
ax.plot(scale,value,color=color,linewidth=LW,label='result')
ax.legend()


plt.savefig('Xi_'+str(J)+'.png',bbox_inches='tight')