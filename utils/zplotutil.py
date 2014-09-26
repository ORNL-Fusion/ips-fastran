#!/usr/bin/env Python

#######################################################################
# Last modified : Mar 2013
# plot utilities
# JM

from numpy import *
from scipy.interpolate import interp1d

import matplotlib
from pylab import *
from matplotlib.backends.backend_pdf import PdfPages
    
matplotlib.rcParams['axes.linewidth']=0.5
matplotlib.rcParams['pdf.fonttype']=42

import Namelist

#######################################################################
# simple wrapper

class zplotutil():

    def __init__(self,f_pdf,size_x=8.5,size_y=11.0):
        self.pdf = PdfPages(f_pdf)
        figure(figsize=(size_x,size_y))

    def savefig(self):
        self.pdf.savefig()
        clf()
    
    def close(self):
        self.pdf.close()

#######################################################################
# plot template

def plot_s(pos,x,y,xax,yax
    ,xerr=None,yerr=None
    ,color='k'
    ,markersize=10
    ,nxtick=2,nytick=2,xlab=None,ylab=None
    ,ylog=False,xlog=False
    ,iref=0
    ,iline=False,isym='o'):

    a = subplot(pos)

    if xlog:
        a.set_xscale('log')
    if ylog:
        a.set_yscale('log')

    if xerr!=None and yerr!=None:
        xx = x
        yy = y
        xe = xerr
        ye = yerr
        a.errorbar(xx,yy,xerr=xe,yerr=ye,
           fmt=color+isym, #'o',
           markerfacecolor=color,
           markeredgecolor=color,
           markersize=markersize)
    elif xerr==None and yerr!=None:
        xx = x
        yy = y
        ye = yerr
        a.errorbar(xx,yy,yerr=ye,
           fmt=color+isym, #'o',
           markerfacecolor=color,
           markeredgecolor=color,
           markersize=markersize)
    elif xerr!=None and yerr==None:
        xx = x
        yy = y
        xe = xerr
        a.errorbar(xx,yy,xerr=xe,
           fmt=color+isym, #'o',
           markerfacecolor=color,
           markeredgecolor=color,
           markersize=markersize)
    else: 
        xx = x
        yy = y
        a.plot(xx,yy,color+isym, #'o',
           markerfacecolor=color,
           markeredgecolor=color,
           markersize=markersize)

    if iline=='spl':
       n0 = 100
       x0 = x[0]+(x[-1]-x[0])*arange(100)/(n0-1.0)
       spl = interp1d(x,y,kind='cubic')
       y0 = spl(x0)
       a.plot(x0,y0,color=color)
    elif iline=='line':
       a.plot(x,y,color=color)


    if yax == None:
       yax = [min(y),1.2*max(y),(1.2*max(y)-min(y))/4]
       num = int(float(("%5.3e"%max(y)).split("e")[0]))+1
       exp = int(float(("%5.3e"%max(y)).split("e")[-1]))
       yax = [0.0,num*pow(10.0,exp),pow(10.0,exp)]
    a.axis(xax[0:2]+yax[0:2])

    x_majL = matplotlib.ticker.MultipleLocator(xax[2])
    x_minL = matplotlib.ticker.MultipleLocator(xax[2]/nxtick)
    y_majL = matplotlib.ticker.MultipleLocator(yax[2])
    y_minL = matplotlib.ticker.MultipleLocator(yax[2]/nytick)

    a.xaxis.set_major_locator(x_majL)
    a.xaxis.set_minor_locator(x_minL)
    a.yaxis.set_major_locator(y_majL)
    a.yaxis.set_minor_locator(y_minL)

    if xlab: a.set_xlabel(xlab)
    if ylab: a.set_ylabel(ylab)

    if iref==1:
       a.plot(xax[0:2],yax[0:2],'k--')

if __name__=="__main__":

    pass

#    dat = Namelist.Namelist("dumpdata")
#    li = dat["scan"]["li"]
#    j1 = dat["scan"]["j1"]
#
#    p = zplot("test.pdf",size_x=4.0,size_y=4.0)
# 
#    plot_s(111,j1,li,[0.5,2.5,0.5],[0.0,1.5,0.5],ispl=True,xlab='j1',ylab='li')
#    p.savefig()
#    
#    p.close()



