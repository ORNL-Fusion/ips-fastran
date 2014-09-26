#!/usr/bin/env Python

"""
 ======================================================================
 Last modified : unknown
 Quick fastran plot
 JM
"""

import sys,os,re,glob,shutil,pickle
from numpy import *
import netCDF4

import matplotlib
from pylab import *
    
matplotlib.rcParams['axes.linewidth']=0.5
matplotlib.rcParams['pdf.fonttype']=42
from matplotlib.backends.backend_pdf import PdfPages

#----------------------------------------------------------------------
# single plot template

def plota(pos,x,y,xax,yax,nxtick=2,nytick=2,color=10*['k']):

    a = subplot(pos)
    for k in range(len(x)):
        a.plot(x[k],y[k],color[k])

    a.axis(xax[0:2]+yax[0:2])

    x_majL = matplotlib.ticker.MultipleLocator(xax[2])
    x_minL = matplotlib.ticker.MultipleLocator(xax[2]/nxtick)
    y_majL = matplotlib.ticker.MultipleLocator(yax[2])
    y_minL = matplotlib.ticker.MultipleLocator(yax[2]/nytick)

    a.xaxis.set_major_locator(x_majL)
    a.xaxis.set_minor_locator(x_minL)
    a.yaxis.set_major_locator(y_majL)
    a.yaxis.set_minor_locator(y_minL)


def plot_fastran(title,ncfile,iconv=False):

    # read fastran nc file

    fastran = netCDF4.Dataset(ncfile,'r',format='NETCDF4')

    rho = fastran.variables["rho"][:]
    ne = fastran.variables["ne"][:,:]
    ni = fastran.variables["ni"][:,:]
    te = fastran.variables["te"][:,:]
    ti = fastran.variables["ti"][:,:]
    j_tot = fastran.variables["j_tot"][:,:]
    j_bs = fastran.variables["j_bs"][:,:]
    j_nb = fastran.variables["j_nb"][:,:]
    j_ec = fastran.variables["j_rf"][:,:]
    q = fastran.variables["q"][:,:]
    omega = fastran.variables["omega"][:,:]
    chie = fastran.variables["chie"][:,:]
    chii = fastran.variables["chii"][:,:]
    chieneo = fastran.variables["chieneo"][:,:]
    chiineo = fastran.variables["chiineo"][:,:]

    fastran.close()
    
    # open pdf file

    pdf = PdfPages('fastran.pdf')

    # set figure

    figure(figsize=(8.5,11.0))

    figtext(0.5,0.975,
        title,
        horizontalalignment='center',
        verticalalignment='top')

    # plot fastran profile
     
    x = [ rho,rho ] 
    y = [ ne[-1],ni[-1] ] 
    xax = [0,1.0,0.2]
    yax = [0,10.0,2.0]
    color = ['r','b']
    plota(321,x,y,xax,yax,color=color)

    x = [ rho,rho ] 
    y = [ te[-1],ti[-1] ] 
    xax = [0,1.0,0.2]
    yax = [0,10.0,2.0]
    color = ['r','b']
    plota(322,x,y,xax,yax,color=color)

    x = [ rho,rho,rho ] 
    y = [ j_tot[-1],j_nb[-1], j_ec[-1] ] 
    xax = [0,1.0,0.2]
    yax = [0,1.5,0.25]
    color = ['r','b','k']
    plota(323,x,y,xax,yax,color=color)

    x = [ rho,rho ] 
    y = [ j_bs[-1],j_tot[-1] ] 
    xax = [0,1.0,0.2]
    yax = [0,1.5,0.25]
    color = ['r','b']
    plota(324,x,y,xax,yax,color=color)

    x = [ rho,rho ] 
    y = [ q[-1],q[0] ] 
    xax = [0,1.0,0.2]
    yax = [0,7.0,1.0]
    color = ['r','k']
    plota(325,x,y,xax,yax,color=color)

    x = [ rho ] 
    y = [ 1.0e-3*omega[-1] ] 
    xax = [0,1.0,0.2]
    yax = [0,200.0,50.0]
    color = ['k']
    plota(326,x,y,xax,yax,color=color)

    pdf.savefig()
    clf()

    x = [ rho, rho ] 
    y = [ chie[-1],chii[-1] ] 
    xax = [0,1.0,0.2]
    yax = [0,5.0,1.0]
    color = ['r','b']
    plota(321,x,y,xax,yax,color=color)

    x = [ rho ] 
    y = [ chieneo[-1] ] 
    xax = [0,1.0,0.2]
    yax = [0,1.0,0.2]
    color = ['r']
    plota(323,x,y,xax,yax,color=color)

    x = [ rho ] 
    y = [ chiineo[-1] ] 
    xax = [0,1.0,0.2]
    yax = [0,1.0,0.2]
    color = ['b']
    plota(324,x,y,xax,yax,color=color)

    pdf.savefig()
    clf()

    pdf.close()

def comp_fastran(key,dirs,ymax=1.5):

    # open pdf file

    pdf = PdfPages('fastran.pdf')

    title = key

    colors = ['k','r','b']

    # set figure

    figure(figsize=(8.5,11.0))

    figtext(0.5,0.975,
        title,
        horizontalalignment='center',
        verticalalignment='top')

    for k, dir in enumerate(dirs):

        # read fastran nc file

        ncfile = dir+'/fastran.nc'

        fastran = netCDF4.Dataset(ncfile,'r',format='NETCDF4')

        rho = fastran.variables["rho"][:]
        profile = fastran.variables[key][:,:]
        fastran.close()

        # overlay

        plota(111,[rho],[profile[-1]],[0.0,1.0,0.2],[0.0,ymax,0.2*ymax],color=[colors[k]])

    pdf.savefig()
    clf()

    pdf.close()

def get_dirs(runs):

    listdir = []

    for run in runs:
        dirs = glob.glob(run+"/simulation_results/iteration_*/components/fastran_tr_fastran_solver_*")
        dirs.sort()
        listdir.append(dirs[-1])
    print listdir
    return listdir


if __name__ == "__main__":

    #listdir = get_dirs(['run00','run01','run02'])
    #comp_fastran("j_tot",listdir,ymax=1.5)
    #comp_fastran("q",listdir,ymax=10.0)

    #dir =  'run01/work/plasma_state/'
    #ncfile = dir+'fastran_140000_fastran.nc' 

    #ncfile = 'fastran.nc'
    #title = 'fastran'
    #plot_fastran(title,ncfile,iconv=True)

    outdir = 'run51'
    list_dir = glob.glob(outdir+"/simulation_results/iteration_*/components/fastran_tr_fastran_solver_*")
    list_dir.sort()
    print list_dir

    ncfile = list_dir[-1]+'/fastran.nc'
    title = 'fastran'
    plot_fastran(title,ncfile,iconv=True)

#   comp_fastran("j_tot",list_dir,ymax=1.5)
#   comp_fastran("q",list_dir,ymax=10.0)

#    ncfile = 'fastran.nc'
#    title = 'fastran'
#    plot_fastran(title,ncfile,iconv=True)
