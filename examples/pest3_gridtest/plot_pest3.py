#!/usr/bin/env Python

import os,sys,glob,re
import numpy as np
import matplotlib.pyplot as plt
        
pat_w = re.compile("ideal instabilitie(s)")
pat_q = re.compile("\s*safety factor")
pat_r = re.compile("\s<rs d psi/d r>")
pat_m = re.compile("D_I")
pat_l = re.compile("Lambdas")

if __name__=="__main__":

   shot=161172
   outdir='/cscratch/kimkyungjin/testgrid/%06d/pest3/test'%shot
   delta_21=[]  
   mu_21=[]
 
   for grid in range(100,1700,100):  
      outfile=os.path.join(outdir,'%d/results'%grid)
      
      lines = open(outfile,"r").readlines()
      ideal = 0
      for k, line in enumerate(lines):
          if pat_w.search(line):
              print line
              ideal = 1
      deltap = []
      qlist = []
      rlist = []
      llist = []
      mlist = []
      for k, line in enumerate(lines):
          if pat_q.search(line):
              q = int(line.split("=")[-1])
              qlist.append(q)
              i_surface = int(lines[k-1].split()[0])
#              print q,i_surface
              for line_m in lines[k+1:]:
                  if pat_m.search(line_m):
#                      print line_m
                      mlist.append(float(line_m.split("=")[-1]))
                      break
              for line_r in lines[k+3:]:
                  if pat_r.search(line_r):
#                      print line_r
                      rlist.append( float(line_r.split("=")[-1].split()[0]))
                      break
              for line_l in lines[k+4:]:
                  if pat_l.search(line_l):
#                      print line_l
                      llist.append( float(line_l.split("=")[-1]))
                      break
              for line_D in lines[k+5:]:
                  pat = re.compile("%d %d"%(i_surface,i_surface))
                  if pat.search(line_D):
#                      print line_D
                      deltap.append( float(line_D.split("=")[-1].split()[0]) )
                      break
   
#      print 'Q =', qlist
#      print 'DELTAP = ', deltap
   
      delta_21.append(deltap[qlist.index(2)])
      mu_21.append(mlist[qlist.index(2)])

#   print delta_21, mu_21
   grid=np.arange(100,1700,100)
   
   plt.subplot(211)
   plt.plot(grid,delta_21,'o-')
   plt.ylabel('deltap')
   plt.title('%06d'%shot)   

   plt.subplot(212)
   plt.plot(grid,mu_21,'o-')
   plt.ylabel('mu')
   plt.xlabel('grids')

#   plt.show()
   plt.savefig('%06d.pdf'%shot)


   
