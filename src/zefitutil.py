#!/usr/bin/env python

#-----------------------------------------------------------------------------#
# efit utility
# JM
#-----------------------------------------------------------------------------#

from numpy import *

def line2vec(line,nlen=16):
    ncol = len(line)/nlen
    vec = zeros(ncol)
    for k in range(ncol): 
        tmp = line[k*nlen:(k+1)*nlen]
        vec[k] = float(tmp)
    return vec

def readvec(lines,k0,nread):
    tmp = []
    for k in range(k0,k0+nread):
        vec = line2vec(lines[k])
        for i in range(len(vec)): tmp.append(vec[i])
    return array(tmp),k0+nread

def get_nread(n,ncol=5):
    nr= n/ncol
    if n%ncol>0: nr+=1
    return nr

def readg(fname):
    
    #-------------------------------------
    # read

    f = open(fname,"r")
    lines = f.readlines()
    f.close()

    #-------------------------------------
    # parse

    header = [ line2vec(lines[k]) for k in range(1,5) ]

    nw     = int(lines[0].split()[-2])
    nh     = int(lines[0].split()[-1])

    xdim   = header[0][0]
    ydim   = header[0][1]
    rzero  = header[0][2]
    rgrid  = header[0][3]
    zmid   = header[0][4]

    rmaxis = header[1][0]
    zmaxis = header[1][1]
    ssimag = header[1][2]
    ssibry = header[1][3]
    bcentr = header[1][4]

    cpasma = header[2][0]
    ssimag = header[2][1]
    rmaxis = header[2][3]

    zmaxis = header[3][0]
    ssibry = header[3][2]

    k = 5

    nread = get_nread(nw)
    fpol,k = readvec(lines,k,nread)
    pres,k = readvec(lines,k,nread)
    ffprim,k = readvec(lines,k,nread)
    pprime,k = readvec(lines,k,nread)

    nread = get_nread(nw*nh)
    tmp,k = readvec(lines,k,nread)
    psirz = zeros([nw,nh])
    kk = 0
    for ii in range(nw):
        for jj in range(nh):
            psirz[ii,jj] =tmp[kk]
            kk+=1

    nread = get_nread(nw)
    qpsi,k = readvec(lines,k,nread)

    nbdry = int(lines[k].split()[-2])
    nlim = int(lines[k].split()[-1])
    k+=1

    nread = get_nread(nbdry*2)
    tmp,k = readvec(lines,k,nread)
    rb = zeros(nbdry)
    zb = zeros(nbdry)
    for ii in range(nbdry):
        rb[ii] = tmp[2*ii]
        zb[ii] = tmp[2*ii+1]

    nread = get_nread(nlim*2)
    tmp, k = readvec(lines,k,nread)
    rlim = zeros(nlim)
    zlim = zeros(nlim)
    for ii in range(nlim):
        rlim[ii] = tmp[2*ii]
        zlim[ii] = tmp[2*ii+1]

    #-------------------------------------
    # return dictionary

    geq = {}

    geq["nw"    ]  = nw
    geq["nh"    ]  = nh
    geq["xdim"  ]  = xdim
    geq["ydim"  ]  = ydim
    geq["rzero" ]  = rzero
    geq["rgrid" ]  = rgrid
    geq["zmid"  ]  = zmid
    geq["rmaxis"]  = rmaxis
    geq["zmaxis"]  = zmaxis
    geq["ssimag"]  = ssimag
    geq["ssibry"]  = ssibry
    geq["bcentr"]  = bcentr
    geq["cpasma"]  = cpasma
    geq["ssimag"]  = ssimag
    geq["rmaxis"]  = rmaxis
    geq["zmaxis"]  = zmaxis
    geq["ssibry"]  = ssibry
    geq["fpol"  ]  = fpol
    geq["pres"  ]  = pres
    geq["ffprim"]  = ffprim
    geq["pprime"]  = pprime
    geq["psirz" ]  = psirz
    geq["qpsi"  ]  = qpsi
    geq["nbdry" ]  = nbdry
    geq["nlim"  ]  = nlim
    geq["rbdry" ]  = rb
    geq["zbdry" ]  = zb
    geq["rlim"  ]  = rlim
    geq["zlim"  ]  = zlim

    return geq

def getbdry(geq,nskip=2):

    nbdry = geq['nbdry']
    rb = geq['rbdry']
    zb = geq['zbdry']

    rb0 = []
    zb0 = []
    for i in range(0,nbdry,nskip):
        rb0.append(rb[i])
        zb0.append(zb[i]) 

    nlim = geq['nlim']
    rlim = geq['rlim']
    zlim = geq['zlim']

    return {'nbdry':len(rb0),
            'rbdry':rb0,
            'zbdry':zb0,
            'nlim':nlim,
            'rlim':rlim,
            'zlim':zlim}

#def writeg(geq,fname):
#
#    f=open(fname,"w")
#
#    header = 6*[""]
#    idum = 0
#    for k in range(6): f.write("%8s"%header[i])
#    f.write("%4d"%idum)
#    f.write("%4d"%geq["nw"])
#    f.write("%4d"%geq["nh"])
#
#
#    write (lun,2020) xdim,zdim,rzero,rgrid,zmid
#    write (lun,2020) rmaxis,zmaxis,ssimag,ssibry,bcentr
#    write (lun,2020) cpasma,ssimag,xdum,rmaxis,xdum
#    write (lun,2020) zmaxis,xdum,ssibry,xdum,xdum
#    write (lun,2020) (fpol(i),i=1,nw)
#    write (lun,2020) (pres(i),i=1,nw)
#    write (lun,2020) (ffprim(i),i=1,nw)
#    write (lun,2020) (pprime(i),i=1,nw)
#    write (lun,2020) ((psirz(i,j),i=1,nw),j=1,nh)
#    write (lun,2020) (qpsi(i),i=1,nw)
#    write (lun,2022) nbbbs,limitr
#    write (lun,2020) (rbbbs(i),zbbbs(i),i=1,nbbbs)
#    write (lun,2020) (xlim(i),ylim(i),i=1,limitr)

def reada(fname):
    
    #-------------------------------------
    # read

    f = open(fname,"r")
    lines = f.readlines()
    f.close()

    #-------------------------------------
    # parse
    for k in range(4,7): print lines[k]

    #data = [ line2vec(lines[k]) for k in range(4,5) ]

    r0    = line2vec(lines[4 ][1:])[1]
    b0    = line2vec(lines[4 ][1:])[2]
    a0    = line2vec(lines[5 ][1:])[3]
    ip    = line2vec(lines[5 ][1:])[0]
    betat = line2vec(lines[7 ][1:])[3]
    qaxis = line2vec(lines[17][1:])[2]
    betan = fabs(ip*1.0e-6)/(0.01*a0*abs(b0))

    print r0,a0 
    print b0,ip
    print betat,betan 

#    read (neqdsk,1060) time(jj),jflag(jj),lflag,limloc(jj),
# .                      mco2v,mco2r,qmflag,nlold,nlnew
#    read (neqdsk,1040) tsaisq(jj),rcencm,bcentr(jj),pasmat(jj)          4
#    read (neqdsk,1040) cpasma(jj),rout(jj),zout(jj),aout(jj)            5 
#    read (neqdsk,1040) eout(jj),doutu(jj),doutl(jj),vout(jj)            6
#    read (neqdsk,1040) rcurrt(jj),zcurrt(jj),qsta(jj),betat(jj)         7
#    read (neqdsk,1040) betap(jj),ali(jj),oleft(jj),oright(jj)           8
#    read (neqdsk,1040) otop(jj),obott(jj),qpsib(jj),vertn(jj)           9
#    read (neqdsk,1040) (rco2v(k,jj),k=1,mco2v) 10
#    read (neqdsk,1040) (dco2v(jj,k),k=1,mco2v) 11 
#    read (neqdsk,1040) (rco2r(k,jj),k=1,mco2r) 12
#    read (neqdsk,1040) (dco2r(jj,k),k=1,mco2r) 13 
#    read (neqdsk,1040) shearb(jj),bpolav(jj),s1(jj),s2(jj)      14
#    read (neqdsk,1040) s3(jj),qout(jj),olefs(jj),orighs(jj)      15
#    read (neqdsk,1040) otops(jj),sibdry(jj),areao(jj),wplasm(jj)  16
#    read (neqdsk,1040) terror(jj),elongm(jj),qqmagx(jj),cdflux(jj) 17
#    read (neqdsk,1040) alpha(jj),rttt(jj),psiref(jj),xndnt(jj)
#    read (neqdsk,1040) rseps(1,jj),zseps(1,jj),rseps(2,jj),zseps(2,jj)
#    read (neqdsk,1040) sepexp(jj),obots(jj),btaxp(jj),btaxv(jj)
#    read (neqdsk,1040) aaq1(jj),aaq2(jj),aaq3(jj),seplim(jj)
#    read (neqdsk,1040) rmagx(jj),zmagx(jj),simagx(jj),taumhd(jj)
#    read (neqdsk,1040) betapd(jj),betatd(jj),wplasmd(jj),diamag(jj)
#    read (neqdsk,1040) vloopt(jj),taudia(jj),qmerci(jj),tavem(jj)
#    read (neqdsk,1041) nsilop0,magpri0,nfcoil0,nesum0
#    read (neqdsk,1040) ((csilop(k,jj),k=1,nsilop0),
# .                      (cmpr2(k,jj),k=1,magpri0))
#    read (neqdsk,1040) (ccbrsp(k,jj),k=1,nfcoil0)
#    read (neqdsk,1040) (eccurt(jj,k),k=1,nesum0)
#    read (neqdsk,1040) pbinj(jj),rvsin(jj),zvsin(jj),rvsout(jj)
#    read (neqdsk,1040) zvsout(jj),vsurfa(jj),wpdot(jj),wbdot(jj)
#    read (neqdsk,1040) slantu(jj),slantl(jj),zuperts(jj),chipre
#    read (neqdsk,1040) cjor95(jj),pp95(jj),ssep(jj),yyy2(jj)
#    read (neqdsk,1040) xnnc(jj),cprof,oring(jj),cjor0(jj)
#    read (neqdsk,1040) fexpan(jj),qqmin(jj),chigamt(jj),ssi01(jj)
#    read (neqdsk,1040) fexpvs(jj),sepnose(jj),ssi95(jj),rqqmin(jj)
#    read (neqdsk,1040) cjor99(jj),cj1ave(jj),rmidin(jj),rmidout(jj)
#    read (neqdsk,1040) psurfa(jj),peak(jj),dminux(jj),dminlx(jj)
#    read (neqdsk,1040) dolubaf(jj),dolubafm(jj),diludom(jj),),diludomm(jj)
#    read (neqdsk,1040) ratsol(jj),rvsiu(jj),zvsiu(jj),rvsid(jj)
#    read (neqdsk,1040) zvsid(jj),rvsou(jj),zvsou(jj),rvsod(jj)
#    read (neqdsk,1040) zvsod(jj),condno,psin32(jj),psin21(jj)
#    read (neqdsk,1040) rq32in(jj),rq21top(jj),chilibt,ali3(jj)
#    read (neqdsk,1040) xbetapr,tflux(jj),xdum,xdum


########################################################################
#  Check

if __name__ == "__main__":

    readg('g000000.00000')





