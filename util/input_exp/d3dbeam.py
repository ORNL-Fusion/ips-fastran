#! /usr/bin/env python

"""
from Tom Osborn's d3d_beams.py
"""

import numpy

###############################################################################################
#
#  THE FOLLOWING ROUTINES COMPUTE THE ION SPECIES DISTRIBUTION
#  THESE ARE A TRANSLATION OF CHUCK GREENFIELD'S IDL ROUTINE d3d_beam_species_mix.pro TO PYTHON
#  WHICH IN TURN WERE ORIGINALLY BASED ON FORTRAN ROUTINES BY JINCHOON KIM.
#
#       Jinchoon's comments:
#               PROGRAM GASEL3
#               PROGRAM AUTHOR J. KIM
#               GENERAL ATOMIC COMPANY, 1980.
#               This program calculates the fractions of all the
#               possible species of ions and neutrals as a funsction
#               of gas cell line density at any energy between 10 keV
#               and 200 keV per amu. The users are asked to input the
#               energy per amu and the species fractions.
###############################################################################################

def beam_species_mix(
    # Compute the beam ion energy species mix.
    # Values are given both before and after the neutralizer
    zkev,                  # beam voltage
    abeam         = 2.0,   # atomic mass of beam ion
    xzbeam        = 1.0,   # z of beam ion
    noneutralizer = False, # if true don't calculate effect of neutralizer
    nl            = 1.0,   # neutralizer density
    # If noneutralizer = false returns a dictionary where 'before'
    # is before the neutralizer and 'after' is after neutralizer
    # and values are [full energy, half energy, third energy]. If
    # noneutralizer = true only an array is returned with values
    # before the neutralizer.
    ):
    """
    d3d_beam_species_mix(
    # Compute the beam ion energy species mix.
    # Values are given both before and after the neutralizer
    zkev,                  # beam voltage
    abeam         = 2.0,   # atomic mass of beam ion
    xzbeam        = 1.0,   # z of beam ion
    noneutralizer = false, # if true don't calculate effect of neutralizer
    nl            = 1.0,   # neutralizer density
    # If noneutralizer = false returns a dictionary where 'before'
    # is before the neutralizer and 'after' is after neutralizer
    # and values are [full energy, half energy, third energy]. If
    # noneutralizer = true only an array is returned with values
    # before the neutralizer.
    )
    """
    # Check for helium added by CMG 5/15/96
    if (xzbeam > 1.0):
        return numpy.array([1.0,0.0,0.0])
    # Use fit for species mix at accelerator from Steve Riggs' memo 3/3/93
    smixin=numpy.zeros(3,'d')
    smixin[0]=(-7.737629E-5*zkev+1.374258e-2)*zkev+0.1854158
    smixin[1]=( 4.049498e-6*zkev-3.088367e-4)*zkev+0.1544813
    smixin[2]=1.0-smixin[0]-smixin[1]
    # Use Jinchoon Kim's neutralizer model to get species mix into tokamak
    if noneutralizer:return smixin
    return { 'after' :_neutralizer(smixin,zkev,abeam,nl),
             'before': smixin                             }
###############################################################################################
def _xsect(
    # GET BEAM NEUTRALIZATION CROSS-SECTIONS
    eb # BEAM EFFECTIVE ENERGY IN KEV = BEAM VOLTAGE/ATOMIC MASS
    ):
    en=numpy.array([10.,20.,30.,40.,60.,80.,100.,120.,150.,200.])
    s=numpy.array([8.50,6.0,4.0,2.6,1.1,.56,.28,.15,.07,.018,9.0,8.5,7.3,6.0,
                     4.0,2.6,1.7,1.1,.65,.28,9.0,8.95,8.5,7.5,6.0,4.2,3.6,2.6, 
                     1.7,0.9,0.92,1.3,1.6,1.6,1.5,1.2,1.1,0.95,0.87,0.7,0.82, 
                     0.92,1.15,1.3,1.6,1.6,1.5,1.4,1.25,1.1,0.77,0.85,0.92,1.1, 
                     1.3,1.55,1.6,1.6,1.5,1.45,4.6,3.6,2.8,2.1,1.2,0.7,0.45,0.3, 
                     0.15,0.08,5.0,4.2,3.6,2.9,2.1,1.5,1.0,0.7,0.44,0.22,0.67, 
                     1.2,1.5,1.7,1.8,2.0,2.0,2.0,2.0,2.1,0.55,0.95,1.2,1.4,1.7, 
                     1.8,1.9,2.0,2.0,2.0,0.38,0.43,0.60,0.67,0.67,0.49,0.30,0.3, 
                     0.25,0.25,0.25,0.4,0.43,0.55,0.67,0.69,0.65,0.49,0.3,0.3, 
                     2.4,1.98,1.63,1.28,1.0,0.85,0.75,0.67,0.6,0.51,2.6,2.3, 
                     1.98,1.67,1.28,1.1,0.9,0.85,0.75,0.82,2.4,2.1,1.9,1.9,2.2, 
                     2.4,2.3,2.25,2.1,1.8,2.4,2.3,2.1,2.0,1.9,2.15,2.3,2.4,2.3, 
                     2.2,7.7,8.5,8.2,7.7,6.1,4.6,3.6,2.5,2.0,1.2,6.7,8.1,8.5, 
                     8.2,7.7,6.6,5.6,4.6,3.6,2.2,1.1,1.25,1.3,1.25,1.15,1.1,0.99, 
                     0.92,0.87,0.79,3.4,4.1,4.3,4.1,3.15,2.4,1.9,1.6,1.2,0.7, 
                     0.86,1.65,1.8,1.95,2.15,2.3,2.3,2.3,2.4,2.4,6.2,8.8,9.6, 
                     9.6,9.0,7.9,6.7,5.2,4.5,3.6,0.047,0.1,0.062,0.025,2.2E-3, 
                     4.6E-4,1.0E-4,3.0E-5,1.0E-5,1.5E-6,.22,.19,.12,.095,.053, 
                     0.029,0.021,8.E-3,4.5E-3,1.4E-3,.4,.44,.44,.43,.4,.34,.3, 
                     0.25,0.22,0.17,10.0,8.7,7.4,6.7,5.4,4.7,4.0,3.6,3.2,2.6])
    s.shape=(26,10)
    if eb < en[0]: 
        print("EBEAM=%10.2f keV, < MIN IN CROSSECTION TABLE=%10.2f, USING EBEAM=%10.2f" %(eb,en[0],en[0]))
        eb = en[0]
    if eb > en[-1]: 
        print("EBEAM=%10.2f keV, > MAX IN CROSSECTION TABLE=%10.2f, USING EBEAM=%10.2f" %(eb,en[-1],en[-1]))
        eb = en[-1]
    j2 = numpy.searchsorted( en, eb )
    j1 = j2 -1 
    xs = numpy.ravel( s[:,j1] + (eb-en[j1])/(en[j2]-en[j1])*(s[:,j2]-s[:,j1]) )
    return xs

###############################################################################################
def _neutralizer(
    # COMPUTE EFFECT OF NEUTRALIZER ON BEAM SPECIES MIX
    smixin, # SPECIES MIX BEFORE NEUTRALIZER
    zkev,   # BEAM ENERGY IN KEV
    abeam,  # BEAM ATOMIC MASS
    nl,     # NEUTRALIZER DENSITY
    ):
    p1=smixin[0]
    p2=smixin[1]
    p3=smixin[2]
    eb=zkev/abeam
    xs = _xsect(eb)
    s1=xs[0]
    s2=xs[1]
    s3=xs[2]
    s4=xs[3]
    s5=xs[4]
    s6=xs[5]
    s7=xs[6]
    s8=xs[7]
    s9=xs[8]
    s10=xs[9]
    s11=xs[10]
    s12=xs[11]
    s13=xs[12]
    s14=xs[13]
    s15=xs[14]
    s16=xs[15]
    s17=xs[16]
    s18=xs[17]
    s19=xs[18]
    s20=xs[19]
    s21=xs[20]
    s22=xs[21]
    e1=xs[22]
    e2=xs[23]
    e3=xs[24]
    e4=xs[25]
    t1=0.333*(s21+s22)+0.667*(s19+s20)
    t2=s8+0.5*(s16+s18)
    t3=s10+0.5*(s12+s14)
    t4=s7+0.5*(s15+s17)
    t5=s9+0.5*(s11+s13)
    t6=s4+e2+e4
    t7=s1+e1+e3
    r1=0.5*(t2+t3)+0.5*numpy.sqrt((t2-t3)**2+4*s8*s10)
    r2=0.5*(t2+t3)-0.5*numpy.sqrt((t2-t3)**2+4*s8*s10)
    q1=0.5*(t4+t5)+0.5*numpy.sqrt((t4-t5)**2+4*s7*s9)
    q2=0.5*(t4+t5)-0.5*numpy.sqrt((t4-t5)**2+4*s7*s9)
    u1=.5*(t6+t7)+.5*numpy.sqrt((t6-t7)**2+4*(s1-e4)*(s4-e3))
    u2=.5*(t6+t7)-.5*numpy.sqrt((t6-t7)**2+4*(s1-e4)*(s4-e3))
    rr=1./(r1-r2)
    tr1=1./(t1-r1)
    tr2=1./(t1-r2)
    ss1=s1+s4
    ss2=s2+s5
    ss3=s3+s6
    h3=(t6*e3+e4*(s4-e3))/(u1*u2)
    h1=(u2*h3-u2+s1+e1)/(u1-u2)
    h2=(-u1*h3+u1-s1-e1)/(u1-u2)
    w3=(t7*e4+e3*(s1-e4))/(u1*u2)
    w1=(u2*w3-s1)/(u1-u2)
    w2=(-u1*w3+s1)/(u1-u2)
    qq=(s2+s5-q1)*(q2-q1)
    yy=(s2+s5-q2)*(q2-q1)
    qqq=s7*(s13-2*s2)
    g=(q2-t4)*(s17-2*s2)/qq+qqq/qq
    gg=(t4-q1)*(s17-2*s2)/yy-qqq/yy
    ta=s14-2*s3
    tb=s18-2*s3
    a=rr*tr1*((r1-t2)*s20-s19*s8)
    aa=rr*tr2*(s19*s8-(r2-t2)*s20)
    aaa=tr1*tr2*(s19*s8-(t1-t2)*s20)
    b=rr*tr1*((r1-t3)*s19-s20*s10)
    bb=rr*tr2*(s20*s10-(r2-t3)*s19)
    bbb=tr1*tr2*(s20*s10-s19*(t1-t3))
    d=(ta*a+tb*b)/(ss3-r1)
    dd=(ta*aa+tb*bb)/(ss3-r2)
    ddd=(ta*aaa+tb*bbb+s22-3*s3)/(ss3-t1)
    x=nl
    ef1=numpy.exp(-ss1*x)
    ef2=numpy.exp(-ss2*x)
    ef3=numpy.exp(-ss3*x)
    fq1=numpy.exp(-q1*x)
    fq2=numpy.exp(-q2*x)
    ft1=numpy.exp(-t1*x)
    fr1=numpy.exp(-r1*x)
    fr2=numpy.exp(-r2*x)
    af1n=w1*numpy.exp(-u1*x)+w2*numpy.exp(-u2*x)+w3
    af1c=h1*numpy.exp(-u1*x)+h2*numpy.exp(-u2*x)+h3
    afnc=1-af1n-af1c
    df2n=(fq1-fq2)*s7/(q2-q1)
    df2c=(t4-q2)*fq1/(q1-q2)+(q1-t4)*fq2/(q1-q2)
    df1n=(2*s2/ss2-ef2*(g+gg+2*s2/ss2)+g*fq1+gg*fq2)*.5
    df1c=1-df1n-df2n-df2c
    tf3c=ft1
    tf2n=(a*fr1+aa*fr2+aaa*ft1)*.667
    tf2c=(b*fr1+bb*fr2+bbb*ft1)*.667
    tf1n=(3*s3/ss3-(3*s3/ss3+d+dd+ddd)*ef3+d*fr1+dd*fr2+ddd*ft1)*0.333
    tf1c=1-tf1n-tf2n-tf2c-tf3c
    pa1n=p1*af1n
    pa1c=p1*af1c
    panc=p1*afnc
    pd1n=p2*df1n
    pd1c=p2*df1c
    pd2n=p2*df2n
    pd2c=p2*df2c
    pt1n=p3*tf1n
    pt1c=p3*tf1c
    pt2n=p3*tf2n
    pt2c=p3*tf2c
    pt3c=p3*tf3c
    ptn=pa1n+pd1n+pd2n+pt1n+pt2n
    af1n=100*af1n
    af1c=100*af1c
    annc=100*afnc
    df1n=100*df1n
    df1c=100*df1c
    df2n=100*df2n
    df2c=100*df2c
    tf1n=100*tf1n
    tf1c=100*tf1c
    tf2n=100*tf2n
    tfc2=100*tf2c
    tfc3=100*tf3c
    pa1n=100*pa1n
    pa1c=100*pa1c
    panc=100*panc
    pd1n=100*pd1n
    pd1c=100*pd1c
    pd2n=100*pd2n
    pd2c=100*pd2c
    pt1n=100*pt1n
    pt1c=100*pt1c
    pt2n=100*pt2n
    ptc2=100*pt2c
    ptc3=100*pt3c
    ptnhalf=pd1n + pd2n
    ptnthird=pt1n + pt2n
    ptn=100*ptn
    full=pa1n
    half=2*(pd1n+pd2n)
    third=3*(pt1n+pt2n)
    total=full+half+third
    return numpy.array([full,half,third])/total
###############################################################################################

