import numpy
from numpy import *
import math
import matplotlib.pyplot as plt

# stole this from http://scatterlib.wikidot.com/mie: credit to Herbert Kaiser (University of Konstanz, Germany)
# and to Bohren and Huffman 

def bhmie(x,refrel,nang):
# This file is converted from mie.m, see http://atol.ucsd.edu/scatlib/index.htm
# Bohren and Huffman originally published the code in their book on light scattering

# Calculation based on Mie scattering theory  
# input:
#      x      - size parameter = k*radius = 2pi/lambda * radius   
#                   (lambda is the wavelength in the medium around the scatterers)
#      refrel - refraction index (n in complex form for example:  1.5+0.02*i;
#      nang   - number of angles for S1 and S2 function in range from 0 to pi/2
# output:
#        S1, S2 - funtion which correspond to the (complex) phase functions - defined page 72, better name is the scattering diagram 
#        Qext   - extinction efficiency
#        Qsca   - scattering efficiency 
#        Qback  - backscatter efficiency
#        gsca   - asymmetry parameter


    nmxx=150000
    
    s1_1=numpy.zeros(nang,dtype=numpy.complex128)
    s1_2=numpy.zeros(nang,dtype=numpy.complex128)
    s2_1=numpy.zeros(nang,dtype=numpy.complex128)
    s2_2=numpy.zeros(nang,dtype=numpy.complex128)
    pi=numpy.zeros(nang)  # pi and tau from 4.46 -> pre-define an array  
    tau=numpy.zeros(nang)
    
    if (nang > 1000):
        print ('error: nang > mxnang=1000 in bhmie')
        return
    
    # Require NANG>1 in order to calculate scattering intensities
    if (nang < 2):
        nang = 2
    
    pii = 4.*numpy.arctan(1.)  
    dx = x
     
    drefrl = refrel
    y = x*drefrl
    ymod = abs(y)
    
    
    #    Series expansion terminated after NSTOP terms
    #    Logarithmic derivatives calculated from NMX on down
    
    xstop = x + 4.*math.pow(x,0.3333) + 2.0  # rough rule of thumb is that need ~x terms for convergence 
    
    nmx = max(xstop,ymod) + 15.0
    nmx=numpy.fix(nmx)
     
    # BTD experiment 91/1/15: add one more term to series and compare resu<s
    #      NMX=AMAX1(XSTOP,YMOD)+16
    # test: compute 7001 wavelen>hs between .0001 and 1000 micron
    # for a=1.0micron SiC grain.  When NMX increased by 1, only a single
    # computed number changed (out of 4*7001) and it only changed by 1/8387
    # conclusion: we are indeed retaining enough terms in series!
    
    nstop = int(xstop)
    
    if (nmx > nmxx):
        print ( "error: nmx > nmxx=%f for |m|x=%f" % ( nmxx, ymod) )
        return
    
    dang = .5*pii/ (nang-1)
    

    amu=numpy.arange(0.0,nang,1)
    amu=numpy.cos(amu*dang)

    pi0=numpy.zeros(nang)
    pi1=numpy.ones(nang)
    
    # Logarithmic derivative D(J) calculated by downward recurrence
    # beginning with initial value (0.,0.) at J=NMX
    
    nn = int(nmx)-1 
    d=numpy.zeros(nn+1,dtype=numpy.complex128)  # there was a bug here -> 
    for n in range(0,nn):                       # was casting the complex values to real and removing 
        en = nmx - n                            # the imaginary part
        d[nn-n-1] = (en/y) - (1./ (d[nn-n]+en/y))  
    
    #*** Riccati-Bessel functions with real argument X
    #    calculated by upward recurrence
    
    psi0 = numpy.cos(dx)
    psi1 = numpy.sin(dx) 
    chi0 = -numpy.sin(dx)
    chi1 = numpy.cos(dx)
    xi1 = psi1-chi1*1j
    qsca = 0.
    gsca = 0.
    p = -1
    
    for n in range(0,nstop):
        en = n+1.0
        fn = (2.*en+1.)/(en* (en+1.))
    
    # for given N, PSI  = psi_n        CHI  = chi_n
    #              PSI1 = psi_{n-1}    CHI1 = chi_{n-1}
    #              PSI0 = psi_{n-2}    CHI0 = chi_{n-2}
    # Calculate psi_n and chi_n
        psi = (2.*en-1.)*psi1/dx - psi0
        chi = (2.*en-1.)*chi1/dx - chi0
        xi = psi-chi*1j
    
    #*** Store previous values of AN and BN for use
    #    in computation of g=<cos(theta)>
        if (n > 0):
            an1 = an
            bn1 = bn
    
    #*** Compute AN and BN:
        an = (d[n]/drefrl+en/dx)*psi - psi1
        an = an/ ((d[n]/drefrl+en/dx)*xi-xi1)
        bn = (drefrl*d[n]+en/dx)*psi - psi1
        bn = bn/ ((drefrl*d[n]+en/dx)*xi-xi1)

    #*** Augment sums for Qsca and g=<cos(theta)>
        qsca += (2.*en+1.)* (abs(an)**2+abs(bn)**2)
        gsca += ((2.*en+1.)/ (en* (en+1.)))*( numpy.real(an)* numpy.real(bn)+numpy.imag(an)*numpy.imag(bn))
    
        if (n > 0):
            gsca += ((en-1.)* (en+1.)/en)*( numpy.real(an1)* numpy.real(an)+numpy.imag(an1)*numpy.imag(an)+numpy.real(bn1)* numpy.real(bn)+numpy.imag(bn1)*numpy.imag(bn))
    
    
    #*** Now calculate scattering intensity pattern
    #    First do angles from 0 to 90
        pi=0+pi1    # 0+pi1 because we want a hard copy of the values
        tau=en*amu*pi-(en+1.)*pi0
        s1_1 += fn* (an*pi+bn*tau)
        s2_1 += fn* (an*tau+bn*pi)
    
    #*** Now do angles greater than 90 using PI and TAU from
    #    angles less than 90.
    #    P=1 for N=1,3,...% P=-1 for N=2,4,...
    #   remember that we have to reverse the order of the elements
    #   of the second part of s1 and s2 after the calculation
        p = -p
        s1_2+= fn*p* (an*pi-bn*tau)
        s2_2+= fn*p* (bn*pi-an*tau)

        psi0 = psi1
        psi1 = psi
        chi0 = chi1
        chi1 = chi
        xi1 = psi1-chi1*1j
    
    #*** Compute pi_n for next value of n
    #    For each angle J, compute pi_n+1
    #    from PI = pi_n , PI0 = pi_n-1
        pi1 = ((2.*en+1.)*amu*pi- (en+1.)*pi0)/ en
        pi0 = 0+pi   # 0+pi because we want a hard copy of the values
    
    #*** Have summed sufficient terms.
    #    Now compute QSCA,QEXT,QBACK,and GSCA

    #   we have to reverse the order of the elements of the second part of s1 and s2
    s1=numpy.concatenate((s1_1,s1_2[-2::-1]))
    s2=numpy.concatenate((s2_1,s2_2[-2::-1]))
    gsca = 2.*gsca/qsca
    qsca = (2./ (dx*dx))*qsca
    qext = (4./ (dx*dx))* numpy.real(s1[0])

    # more common definition of the backscattering efficiency,
    # so that the backscattering cross section really
    # has dimension of length squared
    qback = 4*(abs(s1[2*nang-2])/dx)**2    
    #qback = ((abs(s1[2*nang-2])/dx)**2 )/pii  #old form

    return s1,s2,qext,qsca,qback,gsca

lamb = 532*1e-9 # wavelength 
r = 2*1e-6 # aerosol radius 

x = (2*math.pi*r)/lamb


nr = numpy.arange(1,1.5,0.01)
nI = numpy.arange(0.01,0.5,0.01)
nang = 2

refrel = 1.29 + 0.05*1j

S1 = bhmie(x,refrel,nang)[0]
S2 = bhmie(x,refrel,nang)[1]
Qext = bhmie(x,refrel,nang)[2]
Qsca = bhmie(x,refrel,nang)[3]
Qback = bhmie(x,refrel,nang)[4]
Gsca = bhmie(x,refrel,nang)[5]

xtest = numpy.arange(0.1,25,0.1)
Qexttest = numpy.zeros(len(xtest))
Qscatest = numpy.zeros(len(xtest))
Qexttest1 = numpy.zeros(len(xtest))
Qscatest1 = numpy.zeros(len(xtest))
Qexttest2 = numpy.zeros(len(xtest))
Qscatest2 = numpy.zeros(len(xtest))
Qexttest3 = numpy.zeros(len(xtest))
Qscatest3 = numpy.zeros(len(xtest))
Qexttest4 = numpy.zeros(len(xtest))
Qscatest4 = numpy.zeros(len(xtest))



refrelr = [1 + 0.05*1j,1.1 + 0.05*1j,1.2 + 0.05*1j,1.3 + 0.05*1j,1.4 + 0.05*1j]
refreli = [1.29 + 0.0*1j,1.29 + 0.01*1j,1.29 + 0.02*1j,1.29 + 0.03*1j,1.29 + 0.04*1j]

for i in range(len(xtest)):
    Qexttest[i] = bhmie(xtest[i],refreli[0],nang)[2]
    Qscatest[i] = bhmie(xtest[i],refreli[0],nang)[3]
    Qexttest1[i] = bhmie(xtest[i],refreli[1],nang)[2]
    Qscatest1[i] = bhmie(xtest[i],refreli[1],nang)[3]
    Qexttest2[i] = bhmie(xtest[i],refreli[2],nang)[2]
    Qscatest2[i] = bhmie(xtest[i],refreli[2],nang)[3]
    Qexttest3[i] = bhmie(xtest[i],refreli[3],nang)[2]
    Qscatest3[i] = bhmie(xtest[i],refreli[3],nang)[3]
    Qexttest4[i] = bhmie(xtest[i],refreli[4],nang)[2]
    Qscatest4[i] = bhmie(xtest[i],refreli[4],nang)[3]


#for i in range(len(xtest)):
#    Qexttest[i] = bhmie(xtest[i],refrel,nang)[2]
#    Qscatest[i] = bhmie(xtest[i],refrel,nang)[3]


# =============================================================================
fig = plt.figure()
ax = plt.subplot(111)
label1='Qext' + str(refreli[0])
label2='Qsca' + str(refreli[0])
label3='Qext' + str(refreli[1])
label4='Qsca' + str(refreli[1])
label5='Qext' + str(refreli[2])
label6='Qsca' + str(refreli[2])
label7='Qext' + str(refreli[3])
label8='Qsca' + str(refreli[3])
label9='Qext' + str(refreli[4])
label10='Qsca' + str(refreli[4])

ax.plot(xtest, Qexttest,'*', label=label1 )
#ax.plot(xtest, Qscatest,'o', label=label2 )
ax.plot(xtest, Qexttest1,'*', label=label3 )
#ax.plot(xtest, Qscatest1,'o', label=label4 )
ax.plot(xtest, Qexttest2,'*', label=label5 )
#ax.plot(xtest, Qscatest2,'o', label=label6 )
ax.plot(xtest, Qexttest3,'*', label=label7 )
#ax.plot(xtest, Qscatest3,'o', label=label8 )
ax.plot(xtest, Qexttest4,'*', label=label9 )
#ax.plot(xtest, Qscatest4,'o', label=label10 )

plt.xlabel('x')
plt.title('')
namelabel = 'imagnvariationext.png'
ax.legend(prop={'size': 7})
plt.savefig(namelabel, bbox_inches='tight')
plt.show
# =============================================================================
