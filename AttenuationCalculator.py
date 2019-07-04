# script to read in refractive index data from a CSV and calculate the attentuation of the beam 

import csv
import matplotlib.pyplot as plt
import numpy as np 
import math 
n = []
#z = []
T = []
P = [] 
N = [] # number density for Rayleigh scattering calc 
 
noz = 86001 # equal to (no. km x 1000) + 1
res = 0.001

with open('refractiveindex.csv') as csvDataFile:
    csvReader = csv.reader(csvDataFile)
    for row in csvReader:
        n.append(row[1])
        #z.append(row[0])

with open('temperature.csv') as csvDataFile:
    csvReader = csv.reader(csvDataFile)
    for row in csvReader:
        T.append(row[1])
        
with open('pressure.csv') as csvDataFile:
    csvReader = csv.reader(csvDataFile)
    for row in csvReader:
        P.append(row[1])
        
with open('totalnumberdensity.csv') as csvDataFile:
    csvReader = csv.reader(csvDataFile)
    for row in csvReader:
        N.append(row[1])

noz2 = 400001

z = np.arange(0,noz2*res,res)  # instead of importing z!!

del n[0]  # remove headers 
#del z[0]
del T[0]
del P[0]
del N[0]

for i in range(0,len(n)): 
    n[i] = float(n[i])

for j in range(0,len(z)):
    z[j] = float(z[j])
    
for k in range(0,len(T)):
    T[k] = float(T[k])
    
for l in range(0,len(P)):
    P[l] = float(P[l])

for i in range(0, len(N)):
    N[i] = float(N[i])
    
# constants     
T0 = 273.15
P0 = 1013.25
lamb = 532*1e-9 # wavelength in m 
lambm = 0.532 # wavelength in micrometers

delta = 0.031 # air depolarization factor -> in general will be a function of gas composition and hence z

nc1 = np.ones(noz) # for Rayliegh scattering calculation 
nc = np.array(nc1, dtype = complex)
cn = np.zeros(noz)
bm2 = np.zeros(noz)

Cn2 = np.zeros(noz2)
#wbf = np.zeros(noz2) #  for free space calculations 
#Rbf = np.zeros(noz2)
#wb = np.zeros(noz2) # for Gaussian beam in presence of turbulence 
#Rb = np.zeros(noz2)
cob = np.zeros(noz2)
phi = np.zeros(noz2)

# ============================ molecular scattering -> Rayleigh scattering ============================

for i in range(0,noz):
    
    bm2[i] = ((24*math.pi*math.pi*math.pi)/(N[i]*math.pow(lamb,4)))*math.pow(((n[i]*n[i] - 1)/(n[i]*n[i] + 1)),2)*((6 + 3*delta)/(6 - 7*delta))
    cn[i] = bm2[i]*lamb/(4*math.pi)
    nc[i] = n[i] + bm2[i]*1j  # exact formula 

#-------------------------------------- Turbulence calculations ----------------------------------------------

# ============================ SLC-Day model for turbulence ======================================
for i in range(0, 20001): # up to 20km with 1m resolution 
    if z[i] <= 0.019:
        Cn2[i] = 0
    elif z[i] > 0.019 and z[i] <= 0.230:
        Cn2[i] = 4.008*1e-13*math.pow(z[i]*1000, -1.054)
    elif z[i] > 0.230 and z[i] <= 0.850:
        Cn2[i] = 1.3*1e-15
    elif z[i] > 0.850 and z[i] <= 7:
        Cn2[i] = 6.352*1e-7*math.pow(z[i]*1000, -2.966)
    elif z[i] > 7 and z[i] <= 20:
        Cn2[i] = 6.209*1e-16*math.pow(z[i]*1000, -0.6229)

#=============================== Gurvich model for turbulence ==========================================

Cn2G = np.zeros(noz)
# NOTE: need to check region of validity of this model i.e. altitude 
  # strength of the turbulence at 2.5m 
# note: this code will break if the resolution of z changes -> should find a way to 
# improve this later 

def GurCncalc(Cn2G0):
    
    for i in range(3, noz):  # start at 3m, as the model is valid from 2.5m upwards
        if Cn2G0 > 1e-13: # strong turbulence
            if z[i] <= 1: 
                Cn2G[i] = Cn2G0*math.pow((z[i]*1000)/2.5, -4/3)
            else:
                Cn2G[i] = Cn2G[1000]*math.exp(-(z[i]*1000 - 1000)/9000)
        elif Cn2G0 > 6.5*1e-15 and Cn2G0 <= 1e-13:
            if z[i] <= 0.05:
                Cn2G[i] = Cn2G0*math.pow((z[i]*1000)/2.5, -2/3)
            elif z[i] > 0.05 and z[i] <= 1:
                Cn2G[i] = Cn2G[50]*math.pow((z[i]*1000)/50, -4/3)
            else:
                Cn2G[i] = Cn2G[1000]*math.exp(-(z[i]*1000 - 1000)/9000)
        elif Cn2G0 > 4.3*1e-16 and  Cn2G0 <= 6.5*1e-15:
            if z[i] <= 1:
                Cn2G[i] = Cn2G0*math.pow(((z[i]*1000)/2.5), -2/3)
            else: 
                Cn2G[i] = Cn2G[1000]*math.exp(-((z[i]*1000) - 1000)/9000)
        elif Cn2G0 <= 4.3*1e-16:  # weak turbulence
            if z[i] <= 1:
                Cn2G[i] = Cn2G0
            else:
                Cn2G[i] = Cn2G[1000]*math.exp(-((z[i]*1000) - 1000)/9000) 
 
    return(Cn2G)

Cn2G0 = 1.5e-13  
Cn2G = GurCncalc(Cn2G0)
Cn2G1 = GurCncalc(1e-14)
Cn2G2 = GurCncalc(1e-15)
Cn2G3 = GurCncalc(1e-16)
# Gaussian beam waist calculation -> see how width of beam varies with distance 
w0 = 0.5 # beam waist of Gaussian beam 
R0 = 1e10 # radius of curvature of the beam phase front
k = 2*math.pi/lamb # wavenumber
cobs = 1 # spatial coherence properties of the signal-carrying laser beam as it exits the transmitter 
# (ζS = 1 for a coherent beam, ζS > 1 for a partially coherent beam)
# SLC-Day stuff 
offsetSLC = 20000
offset1 = 19
sig12s = np.zeros(offsetSLC + 1)
zrecgcomps = np.zeros(offsetSLC + 1)

zrecgcomp = np.zeros(noz2) 
offsetGur1 = 2
upperlimGur = 20  # defined upper limit of turbulence model - after which turbulence assumed to be negligable 
offsetGur = upperlimGur*1000 

#------------------------------------ Beam size calcs --------------------------------------------------#
# =============================================================================
# widths = [0.1, 0.2, 0.3, 0.4, 0.5]
# # NOTE: w is the beam radius, not the diameter 
# for j in widths:
# # Free space 
#     for i in range(1, noz2): # don't calculate at transmitter, thus start from 1 
#         wbf[i] = j*math.pow((math.pow((R0 - z[i]*1000)/R0, 2) + math.pow((2*z[i]*1000)/(k*j*j), 2)), 0.5)
#         Rbf[i] = (z[i]*1000*(math.pow((R0 - z[i]*1000)/R0, 2) + math.pow((2*z[i]*1000)/(k*j*j), 2)))/(math.pow((R0 - z[i]*1000)/R0, 2)*(1 - math.pow((R0 - z[i]*1000)/R0, 2)) - math.pow((2*z[i]*1000)/(k*j*j), 2) )
#     fig = plt.figure()
#     ax = plt.subplot(111)
#     alabel='Initial Gaussian beam width = '+ str(j)
#     ax.plot(z, wbf, label=alabel )
#     plt.xlabel('Altitude / km')
#     plt.ylabel('Gaussian beam width / m')
#     plt.xlim(0,400)
#     plt.title('Width as function of altitude for free space')
#     ax.legend()
#     namelabel = 'wbfvsw0'+ str(j) + '.png'
#     plt.savefig(namelabel, bbox_inches='tight')
#     plt.show()
# =============================================================================


def freespacebeam(w0,R0,z,k):
    wbf = w0*math.pow((math.pow((R0 - z*1000)/R0, 2) + math.pow((2*z*1000)/(k*w0*w0), 2)), 0.5)
    Rbf = (z*1000*(math.pow((R0 - z*1000)/R0, 2) + math.pow((2*z*1000)/(k*w0*w0), 2)))/(math.pow((R0 - z*1000)/R0, 2)*(1 - math.pow((R0 - z*1000)/R0, 2)) - math.pow((2*z*1000)/(k*w0*w0), 2) )
    return(wbf,Rbf)

def Turbulencebeam(w0,R0,z,Cn2,k,cobs):
    cob = cobs + (2*w0*w0)/(math.pow(0.55*Cn2*k*k*z*1000, -6/5))
    beta2S = 1.52*Cn2*math.pow(z*1000,3)*math.pow(w0, -1/3)
    wb = w0*math.pow((math.pow((R0 - z*1000)/R0, 2) + cob*math.pow((2*z*1000)/(k*w0*w0), 2)), 0.5) + beta2S
    phi = ((R0 - z*1000)/R0)/((2*z*1000)/(k*w0*w0)) - ((2*z*1000)/(k*w0*w0))*(w0*w0)/(math.pow(0.55*Cn2*k*k*z*1000, -6/5))
    Rb = (z*1000*(math.pow((R0 - z*1000)/R0, 2) + cob*math.pow((2*z*1000)/(k*w0*w0), 2)))/(phi*((2*z*1000)/(k*w0*w0)) - cob*math.pow((2*z*1000)/(k*w0*w0), 2) - math.pow((R0 - z*1000)/R0, 2) )
    sig12s = 1.23*Cn2*math.pow(k,7/6)*math.pow(z, 11/6)  # variables for comparison in (79)
    zrecgcomps = math.pow(z/(0.5*k*math.pow(wb,2)),5/6)*sig12s
    return(wb,Rb,sig12s,zrecgcomps)



def freespace(w0,R0,z,k,noz2):    

    wbf = np.zeros(noz2)
    Rbf = np.zeros(noz2)
    def freespacebeam(w0,R0,z1,k):
        wbf = w0*math.pow((math.pow((R0 - z1*1000)/R0, 2) + math.pow((2*z1*1000)/(k*w0*w0), 2)), 0.5)
        Rbf = (z1*1000*(math.pow((R0 - z1*1000)/R0, 2) + math.pow((2*z1*1000)/(k*w0*w0), 2)))/(math.pow((R0 - z1*1000)/R0, 2)*(1 - math.pow((R0 - z1*1000)/R0, 2)) - math.pow((2*z1*1000)/(k*w0*w0), 2) )
        return(wbf,Rbf)
    
    for i in range(1,noz2):
        wbf[i] = freespacebeam(w0,R0,z[i],k)[0]
        Rbf[i] = freespacebeam(w0,R0,z[i],k)[1]
    return(wbf,Rbf)
    
# free space 
wtest = [0.1,0.5,1]

wbf1 = freespace(wtest[0],R0,z,k,noz2)[0]
Rbf1 = freespace(wtest[0],R0,z,k,noz2)[1] 
#wbf2 = freespace(wtest[1],R0,z,k,noz2)[0]
#Rbf2 = freespace(wtest[1],R0,z,k,noz2)[1] 
#wbf3 = freespace(wtest[2],R0,z,k,noz2)[0]
#Rbf3 = freespace(wtest[2],R0,z,k,noz2)[1] 

wbf1[0] = w0
Rbf1[0] = R0 
#wbf2[0] = w0
#Rbf2[0] = R0 
#wbf3[0] = w0
#Rbf3[0] = R0 


def turbulenceSLC(w0,R0,z,Cn2,k,cobs):
    
    wb = np.zeros(noz2)
    Rb = np.zeros(noz2)
    def Turbulencebeam(w0,R0,z,Cn2,k,cobs):
        cob = cobs + (2*w0*w0)/(math.pow(0.55*Cn2*k*k*z*1000, -6/5))
        beta2S = 1.52*Cn2*math.pow(z*1000,3)*math.pow(w0, -1/3)
        wb = w0*math.pow((math.pow((R0 - z*1000)/R0, 2) + cob*math.pow((2*z*1000)/(k*w0*w0), 2)), 0.5) + beta2S
        phi = ((R0 - z*1000)/R0)/((2*z*1000)/(k*w0*w0)) - ((2*z*1000)/(k*w0*w0))*(w0*w0)/(math.pow(0.55*Cn2*k*k*z*1000, -6/5))
        Rb = (z*1000*(math.pow((R0 - z*1000)/R0, 2) + cob*math.pow((2*z*1000)/(k*w0*w0), 2)))/(phi*((2*z*1000)/(k*w0*w0)) - cob*math.pow((2*z*1000)/(k*w0*w0), 2) - math.pow((R0 - z*1000)/R0, 2) )
        sig12s = 1.23*Cn2*math.pow(k,7/6)*math.pow(z, 11/6)  # variables for comparison in (79)
        zrecgcomps = math.pow(z/(0.5*k*math.pow(wb,2)),5/6)*sig12s
        
        return(wb,Rb,sig12s,zrecgcomps)
        
    for i in range(1, noz2): # don't calculate at transmitter, thus start from 1 
        if z[i] <= 0.019: # Cn2 = 0 here
            wb[i] = freespacebeam(w0,R0,z[i],k)[0]
            Rb[i] = freespacebeam(w0,R0,z[i],k)[1]
        elif z[i] > 0.019 and z[i] <= 20: 
            wb[i] = Turbulencebeam(wb[offset1],Rb[offset1],z[i]-0.019,Cn2[i],k,cobs)[0]
            Rb[i] = Turbulencebeam(wb[offset1],Rb[offset1],z[i]-0.019,Cn2[i],k,cobs)[1]
            sig12s[i] = Turbulencebeam(wb[offset1],Rb[offset1],z[i]-0.019,Cn2[i],k,cobs)[2]
            zrecgcomps[i] = Turbulencebeam(wb[offset1],Rb[offset1],z[i]-0.019,Cn2[i],k,cobs)[3]
        else:
            wb[i] = freespacebeam(wb[offsetSLC],Rb[offsetSLC],z[i]-20,k)[0]
            Rb[i] = freespacebeam(wb[offsetSLC],Rb[offsetSLC],z[i]-20,k)[1]

    return(wb,Rb,sig12s,zrecgcomps)

cobsS = [1,100,500]

# Turbulence - SLC DAY model 
wb = turbulenceSLC(w0,R0,z,Cn2,k,cobsS[0])[0]
Rb = turbulenceSLC(w0,R0,z,Cn2,k,cobsS[0])[1]
#wb1 = turbulenceSLC(w0,R0,z,Cn2,k,cobsS[1])[0]
#Rb1 = turbulenceSLC(w0,R0,z,Cn2,k,cobsS[1])[1]
#wb2 = turbulenceSLC(w0,R0,z,Cn2,k,cobsS[2])[0]
#Rb2 = turbulenceSLC(w0,R0,z,Cn2,k,cobsS[2])[1]


wb[0] = w0
Rb[0] = R0
#wb1[0] = w0
#Rb1[0] = R0
#wb2[0] = w0
#Rb2[0] = R0

# Gurvich stuff


def turbulenceGUR(w0,R0,z,Cn2G,k,cobs):
    
    wbg = np.zeros(noz2)
    Rbg = np.zeros(noz2)
    sig12g = np.zeros(noz2) 
    
    def Turbulencebeam(w0,R0,z,Cn2G,k,cobs):
        cob = cobs + (2*w0*w0)/(math.pow(0.55*Cn2*k*k*z*1000, -6/5))
        beta2S = 1.52*Cn2*math.pow(z*1000,3)*math.pow(w0, -1/3)
        wb = w0*math.pow((math.pow((R0 - z*1000)/R0, 2) + cob*math.pow((2*z*1000)/(k*w0*w0), 2)), 0.5) + beta2S
        phi = ((R0 - z*1000)/R0)/((2*z*1000)/(k*w0*w0)) - ((2*z*1000)/(k*w0*w0))*(w0*w0)/(math.pow(0.55*Cn2*k*k*z*1000, -6/5))
        Rb = (z*1000*(math.pow((R0 - z*1000)/R0, 2) + cob*math.pow((2*z*1000)/(k*w0*w0), 2)))/(phi*((2*z*1000)/(k*w0*w0)) - cob*math.pow((2*z*1000)/(k*w0*w0), 2) - math.pow((R0 - z*1000)/R0, 2) )
        sig12g = 1.23*Cn2*math.pow(k,7/6)*math.pow(z, 11/6)  # variables for comparison in (79)
        zrecgcomp = math.pow(z/(0.5*k*math.pow(wb,2)),5/6)*sig12s
        
        return(wb,Rb,sig12g,zrecgcomp)
    
    for i in range(1, noz2): # don't calculate at transmitter, thus start from 1 
        if z[i] < 0.003: # Gurvich model only defined from 2.5m upwards 
            wbg[i] = freespacebeam(w0,R0,z[i],k)[0]
            Rbg[i] = freespacebeam(w0,R0,z[i],k)[1]
        elif z[i] >= 0.003 and z[i] <= upperlimGur:
            wbg[i] = Turbulencebeam(wbg[offsetGur1],Rbg[offsetGur1],z[i]-0.002,Cn2G[i],k,cobs)[0]
            Rbg[i] = Turbulencebeam(wbg[offsetGur1],Rbg[offsetGur1],z[i]-0.002,Cn2G[i],k,cobs)[1]
            sig12g[i] = Turbulencebeam(wbg[offsetGur1],Rbg[offsetGur1],z[i]-0.002,Cn2G[i],k,cobs)[2]
            zrecgcomp[i] = Turbulencebeam(wbg[offsetGur1],Rbg[offsetGur1],z[i]-0.002,Cn2G[i],k,cobs)[3]
        else:
            wbg[i] = freespacebeam(wbg[offsetGur],Rbg[offsetGur],z[i]-20,k)[0]
            Rbg[i] = freespacebeam(wbg[offsetGur],Rbg[offsetGur],z[i]-20,k)[1]


cobsG = [1,100,500]


wbg = turbulenceSLC(w0,R0,z,Cn2G,k,cobsG[2])[0]
Rbg = turbulenceSLC(w0,R0,z,Cn2G,k,cobsG[2])[1]

wbg1 = turbulenceSLC(w0,R0,z,Cn2G1,k,cobsG[2])[0]
Rbg1 = turbulenceSLC(w0,R0,z,Cn2G1,k,cobsG[2])[1]

wbg2 = turbulenceSLC(w0,R0,z,Cn2G2,k,cobsG[2])[0]
Rbg2 = turbulenceSLC(w0,R0,z,Cn2G2,k,cobsG[2])[1]

wbg3 = turbulenceSLC(w0,R0,z,Cn2G3,k,cobsG[2])[0]
Rbg3 = turbulenceSLC(w0,R0,z,Cn2G3,k,cobsG[2])[1]

# add in initial beam size and radius of curvature 

wbg[0] = w0 
Rbg[0] = R0 
wbg1[0] = w0 
Rbg1[0] = R0 
wbg2[0] = w0 
Rbg2[0] = R0 

# Gurvich model turbulence 
# =============================================================================
#for i in range(1, noz2): # don't calculate at transmitter, thus start from 1 
#    if z[i] < 0.003: # Gurvich model only defined from 2.5m upwards 
#        wbg[i] = freespacebeam(w0,R0,z[i],k)[0]
#        Rbg[i] = freespacebeam(w0,R0,z[i],k)[1]
#    elif z[i] >= 0.003 and z[i] <= upperlimGur:
#         wbg[i] = Turbulencebeam(wbg[offsetGur1],Rbg[offsetGur1],z[i]-0.002,Cn2G[i],k,cobs)[0]
#         Rbg[i] = Turbulencebeam(wbg[offsetGur1],Rbg[offsetGur1],z[i]-0.002,Cn2G[i],k,cobs)[1]
#         sig12g[i] = Turbulencebeam(wbg[offsetGur1],Rbg[offsetGur1],z[i]-0.002,Cn2G[i],k,cobs)[2]
#         zrecgcomp[i] = Turbulencebeam(wbg[offsetGur1],Rbg[offsetGur1],z[i]-0.002,Cn2G[i],k,cobs)[3]
#    else:
#         wbg[i] = freespacebeam(wbg[offsetGur],Rbg[offsetGur],z[i]-20,k)[0]
#         Rbg[i] = freespacebeam(wbg[offsetGur],Rbg[offsetGur],z[i]-20,k)[1]

# =============================================================================
    
#============================ test - compare to examples in paper!=========================================         
# =============================================================================
# lambt = 0.785*1e-6
# w0t = 0.025
# R0t = 100   ## need to determine how this parameter effects the answer!!!!
# kt = 2*math.pi/lambt
# 
# wbt = np.zeros(2000)
# Rbt = np.zeros(2000)
# cobt = np.zeros(2000)
# phit = np.zeros(2000)
# cobst = 500
# 
# wbtf = np.zeros(2000)
# Rbtf = np.zeros(2000)
# 
# 
# Cn2t = np.ones(2000)*1e-14
# 
# wbt[0] = w0t
# wbtf[0] = w0t
# 
# for i in range(1, 2000): # don't calculate at transmitter, thus start from 1 
#     if z[i] <= 0.019: # Cn2 = 0 here
#         wbt[i] = w0t*math.pow((math.pow((R0t - z[i]*1000)/R0t, 2) + math.pow((2*z[i]*1000)/(kt*w0t*w0t), 2)), 0.5)
#         Rbt[i] = (z[i]*1000*(math.pow((R0t - z[i]*1000)/R0t, 2) + math.pow((2*z[i]*1000)/(kt*w0t*w0t), 2)))/(math.pow((R0t - z[i]*1000)/R0t, 2)*(1 - math.pow((R0t - z[i]*1000)/R0t, 2)) - math.pow((2*z[i]*1000)/(kt*w0t*w0t), 2) )
#     else:
#         cobt[i] = cobst + (2*w0t*w0t)/(math.pow(0.55*Cn2t[i]*kt*kt*z[i]*1000, -6/5))
#         wbt[i] = w0t*math.pow((math.pow((R0t - z[i]*1000)/R0t, 2) + cobt[i]*math.pow((2*z[i]*1000)/(kt*w0t*w0t), 2)), 0.5)
#         phit[i] = ((R0t - z[i]*1000)/R0t)/((2*z[i]*1000)/(kt*w0t*w0t)) - ((2*z[i]*1000)/(kt*w0t*w0t))*(w0t*w0t)/(math.pow(0.55*Cn2t[i]*kt*kt*z[i]*1000, -6/5))
#         Rbt[i] = (z[i]*1000*(math.pow((R0t - z[i]*1000)/R0t, 2) + cobt[i]*math.pow((2*z[i]*1000)/(kt*w0t*w0t), 2)))/(phit[i]*((2*z[i]*1000)/(kt*w0t*w0t)) - cobt[i]*math.pow((2*z[i]*1000)/(kt*w0t*w0t), 2) - math.pow((R0t - z[i]*1000)/R0t, 2) )
# for i in range(1, 2000): # don't calculate at transmitter, thus start from 1 
#     wbtf[i] = w0t*math.pow((math.pow((R0t - z[i]*1000)/R0t, 2) + math.pow((2*z[i]*1000)/(kt*w0t*w0t), 2)), 0.5)
#     Rbtf[i] = (z[i]*1000*(math.pow((R0t - z[i]*1000)/R0t, 2) + math.pow((2*z[i]*1000)/(kt*w0t*w0t), 2)))/(math.pow((R0t - z[i]*1000)/R0t, 2)*(1 - math.pow((R0t - z[i]*1000)/R0t, 2)) - math.pow((2*z[i]*1000)/(kt*w0t*w0t), 2) )
#   
# nozt = 2000 # should be equal to (no. of km)/res + 1
# rest = 0.001 # resolution in altitude, 0.001 = 1m precision 

#======================================= optical depth calculation  =====================================
taus = np.zeros(noz)
tau = 0

for i in range(1,noz):
    
    taus[i] = ((bm2[i]+bm2[i-1])/2)*((z[i]-z[i-1])*1000) # use a Riemann sum to calculate an approximate 
    tau = tau + taus[i]                           # value of the optical thickness, tau  
    
I0 = 1000  # initial beam intensity in W/m^2

I = I0*math.exp(-tau)

# ===================================================================================================================== 

# SLC-DAY PLOT

fig = plt.figure()
ax = plt.subplot(111)
#GURlabel='Gurvich, w0 = '+ str(w0) + 'm, R0 = ' + str(R0) + 'm, ζS = ' + str(cobs)
SLClabel='SLC-Day, w0 = '+ str(w0) + 'm, R0 = ' + str(R0) + 'm, ζS = ' + str(cobsS[0])
#SLClabel1='SLC-Day, w0 = '+ str(w0) + 'm, R0 = ' + str(R0) + 'm, ζS = ' + str(cobsS[1])
#SLClabel2='SLC-Day, w0 = '+ str(w0) + 'm, R0 = ' + str(R0) + 'm, ζS = ' + str(cobsS[2])
ax.plot(z, wb, label=SLClabel)
#ax.plot(z, wb1, label=SLClabel1)
#ax.plot(z, wb2, label=SLClabel2)
plt.xlabel('Altitude / km')
plt.ylabel('Gaussian beam width / m')
plt.xlim(0,400)
plt.title('Beam width as function of altitude for SLC-Day model')
ax.legend()
namelabel = 'beamwidthvszSGf'+ str(w0) + 'R0' + str(R0) + 'cobs'+ str(cobs) + '.png'
plt.savefig(namelabel, bbox_inches='tight')
plt.show()

# GURVICH PLOT

fig = plt.figure()
ax = plt.subplot(111)
glabel='w0 = '+ str(w0) + 'm, R0 = ' + str(R0)  + 'm, \n ζS = ' + str(cobsG[2]) + ', Cn2G0 = 1.5e-13' #+ str(Cn2G0)
glabel1='w0 = '+ str(w0) + 'm, R0 = ' + str(R0)  + 'm, \n ζS = ' + str(cobsG[2]) + ', Cn2G0 = 1e-14' #+ str(Cn2G0)
glabel2='w0 = '+ str(w0) + 'm, R0 = ' + str(R0)  + 'm, \n ζS = ' + str(cobsG[2]) + ', Cn2G0 = 1e-15' #+ str(Cn2G0)
glabel3='w0 = '+ str(w0) + 'm, R0 = ' + str(R0)  + 'm, \n ζS = ' + str(cobsG[2]) + ', Cn2G0 = 1e-16' #+ str(Cn2G0)
ax.plot(z, wbg, label=glabel )
ax.plot(z, wbg1, label=glabel1 )
ax.plot(z, wbg2, label=glabel2 )
ax.plot(z, wbg3, label=glabel2 )
plt.xlabel('Altitude / km')
plt.ylabel('Gaussian beam width / m')
plt.xlim(0,400)
plt.ylim(0.5,1.5)
plt.title('Width as function of altitude for Gurvich model')
ax.legend()
namelabel = 'wbGURvsw0'+ str(w0) + 'R0' + str(R0) + 'cobs'+ str(cobs) + 'Cn2G0' + str(Cn2G0) + '.png'
plt.savefig(namelabel, bbox_inches='tight')
plt.show()

# FREE SPACE PLOT

fig = plt.figure()
ax = plt.subplot(111)
flabel1='w0 = '+ str(wtest[0]) + 'm, R0 = ' + str(R0) + 'm' #, ζS = ' + str(cobs) 
#flabel2='w0 = '+ str(wtest[1]) + 'm, R0 = ' + str(R0) + 'm' #, ζS = ' + str(cobs)
#flabel3='w0 = '+ str(wtest[2]) + 'm, R0 = ' + str(R0) + 'm' #, ζS = ' + str(cobs)
ax.plot(z, wbf1, label=flabel1 )
#ax.plot(z, wbf2, label=flabel2 )
#ax.plot(z, wbf3, label=flabel3 )
plt.xlabel('Altitude / km')
plt.ylabel('Gaussian beam width / m')
plt.xlim(0,400)
plt.title('Width as function of altitude for free space')
ax.legend()
namelabel = 'wbfvsw0'+ str(w0) + 'R0' + str(R0) + 'cobs'+ str(cobs) + '.png'
plt.savefig(namelabel, bbox_inches='tight')
plt.show()


# ============================= print out receiver beam size =============================
print("initial beam width",w0)  

print("beam width at receiver for no turbulence w0 = 0.1",wbf1[400000])

#print("beam width at receiver for no turbulence w0 = 0.5",wbf2[400000])

#print("beam width at receiver for no turbulence w0 = 1",wbf3[400000])

print("beam width at receiver for SLC-Day model",wb[400000])  

print("beam width at receiver for Gurvich model",wbg[400000])


#print("sig12g value at 20km for Gurvich model",math.pow(sig12g[20000],0.5))

print("zrecgcomp value at 20km for Gurvich model",math.pow(zrecgcomp[20000],0.5))

print("sig12s value at 20km for SLC-Day model", sig12s[20000])

print("zrecgcomps value at 20km for SLC-Day model",zrecgcomps[20000])

