# efficient version of the analytical model of beam propagation in a turbulent medium 

import csv 
import matplotlib.pyplot as plt
import numpy as np 
import math 
import pandas as pd


n = []

with open('refractiveindex.csv') as csvDataFile:
    csvReader = csv.reader(csvDataFile)
    for row in csvReader:
        n.append(row[1])

del n[0]  # remove headers 


for i in range(0,len(n)): 
    n[i] = float(n[i])
    
    
lamb = 630*1e-9 # wavelength in m 
alt = 35786
res = 0.001
noz2 = alt/res + 1  # = (no. of km)*1000 + 1 
noz2 = int(noz2)

z = np.arange(0,noz2*res,res)  # instead of importing z!!


alt = 2000 # altitude in km 

Cn2 = np.zeros(noz2)
cob = np.zeros(noz2)
phi = np.zeros(noz2)


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

noT = 20001

Cn2G = np.zeros(noT)   # defined up to 20km 

def GurCncalc(Cn2G0):
    
    for i in range(3, noT):  # start at 3m, as the model is valid from 2.5m upwards
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
        
Cn2G0 = 1.5e-13 # in strong turbulence regime  

# Gaussian beam waist calculation -> see how width of beam varies with distance 
w0 = 0.5 # beam waist of Gaussian beam 
R0 = 1e10 # radius of curvature of the beam phase front - large = collimated beam 
k = 2*math.pi/lamb # wavenumber
cobs = 100 # spatial coherence properties of the signal-carrying laser beam as it exits the transmitter 
# (ζS = 1 for a coherent beam, ζS > 1 for a partially coherent beam)
# SLC-Day stuff 
offsetSLC = 20000
offset1 = 19
sig12s = np.zeros(offsetSLC + 1)
zrecgcomps = np.zeros(offsetSLC + 1)

zrecgcomp = np.zeros(noz2) 
offsetGur1 = 2
upperlimGur = 20  # defined upper altitude limit of turbulence model in km - after which turbulence assumed to be negligable 
offsetGur = upperlimGur*1000 


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
    
    

wbfsat = freespacebeam(w0,R0,alt,k)[0]
rbfsat = freespacebeam(w0,R0,alt,k)[1]

ng = np.ones(len(z))   # real part of refractive index for Gaussian bit 
for i in range(0,len(n)):
    ng[i] = n[i]
 
x = 0    
    
def turbulenceSLC(w0,R0,z,Cn2,k,cobs):
    
    wb = np.zeros(noT)
    Rb = np.zeros(noT)
    
    for i in range(1, noT): # don't calculate at transmitter, thus start from 1 
        if z[i]/ng[i] <= 0.019: # Cn2 = 0 here
            wb[i] = freespacebeam(w0,R0,z[i],k)[0]
            Rb[i] = freespacebeam(w0,R0,z[i],k)[1]
        
        else:
            wb[i] = Turbulencebeam(wb[offset1],Rb[offset1],z[i]-0.019,Cn2[i],k,cobs)[0]
            Rb[i] = Turbulencebeam(wb[offset1],Rb[offset1],z[i]-0.019,Cn2[i],k,cobs)[1]
            sig12s[i] = Turbulencebeam(wb[offset1],Rb[offset1],z[i]-0.019,Cn2[i],k,cobs)[2]
            zrecgcomps[i] = Turbulencebeam(wb[offset1],Rb[offset1],z[i]-0.019,Cn2[i],k,cobs)[3]
        
    return(wb,Rb,sig12s,zrecgcomps)

optz = np.multiply(z,ng)  # replace z with optical distance, Re(n)*z to account for non-freespace 


wslc = turbulenceSLC(w0,R0,optz,Cn2,k,cobs)[0][noT-1]
rslc = turbulenceSLC(w0,R0,optz,Cn2,k,cobs)[1][noT-1]


wbsat = freespacebeam(wslc, rslc, alt-20, k)[0]    
rbsat = freespacebeam(wslc, rslc, alt-20, k)[1] 
  


Cn2G = GurCncalc(Cn2G0)

# Gurvich stuff

def turbulenceGUR(w0,R0,z,Cn2G,k,cobs):
    
    wbg = np.zeros(noT)
    Rbg = np.zeros(noT)
    sig12g = np.zeros(noT) 
    
    for i in range(1, noT): # don't calculate at transmitter, thus start from 1 
        if z[i]/ng[i] < 0.003: # Gurvich model only defined from 2.5m upwards 
            wbg[i] = freespacebeam(w0,R0,z[i],k)[0]
            Rbg[i] = freespacebeam(w0,R0,z[i],k)[1]
        else:
            wbg[i] = Turbulencebeam(wbg[offsetGur1],Rbg[offsetGur1],z[i]-0.002,Cn2G[i],k,cobs)[0]
            Rbg[i] = Turbulencebeam(wbg[offsetGur1],Rbg[offsetGur1],z[i]-0.002,Cn2G[i],k,cobs)[1]
            sig12g[i] = Turbulencebeam(wbg[offsetGur1],Rbg[offsetGur1],z[i]-0.002,Cn2G[i],k,cobs)[2]
            zrecgcomp[i] = Turbulencebeam(wbg[offsetGur1],Rbg[offsetGur1],z[i]-0.002,Cn2G[i],k,cobs)[3]
    return(wbg,Rbg,sig12g)   

cobsG = 1

wgur = turbulenceGUR(w0,R0,optz,Cn2G,k,cobs)[0][noT-1]
Rgur = turbulenceGUR(w0,R0,optz,Cn2G,k,cobs)[1][noT-1]

wbgsat = freespacebeam(wgur, Rgur, alt-20, k)[0]  
rbgsat = freespacebeam(wgur, Rgur, alt-20, k)[1]  

lambtest = np.arange(350*1e-9,1550*1e-9,50*1e-9)
ktest = np.zeros(len(lambtest))
w0test = np.arange(0.1,1.1,0.1)
wbftest = np.zeros(len(lambtest))
wbStest = np.zeros(len(lambtest))
wbGtest = np.zeros(len(lambtest))

wbftestr = np.zeros(len(w0test))
wbStestr = np.zeros(len(w0test))
wbGtestr = np.zeros(len(w0test))

kfixed = (2*math.pi)/(800*1e-9)

for i in range(len(lambtest)):
    
    ktest[i] = (2*math.pi)/lambtest[i]
    wbftest[i] = freespacebeam(w0,R0,alt,ktest[i])[0]
    wbStest[i] = freespacebeam(turbulenceSLC(w0,R0,optz,Cn2,ktest[i],cobs)[0][noT-1],turbulenceSLC(w0,R0,optz,Cn2,ktest[i],cobs)[1][noT-1],alt,ktest[i])[0]
    wbGtest[i] = freespacebeam(turbulenceGUR(w0,R0,optz,Cn2G,ktest[i],cobs)[0][noT-1],turbulenceGUR(w0,R0,optz,Cn2G,ktest[i],cobs)[1][noT-1],alt,ktest[i])[0]
    
for i in range(len(w0test)):    
    wbftestr[i] = freespacebeam(w0test[i],R0,alt,kfixed)[0]
    wbStestr[i] = freespacebeam(turbulenceSLC(w0test[i],R0,optz,Cn2,kfixed,cobs)[0][noT-1],turbulenceSLC(w0test[i],R0,optz,Cn2,kfixed,cobs)[1][noT-1],alt,kfixed)[0]
    wbGtestr[i] = freespacebeam(turbulenceGUR(w0test[i],R0,optz,Cn2G,kfixed,cobs)[0][noT-1],turbulenceGUR(w0test[i],R0,optz,Cn2G,kfixed,cobs)[1][noT-1],alt,kfixed)[0]
    
    
print('the value of the beam width for no turbulence is ' + str(wbfsat))      
 
print('the value of the beam width for SLC-Day is ' + str(wbsat))  

print('the value of the beam width for Gurvich is ' + str(wbgsat))  

# diffraction limit 

rlens = 10

rspot = (0.61*alt*1000*lamb)/rlens


print('the diffraction limited spot size with a 10m optical aperture is ' + str(2*rspot) )


df = pd.DataFrame({"beam width as func of wavelength" : wbftest})
df.to_csv("beamwidthGEOwavelengthnoturb.csv", index=False)

dfr = pd.DataFrame({"beam width as func of initial beam width" : wbftestr})
dfr.to_csv("beamwidthGEOradiusnoturb.csv", index=False)

ds = pd.DataFrame({"beam width as func of wavelength" : wbStest})
ds.to_csv("beamwidthGEOwavelengthSLC.csv", index=False)
        
dsr = pd.DataFrame({"beam width as func of initial beam width " : wbStestr})
dsr.to_csv("beamwidthGEOradiusSLC.csv", index=False)

dg = pd.DataFrame({"beam width as func of wavelength" : wbGtest})
dg.to_csv("beamwidthGEOwavelengthGUR.csv", index=False)

dgr = pd.DataFrame({"beam width as func of initial beam width " : wbGtestr})
dgr.to_csv("beamwidthGEOradiusGUR.csv", index=False)


# no turb plot 









    