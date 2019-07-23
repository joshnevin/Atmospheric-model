# Atmospheric model, based on Optical Refractive Index of Air: Dependence on
# Pressure, Temperature and Composition - J C Owens 
# and on the US International Standard Atmosphere 1976, NASA 

import math 
import array 
import numpy 
import matplotlib.pyplot as plt
import scipy.integrate as integrate

import pandas as pd

lamdav = 6300  # vacuum wavelength of Lumi laser in Angstrom 

sig = 1/lamdav # free-space wavenumber 

P1 = 1013.25 # partial pressure of dry, CO2-free air in mb 
P2 = 13.33 # partial pressure of water vapour in mb 
P3 = 1013.25 # partial pressure of CO2 in mb 
T1 = 288.16 # average temp of dry, CO2-free air component 
T2 = 293.16 # average temp of water vapour component 
T3 = 288.16 # average temp of CO2 component 


# --------------------- dry, CO2 free air component -------------------------------

# refractivity -> should be x10^8 
r1 = 8340.78 + 2405640/(130 - math.pow(sig,2)) + 15994/(38.9 - math.pow(sig,2)) 

# parital density for spec. refraction calc
rho10 = 348.328*(P1/T1)*(1 + P1*(57.90*1e-8 - (9.4581*1e-4)/T1 + (0.25844)/math.pow(T1,2)))

n1 = 1 + r1/1e8 

R1 = ((math.pow(n1,2) - 1)/(math.pow(n1,2) + 2))*(1/rho10) # specific refraction 

# ---------------------------- water vapour ---------------------------------------
# refractivity -> should be x10^8 
r2 = 295.235 + 2.6422*math.pow(sig,2) - 0.032380*math.pow(sig,4) + 0.004028*math.pow(sig,6)

# parital density for spec. refraction calc
rho20 = 216.582 *(P2/T2)*(1 + P2*(1 + (3.7*1e-4)*P2)*(-2.37321*1e-3 + (2.23366)/(T2) - (710.792)/(math.pow(T2,2)) + (7.75141*1e4)/(math.pow(T2,3)) ))

n2 = 1 + r2/1e8 

R2 = ((math.pow(n2,2) - 1)/(math.pow(n2,2) + 2))*(1/rho20) # specific refraction 


# -------------------------------- CO2 -------------------------------------------

# refractivity -> should be x10^8
r3 = 22822.1 + 117.8*(math.pow(sig,2)) + 2406030/(130 - math.pow(sig,2)) + 15997/(38.9 - math.pow(sig,2))

# parital density for spec. refraction calc
rho30 = 529.37*(P3/T3)

n3 = 1 + r3/1e8  

R3 = ((math.pow(n3,2) - 1)/(math.pow(n3,2) + 2))*(1/rho30) # specific refraction 

# ---------------------------------------------------------------------------------

# now get the values of P and T from models to use to calculate n 
# linear temperature variation assumed for the 7 layers up to 86km 
# each section of the atmosphere is denoted by a subscript b 
#Lm = [-6.5, 0, 1, 2.8, 0, -2.8, -2]  # linear constants 
#zm = [0, 11, 20, 32, 47, 51, 71, 86]  # height of atmospheric layers 
# for z > 86km, more complex model is used, in terms of kinetic temp, Tk, rather than molecular temp, Tm

# ---------------------------------------------------------------------------------
### set to 7407 = 7.407km to calculate n within its region of validity 
noz = 86001    #7407 # should be equal to (no. of km)/res + 1
noP = noz
res = 0.001 # resolution in altitude, 0.001 = 1m precision 

T = numpy.zeros(noz)
P = numpy.zeros(noz)
#z = numpy.linspace(0,86,861,True)
z = numpy.arange(0,noz*res,res)
nN2 = numpy.zeros(noP) # array to store number densities for N2  # see page 12 of ISA paper
nO = numpy.zeros(noP) # array to store number densities for O
nO2 = numpy.zeros(noP) # array to store number densities for O2
nAr = numpy.zeros(noP) # array to store number densities for Ar
nHe = numpy.zeros(noP) # array to store number densities for He
nH = numpy.zeros(noP) # array to store number densities for H

Lap = numpy.zeros(noP) # array to store lapse rate for He 

nH = numpy.zeros(noP) 

tauH = numpy.zeros(noP)


# ---------------------------------------------------------------------------------------------

icc = 0 # main counter for the loop 

iccP = 0 # counter for the pressure calcs for z > 86Km 

iccL = 0 # counter for lapse rate calcs at high altitudes 

icct = 0 # counter for pressure calculations above 120km

# ------------------------ define constants for T and P calcs --------------------------------


g0 = 9.80665
R = 8.31432
M0 = 28.9644 
r0 = 6.356766*1e6

# number densities at 86km 
n0N2 = 1.129794*1e20 # H2
n0O = 8.6*1e16
n0O2 = 3.030898426*1e19
n0Ar = 1.35140022*1e18
n0He = 7.5817*1e14

T7 = 186.8673

Na = 6.022169*1e26 # Avagadro 

# molar masses 
MN2 = 28.0134
MO = 15.9994
MO2 = 31.9988
MAr = 39.948
MHe = 4.0026 
MH2 = 2.01594
MH = MH2/2

# a constants 
aO = 6.986*1e20
aO2 = 4.863*1e20
aAr = 4.487*1e20
aHe = 1.7*1e21
aH = 3.305*1e21

# b constants 
bO = 0.75
bO2 = 0.75
bAr = 0.87
bHe = 0.691
bH = 0.5

# flux term constants 
# Qi
QO = -5.809644*1e-4
QO2 = 1.366212*1e-4
QAr = 9.434079*1e-5
QHe = -2.457369*1e-4

# qi
qO = -3.416248*1e-3 # only applies for z from 86 to 97km, for z > 97km, qO = 0 
# qi = 0 for all other species for z > 86km 

# Ui 
UO = 56.90311
UO2 = 86
UAr = 86
UHe = 86

# ui
uO = 97
# ui = 0 for all other species 

# Wi 
WO = 2.70624*1e-5
WO2 = 8.33333*1e-5
WAr = 8.33333*1e-5
WHe = 6.666667*1e-4

# wi 
wO = 5.008765*1e-4
# wi = 0 for all other species 

# alpha constants -> equal to 0 for all species except He and H 
alHe = -0.4
alH = -0.25

exp = (g0*M0)/R

# for 86 to 95km 
K1 = 1.2*1e2
# for 95 to 115km 
#K2 = 1.2*1e2*math.exp(1- 400/(400 - math.pow((z[1] - 95),2))) 
# for 115 to 1000km 
K3 = 0

A = -76.3232
a = -19.9429
Z8 = 91
TC = 263.1905


lamb = 0.01875
Z10 = 120
Z12 = 1000 


# ----------------------------- function definitions --------------------------------- # 

def gcalc(z):
    g = (g0*r0*r0)/(r0+z)*(r0+z) # height-dependent gravitational field strength 
    return g

def Dcalculator(a,sign,b,T):
    D = (a/sign)*math.pow((T/273.15),b)
    return D

def func4n(T,D,K,Mi,M,L,g,ali):  # function to calculate f(z) for the calculation of ni for each gas species 
    fz = ((g)/(R*T))*(D/(D+K))*(Mi + (M*K)/(D) + (ali*R*L)/(g))
    return fz

def func4vibit(Qi, z, Ui, Wi, qi, ui, wi):
    vibit = Qi*math.pow((z-Ui),2)*math.exp(-Wi*math.pow((z-Ui),3)) + qi*math.pow((ui-z),2)*math.exp(-wi*math.pow((ui-z),3))
    return vibit

def geomet2geopot(z):  # function for converting geometric height, z, to geopotential height, H 
    H = (r0*z)/(r0 + z)
    return H

def Tm2Tk(Tm,M):
    T = (M/M0)*Tm
    return T

# NOTE that the temperature values calculated here are molecular temp, which should 
# be converted into kinetic temp for use in the calculation of the refractive index n 
# for z < 80km, T = Tm, but this is not true from 80 to 86km -> needs correcting     

# =============================================================================

# =============================================================================
Pb = numpy.zeros(noz)
Pb = 1013.25  # pressure at sea-level in mb 
H = numpy.arange(0,noz*res,res)
 
N = numpy.ones(noz)
  
for ih in range(0,noz):
    H[ih] = geomet2geopot(z[ih])

while icc < noz:  
     if H[icc] >= 0 and H[icc] <= 11: # b = 0 
         T[icc] = 288.15 - 6.5*H[icc]
         P[icc] = 1013.25*math.pow((288.15/(288.15 - 6.5*H[icc])),(exp/-6.5))
         N[icc] = (Na*P[icc])/(R*T[icc]*10)
         icc = icc + 1 
     
      # may need a conversion factor of 1/10 -> look at definition of R in this paper!!!  
        
     elif H[icc] >11 and H[icc] <= 20: # b = 1
          T[icc] = T[numpy.where(z==11)]
          P[icc] = P[numpy.where(z==11)]*numpy.exp((-g0*M0*(H[icc]-11))/(R*T[numpy.where(z==11)]))
          N[icc] = (Na*P[icc])/(R*T[icc]*10)
          icc = icc + 1 
          
        
     elif H[icc] > 20 and H[icc] <= 32: # b = 2
          T[icc] = T[numpy.where(z==20)] + (H[icc]-20)
          P[icc] = P[numpy.where(z==20)]*math.pow((T[numpy.where(z==20)]/(T[numpy.where(z==20)] + (H[icc]-20))),(exp/1))
          N[icc] = (Na*P[icc])/(R*T[icc]*10)
          icc = icc + 1 
      
     elif H[icc] > 32 and H[icc] <= 47: # b = 3
          T[icc] = T[numpy.where(z==32)] + 2.8*(H[icc]-32)
          P[icc] = P[numpy.where(z==32)]*math.pow((T[numpy.where(z==32)]/(T[numpy.where(z==32)] + 2.8*(H[icc]-32))),(exp/2.8))
          N[icc] = (Na*P[icc])/(R*T[icc]*10)
          icc = icc + 1 
 
          
     elif H[icc] > 47 and H[icc] <= 51: # b = 4
          T[icc] = T[numpy.where(z==47)] 
          P[icc] = P[numpy.where(z==47)]*numpy.exp((-g0*M0*(H[icc]-47))/(R*T[numpy.where(z==47)]))
          N[icc] = (Na*P[icc])/(R*T[icc]*10)
          icc = icc + 1
         
     elif H[icc] > 51 and H[icc] <= 71: # b = 5
          T[icc] = T[numpy.where(z==51)]  - 2.8*(H[icc]-51)
          P[icc] = P[numpy.where(z==51)]*math.pow((T[numpy.where(z==51)]/(T[numpy.where(z==51)] - 2.8*(H[icc]-51))),(exp/-2.8))
          N[icc] = (Na*P[icc])/(R*T[icc]*10)
          icc = icc + 1 
     
              
     elif H[icc] > 71 and H[icc] <= 79.5: # b = 6 -> need to modify this, as for z > 86km, Tm is not equal to Tk 
          T[icc] = T[numpy.where(z==71)]  - 2*(H[icc]-71)
          P[icc] = P[numpy.where(z==71)]*math.pow((T[numpy.where(z==71)]/(T[numpy.where(z==71)] - 2*(H[icc]-71))),(exp/-2))
          N[icc] = (Na*P[icc])/(R*T[icc]*10)
          icc = icc + 1 
    
     elif H[icc] > 79.5 and H[icc] <= 80: # b = 6 -> need to modify this, as for z > 86km, Tm is not equal to Tk 
          T[icc] = T[numpy.where(z==79.5)]  - 2*(H[icc]-79.5)
          T[icc] = T[icc]*0.999996
          P[icc] = P[numpy.where(z==79.5)]*math.pow((T[numpy.where(z==79.5)]/(T[numpy.where(z==79.5)] - 2*(H[icc]-79.5))),(exp/-2))
          N[icc] = (Na*P[icc])/(R*T[icc]*10)
          icc = icc + 1 
     elif H[icc] > 80 and H[icc] <= 80.5: # b = 6 -> need to modify this, as for z > 86km, Tm is not equal to Tk 
          T[icc] = T[numpy.where(z==80)]  - 2*(H[icc]-80)
          T[icc] = T[icc]*0.999988
          P[icc] = P[numpy.where(z==80)]*math.pow((T[numpy.where(z==80)]/(T[numpy.where(z==80)] - 2*(H[icc]-80))),(exp/-2))
          N[icc] = (Na*P[icc])/(R*T[icc]*10)
          icc = icc + 1
     elif H[icc] > 80.5 and H[icc] <= 81: # b = 6 -> need to modify this, as for z > 86km, Tm is not equal to Tk 
          T[icc] = T[numpy.where(z==80.5)]  - 2*(H[icc]-80.5)
          T[icc] = T[icc]*0.999969
          P[icc] = P[numpy.where(z==80.5)]*math.pow((T[numpy.where(z==80.5)]/(T[numpy.where(z==80.5)] - 2*(H[icc]-80.5))),(exp/-2))
          N[icc] = (Na*P[icc])/(R*T[icc]*0.999969*10)
          icc = icc + 1
          
     elif H[icc] > 81 and H[icc] <= 81.5: # b = 6 -> need to modify this, as for z > 86km, Tm is not equal to Tk 
          T[icc] = T[numpy.where(z==81)]  - 2*(H[icc]-81)
          T[icc] = T[icc]*0.999938
          P[icc] = P[numpy.where(z==81)]*math.pow((T[numpy.where(z==81)]/(T[numpy.where(z==81)] - 2*(H[icc]-81))),(exp/-2))
          N[icc] = (Na*P[icc])/(R*T[icc]*0.999938*10)
          icc = icc + 1
          
     elif H[icc] > 81.5 and H[icc] <= 82: # b = 6 -> need to modify this, as for z > 86km, Tm is not equal to Tk 
          T[icc] = T[numpy.where(z==81.5)]  - 2*(H[icc]-81.5)
          T[icc] = T[icc]*0.999904
          P[icc] = P[numpy.where(z==81.5)]*math.pow((T[numpy.where(z==81.5)]/(T[numpy.where(z==81.5)] - 2*(H[icc]-81.5))),(exp/-2))
          N[icc] = (Na*P[icc])/(R*T[icc]*0.999904*10)
          icc = icc + 1
     elif H[icc] > 82 and H[icc] <= 82.5: # b = 6 -> need to modify this, as for z > 86km, Tm is not equal to Tk 
          T[icc] = T[numpy.where(z==82)]  - 2*(H[icc]-82)
          T[icc] = T[icc]*0.999864
          P[icc] = P[numpy.where(z==82)]*math.pow((T[numpy.where(z==82)]/(T[numpy.where(z==82)] - 2*(H[icc]-82))),(exp/-2))
          N[icc] = (Na*P[icc])/(R*T[icc]*0.999864*10)
          icc = icc + 1
     elif H[icc] > 82.5 and H[icc] <= 83: # b = 6 -> need to modify this, as for z > 86km, Tm is not equal to Tk 
          T[icc] = T[numpy.where(z==82.5)]  - 2*(H[icc]-82.5)
          T[icc] = T[icc]*0.999822
          P[icc] = P[numpy.where(z==82.5)]*math.pow((T[numpy.where(z==82.5)]/(T[numpy.where(z==82.5)] - 2*(H[icc]-82.5))),(exp/-2))
          N[icc] = (Na*P[icc])/(R*T[icc]*0.999822*10)
          icc = icc + 1
          
          
     elif H[icc] > 83 and H[icc] <= 83.5: # b = 6 -> need to modify this, as for z > 86km, Tm is not equal to Tk 
          T[icc] = T[numpy.where(z==83)]  - 2*(H[icc]-83)
          T[icc] = T[icc]*0.999778
          P[icc] = P[numpy.where(z==83)]*math.pow((T[numpy.where(z==83)]/(T[numpy.where(z==83)] - 2*(H[icc]-83))),(exp/-2))
          N[icc] = (Na*P[icc])/(R*T[icc]*0.999778*10)
          icc = icc + 1
     
     elif H[icc] > 83.5 and H[icc] <= 84: # b = 6 -> need to modify this, as for z > 86km, Tm is not equal to Tk 
          T[icc] = T[numpy.where(z==83.5)]  - 2*(H[icc]-83.5)
          T[icc] = T[icc]*0.999731
          P[icc] = P[numpy.where(z==83.5)]*math.pow((T[numpy.where(z==83.5)]/(T[numpy.where(z==83.5)] - 2*(H[icc]-83.5))),(exp/-2))
          N[icc] = (Na*P[icc])/(R*T[icc]*0.999731*10)
          icc = icc + 1
     elif H[icc] > 84 and H[icc] <= 86: # b = 6 -> need to modify this, as for z > 86km, Tm is not equal to Tk 
          T[icc] = T[numpy.where(z==84)]  - 2*(H[icc]-84)
          T[icc] = T[icc]*0.999681
          P[icc] = P[numpy.where(z==84)]*math.pow((T[numpy.where(z==84)]/(T[numpy.where(z==84)] - 2*(H[icc]-84))),(exp/-2))
          N[icc] = (Na*P[icc])/(R*T[icc]*0.999681*10)
          icc = icc + 1 
          
          
     elif z[icc] > 86 and z[icc] <= 91:
          T[icc] = T7
          nN2[iccP] = n0N2*(T7/T[icc])*math.exp(-integrate.quad(lambda z: (g0*r0*r0*M0)/((r0+z)*(r0+z)*R*T[icc]) ,86 ,z[icc])[0])
          nO[iccP] = n0O*(T7/T[icc])*math.exp(-integrate.quad(lambda z: ((g0*r0*r0)/((r0+z)*(r0+z)*R*T[icc]))*(Dcalculator(aO,nN2[iccP],bO,T[icc])/(Dcalculator(aO,nN2[iccP],bO,T[icc]) + K1))*(MO + (MN2*K1)/Dcalculator(aO,nN2[iccP],bO,T[icc])) + QO*math.pow((z-UO),2)*math.exp(-WO*math.pow((z-UO),3)) + qO*math.pow((uO-z),2)*math.exp(-wO*math.pow((uO-z),3)) ,86 ,z[icc])[0])
          nO2[iccP] = n0O2*(T7/T[icc])*math.exp(-integrate.quad(lambda z: ((g0*r0*r0)/((r0+z)*(r0+z)*R*T[icc]))*(Dcalculator(aO2,nN2[iccP],bO2,T[icc])/(Dcalculator(aO2,nN2[iccP],bO2,T[icc]) + K1))*(MO2 + (MN2*K1)/Dcalculator(aO2,nN2[iccP],bO2,T[icc]) ) + QO2*(z-UO2)*(z-UO2)*math.exp(-WO2*math.pow((z-UO2),3)) , 86, z[icc])[0] )
          nAr[iccP] = n0Ar*(T7/T[icc])*math.exp(-integrate.quad(lambda z: ((g0*r0*r0)/((r0+z)*(r0+z)*R*T[icc]))*(Dcalculator(aAr,nN2[iccP]+nO[iccP]+nO2[iccP],bAr,T[icc])/(Dcalculator(aAr,nN2[iccP]+nO[iccP]+nO2[iccP],bAr,T[icc]) + K1))*(MAr + (((MN2*nN2[iccP] + MO*nO[iccP] + MO2*nO2[iccP])/(nN2[iccP] + nO[iccP] + nO2[iccP]))*K1)/Dcalculator(aAr,nN2[iccP]+nO[iccP]+nO2[iccP],bAr,T[icc]) ) + QAr*(z-UAr)*(z-UAr)*math.exp(-WAr*math.pow((z-UAr),3)) , 86, z[icc])[0] )
          nHe[iccP] = n0He*(T7/T[icc])*math.exp(-integrate.quad(lambda z: ((g0*r0*r0)/((r0+z)*(r0+z)*R*T[icc]))*(Dcalculator(aHe,nN2[iccP]+nO[iccP]+nO2[iccP],bHe,T[icc])/(Dcalculator(aHe,nN2[iccP]+nO[iccP]+nO2[iccP],bHe,T[icc]) + K1))*(MHe + (((MN2*nN2[iccP] + MO*nO[iccP] + MO2*nO2[iccP])/(nN2[iccP] + nO[iccP] + nO2[iccP]))*K1)/Dcalculator(aHe,nN2[iccP]+nO[iccP]+nO2[iccP],bHe,T[icc]) ) + QHe*(z-UHe)*(z-UHe)*math.exp(-WHe*math.pow((z-UHe),3)) , 86, z[icc])[0] )
          
          P[icc] = ((nN2[iccP] + nO[iccP] + nO2[iccP] + nAr[iccP] + nHe[iccP])*1.380622*1e-23*T[icc])/100 # divide by 100 to get to mb
           
          icc = icc + 1 
          iccP = iccP + 1 
      
     elif z[icc] > 91 and z[icc] <= 95:
          T[icc] = 263.1905 + A*math.pow((1-math.pow(((z[icc]-Z8)/a),2)),0.5)
          
          Lap[iccL] = (-A/a)*((z[icc]-Z8)/a)*math.pow((1-math.pow(((z[icc]-Z8)/a),2)),-0.5)
          
          nN2[iccP] = n0N2*(T7/T[icc])*math.exp(-integrate.quad(lambda z: (g0*r0*r0*M0)/((r0+z)*(r0+z)*R*T[icc]) ,86 ,z[icc])[0])
          nO[iccP] = n0O*(T7/T[icc])*math.exp(-integrate.quad(lambda z: ((g0*r0*r0)/((r0+z)*(r0+z)*R*T[icc]))*(Dcalculator(aO,nN2[iccP],bO,T[icc])/(Dcalculator(aO,nN2[iccP],bO,T[icc]) + K1))*(MO + (MN2*K1)/Dcalculator(aO,nN2[iccP],bO,T[icc])) + QO*math.pow((z-UO),2)*math.exp(-WO*math.pow((z-UO),3)) + qO*math.pow((uO-z),2)*math.exp(-wO*math.pow((uO-z),3)) ,86 ,z[icc])[0])
          nO2[iccP] = n0O2*(T7/T[icc])*math.exp(-integrate.quad(lambda z: ((g0*r0*r0)/((r0+z)*(r0+z)*R*T[icc]))*(Dcalculator(aO2,nN2[iccP],bO2,T[icc])/(Dcalculator(aO2,nN2[iccP],bO2,T[icc]) + K1))*(MO2 + (MN2*K1)/Dcalculator(aO2,nN2[iccP],bO2,T[icc]) ) + QO2*(z-UO2)*(z-UO2)*math.exp(-WO2*math.pow((z-UO2),3)) , 86, z[icc])[0] )
          nAr[iccP] = n0Ar*(T7/T[icc])*math.exp(-integrate.quad(lambda z: ((g0*r0*r0)/((r0+z)*(r0+z)*R*T[icc]))*(Dcalculator(aAr,nN2[iccP]+nO[iccP]+nO2[iccP],bAr,T[icc])/(Dcalculator(aAr,nN2[iccP]+nO[iccP]+nO2[iccP],bAr,T[icc]) + K1))*(MAr + (((MN2*nN2[iccP] + MO*nO[iccP] + MO2*nO2[iccP])/(nN2[iccP] + nO[iccP] + nO2[iccP]))*K1)/Dcalculator(aAr,nN2[iccP]+nO[iccP]+nO2[iccP],bAr,T[icc]) ) + QAr*(z-UAr)*(z-UAr)*math.exp(-WAr*math.pow((z-UAr),3)) , 86, z[icc])[0] )
          nHe[iccP] = n0He*(T7/T[icc])*math.exp(-integrate.quad(lambda z: ((g0*r0*r0)/((r0+z)*(r0+z)*R*T[icc]))*(Dcalculator(aHe,nN2[iccP]+nO[iccP]+nO2[iccP],bHe,T[icc])/(Dcalculator(aHe,nN2[iccP]+nO[iccP]+nO2[iccP],bHe,T[icc]) + K1))*(MHe + (((MN2*nN2[iccP] + MO*nO[iccP] + MO2*nO2[iccP])/(nN2[iccP] + nO[iccP] + nO2[iccP]))*K1)/Dcalculator(aHe,nN2[iccP]+nO[iccP]+nO2[iccP],bHe,T[icc]) + (alHe*R*Lap[iccL]*(r0+z)*(r0+z))/(g0*r0*r0) ) + QHe*(z-UHe)*(z-UHe)*math.exp(-WHe*math.pow((z-UHe),3)) , 86, z[icc])[0] )
          
          P[icc] = ((nN2[iccP] + nO[iccP] + nO2[iccP] + nAr[iccP] + nHe[iccP])*1.380622*1e-23*T[icc])/100 # divide by 100 to get to mb
          
          icc = icc + 1 
          iccP = iccP + 1 
          iccL = iccL + 1 
          
     elif z[icc] > 95 and z[icc] <= 97:
          T[icc] = 263.1905 + A*math.pow((1-math.pow(((z[icc]-Z8)/a),2)),0.5)
          
          Lap[iccL] = (-A/a)*((z[icc]-Z8)/a)*math.pow((1-math.pow(((z[icc]-Z8)/a),2)),-0.5)
          
          nN2[iccP] = n0N2*(T7/T[icc])*math.exp(-integrate.quad(lambda z: (g0*r0*r0*M0)/((r0+z)*(r0+z)*R*T[icc]) ,86 ,z[icc])[0])
          nO[iccP] = n0O*(T7/T[icc])*math.exp(-integrate.quad(lambda z: ((g0*r0*r0)/((r0+z)*(r0+z)*R*T[icc]))*(Dcalculator(aO,nN2[iccP],bO,T[icc])/(Dcalculator(aO,nN2[iccP],bO,T[icc]) + 1.2*1e2*math.exp(1- 400/(400 - math.pow((z - 95),2)))))*(MO + (MN2*1.2*1e2*math.exp(1- 400/(400 - math.pow((z - 95),2))))/Dcalculator(aO,nN2[iccP],bO,T[icc])) + QO*math.pow((z-UO),2)*math.exp(-WO*math.pow((z-UO),3)) + qO*math.pow((uO-z),2)*math.exp(-wO*math.pow((uO-z),3)) ,86 ,z[icc])[0])
          nO2[iccP] = n0O2*(T7/T[icc])*math.exp(-integrate.quad(lambda z: ((g0*r0*r0)/((r0+z)*(r0+z)*R*T[icc]))*(Dcalculator(aO2,nN2[iccP],bO2,T[icc])/(Dcalculator(aO2,nN2[iccP],bO2,T[icc]) + 1.2*1e2*math.exp(1- 400/(400 - math.pow((z - 95),2)))))*(MO2 + (MN2*1.2*1e2*math.exp(1- 400/(400 - math.pow((z - 95),2))))/Dcalculator(aO2,nN2[iccP],bO2,T[icc]) ) + QO2*(z-UO2)*(z-UO2)*math.exp(-WO2*math.pow((z-UO2),3)) , 86, z[icc])[0] )
          nAr[iccP] = n0Ar*(T7/T[icc])*math.exp(-integrate.quad(lambda z: ((g0*r0*r0)/((r0+z)*(r0+z)*R*T[icc]))*(Dcalculator(aAr,nN2[iccP]+nO[iccP]+nO2[iccP],bAr,T[icc])/(Dcalculator(aAr,nN2[iccP]+nO[iccP]+nO2[iccP],bAr,T[icc]) + 1.2*1e2*math.exp(1- 400/(400 - math.pow((z - 95),2)))))*(MAr + (((MN2*nN2[iccP] + MO*nO[iccP] + MO2*nO2[iccP])/(nN2[iccP] + nO[iccP] + nO2[iccP]))*1.2*1e2*math.exp(1- 400/(400 - math.pow((z - 95),2))))/Dcalculator(aAr,nN2[iccP]+nO[iccP]+nO2[iccP],bAr,T[icc]) ) + QAr*(z-UAr)*(z-UAr)*math.exp(-WAr*math.pow((z-UAr),3)) , 86, z[icc])[0] )
          nHe[iccP] = n0He*(T7/T[icc])*math.exp(-integrate.quad(lambda z: ((g0*r0*r0)/((r0+z)*(r0+z)*R*T[icc]))*(Dcalculator(aHe,nN2[iccP]+nO[iccP]+nO2[iccP],bHe,T[icc])/(Dcalculator(aHe,nN2[iccP]+nO[iccP]+nO2[iccP],bHe,T[icc]) + 1.2*1e2*math.exp(1- 400/(400 - math.pow((z - 95),2)))))*(MHe + (((MN2*nN2[iccP] + MO*nO[iccP] + MO2*nO2[iccP])/(nN2[iccP] + nO[iccP] + nO2[iccP]))*1.2*1e2*math.exp(1- 400/(400 - math.pow((z - 95),2))))/Dcalculator(aHe,nN2[iccP]+nO[iccP]+nO2[iccP],bHe,T[icc]) + (alHe*R*Lap[iccL]*(r0+z)*(r0+z))/(g0*r0*r0) ) + QHe*(z-UHe)*(z-UHe)*math.exp(-WHe*math.pow((z-UHe),3)) , 86, z[icc])[0] )
          
          P[icc] = ((nN2[iccP] + nO[iccP] + nO2[iccP] + nAr[iccP] + nHe[iccP])*1.380622*1e-23*T[icc])/100 # divide by 100 to get to mb
                  
          
          icc = icc + 1 
          iccP = iccP + 1 
          iccL = iccL + 1 
      
     elif z[icc] > 97 and z[icc] <= 110:
          T[icc] = 263.1905 + A*math.pow((1-math.pow(((z[icc]-Z8)/a),2)),0.5)
          
          Lap[iccL] = (-A/a)*((z[icc]-Z8)/a)*math.pow((1-math.pow(((z[icc]-Z8)/a),2)),-0.5)
          
          nN2[iccP] = n0N2*(T7/T[icc])*math.exp(-integrate.quad(lambda z: (g0*r0*r0*M0)/((r0+z)*(r0+z)*R*T[icc]) ,86 ,z[icc])[0])
          nO[iccP] = n0O*(T7/T[icc])*math.exp(-integrate.quad(lambda z: ((g0*r0*r0)/((r0+z)*(r0+z)*R*T[icc]))*(Dcalculator(aO,nN2[iccP],bO,T[icc])/(Dcalculator(aO,nN2[iccP],bO,T[icc]) + 1.2*1e2*math.exp(1- 400/(400 - math.pow((z - 95),2)))))*(MO + (MN2*1.2*1e2*math.exp(1- 400/(400 - math.pow((z - 95),2))))/Dcalculator(aO,nN2[iccP],bO,T[icc])) + QO*math.pow((z-UO),2)*math.exp(-WO*math.pow((z-UO),3)) ,86 ,z[icc])[0])
          nO2[iccP] = n0O2*(T7/T[icc])*math.exp(-integrate.quad(lambda z: ((g0*r0*r0)/((r0+z)*(r0+z)*R*T[icc]))*(Dcalculator(aO2,nN2[iccP],bO2,T[icc])/(Dcalculator(aO2,nN2[iccP],bO2,T[icc]) + 1.2*1e2*math.exp(1- 400/(400 - math.pow((z - 95),2)))))*(MO2 + (MN2*1.2*1e2*math.exp(1- 400/(400 - math.pow((z - 95),2))))/Dcalculator(aO2,nN2[iccP],bO2,T[icc]) ) + QO2*(z-UO2)*(z-UO2)*math.exp(-WO2*math.pow((z-UO2),3)) , 86, z[icc])[0] )
          nAr[iccP] = n0Ar*(T7/T[icc])*math.exp(-integrate.quad(lambda z: ((g0*r0*r0)/((r0+z)*(r0+z)*R*T[icc]))*(Dcalculator(aAr,nN2[iccP]+nO[iccP]+nO2[iccP],bAr,T[icc])/(Dcalculator(aAr,nN2[iccP]+nO[iccP]+nO2[iccP],bAr,T[icc]) + 1.2*1e2*math.exp(1- 400/(400 - math.pow((z - 95),2)))))*(MAr + (((MN2*nN2[iccP] + MO*nO[iccP] + MO2*nO2[iccP])/(nN2[iccP] + nO[iccP] + nO2[iccP]))*1.2*1e2*math.exp(1- 400/(400 - math.pow((z - 95),2))))/Dcalculator(aAr,nN2[iccP]+nO[iccP]+nO2[iccP],bAr,T[icc]) ) + QAr*(z-UAr)*(z-UAr)*math.exp(-WAr*math.pow((z-UAr),3)) , 86, z[icc])[0] )
          nHe[iccP] = n0He*(T7/T[icc])*math.exp(-integrate.quad(lambda z: ((g0*r0*r0)/((r0+z)*(r0+z)*R*T[icc]))*(Dcalculator(aHe,nN2[iccP]+nO[iccP]+nO2[iccP],bHe,T[icc])/(Dcalculator(aHe,nN2[iccP]+nO[iccP]+nO2[iccP],bHe,T[icc]) + 1.2*1e2*math.exp(1- 400/(400 - math.pow((z - 95),2)))))*(MHe + (((MN2*nN2[iccP] + MO*nO[iccP] + MO2*nO2[iccP])/(nN2[iccP] + nO[iccP] + nO2[iccP]))*1.2*1e2*math.exp(1- 400/(400 - math.pow((z - 95),2))))/Dcalculator(aHe,nN2[iccP]+nO[iccP]+nO2[iccP],bHe,T[icc]) + (alHe*R*Lap[iccL]*(r0+z)*(r0+z))/(g0*r0*r0) ) + QHe*(z-UHe)*(z-UHe)*math.exp(-WHe*math.pow((z-UHe),3)) , 86, z[icc])[0] )
          
          P[icc] = ((nN2[iccP] + nO[iccP] + nO2[iccP] + nAr[iccP] + nHe[iccP])*1.380622*1e-23*T[icc])/100 # divide by 100 to get to mb
                  
           
          icc = icc + 1 
          iccP = iccP + 1 
          iccL = iccL + 1 
     elif z[icc] > 110 and z[icc] <= 115:
          T[icc] = 240 + 12*(z[icc] - 110)
          
          Lap[iccL] = 12
          
          nN2[iccP] = n0N2*(T7/T[icc])*math.exp(-integrate.quad(lambda z: (g0*r0*r0*M0)/((r0+z)*(r0+z)*R*T[icc]) ,86 ,z[icc])[0])
          nO[iccP] = n0O*(T7/T[icc])*math.exp(-integrate.quad(lambda z: ((g0*r0*r0)/((r0+z)*(r0+z)*R*T[icc]))*(Dcalculator(aO,nN2[iccP],bO,T[icc])/(Dcalculator(aO,nN2[iccP],bO,T[icc]) + 1.2*1e2*math.exp(1- 400/(400 - math.pow((z - 95),2)))))*(MO + (MN2*1.2*1e2*math.exp(1- 400/(400 - math.pow((z - 95),2))))/Dcalculator(aO,nN2[iccP],bO,T[icc])) + QO*math.pow((z-UO),2)*math.exp(-WO*math.pow((z-UO),3)) ,86 ,z[icc])[0])
          nO2[iccP] = n0O2*(T7/T[icc])*math.exp(-integrate.quad(lambda z: ((g0*r0*r0)/((r0+z)*(r0+z)*R*T[icc]))*(Dcalculator(aO2,nN2[iccP],bO2,T[icc])/(Dcalculator(aO2,nN2[iccP],bO2,T[icc]) + 1.2*1e2*math.exp(1- 400/(400 - math.pow((z - 95),2)))))*(MO2 + (MN2*1.2*1e2*math.exp(1- 400/(400 - math.pow((z - 95),2))))/Dcalculator(aO2,nN2[iccP],bO2,T[icc]) ) + QO2*(z-UO2)*(z-UO2)*math.exp(-WO2*math.pow((z-UO2),3)) , 86, z[icc])[0] )
          nAr[iccP] = n0Ar*(T7/T[icc])*math.exp(-integrate.quad(lambda z: ((g0*r0*r0)/((r0+z)*(r0+z)*R*T[icc]))*(Dcalculator(aAr,nN2[iccP]+nO[iccP]+nO2[iccP],bAr,T[icc])/(Dcalculator(aAr,nN2[iccP]+nO[iccP]+nO2[iccP],bAr,T[icc]) + 1.2*1e2*math.exp(1- 400/(400 - math.pow((z - 95),2)))))*(MAr + (((MN2*nN2[iccP] + MO*nO[iccP] + MO2*nO2[iccP])/(nN2[iccP] + nO[iccP] + nO2[iccP]))*1.2*1e2*math.exp(1- 400/(400 - math.pow((z - 95),2))))/Dcalculator(aAr,nN2[iccP]+nO[iccP]+nO2[iccP],bAr,T[icc]) ) + QAr*(z-UAr)*(z-UAr)*math.exp(-WAr*math.pow((z-UAr),3)) , 86, z[icc])[0] )
          nHe[iccP] = n0He*(T7/T[icc])*math.exp(-integrate.quad(lambda z: ((g0*r0*r0)/((r0+z)*(r0+z)*R*T[icc]))*(Dcalculator(aHe,nN2[iccP]+nO[iccP]+nO2[iccP],bHe,T[icc])/(Dcalculator(aHe,nN2[iccP]+nO[iccP]+nO2[iccP],bHe,T[icc]) + 1.2*1e2*math.exp(1- 400/(400 - math.pow((z - 95),2)))))*(MHe + (((MN2*nN2[iccP] + MO*nO[iccP] + MO2*nO2[iccP])/(nN2[iccP] + nO[iccP] + nO2[iccP]))*1.2*1e2*math.exp(1- 400/(400 - math.pow((z - 95),2))))/Dcalculator(aHe,nN2[iccP]+nO[iccP]+nO2[iccP],bHe,T[icc]) + (alHe*R*Lap[iccL]*(r0+z)*(r0+z))/(g0*r0*r0) ) + QHe*(z-UHe)*(z-UHe)*math.exp(-WHe*math.pow((z-UHe),3)) , 86, z[icc])[0] )
          
          P[icc] = ((nN2[iccP] + nO[iccP] + nO2[iccP] + nAr[iccP] + nHe[iccP])*1.380622*1e-23*T[icc])/100 # divide by 100 to get to mb
                  
          
          icc = icc + 1 
          iccP = iccP + 1 
          iccL = iccL + 1   
          
     elif z[icc] > 115 and z[icc] <= 120:
          T[icc] = 240 + 12*(z[icc] - 110)
          
          Lap[iccL] = 12
          
          nN2[iccP] = n0N2*(T7/T[icc])*math.exp(-integrate.quad(lambda z: (g0*r0*r0*M0)/((r0+z)*(r0+z)*R*T[icc]) ,86 ,z[icc])[0])
          nO[iccP] = n0O*(T7/T[icc])*math.exp(-integrate.quad(lambda z: ((g0*r0*r0)/((r0+z)*(r0+z)*R*T[icc]))*(Dcalculator(aO,nN2[iccP],bO,T[icc])/(Dcalculator(aO,nN2[iccP],bO,T[icc]) ))*MO + QO*math.pow((z-UO),2)*math.exp(-WO*math.pow((z-UO),3)) ,86 ,z[icc])[0])
          nO2[iccP] = n0O2*(T7/T[icc])*math.exp(-integrate.quad(lambda z: ((g0*r0*r0)/((r0+z)*(r0+z)*R*T[icc]))*(Dcalculator(aO2,nN2[iccP],bO2,T[icc])/(Dcalculator(aO2,nN2[iccP],bO2,T[icc])))*MO2  + QO2*(z-UO2)*(z-UO2)*math.exp(-WO2*math.pow((z-UO2),3)) , 86, z[icc])[0] )
          nAr[iccP] = n0Ar*(T7/T[icc])*math.exp(-integrate.quad(lambda z: ((g0*r0*r0)/((r0+z)*(r0+z)*R*T[icc]))*(Dcalculator(aAr,nN2[iccP]+nO[iccP]+nO2[iccP],bAr,T[icc])/(Dcalculator(aAr,nN2[iccP]+nO[iccP]+nO2[iccP],bAr,T[icc]) ))*MAr + QAr*(z-UAr)*(z-UAr)*math.exp(-WAr*math.pow((z-UAr),3)) , 86, z[icc])[0] )
          nHe[iccP] = n0He*(T7/T[icc])*math.exp(-integrate.quad(lambda z: ((g0*r0*r0)/((r0+z)*(r0+z)*R*T[icc]))*(Dcalculator(aHe,nN2[iccP]+nO[iccP]+nO2[iccP],bHe,T[icc])/(Dcalculator(aHe,nN2[iccP]+nO[iccP]+nO2[iccP],bHe,T[icc]) ))*(MHe + (alHe*R*Lap[iccL]*(r0+z)*(r0+z))/(g0*r0*r0) ) + QHe*(z-UHe)*(z-UHe)*math.exp(-WHe*math.pow((z-UHe),3)) , 86, z[icc])[0] )
          
          P[icc] = ((nN2[iccP] + nO[iccP] + nO2[iccP] + nAr[iccP] + nHe[iccP])*1.380622*1e-23*T[icc])/100 # divide by 100 to get to mb
                  
           #K2 = 1.2*1e2*math.exp(1- 400/(400 - math.pow((z[1] - 95),2)))  
          icc = icc + 1 
          iccP = iccP + 1 
          iccL = iccL + 1 
          
     elif z[icc] > 120 and z[icc] <= 150: # z>120km but no H conc
          T[icc] = 1000 - (1000 - 360)*math.exp(-0.01875*((z[icc] - 120)*(r0 + 120)/(r0+z[icc])))
          
          Lap[iccL] = 0.01875*(1000 - 360)*math.pow(((r0 + 120)/(r0 + z[icc])), 2)*math.exp(-0.01875*((z[icc] - 120)*(r0 + 120)/(r0+z[icc])))
          
          nN2[iccP] = n0N2*(T7/T[icc])*math.exp(-integrate.quad(lambda z: (g0*r0*r0*M0)/((r0+z)*(r0+z)*R*T[icc]) ,86 ,z[icc])[0])
          nO[iccP] = n0O*(T7/T[icc])*math.exp(-integrate.quad(lambda z: ((g0*r0*r0)/((r0+z)*(r0+z)*R*T[icc]))*(Dcalculator(aO,nN2[iccP],bO,T[icc])/(Dcalculator(aO,nN2[iccP],bO,T[icc]) ))*MO + QO*math.pow((z-UO),2)*math.exp(-WO*math.pow((z-UO),3)) ,86 ,z[icc])[0])
          nO2[iccP] = n0O2*(T7/T[icc])*math.exp(-integrate.quad(lambda z: ((g0*r0*r0)/((r0+z)*(r0+z)*R*T[icc]))*(Dcalculator(aO2,nN2[iccP],bO2,T[icc])/(Dcalculator(aO2,nN2[iccP],bO2,T[icc])))*MO2  + QO2*(z-UO2)*(z-UO2)*math.exp(-WO2*math.pow((z-UO2),3)) , 86, z[icc])[0] )
          nAr[iccP] = n0Ar*(T7/T[icc])*math.exp(-integrate.quad(lambda z: ((g0*r0*r0)/((r0+z)*(r0+z)*R*T[icc]))*(Dcalculator(aAr,nN2[iccP]+nO[iccP]+nO2[iccP],bAr,T[icc])/(Dcalculator(aAr,nN2[iccP]+nO[iccP]+nO2[iccP],bAr,T[icc]) ))*MAr + QAr*(z-UAr)*(z-UAr)*math.exp(-WAr*math.pow((z-UAr),3)) , 86, z[icc])[0] )
          nHe[iccP] = n0He*(T7/T[icc])*math.exp(-integrate.quad(lambda z: ((g0*r0*r0)/((r0+z)*(r0+z)*R*T[icc]))*(Dcalculator(aHe,nN2[iccP]+nO[iccP]+nO2[iccP],bHe,T[icc])/(Dcalculator(aHe,nN2[iccP]+nO[iccP]+nO2[iccP],bHe,T[icc]) ))*(MHe + (alHe*R*Lap[iccL]*(r0+z)*(r0+z))/(g0*r0*r0) ) + QHe*(z-UHe)*(z-UHe)*math.exp(-WHe*math.pow((z-UHe),3)) , 86, z[icc])[0] )
          
          P[icc] = ((nN2[iccP] + nO[iccP] + nO2[iccP] + nAr[iccP] + nHe[iccP])*1.380622*1e-23*T[icc])/100 # divide by 100 to get to mb
                   
          icc = icc + 1 
          iccP = iccP + 1 
          iccL = iccL + 1 
         
     elif z[icc] > 150 and z[icc] <= 500: 
          T[icc] = 1000 - (1000 - 360)*math.exp(-0.01875*((z[icc] - 120)*(r0 + 120)/(r0+z[icc])))
          
          Lap[iccL] = 0.01875*(1000 - 360)*math.pow(((r0 + 120)/(r0 + z[icc])), 2)*math.exp(-0.01875*((z[icc] - 120)*(r0 + 120)/(r0+z[icc])))
          
          nN2[iccP] = n0N2*(T7/T[icc])*math.exp(-integrate.quad(lambda z: (g0*r0*r0*M0)/((r0+z)*(r0+z)*R*T[icc]) ,86 ,z[icc])[0])
          nO[iccP] = n0O*(T7/T[icc])*math.exp(-integrate.quad(lambda z: ((g0*r0*r0)/((r0+z)*(r0+z)*R*T[icc]))*(Dcalculator(aO,nN2[iccP],bO,T[icc])/(Dcalculator(aO,nN2[iccP],bO,T[icc]) ))*MO + QO*math.pow((z-UO),2)*math.exp(-WO*math.pow((z-UO),3)) ,86 ,z[icc])[0])
          nO2[iccP] = n0O2*(T7/T[icc])*math.exp(-integrate.quad(lambda z: ((g0*r0*r0)/((r0+z)*(r0+z)*R*T[icc]))*(Dcalculator(aO2,nN2[iccP],bO2,T[icc])/(Dcalculator(aO2,nN2[iccP],bO2,T[icc])))*MO2  + QO2*(z-UO2)*(z-UO2)*math.exp(-WO2*math.pow((z-UO2),3)) , 86, z[icc])[0] )
          nAr[iccP] = n0Ar*(T7/T[icc])*math.exp(-integrate.quad(lambda z: ((g0*r0*r0)/((r0+z)*(r0+z)*R*T[icc]))*(Dcalculator(aAr,nN2[iccP]+nO[iccP]+nO2[iccP],bAr,T[icc])/(Dcalculator(aAr,nN2[iccP]+nO[iccP]+nO2[iccP],bAr,T[icc]) ))*MAr + QAr*(z-UAr)*(z-UAr)*math.exp(-WAr*math.pow((z-UAr),3)) , 86, z[icc])[0] )
          nHe[iccP] = n0He*(T7/T[icc])*math.exp(-integrate.quad(lambda z: ((g0*r0*r0)/((r0+z)*(r0+z)*R*T[icc]))*(Dcalculator(aHe,nN2[iccP]+nO[iccP]+nO2[iccP],bHe,T[icc])/(Dcalculator(aHe,nN2[iccP]+nO[iccP]+nO2[iccP],bHe,T[icc]) ))*(MHe + (alHe*R*Lap[iccL]*(r0+z)*(r0+z))/(g0*r0*r0) ) + QHe*(z-UHe)*(z-UHe)*math.exp(-WHe*math.pow((z-UHe),3)) , 86, z[icc])[0] )
          
          tauH[icct] = integrate.quad(lambda z: (g0*r0*r0*MH)/((r0+z)*(r0+z)*R*T[icc]), 500, z[icc])[0]
          
          nH[iccP] = 8*1e10 - integrate.quad(lambda z: (7.2*1e11/Dcalculator(aH,nN2[iccP]+nO[iccP]+nO2[iccP]+nAr[icc]+nHe[icc],bH,T[icc]))*math.pow((T[icc]/999.2356),(1 + alH))*math.exp(tauH[icct]) , 500, z[icc])[0]*math.pow((999.2356/T[icc]),(1+alH))*math.exp(-tauH[icct])
          
          #nH[iccP] = 8*1e10 - integrate.quad( (7.2*1e11/Dcalculator(aH,nN2[iccP]+nO[iccP]+nO2[iccP]+nAr[icc]+nHe[icc],bH,T[icc])) ,500,z[icc])[0]
          
          P[icc] = ((nN2[iccP] + nO[iccP] + nO2[iccP] + nAr[iccP] + nHe[iccP] + nH[iccP])*1.380622*1e-23*T[icc])/100 # divide by 100 to get to mb
                  
           #K2 = 1.2*1e2*math.exp(1- 400/(400 - math.pow((z[1] - 95),2)))  
          icc = icc + 1 
          iccP = iccP + 1 
          iccL = iccL + 1      
          icct = icct + 1 
     
     elif z[icc] > 500 and z[icc] <= 1000: # H in diffusive equillibruim above 500km 
          T[icc] = 1000 - (1000 - 360)*math.exp(-0.01875*((z[icc] - 120)*(r0 + 120)/(r0+z[icc])))
          
          Lap[iccL] = 0.01875*(1000 - 360)*math.pow(((r0 + 120)/(r0 + z[icc])), 2)*math.exp(-0.01875*((z[icc] - 120)*(r0 + 120)/(r0+z[icc])))
          
          nN2[iccP] = n0N2*(T7/T[icc])*math.exp(-integrate.quad(lambda z: (g0*r0*r0*M0)/((r0+z)*(r0+z)*R*T[icc]) ,86 ,z[icc])[0])
          nO[iccP] = n0O*(T7/T[icc])*math.exp(-integrate.quad(lambda z: ((g0*r0*r0)/((r0+z)*(r0+z)*R*T[icc]))*(Dcalculator(aO,nN2[iccP],bO,T[icc])/(Dcalculator(aO,nN2[iccP],bO,T[icc]) ))*MO + QO*math.pow((z-UO),2)*math.exp(-WO*math.pow((z-UO),3)) ,86 ,z[icc])[0])
          nO2[iccP] = n0O2*(T7/T[icc])*math.exp(-integrate.quad(lambda z: ((g0*r0*r0)/((r0+z)*(r0+z)*R*T[icc]))*(Dcalculator(aO2,nN2[iccP],bO2,T[icc])/(Dcalculator(aO2,nN2[iccP],bO2,T[icc])))*MO2  + QO2*(z-UO2)*(z-UO2)*math.exp(-WO2*math.pow((z-UO2),3)) , 86, z[icc])[0] )
          nAr[iccP] = n0Ar*(T7/T[icc])*math.exp(-integrate.quad(lambda z: ((g0*r0*r0)/((r0+z)*(r0+z)*R*T[icc]))*(Dcalculator(aAr,nN2[iccP]+nO[iccP]+nO2[iccP],bAr,T[icc])/(Dcalculator(aAr,nN2[iccP]+nO[iccP]+nO2[iccP],bAr,T[icc]) ))*MAr + QAr*(z-UAr)*(z-UAr)*math.exp(-WAr*math.pow((z-UAr),3)) , 86, z[icc])[0] )
          nHe[iccP] = n0He*(T7/T[icc])*math.exp(-integrate.quad(lambda z: ((g0*r0*r0)/((r0+z)*(r0+z)*R*T[icc]))*(Dcalculator(aHe,nN2[iccP]+nO[iccP]+nO2[iccP],bHe,T[icc])/(Dcalculator(aHe,nN2[iccP]+nO[iccP]+nO2[iccP],bHe,T[icc]) ))*(MHe + (alHe*R*Lap[iccL]*(r0+z)*(r0+z))/(g0*r0*r0) ) + QHe*(z-UHe)*(z-UHe)*math.exp(-WHe*math.pow((z-UHe),3)) , 86, z[icc])[0] )
          nH[iccP] = 8*1e10*math.pow( 999.2356/T[icc] ,(1+alH))
          
          P[icc] = ((nN2[iccP] + nO[iccP] + nO2[iccP] + nAr[iccP] + nHe[iccP] + nH[iccP])*1.380622*1e-23*T[icc])/100 # divide by 100 to get to mb
                  
          icc = icc + 1 
          iccP = iccP + 1 
          iccL = iccL + 1      
          icct = icct + 1 
          
          
# number density of H is negligable compared to other species until 150Km 

# refractive index calc 
# general values of density -> depend on P and T 
# =============================================================================

# check refractive index calc -> doesn't agree with refractiveindex.info 
# currently a factor of 3 too large at sea level - look at JC Owens paper again           

ic = 0


 
n = numpy.zeros(len(P)) # pre-allocate an array to overwrite 
# =============================================================================
# for i, j in zip(range(7407), range(7407)):  # simulataneously loop over i,j 
#     
#     if i > 100:
#     
#         rho1 = 348.328*(i/j)*(1 + i*(57.90*1e-8 - (9.4581*1e-4)/j + (0.25844)/math.pow(j,2)))
#         #rho2 = 216.582 *(i/j)*(1 + i*(1 + (3.7*1e-4)*i)*(-2.37321*1e-3 + (2.23366)/(j) - (710.792/(math.pow(j,2))) + (7.75141*1e4)/(math.pow(j,3)) ))
#         #rho3 = 529.37*(i/j)
#  
#         A = R1*rho1  # equation (4) 
#  
#         B = (2*A + 1)/(1 - A) 
#  
#         n[ic] = math.sqrt(B)  
#  
#         ic = ic + 1
#  
#     elif i > 17:
#         rho1 = 348.328*(i/j)*(1 + i*(57.90*1e-8 - (9.4581*1e-4)/j + (0.25844)/math.pow(j,2)))
#         rho2 = 216.582 *(i/j)*(1 + i*(1 + (3.7*1e-4)*i)*(-2.37321*1e-3 + (2.23366)/(j) - (710.792/(math.pow(j,2))) + (7.75141*1e4)/(math.pow(j,3)) ))
#     
#         A = R1*rho1 + R2*rho2  # equation (4) 
#  
#         B = (2*A + 1)/(1 - A) 
#  
#         n[ic] = math.sqrt(B)  
#  
#         ic = ic + 1
#     
#     else:
#         rho1 = 348.328*(i/j)*(1 + i*(57.90*1e-8 - (9.4581*1e-4)/j + (0.25844)/math.pow(j,2)))
#         rho2 = 216.582 *(i/j)*(1 + i*(1 + (3.7*1e-4)*i)*(-2.37321*1e-3 + (2.23366)/(j) - (710.792/(math.pow(j,2))) + (7.75141*1e4)/(math.pow(j,3)) ))
#         rho3 = 529.37*(i/j)
#  
#         A = R1*rho1 + R2*rho2 + R3*rho3 # equation (4) 
#  
#         B = (2*A + 1)/(1 - A) 
#  
#         n[ic] = math.sqrt(B)  
#  
#         ic = ic + 1
# =============================================================================
 
for i, j in zip(P, T):
    rho1 = 348.328*(i/j)*(1 + i*(57.90*1e-8 - (9.4581*1e-4)/j + (0.25844)/math.pow(j,2)))
    A = R1*rho1  # equation (4)   
    B = (2*A + 1)/(1 - A) 
  
    n[ic] = math.sqrt(B)  
  
    ic = ic + 1
   
 
nminusone = numpy.zeros(noz)
for il in range(0,noz - 1):
     nminusone[il] = n[il] - 1     
    
 
plt.semilogy(z,nminusone,'*')   
plt.title('log(refractive index - 1) vs height')
plt.xlim(0, 7.407) 
plt.show()
# =============================================================================

plt.plot(z,n,'*')
plt.title('refractive index vs height')
plt.xlabel('Altitude / km')
plt.ylabel('Refractive index')
plt.xlim(0, 7.407) 
#plt.savefig('nvsz.png', bbox_inches='tight')
plt.show()


fig, ax = plt.subplots()
plt.plot(T,z,'*')
ax.set_aspect(1.5)
plt.title('temperature vs height')
plt.xlabel('Temperature / K')
plt.ylabel('Altitude / km')
plt.savefig('Tvsz.png', bbox_inches='tight')
plt.show()

axes = plt.gca()
#axes.set_xlim([100,1000])
plt.semilogy(z,P,'*')
#axes.set_ylim([0,86])
plt.title('log pressure vs height')
plt.xlabel('Altitude / km')
plt.ylabel('Pressure / mb ')
plt.savefig('logPvsz.png', bbox_inches='tight')
plt.show()

axes = plt.gca()
#axes.set_xlim([100,1000])
plt.plot(z,P,'*')
#axes.set_ylim([0,86])
plt.title('Linear pressure vs height')
plt.xlabel('Altitude / km')
plt.ylabel('Pressure / mb ')
plt.savefig('linearPvsz.png', bbox_inches='tight')
plt.show()




# write refractive index model to a CSV for processing 
# =============================================================================
#df = pd.DataFrame({"altitude" : z, "refractive index" : n})
#df.to_csv("refractiveindex.csv", index=False)

#dfp = pd.DataFrame({"altitude" : z, "pressure" : P})
#dfp.to_csv("pressure.csv", index=False)

#dfT = pd.DataFrame({"altitude" : z, "temperature" : T})
#dfT.to_csv("temperature.csv", index=False)

# =============================================================================
# dfnN2 = pd.DataFrame({"altitude" : z, "no density of N2" : nN2})
# dfnN2.to_csv("temperature.csv", index=False)
# 
# dfnO = pd.DataFrame({"altitude" : z, "no density of O" : nO})
# dfnO.to_csv("temperature.csv", index=False)
# 
# dfnO2 = pd.DataFrame({"altitude" : z, "no density of O2" : nO2})
# dfnO2.to_csv("temperature.csv", index=False)
# 
# dfnHe = pd.DataFrame({"altitude" : z, "no density of He" : nHe})
# dfnHe.to_csv("temperature.csv", index=False)
# 
# dfnH = pd.DataFrame({"altitude" : z, "no density of H" : nH})
# dfnH.to_csv("temperature.csv", index=False)
# 
# dfnAr = pd.DataFrame({"altitude" : z, "no density of Ar" : nAr})
# dfnAr.to_csv("temperature.csv", index=False)
# =============================================================================

dfN = pd.DataFrame({"altitude" : z, "total no density for z below 86km" : N})
dfN.to_csv("totalnumberdensity.csv", index=False)

# =============================================================================
