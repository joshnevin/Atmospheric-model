import csv
import matplotlib.pyplot as plt
import numpy as np 
import math 



n1 = []
N = [] # number density for Rayleigh scattering calc 

with open('refractiveindex.csv') as csvDataFile:
    csvReader = csv.reader(csvDataFile)
    for row in csvReader:
        n1.append(row[1])
        #z.append(row[0])
        
        
with open('totalnumberdensity.csv') as csvDataFile:
    csvReader = csv.reader(csvDataFile)
    for row in csvReader:
        N.append(row[1])
        
del n1[0]  # remove headers 
del N[0]

n = np.ones(30000)

res = 0.001
noz = len(n)

z = np.arange(0,noz*res,res)  # instead of importing z!!

for i in range(len(n1)):
    n1[i] = float(n1[i])


for i in range(len(n1)):
    n[i] = n1[i]


for i in range(0, len(N)):
    N[i] = float(N[i])
    
       
bm2 = np.zeros(noz)
cn = np.zeros(noz)
nc = np.zeros(noz,dtype = complex)

# constants 
lamb = 630*1e-9
delta = 0.031

for i in range(0,noz):
    
    bm2[i] = ((24*math.pi*math.pi*math.pi)/(N[i]*math.pow(lamb,4)))*math.pow(((n[i]*n[i] - 1)/(n[i]*n[i] + 1)),2)*((6 + 3*delta)/(6 - 7*delta))
    cn[i] = bm2[i]*lamb/(4*math.pi)
    nc[i] = n[i] + bm2[i]*1j  # exact formula 
 
    
taus = np.zeros(noz)
tau = 0

for i in range(1,noz):
    
    taus[i] = ((bm2[i]+bm2[i-1])/2)*((z[i]-z[i-1])*1000) # use a Riemann sum to calculate an approximate 
    tau = tau + taus[i]                           # value of the optical thickness, tau  
    
I0 = 1  # initial beam intensity in W/m^2

I = I0*math.exp(-tau)


### NOTE: if n=0 then the attenuation coefficient is 0,
### thus the optical depth for z > 7.4km is just the same as calculating up to 7.4km 

# ============================================ test ==============================================


plt.plot(z,bm2,'*')
plt.xlabel('altitude')
plt.ylabel('attenuation coefficient')
plt.show()

plt.plot(z,nc,'*')
plt.xlabel('altitude')
plt.ylabel('complex refractive index')
plt.show()    
 
plt.semilogy(z,bm2,'*')
plt.xlabel('altitude')
plt.ylabel('attenuation coefficient')
plt.show()  


print('the optical depth is ' + str(tau))
print('the normalised intensity after 7.407km is ' + str(I))

