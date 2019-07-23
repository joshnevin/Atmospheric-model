## script for quick estimates of beam size for BeamPROP testing 
import numpy as np
import math
import matplotlib.pyplot as plt

alt = 2000 # altitude in km 
res = 0.001 # resolution in km 

noz2 = alt/res + 1 
noz2 = int(noz2)

z = np.arange(0,noz2*res,res)  # instead of importing z!!
xlim = 5
ylim = 5 
xyres = 0.1


x = np.arange(-xlim, xlim+xyres, xyres)
y = np.arange(-xlim, xlim+xyres, xyres)


lam = 800*1e-9 # wavelength in nm 
k = (2*math.pi)/lam

n = np.ones(len(z))*1.1

zopt = np.multiply(z,n)

def freespacebeam(w0,R0,z,k):
    wbf = w0*math.pow((math.pow((R0 - z*1000)/R0, 2) + math.pow((2*z*1000)/(k*w0*w0), 2)), 0.5)
    Rbf = (z*1000*(math.pow((R0 - z*1000)/R0, 2) + math.pow((2*z*1000)/(k*w0*w0), 2)))/(math.pow((R0 - z*1000)/R0, 2)*(1 - math.pow((R0 - z*1000)/R0, 2)) - math.pow((2*z*1000)/(k*w0*w0), 2) )
    return(wbf,Rbf)
    
    
def freespace(w0,R0,z,k,noz2):    
    wbf = np.zeros(noz2)
    Rbf = np.zeros(noz2)
    #def freespacebeam(w0,R0,z1,k):
    #    wbf = w0*math.pow((math.pow((R0 - z1*1000)/R0, 2) + math.pow((2*z1*1000)/(k*w0*w0), 2)), 0.5)
    #    Rbf = (z1*1000*(math.pow((R0 - z1*1000)/R0, 2) + math.pow((2*z1*1000)/(k*w0*w0), 2)))/(math.pow((R0 - z1*1000)/R0, 2)*(1 - math.pow((R0 - z1*1000)/R0, 2)) - math.pow((2*z1*1000)/(k*w0*w0), 2) )
    #    return(wbf,Rbf)
    
    for i in range(1,noz2):
        wbf[i] = freespacebeam(w0,R0,z[i],k)[0]
        Rbf[i] = freespacebeam(w0,R0,z[i],k)[1]
    return(wbf,Rbf)
    
w0 = 0.5
R0 = 1e10

wbf = freespace(w0,R0,zopt,k,noz2)[0]
Rbf = freespace(w0,R0,zopt,k,noz2)[1] 

wbf[0] = w0

# power of Gaussian beams 

zmax = len(z) 
xmax = len(x) 

# function for calculating power at all z values along path 
# =============================================================================
# def powercalc(wbf,x,y):
#     
#     Ig = np.zeros((zmax,xmax))
#     
#     for i in range(zmax):
#     
#         for j in range(xmax):
#        
#             Ig[i][j] = ((w0*w0)/(wbf[i]*wbf[i]))*math.exp((-2*(x[j]*x[j]+y[j]*y[j]))/(wbf[i]*wbf[i]))
# 
#     return(Ig)
# =============================================================================

# calculate power only at last z value 

Ig2 = np.zeros(xmax)

for i in range(xmax):

        Ig2[i] = ((w0*w0)/(wbf[zmax-1]*wbf[zmax-1]))*math.exp(-2*(x[i]*x[i]+y[i]*y[i])/(wbf[zmax-1]*wbf[zmax-1]))

print(wbf[len(wbf)-1])

plt.plot(z,wbf,'*') 
plt.show()

#plt.plot(x,Ig[len(z)-1],'*')
#plt.show()

plt.plot(x,Ig2,'*')
plt.show()