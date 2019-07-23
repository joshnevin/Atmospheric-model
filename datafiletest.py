### test script to see if I understand the RSoft standard data file format 

import numpy as np


# 1D data file for f(x,y,z) = x, for x = [0,1] and grid of 0.1

x = [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]

Nxs = len(x)

X0s = 0

Xns = 1

Zposs = 0

with open('1dsimpleind.txt', 'w') as f:
    f.write('/rn,a,b/nx0\n')
    f.write(str(Nxs) + ' ' + str(X0s) + ' ' + str(Xns) + ' ' + str(Zposs) +  ' ' + 'OUTPUT_REAL' + '\n')
    for i in range(0, len(x)):
        f.write("%.2f" % x[i] + '\n')

# 2D data file for f(x,y,z) = x + y 

x = [-1,-0.5,0,0.5,1]

y = [-1,-0.5,0,0.5,1]

z = [0,1,2,3,4,5]

Nx2 = 5
Ny2 = 5
Nz2 = 6

X02 = -1
Y02 = -1
Z02 = 0
Xn2 = 1
Yn2 = 1
Zn2 = 1
Zpos2 = 0 

ic = 0 

iz = 0

xm = np.zeros(len(x))

for i in range(len(x)):  # for n0 = 1: need to subtract n0 from all real values 
    xm[i] = x[i] - 1  # try abs( abs(x[i]) - 1 ) 
    

with open('3dsimpleind.txt', 'w') as f:
    f.write('/rn,a,b/nx0\n')
    f.write('/r,qa,qb\n')
    f.write('/r\n')
    f.write(str(Nx2) + ' ' + str(X02) + ' ' + str(Xn2) + ' Z_DEPENDENT ' + 'OUTPUT_REAL_IMAG_3D' +  '\n' )
    f.write(str(Ny2) + ' ' + str(Y02) + ' ' + str(Yn2) + '\n')
    f.write(str(Nz2) + ' ' + str(Z02) + ' ' + str(Zn2) + '\n')
    while iz < len(z):
        
        
        for i in range(Nx2):
            while ic < Ny2:
                if ic < Ny2 - 1:
                    f.write("%.2f" %  xm[i] + '    ' + "%.2f" % x[i] + '    '  )
                else:
                    f.write("%.2f" %  xm[i] + '    ' + "%.2f" % x[i] + '\n')
                ic = ic + 1
            ic = 0
        iz = iz + 1 


izt = 0
izt2 = 0 

nzrealt = [1.0,2.0,3.0,4.0,5.0,6.0]
nzimagt = [1.0,2.0,3.0,4.0,5.0,6.0]

for i in range(len(nzrealt)):
    nzrealt[i] = nzrealt[i] - 1

with open('3dzdepindtestv2.txt', 'w') as f:
    f.write('/rn,a,b/nx0\n')
    f.write('/r,qa,qb\n')
    f.write('/r\n')
    f.write(str(Nx2) + ' ' + str(X02) + ' ' + str(Xn2) + ' Z_DEPENDENT ' + 'OUTPUT_REAL_IMAG_3D' +  '\n' )
    f.write(str(Ny2) + ' ' + str(Y02) + ' ' + str(Yn2) + '\n')
    f.write(str(Nz2) + ' ' + str(Z02) + ' ' + str(Zn2) + '\n')
    while izt2 < len(z):
        
        
        for i in range(Nx2):
            while ic < Ny2:
                if ic < Ny2 - 1:
                    f.write("%.2f" %  nzrealt[izt2] + '    ' + "%.2f" % nzimagt[izt2] + '    '  )
                else:
                    f.write("%.2f" %  nzrealt[izt2] + '    ' + "%.2f" % nzimagt[izt2] + '\n')
                ic = ic + 1
            ic = 0
        izt2 = izt2 + 1 


with open('3dzdepindtest.txt', 'w') as f:
    f.write('/rn,a,b/nx0\n')
    f.write('/r,qa,qb\n')
    f.write('/r\n')
    f.write(str(Nx2) + ' ' + str(X02) + ' ' + str(Xn2) + ' Z_DEPENDENT ' + 'OUTPUT_REAL_IMAG_3D' +  '\n' )
    f.write(str(Ny2) + ' ' + str(Y02) + ' ' + str(Yn2) + '\n')
    f.write(str(Nz2) + ' ' + str(Z02) + ' ' + str(Zn2) + '\n')
    while izt < len(z):
        
        
        for i in range(Nx2):
            while ic < Ny2:
                if ic < Ny2 - 1:
                    f.write("%.2f" %  (izt-1)   + '    ' + "%.2f" % izt + '    '  )
                else:
                    f.write("%.2f" %  (izt-1)  + '    ' + "%.2f" % izt + '\n')
                ic = ic + 1
            ic = 0
        izt = izt + 1 


izt3 = 0

nzrealt3 = [6.0,5.0,4.0,3.0,2.0,1.0]
nzimagt3 = [6.0,5.0,4.0,3.0,2.0,1.0]


for i in range(len(nzrealt3)):
    nzrealt3[i] = nzrealt3[i] - 1

with open('3dzdepindtestv3.txt', 'w') as f:
    f.write('/rn,a,b/nx0\n')
    f.write('/r,qa,qb\n')
    f.write('/r\n')
    f.write(str(Nx2) + ' ' + str(X02) + ' ' + str(Xn2) + ' Z_DEPENDENT ' + 'OUTPUT_REAL_IMAG_3D' +  '\n' )
    f.write(str(Ny2) + ' ' + str(Y02) + ' ' + str(Yn2) + '\n')
    f.write(str(Nz2) + ' ' + str(Z02) + ' ' + str(Zn2) + '\n')
    while izt3 < len(z):
        
        
        for i in range(Nx2):
            while ic < Ny2:
                if ic < Ny2 - 1:
                    f.write("%.2f" %  nzrealt3[izt3] + '    ' + "%.2f" % nzimagt3[izt3] + '    '  )
                else:
                    f.write("%.2f" %  nzrealt3[izt3] + '    ' + "%.2f" % nzimagt3[izt3] + '\n')
                ic = ic + 1
            ic = 0
        izt3 = izt3 + 1 




# try and recreate 3D data file example 
        
def f3d(x,y,z): 
    f3d = x + y + z
    return f3d 

test3dy0z1 = np.zeros(3)
test3dy0z2 = np.zeros(3)
test3dy0z3 = np.zeros(3)


test3dy1z1 = np.zeros(3)
test3dy1z2 = np.zeros(3)
test3dy1z3 = np.zeros(3)

test3dy2z1 = np.zeros(3)
test3dy2z2 = np.zeros(3)
test3dy2z3 = np.zeros(3)

x = [-1,0,1]

y = [-1,0,1]

z = [0,0.5,1]


for i in range(0,3): # find first column 

    test3dy0z1[i] = f3d(x[i],y[0],z[0])  
    test3dy0z2[i] = f3d(x[i],y[0],z[1])
    test3dy0z3[i] = f3d(x[i],y[0],z[2])

for i in range(0,3): # find second column 

    test3dy1z1[i] = f3d(x[i],y[1],z[0])  
    test3dy1z2[i] = f3d(x[i],y[1],z[1])
    test3dy1z3[i] = f3d(x[i],y[1],z[2])
    
for i in range(0,3): # find thrid column 

    test3dy2z1[i] = f3d(x[i],y[2],z[0])  
    test3dy2z2[i] = f3d(x[i],y[2],z[1])
    test3dy2z3[i] = f3d(x[i],y[2],z[2])
    

# real component = x + y + z
# imag component = x

with open('testindexfile3D.txt', 'w') as f:
    for i in range(0,3):
        f.write("%.3f" % test3dy0z1[i] + '    ' +  "%.3f" % x[i] + '    ' +  "%.3f" %  test3dy1z1[i] + '    '  +  "%.3f" % x[i] + '    ' +  "%.3f" %  test3dy2z1[i] +  '    ' + "%.3f" % x[i]  + '\n')
    for i in range(0,3):
        f.write("%.3f" % test3dy0z2[i] + '    ' +  "%.3f" % x[i] + '    ' +  "%.3f" %  test3dy1z2[i] + '    '  +  "%.3f" % x[i] + '    ' +  "%.3f" %  test3dy2z2[i] +  '    ' + "%.3f" % x[i] + '\n')
    for i in range(0,3):
        f.write("%.3f" % test3dy0z3[i] + '    ' +  "%.3f" % x[i] + '    ' +  "%.3f" %  test3dy1z3[i] + '    '  +  "%.3f" % x[i] + '    ' +  "%.3f" %  test3dy2z3[i] +  '    ' + "%.3f" % x[i] + '\n')
    #f.write("%.3f" % test3dy0z2[i] + '    ' +  "%.3f" %  test3dy1z2[i] + '    '  +  "%.3f" %  test3dy2z2[i] + '\n')
        #f.write("%.3f" % test3dy0z3[i] + '    ' +  "%.3f" %  test3dy1z3[i] + '    '  +  "%.3f" %  test3dy2z3[i] + '\n')

## example 6 - imaginary values 
        
real = [1,1.25,1.5,1.75,2]
imag = [0,0.5,1,1.5,2]

Nx6 = len(real)
Ny6 = len(imag)

X06 = -5
Xn6 = 5
Y06 = -5
Yn6 = 5

Zpos = 0 # numerical value is irrelevant -> this is set to Z_DEPENDENT keyword for 3D data 


# fex6 = x1 + y1*1j
with open('testindexfileexample6v2.txt', 'w') as f:
    f.write('/rn,a,b/nx0\n')
    f.write('/r,qa,qb\n')
    f.write('/r\n')
    f.write(str(Nx6) + ' ' +  str(X06) + ' ' + str(Xn6) + ' ' + str(Zpos)  + ' ' + 'OUTPUT_REAL_IMAG_3D\n' )
    f.write(str(Ny6) + ' ' + str(Y06) + ' ' + str(Yn6) + '\n')
    for i in range(len(real)):
        f.write("%.3f" % real[i] + '    ' + "%.3f" % imag[0] + '    ' + "%.3f" % real[i] + '    ' + "%.3f" % imag[1] + '    ' +  "%.3f" % real[i] + '    ' + "%.3f" % imag[2] + '    ' + "%.3f" % real[i] + '    ' + "%.3f" % imag[3] + '    ' + "%.3f" % real[i] + '    ' + "%.3f" % imag[4] + '\n' ) 
       
       

#========================================== 3D datafile generator ==================================================================

Nx = 5 
Ny = 5 

X0 = -1 # for x range and y range as described above 
Y0 = -1 
Z0 = 0

Xn = 1   # in normalised coordinates, x' = 2x/w, y' = 2y/h and z' = z/l
Yn = 1 # so -1 to 1 in x and y means -w/2 to w/2 etc and 0 to 1 means from 0 to l in z 
Zn = 500


#z = [0,1,2,3,4,5,6,7,8,9,10]
nzreal = [1.1,1.09,1.08,1.07,1.06,1.05,1.04,1.03,1.02,1.01,1]
nzimag = [0.1,0.09,0.08,0.07,0.06,0.05,0.04,0.03,0.02,0.01,0]

for i in range(len(nzreal)):
    nzreal[i] = nzreal[i] - 1 

Nz = len(nzreal)


j =0 #counter for while loop 
l = 0

#with open('testindexfilecomplex.txt', 'w') as f:
#    for i in range(len(nzreal)):
#        while j < xres:  
#            f.write("%.3f" % nzreal[i] + '    ' + "%.3f" % nzimag[i] + '    ' + "%.3f" % nzreal[i] + '    ' + "%.3f" % nzimag[i] + '    ' +  "%.3f" % nzreal[i] + '    ' + "%.3f" % nzimag[i] + '\n' ) 
#            j = j + 1 
#        j = 0


with open('3dzdependentindex.txt', 'w') as f:
    f.write('/rn,a,b/nx0\n')
    f.write('/r,qa,qb\n')
    f.write('/r\n')
    f.write(str(Nx) + ' ' +  str(X0) + ' ' + str(Xn) + ' ' + 'Z_DEPENDENT' + ' ' + 'OUTPUT_REAL_IMAG_3D\n' )
    f.write(str(Ny) + ' ' + str(Y0) + ' ' + str(Yn) + '\n')
    f.write(str(Nz) + ' ' + str(Z0) + ' ' + str(Zn) + '\n')
    
    for i in range(Nz):
        while j < Nx:  
            while l < Ny:
                if(l < Ny - 1):
                    f.write("%.2f" % nzreal[i] + '    ' + "%.2f" % nzimag[i] + '    '  ) 
                else:
                    f.write("%.2f" % nzreal[i] + '    ' + "%.2f" % nzimag[i] + '\n' ) 
                l = l + 1
            l = 0    
            j = j + 1 
        j = 0







    