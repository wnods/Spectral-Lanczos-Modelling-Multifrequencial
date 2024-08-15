'''

'''

import numpy as np
import matplotlib.pyplot as plt
import math as ma

ndata = 81

x  = np.zeros(ndata,np.float)
rHy = np.zeros(ndata,np.float)
iHy = np.zeros(ndata,np.float)

rHya = np.zeros(ndata,np.float)
iHya = np.zeros(ndata,np.float)

rEx = np.zeros(ndata,np.float)
iEx = np.zeros(ndata,np.float)
rEz = np.zeros(ndata,np.float)
iEz = np.zeros(ndata,np.float)

rho = np.zeros(ndata,np.float)

xi = -3000.0
xf =  3000.0
dx = (xf-xi)/(ndata-1)
for i in range(ndata):
    x[i] = xi + i * dx


lw = 2.0
MEDIUM_SIZE = 13
filefig = 'campoHy.png'

cH = open('Hyfield.dat','r')

cExEz = open('ExEzfield.dat','r')

rhoap = open('Res_apar.dat','r')

HyExa = open('perfilagem.dat','r')

m = -1
for i in cH:
	m += 1
	rw = i.split()
	rHy[m]  = np.float(rw[0])
	iHy[m] = np.float(rw[1])

m = -1
for i in cExEz:
	m += 1
	rw = i.split()
	rEx[m]  = np.float(rw[0])
	iEx[m] = np.float(rw[1])
	rEz[m]  = np.float(rw[2])
	iEz[m] = np.float(rw[3])


m = -1
for i in rhoap:
	m += 1
	rw = i.split()
	rho[m]  = np.float(rw[0])

m = -1
for i in HyExa:
	m += 1
	rw = i.split()
	rHya[m]  = np.float(rw[0])
	iHya[m]  = np.float(rw[1])

plt.rc('font', size=MEDIUM_SIZE)

plt.figure(figsize=(10.5,6.5))

plt.subplot(2,1,1)

plt.plot(x,rHy,'or',label= r'$ re(2D) $',linewidth=lw,mfc='None')
plt.plot(x,iHy,'sb',label= r'$ im(2D) $',linewidth=lw,mfc='None')
plt.plot(x,rHya,'-k',label= r'$ re(analytic-1D) $',linewidth=lw,mfc='None')
plt.plot(x,iHya,'-g',label= r'$ im(analytic-1D) $',linewidth=lw,mfc='None')
plt.legend(loc='best')
#plt.legend(loc='upper right')
plt.ylabel(r'$ H_y^s \ (A/m) $')
#plt.xlabel('x (m)')
#plt.xlim([xi,xf])
plt.grid()

plt.subplot(2,1,2)

plt.plot(x,100*abs((rHya-rHy)/rHya),'-k',label= r'$ re $',linewidth=lw,mfc='None')
plt.plot(x,100*abs((iHya-iHy)/iHya),'-g',label= r'$ im $',linewidth=lw,mfc='None')
plt.legend(loc='best')
#plt.legend(loc='upper right')
plt.ylabel(r' Rel. Dif. $ \% $')
plt.xlabel('x (m)')
#plt.xlim([xi,xf])
plt.grid()

#plt.figure()
#plt.subplot(3,1,2)

#plt.plot(x,rEx,'-or',label= r'$ re $',linewidth=lw,mfc='None')
#plt.plot(x,iEx,'-ob',label= r'$ im $',linewidth=lw,mfc='None')
#plt.legend(loc='upper right')
#plt.ylabel(r'$ E_x \ (V/m) $')
#plt.xlabel('x (m)')
#plt.xlim([xi,xf])
#plt.ylim([0.02,0.08])
#plt.grid()

#plt.figure()
#plt.subplot(3,1,3)

#plt.plot(x,rEz,'-or',label= r'$ real. $',linewidth=lw,mfc='None')
#plt.plot(x,iEz,'-ob',label= r'$ imag. $',linewidth=lw,mfc='None')
#plt.legend(loc='upper right')
#plt.ylabel(r'$ E_z \ (V/m) $')
#plt.xlabel('x (m)')
#plt.xlim([xi,xf])
#plt.grid()

#plt.figure()

#plt.plot(x,rho,'-sb',label= r'$ \rho_a $',linewidth=lw,mfc='None')
#plt.legend(loc='upper right')
#plt.ylabel(r'$ \rho_a \ (\Omega m) $')
#plt.xlabel('x (m)')
#plt.yscale('log')
#plt.xlim([xi,xf])
#plt.grid()

plt.savefig(filefig)
plt.show()
