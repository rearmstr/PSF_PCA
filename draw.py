import pyfits
import numpy as np
import matplotlib.pyplot as plt


f=pyfits.open('test.fits')
coeff=f['vec'].data.transpose()
data=f['odata'].data.transpose()
tdata=f['tdata'].data.transpose()
rdata=f['reco'].data.transpose()
x=[]
n=len(coeff[0])
for i in range(0,n):
    x.append(2*np.pi/n*(i+0.5))
f=plt.figure(figsize=(14,8))
ax1=f.add_subplot(121)
ax2=f.add_subplot(122)

nn=3
for i in range(0,nn):
    ax1.plot(x,coeff[i],'o-')

nn=5
for i in range(0,nn):

    ax2.plot(x,data[i],'o',0,5)
    ax2.plot(x,tdata[i],'-',color='black')
    ax2.plot(x,rdata[i],'-',color='red')
plt.ylim([-250,250])
plt.xlim([0,6.28])
f.show()
