import numpy as np
import matplotlib.pyplot as plt

savePlots = 0
showPlots = 1
h=16/1000
psi0 = np.genfromtxt("../data/harmE0.dat",unpack=True)
psi1 = np.genfromtxt("../data/harmE1.dat",unpack=True)
psi2 = np.genfromtxt("../data/harmE2.dat",unpack=True)
psi3 = np.genfromtxt("../data/harmE3.dat",unpack=True)
psim8 = np.genfromtxt("../data/harmE0m8.dat",unpack=True)
x1=np.linspace(-8+h,-8+501*h,501)
x2 = np.linspace(8-h,8-501*h,501)
xx = np.linspace(-4,4,1000)
yy = xx*xx/2

plt.figure("Eigenstates")
plt.plot(x1,psi0[:501]/np.max(psi0),'b-')
plt.plot(x2,psi0[501:]/np.max(psi0),'b-')

plt.plot(x1,2+psi1[:501]/np.max(psi1),'r-')
plt.plot(x2,2+psi1[501:]/np.max(psi1),'r-')

plt.plot(x1,4+psi2[:501]/np.max(psi2),'g-')
plt.plot(x2,4+psi2[501:]/np.max(psi2),'g-')

plt.plot(x1,6+psi3[:501]/np.max(psi3),'y-')
plt.plot(x2,6+psi3[501:]/np.max(psi3),'y-')

plt.plot(xx,yy)

plt.figure("comparison")
#plt.semilogy(x1,psi0[:501],'b-')
#plt.semilogy(x2,psi0[501:],'b-')
#plt.plot(x1,psi0[:501]/psim8[:501],'r')
#plt.plot(x2,psi0[501:]/psim8[501:],'r')

#plt.semilogy(x1,psim8[:501],'r-')
#plt.semilogy(x2,psim8[501:],'r-')


psi0 = np.genfromtxt("../data/harm2E0.dat",unpack=True)
psi1 = np.genfromtxt("../data/harm2E1.dat",unpack=True)
psi2 = np.genfromtxt("../data/harm2E2.dat",unpack=True)
psi3 = np.genfromtxt("../data/harm2E3.dat",unpack=True)

plt.figure("EigenStates")
plt.plot(np.linspace(-8,8,1001),psi0/np.max(psi0))
plt.plot(np.linspace(-8,8,1001),2+psi1/np.max(psi1))
plt.plot(np.linspace(-8,8,1001),4+psi2/np.max(psi2))
plt.plot(np.linspace(-8,8,1001),6+psi3/np.max(psi3))
plt.plot(xx,yy)
plt.ylim(0,8)

psiLeft = np.genfromtxt("../data/harmleftE0.dat")
plt.figure("Blowup")
plt.plot(np.linspace(-8,8,1001),psiLeft,'r')
plt.ylim(-10000000,5000000)

if showPlots:
    plt.show()
