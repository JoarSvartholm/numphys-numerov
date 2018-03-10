import numpy as np
import matplotlib.pyplot as plt

savePlots = 0
showPlots = 1
h=16/1000

psi0 = np.genfromtxt("../data/harm2E0.dat",unpack=True)
psi1 = np.genfromtxt("../data/harm2E1.dat",unpack=True)
psi2 = np.genfromtxt("../data/harm2E2.dat",unpack=True)
psi3 = np.genfromtxt("../data/harm2E3.dat",unpack=True)
psim8 = np.genfromtxt("../data/harm2E0m8.dat")
x = np.linspace(-8,8,1001)
xx = np.linspace(-4,4,1000)
yy = xx*xx/2
exx = -0.00019999999*np.exp(x*x/2)

plt.figure("EigenStates")
plt.plot(xx,yy,label=r"$\mathcal{V}(x)$")
plt.plot(x,psi0/np.max(2*psi0),label=r"$\Psi_0$")
plt.plot(x,1+psi1/np.max(2*psi1),label=r"$\Psi_1$")
plt.plot(x,2+psi2/np.max(2*psi2),label=r"$\Psi_2$")
plt.plot(x,3+psi3/np.max(2*psi3),label=r"$\Psi_3$")
plt.xlabel("x [a.u]")
plt.yticks((0,1,2,3),("E0","E1","E2","E3"))
plt.ylim(0,4)
plt.legend(loc="upper right")
if savePlots:
    plt.savefig("../figs/harm-Eigenstates.pdf")

plt.figure("comparison")
plt.plot(x,psi0,label=r"$Y_1 = 10^{-12}$")
plt.plot(x,psim8,'--',label=r"$Y_1 = 10^{-8}$")
plt.legend()
plt.xlabel("x [a.u]")
if savePlots:
    plt.savefig("../figs/harm-comparison.pdf")

psiLeft = np.genfromtxt("../data/harmleftE0.dat")
plt.figure("Blowup")
plt.plot(np.linspace(-8,8,1001),psiLeft,'r',label=r"$\Psi_{left}$")
plt.plot(x,exx,'--',label=r"$-0.00019999999 \cdot e^{x^2/2}$")
plt.xlabel("x [a.u]")
plt.legend()
plt.ylim(-5000000,5000000)
if savePlots:
    plt.savefig("../figs/harm-blowup.pdf")

if showPlots:
    plt.show()
