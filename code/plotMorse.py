import numpy as np
import matplotlib.pyplot as plt

savePlots = 1
showPlots = 1

EB = 7.392*10**-19
a = 1.812*10**10
Re = 1.275*10**-10


xx = np.linspace(0.5*10**-10,3*10**-10,10000)
vv = EB*(1-np.exp(-a*(xx-Re)))**2
psi0 = np.genfromtxt("../data/morseE0.dat")
psi1 = np.genfromtxt("../data/morseE1.dat")
psi2 = np.genfromtxt("../data/morseE2.dat")
psi3 = np.genfromtxt("../data/morseE3.dat")
psi4 = np.genfromtxt("../data/morseE4.dat")
psi5 = np.genfromtxt("../data/morseE5.dat")

plt.figure("morse")
plt.plot(xx,vv/(0.2*vv[-1]),label=r"$\mathcal{V}(x)$")

plt.plot(np.linspace(0.5*10**-10,3*10**-10,1001),psi0/np.max(2*psi0),label=r"$\Psi_0$")
plt.plot(np.linspace(0.5*10**-10,3*10**-10,1001),1+psi1/np.max(2*psi1),label=r"$\Psi_1$")
plt.plot(np.linspace(0.5*10**-10,3*10**-10,1001),2+psi2/np.max(2*psi2),label=r"$\Psi_2$")
plt.plot(np.linspace(0.5*10**-10,3*10**-10,1001),3+psi3/np.max(2*psi3),label=r"$\Psi_3$")
plt.plot(np.linspace(0.5*10**-10,3*10**-10,1001),4+psi4/np.max(2*psi4),label=r"$\Psi_4$")
plt.plot(np.linspace(0.5*10**-10,3*10**-10,1001),5+psi5/np.max(2*psi5),label=r"$\Psi_5$")
plt.xlabel("R [m]")
plt.yticks((0,1,2,3,4,5),("E0","E1","E2","E3","E4","E5"))
plt.ylim(0,6)
plt.legend(loc="upper right")
if savePlots:
    plt.savefig("../figs/morse-Eigenstates.pdf")

if showPlots:
    plt.show()

print(len(psi0))
print(psi0[651])
print(psi0[650])
