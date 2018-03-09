import numpy as np
import matplotlib.pyplot as plt

savePlots = 0
showPlots = 1

EB = 7.392*10**-19
a = 1.812*10**10
Re = 1.275*10**-10


xx = np.linspace(0,3.5*10**-10,10000)
vv = EB*(1-np.exp(-a*(xx-Re)))**2
psi0 = np.genfromtxt("../data/morsetest.dat")

plt.figure("morse")
plt.plot(xx,vv)
#plt.xlim(0.5e-10,3.5e-10)

plt.figure("morse0")
plt.plot(psi0,'-')
plt.show()

print(len(psi0))
print(psi0[499])
print(psi0[500])
