import numpy as np
import matplotlib.pyplot as plt

h=16/1000
psi = np.genfromtxt("../data/harm1.dat",unpack=True)
x1=np.linspace(-8+h,-8+501*h,501)
x2 = np.linspace(8-h,8-501*h,501)

plt.figure("psi")
plt.plot(x1,psi[:501],'.-')
plt.plot(x2,psi[501:],'.-')
plt.show()
