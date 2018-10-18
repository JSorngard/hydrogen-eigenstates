import matplotlib.pyplot as plt
import numpy as np

x,y=np.loadtxt("wavefunc.txt",delimiter=",",unpack=True)
plt.plot(x,y,label="Wavefunction")
plt.xlabel("x")
plt.legend()
plt.show()
