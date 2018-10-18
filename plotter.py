import matplotlib.pyplot as plt
import numpy as np

xd,yd=np.loadtxt("/home/johan/Documents/Assignment 5/wavefunc.txt",delimiter=",",unpack=True)
#x,y=np.loadtxt("/home/johan/Documents/Assignment 5/spl.txt",delimiter=",",unpack=True)
#plt.plot(x,y,label="spline")
plt.plot(xd,yd,label="Wavefunction")
plt.xlabel("x")
#plt.legend()
plt.show()
