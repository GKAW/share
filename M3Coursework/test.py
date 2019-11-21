import os
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate
import pandas as pd


xlist = [0,1,2,3,4]
ylist = [5,4,3,2,1]
interplist = np.linspace(0,4,10)
print(interplist)

ynew = np.interp(interplist,xlist,ylist)
print(ynew)

plt.plot(interplist,ynew,'o',alpha=0.8)
plt.plot(xlist,ylist,'o',alpha=0.8)
plt.show()
