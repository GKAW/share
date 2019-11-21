import scipy.integrate
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
print("all imports complete.")

def dydt(t,y,k0,k1):
    m,n = y # looks like puts A and B into vector y
    dmdt = k0 - m * k1
    dndt = k0 - m * k1
    return(dmdt,dndt)

dydt_withks = lambda t,y: dydt(t,y,0.2,0.01)
sol = scipy.integrate.solve_ivp(dydt_withks, t_span=(0,1200), y0=(1,1), method="RK45", rtol=1e-6)


plt.plot(sol.t,sol.y[0],'b',linewidth="2") #for A
plt.title("Chemical Reactions (ex1)")
plt.xlabel("time (arbitrary)")
plt.ylabel("Relative Concentration")
plt.grid(True)
plt.show()

#done.
print(sol.y[1])
