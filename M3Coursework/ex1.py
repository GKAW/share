import sys
import numpy as np

#===DEFINE VARIABLES===#
k0 = 0.2  #s^-1
k1 = 0.01 #s^-1
mstar = k0/k1
m0 = 0 # initial condition
m = 5 # init m
#===DEFINE EQUATIONS===#
dmdt = k0 - m * k1

def integrate(eqn,var,y0,inc,time):
    for i in np.arange(start=inc, stop=(time+inc), step=inc):
        try:
            step += 1 # step counter
            print('step: ',step)
        except Exception:
            step = 0 # create step counter
            print('initialised step counter')
        try:
            if var == -5:
                print(time/inc+1)
                print(int(time/inc+1))
                var=np.zeros(int(time/inc+1))
        except Exception:
            if step==1:
                print(var)
        print(i)

integrate(dmdt,m,0,0.1,1)
