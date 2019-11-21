import os
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate
import pandas as pd
import math


#===INITIALISE===
pltno = 100 # how many gillespie do you want?
endTime = 12000 # When should Gillespie stop? 25k seconds is ca 20 cycles (24000)

cutoff = 10

k0 = 0.2 #s^-1
k1 = 0.02 #s^-1
k0store = []
k1store = []
k0instore = []
k1instore = []
#k0 = np.random.uniform(0,0.2)
#k1 = np.random.uniform(0,0.02)
mstar = k0/k1 # from ODE: k0/k1

measuredVar = 10
measuredMean = 10
sqdmean = []
sqdvar = []

m0 = 20 # initial condition
t0 = 0 # intitial time

mstore = []
mmegastore = [m0]
mmegavarstore = []
mmegaavgstore = []
mmegafanostore = []
tstore = []
tmegastore = [t0]
cycle = 1200
def gil(t0,tmax,y0,split):
    # get stuff in the correct data format
    t0 = float(t0)
    tmax = float(tmax)
    y0 = float(y0)
    t = t0
    m = y0 # initialise substance
    nthCycle = 0

    # main loop
    while t < tmax:
        if split == 1:
            newt = t - (nthCycle*cycle)
            if newt >= cycle:
                m = np.random.binomial(m,0.5) # to account for 50/50 binomial distribution at division cell cycle
                nthCycle += 1

        r0 = k0
        r1 = k1*m

        t0step = np.random.exponential(1/r0)
        if r1 == 0:
            t1step = math.inf
        else:
            t1step = np.random.exponential(1/r1)

        if t0step < t1step:
            m += 1
            t += t0step
        else:
            m -= 1
            t += t1step

        # STORE DATA
        global mstore
        global tstore
        mstore.append(m)
        tstore.append(t)


def megagil(n,t0,tmax,y0):
    # initialise correct size of store nested list
    global mmegastore
    global tmegastore
    global mmegavarstore
    global mmegaavgstore
    global mmegafanostore

    mmegastore = [m0] * n
    tmegastore = [None] * n
    mmegavarstore = [None] * n
    mmegaavgstore = [None] * n
    mmegafanostore = [None] * n


    for i in range(n):
        '''
        #print("\n \n \n \n \n \n \n \n \n \n \n") #if you want to clear the console in a hacky way.
        megagilprog = (float(i+1)/float(n))*100
        print("Progress: "+str(megagilprog)+"%")
        if megagilprog > 99.9999:
            print("Gillespie completed. Now generating stats and plots...")
        '''
        global k0
        global k1
        global k0store
        global k1store
        #randomise k rates for bayesian ABC
        k0 = np.random.uniform(0,0.2)
        k1 = np.random.uniform(0,0.02)

        split = 0

        gil(t0,tmax,y0,split)
        global mstore
        mmegastore[i] = mstore
        mstore=[]
        global tstore
        tmegastore[i] = tstore
        tstore=[]

        mmegavarstore[i] = np.var(mmegastore[i])
        mmegaavgstore[i] = np.average(mmegastore[i])
        mmegafanostore[i] = (mmegavarstore[i]/mmegaavgstore[i])

        #collect stuff for sorting to successful and unsuccessful experiment
        sqdmeancurrent = (measuredMean-mmegaavgstore[i])**2
        sqdvarcurrent = (measuredVar-mmegavarstore[i])**2
        sqdmean.append(sqdmeancurrent)
        sqdvar.append(sqdvarcurrent)

        k0store.append(k0)
        k1store.append(k1)

        if sqdmeancurrent+sqdvarcurrent < cutoff:
            global k0instore
            global k1instore
            k0instore.append(k0)
            k1instore.append(k1)

# HOW WOULD AVERAGING WORK? MAYBE NOT THE SAME AMOUNT OF TIME POINTS. ONLY AVERAGE FIRST 100s??
# to average: unfold megagil listed list. create a list per iteration point. average out these points. save stats about those too.

def stats(llist):
    global avglist
    global stduplist
    global stddownlist
    global maxlist
    global minlist
    avglist=[]
    stduplist=[]
    stddownlist=[]
    maxlist=[]
    minlist=[]
    actlist = [list(i) for i in zip(*llist)] # unfold the llist (with * operator), then zip it. last timesteps if lists are different length
    #(loss of data though... maybe better to interpolate based on a fixed value in time?) Would do this if not just coursework probs. for now too much effort.

    for i in range(len(actlist)):
        """
        statsprog = (float(i+1)/float(len(actlist)))*100
        print("\n \n \n \n \n \n \n \n \n \n \n") #if you want to clear the console in a hacky way.
        print("Progress: "+str(statsprog)+"%")
        if statsprog > 99.9999:
            print("Stats calc completed.")
        """
        avglist.append(np.average(actlist[i]))
        stduplist.append(np.average(actlist[i])+np.std(actlist[i]))
        stddownlist.append(np.average(actlist[i])-np.std(actlist[i]))
        maxlist.append(np.max(actlist[i]))
        minlist.append(np.min(actlist[i]))


#===ODE for overlay===
def dydt(t,y,k0,k1):
    m,n = y # looks like puts A and B into vector y
    dmdt = k0 - m * k1
    dndt = k0 - m * k1
    return(dmdt,dndt)

dydt_withks = lambda t,y: dydt(t,y,k0,k1)
sol = scipy.integrate.solve_ivp(dydt_withks, t_span=(0,endTime), y0=(m0,m0), method="RK45", rtol=1e-6)


#===EXECUTE===
#try:
megagil(pltno,t0,endTime,m0) # returns two nested lists (mmegastore and tmegastore) in the format mmegastore[i][j] where i is per Gillepsie simulation and j is the timepoint/mRNA number
#except Exception:
#    print('Could not compute Gillepsie.')

#try:
stats(mmegastore)
mavgstore = avglist
mstdupstore = stduplist
mstddownstore = stddownlist
mmaxstore = maxlist
mminstore = minlist

flat_m = [item for sublist in mmegastore for item in sublist]

mtotavg = np.average(flat_m)
mtotvar = np.var(flat_m)
mtotfano = (mtotvar/mtotavg)

stats(tmegastore)
tavgstore = avglist
tstdupstore = stduplist
tstddownstore = stddownlist
tmaxstore = maxlist
tminstore = minlist


print(mmegaavgstore)
print(mmegavarstore)
print(mmegafanostore)
print(mtotavg)
print(mtotvar)
print(mtotfano)

#except Exception:
#    print("Stats could not be generated from Gillespie data.")


# needs patching: include t0 and m0


#===PLOT===
#for i in range(pltno):
#    plt.plot(tmegastore[i],mmegastore[i])
#mstar = [k0/k1] * 2
mstar = [mstar] * 2
try:
    #plt.figure()
    #plt.plot([0,int(tavgstore[-1])],mstar,linewidth=1.5,color='r') # analytical steady state

    #plt.plot(sol.t,sol.y[0],'k--',linewidth=1.5) #plot the ODE

    #plt.plot(tavgstore,mavgstore,linewidth=1.25,color='b') # average of all Gillespie sims

    #plt.fill_between(tavgstore,mstdupstore,mstddownstore,facecolor='blue',alpha=0.5)
    #plt.fill_between(tavgstore,mmaxstore,mminstore,facecolor='blue',alpha=0.3)


    #plt.xlabel('time (seconds)')
    #plt.ylabel('mRNA (absolute number)')
    #plt.title('Gillespie simulation of mRNA production')

    plt.figure()
    plt.hist(mmegaavgstore,alpha=0.75,bins=range(20),label="Squared diff. mean")
    plt.hist(mmegavarstore,alpha=0.75,bins=range(30),label="Squared diff. variance")
    plt.xlim((0,30))
    plt.xlabel("Squared difference")
    plt.ylabel("Abundance")
    plt.legend()
    plt.title("Squared difference: Without cell division")

    plt.figure()
    plt.hist(k0instore,alpha=0.75,color='g')
    #plt.hist(k0store,alpha=0.75,color='r')

    plt.figure()
    plt.hist(k1instore,alpha=0.75,color='g')
    #plt.hist(k1store,alpha=0.75,color='r')
    plt.xlim((0,0.02))

    plt.show()

    print(k0store)
    print(len(k0store))
    print(k0instore)
    print(len(k0instore))
    print(k1store)
    print(len(k1store))
    print(k1instore)
    print(len(k1instore))
except Exception:
    print('There was an error generating graphs from data. View raw avg data below:')
    print(['%.2f' % elem for elem in mavgstore]) # also do the redundant formatting here so values are below each other
    #print(tavgstore)
    print(['%.2f' % elem for elem in tavgstore])
