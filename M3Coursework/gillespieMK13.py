import os
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate
import pandas as pd


#===INITIALISE===
pltno = 5 # how many gillespie do you want?
endTime = 12000 # When should Gillespie stop? 25k seconds is ca 20 cycles (24000)

ABCiter = 2 # how many rounds of ABC do you want?
initk0min = 0
initk0max = 0.2
initk1min = 0
initk1max = 0.02



cutoff1 = 25
cutoff2 = 5
cutoff3 = 1

k0 = 0.2  #s^-1
k1 = 0.02 #s^-1
m0 = 0 # initial condition
t0 = 0 # intitial time

k0store = []
k1store = []
k0instore1 = []
k1instore1 = []
k0instore2 = []
k1instore2 = []
k0instore3 = []
k1instore3 = []

measuredVar = 10
measuredMean = 10
sqdmean = []
sqdvar = []

mmegavarstore = []
mmegaavgstore = []
mmegafanostore = []

mnewmegastore = []
mnewstore = []
mstore = []
mmegastore = [m0]
tstore = []
tmegastore = [t0]
cycle = 1200
def gil(t0,tmax,y0):
    # get stuff in the correct data format
    t0 = float(t0)
    tmax = float(tmax)
    y0 = float(y0)
    t = t0
    m = y0 # initialise substance
    nthCycle = 0

    # main loop
    while t < tmax:
        newt = t - (nthCycle*cycle)
        if newt >= cycle:
            m = np.random.binomial(m,0.5)
            nthCycle += 1
            mnewstore.append(m)

        r0 = k0
        r1 = k1*m
        rtot = r0+r1

        p0 = r0/rtot
        p1 = r1/rtot

        tstep = np.random.exponential(1/rtot)
        t += tstep

        rand0 = np.random.uniform()
        rand1 = np.random.uniform()

        if rand0 < p0:
            m += 1
        if rand1 < p1:
            m -= 1
        if m < 0:
            m = 0


        # STORE DATA
        global mstore
        global tstore
        mstore.append(m)
        tstore.append(t)

def megagil(n,t0,tmax,y0):
    # initialise correct size of store nested list
    global mmegastore
    global tmegastore
    global mnewmegastore
    mmegastore = [None] * n
    tmegastore = [None] * n
    mnewmegastore = [None] * n

    global mmegavarstore
    global mmegaavgstore
    global mmegafanostore
    mmegavarstore = [None] * n
    mmegaavgstore = [None] * n
    mmegafanostore = [None] * n

    for i in range(n):
        #print("\n \n \n \n \n \n \n \n \n \n \n") #if you want to clear the console in a hacky way.
        megagilprog = (float(i+1)/float(n))*100
        print("Progress: "+str(megagilprog)+"%")
        if megagilprog > 99.9999:
            print("Gillespie completed. Now generating stats and plots...")

        global k0
        global k1
        global k0store
        global k1store
        #randomise k rates for bayesian ABC
        k0 = np.random.uniform(initk0min,initk0max)
        k1 = np.random.uniform(initk1min,initk1max)

        gil(t0,tmax,y0)
        global mstore
        mmegastore[i] = mstore
        mstore=[]
        global tstore
        tmegastore[i] = tstore
        tstore=[]
        global mnewstore
        mnewmegastore[i] = mnewstore
        mnewstore=[]

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

        if sqdmeancurrent+sqdvarcurrent < cutoff1:
            global k0instore1
            global k1instore1
            k0instore1.append(k0)
            k1instore1.append(k1)

        if sqdmeancurrent+sqdvarcurrent < cutoff2:
            global k0instore2
            global k1instore2
            k0instore2.append(k0)
            k1instore2.append(k1)

        if sqdmeancurrent+sqdvarcurrent < cutoff3:
            global k0instore3
            global k1instore3
            k0instore3.append(k0)
            k1instore3.append(k1)

# HOW WOULD AVERAGING WORK? MAYBE NOT THE SAME AMOUNT OF TIME POINTS. ONLY AVERAGE FIRST 100s??
# to average: unfold megagil listed list. create a list per iteration point. average out these points. save stats about those too.
timevec = []
linter = []
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
    actlist = llist
    actlist = [list(i) for i in zip(*llist)] # unfold the llist (with * operator), then zip it. last timesteps if lists are different length

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


def prestats(tllist,mllist):
    global timevec
    timevec = np.linspace(0,endTime,(endTime*10))
    for i in range(len(tllist)):
        actinter = np.interp(list(timevec),list(tllist[i]),list(mllist[i]))
        global linter
        linter.append(actinter)


def statgil():
    #===EXECUTE===
    #try:
    megagil(pltno,0,endTime,0) # returns two nested lists (mmegastore and tmegastore) in the format mmegastore[i][j] where i is per Gillepsie simulation and j is the timepoint/mRNA number
    #except Exception:
    #    print('Could not compute Gillepsie.')

    #try:
    prestats(tmegastore, mmegastore)

    stats(linter)
    global mavgstore
    global mstdupstore
    global mstddownstore
    global mmaxstore
    global minter
    global mminstore
    mavgstore = avglist
    mstdupstore = stduplist
    mstddownstore = stddownlist
    mmaxstore = maxlist
    mminstore = minlist
    minter = linter

    global flat_m
    flat_m = [item for sublist in mmegastore for item in sublist]

    global mtotavg
    global mtotvar
    global mtotfano
    mtotavg = np.average(flat_m)
    mtotvar = np.var(flat_m)
    mtotfano = (mtotvar/mtotavg)

    print('mean mRNA:',mtotavg)
    print('var mRNA:',mtotvar)
    print('fano mRNA:',mtotfano)

    global flat_newm
    flat_newm = [item for sublist in mnewmegastore for item in sublist]

    global mnewtotavg
    global mnewtotvar
    global mnewtotfano
    mnewtotavg = np.average(flat_newm)
    mnewtotvar = np.var(flat_newm)
    mnewtotfano = (mnewtotvar/mnewtotavg)

    print('mean new mRNA:',mnewtotavg)
    print('var new mRNA:',mnewtotvar)
    print('fano new mRNA:',mnewtotfano)

#===start parameter search===
def ABC():
    statgil()
    for i in range(ABCiter):
        initk0min = min(k0instore1)
        initk0max = max(k0instore1)
        initk1min = min(k1instore1)
        initk1max = max(k1instore1)

        global cutoff1
        cutoff1 = cutoff1/5

        print('===ITERATION '+str(i+1)+'===')
        statgil()


#===ABC EXECUTION===
ABC()





#===ODE for overlay===
def dydt(t,y,k0,k1):
    m,n = y # looks like puts A and B into vector y
    dmdt = k0 - m * k1
    dndt = k0 - m * k1
    return(dmdt,dndt)

dydt_withks = lambda t,y: dydt(t,y,0.2,0.01)
sol = scipy.integrate.solve_ivp(dydt_withks, t_span=(0,endTime), y0=(1,1), method="RK45", rtol=1e-6)




#===PLOT===
#for i in range(pltno):
#    plt.plot(tmegastore[i],mmegastore[i])

mstar = [k0/k1] * 2
#try:
plt.plot([0,int(timevec[-1])],mstar,linewidth=1.5,color='r') # analytical steady state

plt.plot(sol.t,sol.y[0],'k',linewidth=1.5) #plot the ODE
'''
plt.plot(tmegastore[0],mmegastore[0])
plt.plot(tmegastore[1],mmegastore[1])
plt.plot(tmegastore[2],mmegastore[2])
plt.plot(tmegastore[3],mmegastore[3])
plt.plot(tmegastore[4],mmegastore[4])
'''
plt.plot(timevec,mavgstore,linewidth=1.25,color='b') # average of all Gillespie sims

plt.fill_between(timevec,mstdupstore,mstddownstore,facecolor='blue',alpha=0.5)
plt.fill_between(timevec,mmaxstore,mminstore,facecolor='blue',alpha=0.3)


plt.xlabel('time (seconds)')
plt.ylabel('mRNA (absolute number)')
plt.title('Gillespie simulation of mRNA production')

plt.figure()
plt.plot(k0store,k1store,'o')
plt.plot(k0instore1,k1instore1,'o')
plt.plot(k0instore2,k1instore2,'o')
plt.plot(k0instore3,k1instore3,'o')

plt.show()
#except Exception:
#    print('There was an error generating graphs from data. View raw avg data below:')
#    print(['%.2f' % elem for elem in mavgstore]) # also do the redundant formatting here so values are below each other
    #print(tavgstore)
#    print(['%.2f' % elem for elem in tavgstore])
