# -*- coding: utf-8 -*-
"""
Created on Mon Sep 19 17:16:39 2022

@author: araujorg
"""

##############

# In[]--------------------------------------------------------------
#### Define Analysis parameters

# Perform analysis to get equilibrium with lateral loads
tol = 1.0e-8*1; # max tolerance allowed

# Analysis
ops.wipeAnalysis()
ops.system("UmfPack") # create SOE
ops.numberer("RCM") # create DOF number
ops.constraints("Plain") # create constraint handler
ops.integrator('Newmark',  0.5,  0.25)
ops.test('EnergyIncr', tol, 50) # create test
ops.algorithm("KrylovNewton") # create algorithm
ops.analysis("Transient") # create analysis object



# In[]--------------------------------------------------------------
import time

Time_start = time.time()


#### Solve Loop
print('Running Nonlinear Dynamic Analysis...')

# initialize
tol0 = tol*1.0
dt0 = dtA*1.0

# convergence parameters
numConv = 30
timesConv = 0
stop = 0

# Returns current time in the analysis
ok = ops.analyze(1,dtA)

# Analysis time step and duration
tcurr = ops.getTime()


# enter analysis loop
while ok == 0 and np.abs(tcurr) < np.abs(TimeMax) and stop == 0:
    CurrDrift = CurrentDriftRatio(center_nodes)
    if CurrDrift>0.10:
        stop = 1
        ok = -9
        print('\n\n\nLarge drift -> %.3f Collapse'%CurrDrift,', time =',tcurr,'\n\n\n')
        continue


    ok = ops.analyze(1,dtA) 

    while ok != 0 and stop == 0:

        CurrDrift = CurrentDriftRatio(center_nodes)
        if CurrDrift>0.10:
            stop = 1
            ok = -9
            print('\n\n\nLarge drift -> %.3f Collapse'%CurrDrift,', time =',tcurr,'\n\n\n')
            continue


        # reset number of converged steps
        timesConv = 0

        # cut time step
        if np.abs(dtA) > np.abs(dt0/(10)):
            dtA = dtA/2.0
            print('\t\t *** cut dtA, dtA =', dtA,', time =',tcurr)
        
        else:
            # increase tol (if dtA is too small)
            dtA = dt0/(10)
            tol = tol*10.0
            print('\t\t *** increase tol, tol =', tol,', time =',tcurr)
            
        # run a step
        ops.test('EnergyIncr', tol, 100)
        ok = ops.analyze(1,dtA)

        # switch algorithms
        if ok != 0:
            ops.algorithm('BFGS')
            print('\t\t *** trying BFGS',', time =',tcurr)
            ok = ops.analyze(1,dtA)
            ops.algorithm('KrylovNewton')
            


        # convergence failure
        if tol > tol0*100 and ok != 0:
            stop = 1
            print('\t\t *** convergence failure',', time =',tcurr)

    # reset
    if ok == 0:
        timesConv = timesConv+1

        # reset tol
        if timesConv >= numConv and tol/10.0 >= tol0:
            timesConv = 0
            tol = tol/10.0
            ops.test('EnergyIncr', tol, 50)
            print('\t\t *** cut tol',', time =',tcurr)

        # reset dtA
        if timesConv >= numConv and np.abs(dtA)*2.0 <= np.abs(dt0):
            timesConv = 0
            dtA = dtA*2.0
            ops.test('EnergyIncr', tol, 50)
            ops.algorithm("KrylovNewton")
            print('\t\t *** increase dtA',', time =',tcurr)
        
    # current displacement
    tcurr = ops.getTime()
    

print ('ok: %d' %ok)

Time_end = time.time()
print(Time_end - Time_start,'s. for nonlinear THA')