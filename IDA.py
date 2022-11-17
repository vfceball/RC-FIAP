global Loc_span, Loc_heigth, ListNodes, Elements, DataBeamDesing, DataColDesing, Wtotal, num_elems, ListNodesDrift,\
	cIndex, ListNodesBasal, T1m, Wtotal, IMv, Sa_maxv, RDR_maxv, SDR_maxv, nrecs, dCap, T_v

# set pid [getPID]
# set np [getNP]
from datetime import datetime
from ReadRecord import ReadRecord
from getSaT import getSaT
from runNRHA3D import runNRHA3D
from Exci_pattern import Exci_pattern

# exec(open("runNRHA3D.py").read())

from load_PEERNGA_record import load_PEERNGA_record

start_time = datetime.now()

firstInt = float(self.ui.firstInt.text())  # This is the first intensity to run the elastic run (e.g. 0.05g)
incrStep = float(self.ui.incrStep.text())  # This is the increment used during the hunting phase (e.g. 0.10g)
maxRuns = int(self.ui.maxRuns.text())  # This is the maximum number of runs to use (e.g. 20)
dCap = float(self.ui.dCap.text())/100
xi = float(self.ui.Damping.text())/100  # Damping
IMtype = self.ui.comboBoxIM.currentText()
if not os.path.exists("IDA"):
    os.mkdir("IDA")
pflag = 1
# Open an error file that will log the IDA_HTF errors
error_log = open('IDA/IDA_HTF_error_log.txt', 'w')
print('^^^^^^^^ STARTING IDA HTF ^^^^^^^^')
error_log.write('^^^^^^^^ STARTING IDA HTF ^^^^^^^^')
# Get the ground motion set information
nmsfile = self.ui.nmsfile.text()
print('nmsfile =', nmsfile)
with open('IDA_Gmotions/' + nmsfile + '.txt', "r") as archivo:
    eqnms_list = [linea.rstrip() for linea in archivo]
print('eqnms_list =', eqnms_list)
nrecs = len(eqnms_list)
print('nrecs =', nrecs)
Sa_maxv, IMv, RDR_maxv, SDR_maxv = np.zeros([nrecs, maxRuns+1]), np.zeros([nrecs, maxRuns+1]),\
                                      np.zeros([nrecs, maxRuns+1]), np.zeros([nrecs, maxRuns+1])

for i in range(nrecs):
    IM_log = open('IDA/IM_' + str(i) + '.txt', "w")
    EQname = 'IDA_Gmotions/' + eqnms_list[i]
    inFile = EQname + '.at2'
    outFile = EQname + '.g4'
    dt, npts = ReadRecord(inFile, outFile)
    dur = npts * dt

    if IMtype == 'PGA':
        Sa, Sv, Sd, pga, amax, accelg = getSaT(outFile, dt, 0.0, xi, npts)  # Get the PGA of the record in the X direction
        IMgeomean = pga

    elif IMtype == 'Sa(T1)':
        Sa, Sv, Sd, pga, amax, accelg = getSaT(outFile, dt, T1m, xi, npts)  # Get the PGA of the record in the X direction
        IMgeomean = Sa
        print('T1m ={0:.3f} Sa ={1:.3f} amax= {2:.3f} pga = {3:.3f}'.format(T1m, Sa, amax, pga))

    print('IMgeomean = %.3f' % IMgeomean)
    # Set up the initial indices for HTF
    j = 1
    IM = []  # Initialise the list of IM used for printing
    # Sa_max, Sa_T1, RDR_max, SDR_max = [], [], [], []
    Sa_max, RDR_max, SDR_max = [], [], []
    IMlist = []  # This is just a list that will be used in filling
    hFlag = 1  # Hunting flag (1 for when we're hunting)
    tFlag = 0  # Tracing flag (0 at first)
    fFlag = 0  # Filling flag (0 at first)
    # print('j=', j, ' maxRuns = ', maxRuns)
    Tabla_IDA = open('IDA/Tabla_IDA_Record' + str(i) + '.txt', 'w')
    while j <= maxRuns:
        # As long as the hunting flag is 1, meaning we havent reached a collapse
        if hFlag == 1:
            # Determine the intensity to run at during the hunting (Andiamo a cacciare!)
            if j == 1:
                IM = np.append(IM, firstInt)
            else:
                IM = np.append(IM, IM[j - 2] + (j - 1) * incrStep)
            print('IM =', IM)
            # Determine the scale factor that needs to be applied to the record
            # Sa_T1 = np.append(Sa_T1, IM[j - 1])
            sf = IM[j - 1] / IMgeomean * g
            print('sf =', sf)
            run = 'Record' + str(i) + '_Run' + str(j)  # This is a tag that outputs will be labelled with
            log = open('IDA/log_IDA_' + run + '.txt', 'w')
            # The hunting intensity has been determined, now we can analyse
            self.CreateNLM()
            Exci_pattern(dt, outFile, sf)
            if pflag > 0:
                print('Record:{0:2d} Run:{1:2d} IM:{2:.3f}'.format(i, j, IM[j - 1]))
            cIndex, mDT, mDPB, maxVB = runNRHA3D(dt, dur, dCap, log, pflag, Loc_heigth, ListNodesDrift, ListNodesBasal)
            log.close
            j += 1
            print('cIndex =', cIndex, 'j =', j)
            # print('Sa_T1 =', Sa_T1)
            Sa_max = np.append(Sa_max, maxVB / Wtotal)
            RDR_max = np.append(RDR_max, mDT)
            SDR_max = np.append(SDR_max, mDPB)

            # Check the hunted run for collapse
            if cIndex > 0:
                hFlag = 0  # Stop hunting
                tFlag = 1  # Start tracing
                j -= 1  # Reduce by 1 because j was increased at end of hunting and we want to redo that point
                jhunt = j  # The value of j we hunted to
                if jhunt == 2:
                    error_log.write('WARNING:' + run + ' - Collapsed achieved on first increment, reduce increment...\n')
                else:
                    IM_log.write('{0:.3f}\n'.format(IM[j - 2]))  # j-2 because we've already increased j, but need to know if collapsed

            op.wipe
        # Close hunting

        # When the first collapse is reached, we start tracing between last convergence and the first collapse
        if tFlag == 1 and j <= maxRuns:
            # The first phase is to trace the last DeltaIM to get within the resolution
            if j == jhunt:
                firstC = IM[j - 1]  # This is the IM of the hunting collapse
                IM = IM[:-1]  # Remove that value of IM from the array (it's already been appended)
                # Sa_T1 = Sa_T1[:-1]
                # print('Sa_T1 =', Sa_T1)
                Sa_max = Sa_max[:-1]
                RDR_max = RDR_max[:-1]
                SDR_max = SDR_max[:-1]
            if j == 1:
                diff = firstC
            else:
                diff = firstC - IM[j - 2]  # Determine the difference between the hunting's noncollapse and collapse IM
            inctr = 0.20 * diff  # Take 0.2 of the difference
            if inctr < 0.05:
                inctr = 0.025  # Place a lower threshold on the increment so it doesnt start tracing too fine
            if j == 1:
                IMtr = inctr
            else:
                IMtr = IM[j - 2] + inctr  # Calculate new tracing IM, which is previous noncollapse plus increment
            IM = np.append(IM, IMtr)
            # Sa_T1 = np.append(Sa_T1, IM[j - 1])
            sf = IM[j - 1] / IMgeomean * g
            IM_log.write('{0:.3f}\n'.format(IMtr))
            run = 'Record' + str(i) + '_Run' + str(j)  # This is a tag that outputs will be labelled with
            log = open('IDA/log_IDA_' + run + '.txt', 'w')

            # The trace intensity has been determined, now we can analyse
            self.CreateNLM()
            # applying Dynamic Ground motion analysis
            Exci_pattern(dt, outFile, sf)
            if pflag > 0:
                print('Record:{0:2d} Run:{1:2d} IM:{2:.3f}'.format(i, j, IMtr))
            cIndex, mDT, mDPB, maxVB = runNRHA3D(dt, dur, dCap, log, pflag, Loc_heigth, ListNodesDrift, ListNodesBasal)
            log.close
            if cIndex == 0:
                Tabla_IDA.write('{0:.3f} {1:.3f}\n'.format(IMtr, mDT))

            # Check the hunted run for collapse
            if cIndex > 0:
                # Not sure if this is the best way, to just trace back up to collapse again
                tFlag = 0  # Stop tracing
                fFlag = 1  # Start filling
                jtrace = j  # The value of j we traced to
                IMlist = IM  # Get the list of IMs
                if j == jhunt:
                    # This means the first trace collapsed, should reduce the increment
                    error_log.write('WARNING:' + run + ' - First trace for collapse resulted in collapse...\n')
            j += 1
            # print('Sa_T1 =', Sa_T1)
            Sa_max = np.append(Sa_max, maxVB / Wtotal)
            RDR_max = np.append(RDR_max, mDT)
            SDR_max = np.append(SDR_max, mDPB)
            op.wipe
        # Close the tracing
        # When the required resolution is reached, we start filling
        # mDPB_collap = np.max(SDR_max)
        if fFlag == 1 and j <= maxRuns:
            # Reorder the list so we can account for filled runs
            IMlist = np.sort(IMlist)
            # Determine the biggest gap in IM for the hunted runs
            # mDPB = mDPB_collap
            # while mDPB >= mDPB_collap:
            gap = 0.0
            # We go to the end of the list minus 1 because, if not we would be filling between a noncollapsing and a collapsing run,
            # for which we are not sure if that filling run would be a non collapse - In short, does away with collapsing fills
            for ii in range(len(IMlist) - 1):
                temp = IMlist[ii] - IMlist[ii - 1]
                if temp > gap:
                    gap = temp
                    IMfil = IMlist[ii - 1] + gap / 2
            # print('j =', j)
            # print('IMfil =', IMfil)
            sf = IMfil / IMgeomean * g
            IM_log.write('{0:.3f}\n'.format(IMfil))
            run = 'Record' + str(i) + '_Run' + str(j)  # This is a tag that outputs will be labelled with
            log = open('IDA/log_IDA_' + run + '.txt', 'w')

            # The trace intensity has been determined, now we can analyse
            self.CreateNLM()
            Exci_pattern(dt, outFile, sf)

            if pflag > 0:
                print('Record:{0:2d} Run:{1:2d} IM:{2:.3f}'.format(i, j, IMfil))
            cIndex, mDT, mDPB, maxVB = runNRHA3D(dt, dur, dCap, log, pflag, Loc_heigth, ListNodesDrift, ListNodesBasal)
            IM = np.append(IM, IMfil)
            IMlist = np.append(IMlist, IMfil)
            log.close
            if cIndex == 0:
                Tabla_IDA.write('{0:.3f} {1:.3f}\n'.format(IMfil, mDT))
            j += 1
            # print('Sa_T1 =', Sa_T1)
            Sa_max = np.append(Sa_max, maxVB / Wtotal)
            RDR_max = np.append(RDR_max, mDT)
            SDR_max = np.append(SDR_max, mDPB)
            op.wipe
        # Close the filling
        # Wrap it up and finish
        if j == maxRuns and hFlag == 1:
            error_log.write('WARNING:' + run + ' - Collapse not achieved, increase increment or number of runs...\n')
        if j == maxRuns and fFlag == 0:
            error_log.write(
                'WARNING:' + run + ' - No filling, algorithm still tracing for collapse (reduce increment & increase runs)...\n')
        op.wipe
    # Close the maxRuns while loop
    IM_log.close
    Tabla_IDA.close
    ind = np.argsort(IM)
    # print('ind =', ind)
    IM = IM[ind]
    # print('IM =', IM)
    Sa_max = Sa_max[ind]
    RDR_max = RDR_max[ind]
    SDR_max = SDR_max[ind]
    # print('SDR_max', SDR_max)
    # print('Zeros_SDR_max', np.count_nonzero(SDR_max))
    ind = np.max(ind)+1
    IMv[i, 1:ind+1] = IM
    Sa_maxv[i, 1:ind+1] = Sa_max
    RDR_maxv[i, 1:ind+1] = RDR_max
    SDR_maxv[i, 1:ind+1] = SDR_max
print('^^^^^^^^ FINISHED IDA HTF ^^^^^^^^')
# print('IMv =', IMv)
# print('Sa_maxv =', Sa_maxv)
# print('RDR_maxv =', RDR_maxv)
# print('SDR_maxv =', SDR_maxv)

error_log.write('^^^^^^^^ FINISHED IDA HTF ^^^^^^^^\n')
error_log.close

time_elapsed = datetime.now() - start_time
print('Time elapsed (hh:mm:ss.ms) {}'.format(time_elapsed))

