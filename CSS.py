global Loc_span, Loc_heigth, ListNodes, Elements, DataBeamDesing, DataColDesing, Wtotal, num_elems, ListNodesDrift, \
    cIndex, ListNodesBasal, T1m, Wtotal, IM, Sa_max, RDR_max, SDR_max, nrecs, RA_max, EleCol, EleBeam, DataColPhl,\
    VnVu_max

from datetime import datetime
from ReadRecord import ReadRecord
from getSaT import getSaT
from runNRHA_CSS import runNRHA_CSS
from Exci_pattern import Exci_pattern


def Read_Two_Column_File(file_name):
    with open(file_name, 'r') as data:
        x = []
        y = []
        for line in data:
            p = line.split()
            if p:
                x.append(str(p[0]).replace('a/', ''))
                y.append(float(p[1]))
    return x, y


start_time = datetime.now()
xi = float(self.ui.Damping.text()) / 100  # Damping
HL_directory = self.ui.HL_directory.text()
OutputCSSFile = self.ui.OutputCSSFile.text()
IMtype = self.ui.comboBoxIM_CSS.currentText()

print('^^^^^^^^ STARTING CSS ^^^^^^^^')
ind = 1
Sa_max, RDR_max, SDR_max, IM, RA_max, VnVu_max = [], [], [], [], [], []
Full_eqnms_list, Full_sf_list = [], []
EQnamev, Testv, HL_v, EQname_v, controlTimev, Tmaxv = [], [], [], [], [], []
while os.path.exists(HL_directory + '/Hazard_Level_' + str(ind)):
    archivo = HL_directory + '/Hazard_Level_' + str(ind) + '/listofgmotions.txt'
    print('Read: ' + archivo)
    eqnms_list, sf_list = Read_Two_Column_File(archivo)
    # print('eqnms_list =', eqnms_list)
    # print('sf_list =', sf_list)
    nrecs = len(eqnms_list)
    # print('nrecs =', nrecs)
    for i in range(nrecs):
        EQname = HL_directory + '/Hazard_Level_' + str(ind) + '/gmotions/' + eqnms_list[i].replace('.AT2', '')
        Full_eqnms_list.append(EQname)
        HL_v = np.append(HL_v, str(ind))
        EQname_v = np.append(EQname_v, eqnms_list[i].replace('.AT2', ''))
        Full_sf_list.append(sf_list[i])
    ind += 1

for EQname in Full_eqnms_list:
    if not os.path.exists(EQname + '.AT2'):
        print('No exist: ' + EQname + '.AT2')

for (EQname, sf) in zip(Full_eqnms_list, Full_sf_list):
    print('Run: ' + EQname)
    inFile = EQname + '.AT2'
    outFile = EQname + '.G4'
    dt, npts = ReadRecord(inFile, outFile)
    dur = npts * dt
    Sa, Sv, Sd, pga, amax = getSaT(outFile, dt, T1m, xi, npts)  # Get the PGA of the record in the X direction
    IMgeomean = Sa
    IM = np.append(IM, sf * IMgeomean)
    self.CreateNLM()
    Exci_pattern(dt, outFile, sf * g)
    if self.ui.radioButtonYesLC.isChecked():
        LC = 0
    else:
        LC = 1

    cIndex, mDT, mDPB, maxVB, maxRA, maxVnVu, Test, controlTime, Tmax = runNRHA_CSS(dt, dur, Loc_heigth, ListNodesDrift,
                                                                                    ListNodesBasal, EleCol, LC,
                                                                                    DataColPhl, DataColDesing,
                                                                                    ListNodesLC)
    Sa_max = np.append(Sa_max, maxVB / Wtotal)
    RDR_max = np.append(RDR_max, mDT)
    SDR_max = np.append(SDR_max, mDPB)
    EQnamev = np.append(EQnamev, EQname)
    RA_max = np.append(RA_max, maxRA)
    VnVu_max = np.append(VnVu_max, maxVnVu)
    Testv = np.append(Testv, Test)
    controlTimev = np.append(controlTimev, controlTime)
    Tmaxv = np.append(Tmaxv, Tmax)
    op.wipe
time_elapsed = datetime.now() - start_time
print('Time elapsed (hh:mm:ss.ms) {}'.format(time_elapsed))
controlTimev = np.around(controlTimev, decimals=2)
Tmaxv = np.around(Tmaxv, decimals=2)
maxDriftPisoBdg = np.around(SDR_max*100, decimals=2)
Run_Earthquake = np.vstack((HL_v, EQname_v, controlTimev, Tmaxv, maxDriftPisoBdg))
print(Run_Earthquake.T)
np.savetxt('CSS/' + OutputCSSFile + '_Test_Run_Earthquake.txt', Run_Earthquake.T, fmt='%s, %s, %s, %s, %s')
# Run_Earthquake = pd.DataFrame(Run_Earthquake).T
# Run_Earthquake.to_excel(excel_writer="CSS/Test_Run_Earthquake.xlsx")
