global Loc_span, Loc_heigth, ListNodes, Elements, DataBeamDesing, DataColDesing, Wtotal, num_elems, ListNodesDrift, \
    cIndex, ListNodesBasal, T1m, Wtotal, IM, Sa_max, RDR_max, SDR_max, nrecs, RA_max, EleCol, EleBeam, DataColPhl,\
    VnVu_max, HL_directory, list_beams, list_cols, maxPhRot_Colv, DataBeamPhl, maxPhRot_Beamv, maxSDRBdg, maxSDRBdgv,\
    maxPhRot_Colcv, maxPhRot_Beamcv, MedPhRot_Colmv_v, MedPhRot_Beammv_v

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
num_beams = len(EleBeam)
num_cols = len(EleCol)
rangCols01 = np.sort(np.append(np.append(0, np.arange(1, 12*num_cols+1, 12)), np.arange(2, 12*num_cols+1, 12)))
rangCols06 = np.sort(np.append(np.append(0, np.arange(11, 12*num_cols+1, 12)), np.arange(12, 12*num_cols+1, 12)))
rangBeams01 = np.sort(np.append(np.append(0, np.arange(1, 12*num_beams+1, 12)), np.arange(2, 12*num_beams+1, 12)))
rangBeams06 = np.sort(np.append(np.append(0, np.arange(11, 12*num_beams+1, 12)), np.arange(12, 12*num_beams+1, 12)))
rangDriftNode = np.int_(np.append(0, ListNodesDrift[:, 0]+1))

print('num_beams =', num_beams)
print('rangBeams01 =', rangBeams01)
print('rangDriftNode =', rangDriftNode)
# rangCols01.astype(int)
# rangCols06.astype(int)
# rangBeams01.astype(int)
# rangBeams06.astype(int)


print('^^^^^^^^ STARTING CSS ^^^^^^^^')
ind = 1
floors_num = len(Loc_heigth) - 1
Sa_max, RDR_max, SDR_max, IM, RA_max, VnVu_max = [], [], [], [], [], []
Full_eqnms_list, Full_sf_list = [], []
EQnamev, Testv, HL_v, EQname_v, controlTimev, Tmaxv = [], [], [], [], [], []
maxPhRot_Colv = np.zeros(floors_num)
maxPhRot_Beamv = np.zeros(floors_num)
maxPhRot_Colcv = np.zeros(floors_num)
maxPhRot_Beamcv = np.zeros(floors_num)
maxSDRBdgv = np.zeros(floors_num)

MedPhRot_Colmv_v = np.zeros(floors_num)
MedPhRot_Beammv_v = np.zeros(floors_num)

while os.path.exists(HL_directory + '/Hazard_Level_' + str(ind)):
    data_dir = HL_directory + '/Hazard_Level_' + str(ind) + '/Data'
    if not os.path.exists(data_dir):
        os.mkdir(data_dir)
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

    op.wipeAnalysis()

    # Recording of forces and deformations from nonlinear analysis
    # slist_beams = str(list_beams)
    # slist_cols = str(list_cols)
    # sListNodes = str(ListNodes)
    slist_beams = np.array(list_beams)
    slist_cols = np.array(list_cols)
    sListNodes = ListNodes[:, 0]
    # print(type(slist_beams))
    # print(type(slist_cols))
    # print(type(sListNodes))
    # print('sListNodes =', sListNodes)

    if self.ui.checkBoxSaveCSS.isChecked() == True:
        data_dir = EQname.replace('gmotions', 'Data')
        op.recorder('Element', '-file', data_dir + '_DefoSec_Beams.out', '-closeOnWrite',
                    '-time', '-eleRange', int(slist_beams[0]), int(slist_beams[-1]), 'section', 'deformation')
        op.recorder('Element', '-file', data_dir + '_ForceSec_Beams.out', '-closeOnWrite',
                    '-time', '-eleRange', int(slist_beams[0]), int(slist_beams[-1]), 'section', 'force')
        op.recorder('Element', '-file', data_dir + '_DefoSec_Cols.out', '-closeOnWrite',
                    '-time', '-eleRange', int(slist_cols[0]), int(slist_cols[-1]), 'section', 'deformation')
        op.recorder('Element', '-file', data_dir + '_ForceSec_Cols.out', '-closeOnWrite',
                    '-time', '-eleRange', int(slist_cols[0]), int(slist_cols[-1]), 'section', 'force')
        op.recorder('Node', '-file', data_dir + '_HoriNodes.out', '-closeOnWrite',
                    '-time', '-nodeRange', int(sListNodes[0]), int(sListNodes[-1]), '-dof', 1, 'disp')
        op.recorder('Node', '-file', data_dir + '_VertNodes.out', '-closeOnWrite',
                    '-time', '-nodeRange', int(sListNodes[0]), int(sListNodes[-1]), '-dof', 2, 'disp')

    cIndex, mDT, mDPB, maxVB, maxRA, maxVnVu, Test, controlTime, Tmax, maxPhRot_Col, maxPhRot_Beam, maxSDRBdg, \
    maxPhRot_Colc, maxPhRot_Beamc, MedPhRot_Colm_v, MedPhRot_Beamm_v = runNRHA_CSS(dt, dur, Loc_heigth, Loc_span,
                                                                                   ListNodesDrift, ListNodesBasal,
                                                                                   EleCol, EleBeam, LC, DataColPhl,
                                                                                   DataColDesing, ListNodesLC, EQname,
                                                                                   HL_directory, ListNodes, DataBeamPhl,
                                                                                   DataBeamDesing)

    if self.ui.checkBoxSaveCSS.isChecked() == True:
        DefoSec_Beams = np.loadtxt(data_dir + '_DefoSec_Beams.out')
        DefoSec_Beams = np.array(DefoSec_Beams)
        DefoSec01_Beams = DefoSec_Beams[:, rangBeams01]
        DefoSec06_Beams = DefoSec_Beams[:, rangBeams06]
        np.savetxt(data_dir + '_DefoSec01_Beams.out', DefoSec01_Beams)
        np.savetxt(data_dir + '_DefoSec06_Beams.out', DefoSec06_Beams)

        ForceSec_Beams = np.loadtxt(data_dir + '_ForceSec_Beams.out')
        ForceSec_Beams = np.array(ForceSec_Beams)
        ForceSec01_Beams = ForceSec_Beams[:, rangBeams01]
        ForceSec06_Beams = ForceSec_Beams[:, rangBeams06]
        np.savetxt(data_dir + '_ForceSec01_Beams.out', ForceSec01_Beams)
        np.savetxt(data_dir + '_ForceSec06_Beams.out', ForceSec06_Beams)

        DefoSec_Cols = np.loadtxt(data_dir + '_DefoSec_Cols.out')
        DefoSec_Cols = np.array(DefoSec_Cols)
        DefoSec01_Cols = DefoSec_Cols[:, rangCols01]
        DefoSec06_Cols = DefoSec_Cols[:, rangCols06]
        np.savetxt(data_dir + '_DefoSec01_Cols.out', DefoSec01_Cols)
        np.savetxt(data_dir + '_DefoSec06_Cols.out', DefoSec06_Cols)

        ForceSec_Cols = np.loadtxt(data_dir + '_ForceSec_Cols.out')
        ForceSec_Cols = np.array(ForceSec_Cols)
        ForceSec01_Cols = ForceSec_Cols[:, rangCols01]
        ForceSec06_Cols = ForceSec_Cols[:, rangCols06]
        np.savetxt(data_dir + '_ForceSec01_Cols.out', ForceSec01_Cols)
        np.savetxt(data_dir + '_ForceSec06_Cols.out', ForceSec06_Cols)

        HoriNodes = np.loadtxt(data_dir + '_HoriNodes.out')
        HoriNodes = np.array(HoriNodes)
        HoriNodes = HoriNodes[:, rangDriftNode]
        np.savetxt(data_dir + '_HoriNodes.out', HoriNodes)

        VertNodes = np.loadtxt(data_dir + '_VertNodes.out')
        VertNodes = np.array(VertNodes)
        VertNodes = VertNodes[:, rangDriftNode]
        np.savetxt(data_dir + '_VertNodes.out', VertNodes)

    Sa_max = np.append(Sa_max, maxVB / Wtotal)
    RDR_max = np.append(RDR_max, mDT)
    SDR_max = np.append(SDR_max, mDPB)
    EQnamev = np.append(EQnamev, EQname)
    RA_max = np.append(RA_max, maxRA)
    VnVu_max = np.append(VnVu_max, maxVnVu)
    Testv = np.append(Testv, Test)
    controlTimev = np.append(controlTimev, controlTime)
    Tmaxv = np.append(Tmaxv, Tmax)
    maxPhRot_Colv = np.vstack((maxPhRot_Colv, maxPhRot_Col))
    maxPhRot_Beamv = np.vstack((maxPhRot_Beamv, maxPhRot_Beam))
    maxPhRot_Colcv = np.vstack((maxPhRot_Colcv, maxPhRot_Colc))
    maxPhRot_Beamcv = np.vstack((maxPhRot_Beamcv, maxPhRot_Beamc))
    MedPhRot_Colmv_v = np.vstack((MedPhRot_Colmv_v, MedPhRot_Colm_v))
    MedPhRot_Beammv_v = np.vstack((MedPhRot_Beammv_v, MedPhRot_Beamm_v))
    maxSDRBdgv = np.vstack((maxSDRBdgv, maxSDRBdg))
    op.wipe
maxPhRot_Colv = np.delete(maxPhRot_Colv, 0, axis=0)
maxPhRot_Beamv = np.delete(maxPhRot_Beamv, 0, axis=0)
maxPhRot_Colcv = np.delete(maxPhRot_Colcv, 0, axis=0)
maxPhRot_Beamcv = np.delete(maxPhRot_Beamcv, 0, axis=0)
MedPhRot_Colmv_v = np.delete(MedPhRot_Colmv_v, 0, axis=0)
MedPhRot_Beammv_v = np.delete(MedPhRot_Beammv_v, 0, axis=0)
maxSDRBdgv = np.delete(maxSDRBdgv, 0, axis=0)
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
