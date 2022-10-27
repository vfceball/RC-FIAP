def runNRHA_CSS(Dt, Tmax, Loc_heigth, Loc_span, ListNodesDrift, ListNodesBasal, EleCol, EleBeam, LC, DataColPhl,
                DataColDesing, ListNodesLC, EQname, HL_directory, ListNodes, DataBeamPhl, DataBeamDesing):
    # --------------------------------------------------
    # Description of Parameters
    # --------------------------------------------------
    # Dt:		Analysis time step
    # Tmax:	Length of the record (including padding of 0's)
    # Dc:		Deformation capacity (radians or metres)
    # --------------------------------------------------

    # Make the control index available outside the proc
    # Define the initial Analysis parameters
    import openseespy.opensees as op
    import numpy as np  # load the numpy module, calling it np
    print("Starting runNRHA")

    op.constraints('Transformation')
    op.numberer('RCM')
    op.system('BandGeneral')

    testType = 'NormDispIncr'
    tolInit = 1.0e-7  # the initial Tolerance, so it can be referred back to
    iterInit = 50  # the initial Max Number of Iterations
    algorithmType = 'KrylovNewton'  # the algorithm type

    op.test(testType, tolInit, iterInit)  # determine if convergence has been achieved at the end of an iteration step
    op.algorithm(algorithmType)  # use Newton solution algorithm: updates tangent stiffness at every iteration
    NewmarkGamma = 0.5
    NewmarkBeta = 0.25
    op.integrator('Newmark', NewmarkGamma, NewmarkBeta)
    op.analysis('Transient')

    # Set up analysis parameters
    Nsteps = int(Tmax / Dt)  # integer of the number of steps needed
    cIndex = 0  # Initially define the control index (-1 for non-converged, 0 for stable, 1 for global
    # collapse, 2 for local collapse)
    controlTime = 0.0  # Start the controlTime
    ok = 0  # Set the convergence to 0 (initially converged)
    md = 0.0  # Set the initial storey drift
    mflr = 0  # Set the initial storey collapse location

    # Set up the storey drift and acceleration values
    maxVu_Vn = 0.0
    htot = Loc_heigth[-1]
    floors_num = len(Loc_heigth) - 1
    axes_num = len(Loc_span)

    VBasal_v = []
    DriftTecho_v = []
    maxDriftPiso_v = []
    maxRA_v = []
    maxRPD_v = []
    maxVu_Vn_v = []
    PhRot_Col_v = np.zeros(floors_num)
    PhRot_Colm_v = np.zeros((1, floors_num, 2*axes_num))
    PhRot_Colc_v = np.zeros(floors_num)
    PhRot_Beam_v = np.zeros(floors_num)
    PhRot_Beamm_v = np.zeros((1, floors_num, 2*(axes_num-1)))
    PhRot_Beamc_v = np.zeros(floors_num)
    maxSDR_v = np.zeros(floors_num)

    # if self.ui.checkBoxSaveCSS.isChecked() == True:
    #     data_dir = EQname.replace('gmotions', 'Data')
    #     ForceSec01_Beams, DefoSec01_Beams, ForceSec06_Beams, DefoSec06_Beams = np.zeros(2 * nele), np.zeros(2 * nele), \
    #                                                                            np.zeros(2 * nele), np.zeros(2 * nele)
    #     ForceSec01_Cols, DefoSec01_Cols, ForceSec06_Cols, DefoSec06_Cols = np.zeros(2 * nele), np.zeros(2 * nele), \
    #                                                                        np.zeros(2 * nele), np.zeros(2 * nele)

    # Run the actual analysis now
    while cIndex == 0 and controlTime <= Tmax and ok == 0:
        # Runs while the building is stable, time is less
        # than that of the length of the record (plus buffering)
        # and the analysis is still converging

        # Do the analysis
        vForceSec01_Beams, vDefoSec01_Beams, vForceSec06_Beams, vDefoSec06_Beams = [], [], [], []
        ok = op.analyze(1, Dt)
        controlTime = op.getTime()  # Update the control time

        # First changes are to change algorithm to achieve convergence...
        if ok != 0:
            print('~~~ Failed at ', controlTime, ' - Reduced timestep by half......')
            Dtt = 0.5 * Dt
            ok = op.analyze(1, Dtt)
        if ok != 0:
            print('~~~ Failed at ', controlTime, ' - Reduced timestep by quarter......')
            Dtt = 0.25 * Dt
            ok = op.analyze(1, Dtt)
        if ok != 0:
            print('~~~ Failed at ', controlTime, ' - Trying Broyden......')
            op.algorithm('Broyden', 8)
            ok = op.analyze(1, Dt)
            op.algorithm(algorithmType)
        if ok != 0:
            print('~~~ Failed at ', controlTime, ' - Trying Newton with Initial Tangent ......')
            op.algorithm('Newton', '-initial')
            ok = op.analyze(1, Dt)
            op.algorithm(algorithmType)
        if ok != 0:
            print('~~~ Failed at ', controlTime, ' - Trying NewtonWithLineSearch...... ......')
            op.algorithm('NewtonLineSearch', 0.8)
            ok = op.analyze(1, Dt)
            op.algorithm(algorithmType)
        if ok != 0:
            print('~~~ Failed at ', controlTime, ' - Trying Newton with Initial Tangent & relaxed convergence......')
            op.test(testType, tolInit * 0.1, iterInit * 50)
            op.algorithm('Newton', '-initial')
            ok = op.analyze(1, Dt)
            op.test(testType, tolInit, iterInit)
            op.algorithm(algorithmType)
        if ok != 0:
            print('~~~ Failed at ', controlTime, ' - Trying Newton with Initial Tangent & relaxed convergence......')
            op.test(testType, tolInit * 0.1, iterInit * 50)
            op.algorithm('Newton', '-initial')
            ok = op.analyze(1, Dt)
            op.test(testType, tolInit, iterInit)
            op.algorithm(algorithmType)
        if ok != 0:
            print('~~~ Failed at ', controlTime, ' - Trying Newton with Initial Tangent & relaxed convergence......')
            op.test(testType, tolInit * 0.1, iterInit * 50)
            op.algorithm('Newton', '-initial')
            ok = op.analyze(1, Dt)
            op.test(testType, tolInit, iterInit)
            op.algorithm(algorithmType)
        if ok != 0:
            print('~~~ Failed at ', controlTime, ' - Trying NewtonWithLineSearch & relaxed convergence......')
            op.test(testType, tolInit * 0.1, iterInit * 50)
            op.algorithm('NewtonLineSearch', '-initial')
            ok = op.analyze(1, Dt)
            op.test(testType, tolInit, iterInit)
            op.algorithm(algorithmType)
        if ok != 0:
            print('~~~ Failed at ', controlTime,
                  ' - Trying Newton with Initial Tangent, reduced timestep & relaxed convergence......')
            op.test(testType, tolInit * 0.1, iterInit * 50)
            op.algorithm('Newton', '-initial')
            Dtt = 0.5 * Dt
            ok = op.analyze(1, Dtt)
            op.test(testType, tolInit, iterInit)
            op.algorithm(algorithmType)
        if ok != 0:
            print('~~~ Failed at ', controlTime,
                  ' - Trying Newton with Initial Tangent, reduced timestep & relaxed convergence......')
            op.test(testType, tolInit * 0.1, iterInit * 50)
            op.algorithm('Newton', '-initial')
            Dtt = 0.5 * Dt
            ok = op.analyze(1, Dtt)
            op.test(testType, tolInit, iterInit)
            op.algorithm(algorithmType)
        if ok != 0:
            print('~~~ Failed at ', controlTime,
                  ' - Trying NewtonWithLineSearch, reduced timestep & relaxed convergence......')
            op.test(testType, tolInit * 0.1, iterInit * 50)
            op.algorithm('NewtonLineSearch', 0.8)
            Dtt = 0.5 * Dt
            ok = op.analyze(1, Dtt)
            op.test(testType, tolInit, iterInit)
            op.algorithm(algorithmType)
        # Game over......
        if ok != 0:
            print('~~~ Failed at ', controlTime, ' - exit analysis......')
            cIndex = -1

        # Check the actual state of the model with respect to the limits provided

        # Calculation of maximum Drift between floors
        maxDriftPiso = 0.0
        nfloor = 0
        maxSDR = np.zeros(floors_num)
        # print('maxSDR', maxSDR)
        for (nod_ini, nod_end) in zip(ListNodesDrift[:-1, 0], ListNodesDrift[1:, 0]):
            # print('nod_ini ', nod_ini, 'nod_end', nod_end)
            nod_ini = int(nod_ini)
            nod_end = int(nod_end)
            pos_i = op.nodeCoord(nod_ini, 2)
            pos_s = op.nodeCoord(nod_end, 2)
            hpiso = pos_s - pos_i
            desp_i = op.nodeDisp(nod_ini, 1)
            desp_s = op.nodeDisp(nod_end, 1)
            desp_piso = abs(desp_s - desp_i)
            drift_piso = desp_piso / hpiso
            if drift_piso >= maxDriftPiso:
                maxDriftPiso = drift_piso
            maxSDR[nfloor] = drift_piso
            nfloor += 1
            # print('maxSDR', maxSDR)
            # print('nfloor', nfloor)
        maxDriftPiso_v = np.append(maxDriftPiso_v, maxDriftPiso)
        maxSDR_v = np.vstack((maxSDR_v, maxSDR))

        VBasal = 0.
        op.reactions()
        for node in ListNodesBasal:
            # print('ind Basal ', node[0])
            VBasal = VBasal + op.nodeReaction(node[0], 1)
        if LC == 0:
            VBasal = VBasal + op.nodeReaction(int(ListNodesLC[0, 0]), 1)
        ctrlNode = int(ListNodesDrift[-1, 0])
        VBasal_v = np.append(VBasal_v, VBasal)
        DriftTecho = op.nodeDisp(ctrlNode, 1) / htot
        DriftTecho_v = np.append(DriftTecho_v, DriftTecho)

        maxRA = 0.0  # Max Rotation Angle
        maxRPD = 0.0
        PhRot_Col = np.zeros(floors_num)
        PhRot_Colc = np.zeros(floors_num)
        PhRot_Colm = np.zeros((floors_num, 2*axes_num))
        nfloor = 1
        naxe = 1
        for (Ele, DCPhl, DC) in zip(EleCol, DataColPhl, DataColDesing):
            DeforsS1 = np.array(op.eleResponse(Ele.EleTag, 'section', 1, 'deformation'))
            DeforsS6 = np.array(op.eleResponse(Ele.EleTag, 'section', 6, 'deformation'))
            fi_S1, fi_S6 = DeforsS1[1], DeforsS6[1]
            RA = DCPhl.phl1 * max(map(abs, [fi_S1, fi_S6]))
            PD = op.eleResponse(Ele.EleTag, 'plasticDeformation')
            PD1 = abs(PD[1])
            PD2 = abs(PD[2])
            RPD = max(PD1, PD2)
            if RA >= maxRA:
                maxRA = RA
            if RPD >= maxRPD:
                maxRPD = RPD
            ForcesS1 = np.array(op.eleResponse(Ele.EleTag, 'force'))
            # print('ForcesS1', ForcesS1)
            Vu = abs(ForcesS1[0])
            # print('Vu', Vu)
            Vu_Vn = Vu/DC.Vn
            if Vu_Vn >= maxVu_Vn:
                maxVu_Vn = Vu_Vn
            if ListNodes[Ele.Nod_end, 2] == Loc_heigth[nfloor]:
                if PD1 > PhRot_Colm[nfloor-1, 2*naxe-2]:
                    PhRot_Colm[nfloor-1, 2*naxe-2] = PD1
                if PD2 > PhRot_Colm[nfloor-1, 2*naxe-1]:
                    PhRot_Colm[nfloor-1, 2*naxe-1] = PD2
                naxe += 1
                if RA > PhRot_Col[nfloor-1]:
                    PhRot_Col[nfloor-1] = RA
                if RPD > PhRot_Colc[nfloor-1]:
                    PhRot_Colc[nfloor-1] = RPD
            else:
                naxe = 1
                nfloor += 1
                if PD1 > PhRot_Colm[nfloor-1, 2*naxe-2]:
                    PhRot_Colm[nfloor-1, 2*naxe-2] = PD1
                if PD2 > PhRot_Colm[nfloor-1, 2*naxe-1]:
                    PhRot_Colm[nfloor-1, 2*naxe-1] = PD2
                naxe += 1
                if RA > PhRot_Col[nfloor-1]:
                    PhRot_Col[nfloor-1] = RA
                if RPD > PhRot_Colc[nfloor-1]:
                    PhRot_Colc[nfloor-1] = RPD
        maxRA_v = np.append(maxRA_v, maxRA)
        maxRPD_v = np.append(maxRPD_v, maxRPD)
        maxVu_Vn_v = np.append(maxVu_Vn_v, maxVu_Vn)
        PhRot_Col_v = np.vstack((PhRot_Col_v, PhRot_Col))
        PhRot_Colc_v = np.vstack((PhRot_Colc_v, PhRot_Colc))
        PhRot_Colm_v = np.concatenate((PhRot_Colm_v, [PhRot_Colm]), axis=0)
        maxRA = 0.0  # Max Rotation Angle
        PhRot_Beam = np.zeros(floors_num)
        PhRot_Beamc = np.zeros(floors_num)
        PhRot_Beamm = np.zeros((floors_num, 2*(axes_num-1)))
        nfloor = 1
        naxe = 1
        for (Ele, DBPhl, DC) in zip(EleBeam, DataBeamPhl, DataBeamDesing):
            DeforsS1 = np.array(op.eleResponse(Ele.EleTag, 'section', 1, 'deformation'))
            DeforsS6 = np.array(op.eleResponse(Ele.EleTag, 'section', 6, 'deformation'))
            fi_S1, fi_S6 = DeforsS1[1], DeforsS6[1]
            RA1 = DBPhl.phl1 * abs(fi_S1)
            RA6 = DBPhl.phl2 * abs(fi_S6)
            RA = max([RA1, RA6])
            PD = op.eleResponse(Ele.EleTag, 'plasticDeformation')
            PD1 = abs(PD[1])
            PD2 = abs(PD[2])
            RPD = max(PD1, PD2)
            if ListNodes[Ele.Nod_end, 2] == Loc_heigth[nfloor]:
                if PD1 > PhRot_Beamm[nfloor-1, 2*naxe-2]:
                    PhRot_Beamm[nfloor-1, 2*naxe-2] = PD1
                if PD2 > PhRot_Beamm[nfloor-1, 2*naxe-1]:
                    PhRot_Beamm[nfloor-1, 2*naxe-1] = PD2
                naxe += 1
                if RA > PhRot_Beam[nfloor-1]:
                    PhRot_Beam[nfloor - 1] = RA
                if RPD > PhRot_Beamc[nfloor-1]:
                    PhRot_Beamc[nfloor - 1] = RPD
            else:
                naxe = 1
                nfloor += 1
                if PD1 > PhRot_Beamm[nfloor-1, 2*naxe-2]:
                    PhRot_Beamm[nfloor-1, 2*naxe-2] = PD1
                if PD2 > PhRot_Beamm[nfloor-1, 2*naxe-1]:
                    PhRot_Beamm[nfloor-1, 2*naxe-1] = PD2
                naxe += 1
                if RA > PhRot_Beam[nfloor-1]:
                    PhRot_Beam[nfloor-1] = RA
                if RPD > PhRot_Beamc[nfloor-1]:
                    PhRot_Beamc[nfloor-1] = RPD
        PhRot_Beam_v = np.vstack((PhRot_Beam_v, PhRot_Beam))
        PhRot_Beamc_v = np.vstack((PhRot_Beamc_v, PhRot_Beamc))
        PhRot_Beamm_v = np.concatenate((PhRot_Beamm_v, [PhRot_Beamm]), axis=0)
    maxDriftTecho = np.abs(DriftTecho_v).max()
    maxDriftPisoBdg = np.abs(maxDriftPiso_v).max()
    maxVBasal = np.abs(VBasal_v).max()
    maxRABdg = np.abs(maxRA_v).max()
    maxVuVnBdg = np.abs(maxVu_Vn_v).max()
    maxPhRot_Col = np.max(PhRot_Col_v, axis=0)
    maxPhRot_Beam = np.max(PhRot_Beam_v, axis=0)
    maxPhRot_Colc = np.max(PhRot_Colc_v, axis=0)
    # print('max =', maxPhRot_Colc)
    maxPhRot_Beamc = np.max(PhRot_Beamc_v, axis=0)
    maxSDRBdg = np.max(maxSDR_v, axis=0)
    MaxPhRot_Colm_v = np.max(PhRot_Colm_v, axis=0)
    MedPhRot_Colm_v = np.median(MaxPhRot_Colm_v, axis=1)
    # prueba = np.max(MaxPhRot_Colm_v, axis=1)
    # print('max1 =', prueba)
    # print('med =', MedPhRot_Colm_v)
    MaxPhRot_Beamm_v = np.max(PhRot_Beamm_v, axis=0)
    MedPhRot_Beamm_v = np.median(MaxPhRot_Beamm_v, axis=1)
    md = maxDriftPisoBdg
    print('maxDriftTecho={0:.3f} maxDriftPisoBdg={1:.3f} maxVBasal={2:.3f}'.format(maxDriftTecho, maxDriftPisoBdg,
                                                                                   maxVBasal))
    # print('maxPhRot_Col', maxPhRot_Col)

    # Print to the max interstorey drifts

    print('PeakDemand:{0:.5f}'.format(md))  # Print to the max demand
    if cIndex == -1:
        Test = ':::::: ANALYSIS FAILED TO CONVERGE at ', controlTime, ' of ', Tmax, ' :::::', 'Max Drift Piso=', maxDriftPisoBdg
    if cIndex == 0:
        Test = '######## ANALYSIS COMPLETED SUCCESSFULLY #####', 'Max Drift Piso=', maxDriftPisoBdg
    if cIndex == 1:
        print('========== LOCAL STRUCTURE COLLAPSE ==========')
    print(Test)
    Test = np.array([controlTime, Tmax, maxDriftPisoBdg*100])
    Test = np.around(Test, decimals=2)
    # Test = "".join(map(str, Test))

    return cIndex, maxDriftTecho, maxDriftPisoBdg, maxVBasal, maxRABdg, maxVuVnBdg, Test, controlTime, Tmax,\
           maxPhRot_Col, maxPhRot_Beam, maxSDRBdg, maxPhRot_Colc, maxPhRot_Beamc, MedPhRot_Colm_v, MedPhRot_Beamm_v
