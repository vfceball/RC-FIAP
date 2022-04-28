def runNRHA_CSS(Dt, Tmax, Loc_heigth, ListNodesDrift, ListNodesBasal, EleCol, LC, DataColPhl, DataColDesing,
                ListNodesLC):
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
    op.wipeAnalysis()
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
    VBasal_v = []
    DriftTecho_v = []
    maxDriftPiso_v = []
    maxRA_v = []
    maxVu_Vn_v = []

    # Run the actual analysis now
    while cIndex == 0 and controlTime <= Tmax and ok == 0:
        # Runs while the building is stable, time is less
        # than that of the length of the record (plus buffering)
        # and the analysis is still converging

        # Do the analysis
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
        maxDriftPiso_v = np.append(maxDriftPiso_v, maxDriftPiso)
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
        for (Ele, DCPhl, DC) in zip(EleCol, DataColPhl, DataColDesing):
            DeforsS1 = np.array(op.eleResponse(Ele.EleTag, 'section', 1, 'deformation'))
            DeforsS6 = np.array(op.eleResponse(Ele.EleTag, 'section', 6, 'deformation'))
            fi_S1, fi_S6 = DeforsS1[1], DeforsS6[1]
            RA = DCPhl.phl1 * max(map(abs, [fi_S1, fi_S6]))
            if RA >= maxRA:
                maxRA = RA
            ForcesS1 = np.array(op.eleResponse(Ele.EleTag, 'force'))
            # print('ForcesS1', ForcesS1)
            Vu = abs(ForcesS1[0])
            # print('Vu', Vu)
            Vu_Vn = Vu/DC.Vn
            if Vu_Vn >= maxVu_Vn:
                maxVu_Vn = Vu_Vn
        maxRA_v = np.append(maxRA_v, maxRA)
        maxVu_Vn_v = np.append(maxVu_Vn_v, maxVu_Vn)
            # if Element.EleTag <
            # ElemsDeforS1 = np.append(ElemsDeforS1, DeforsS1)

        #     ElemsDeforS1 = np.append(ElemsDeforS1, DeforsS1)
        #     ElemsDeforS6 = np.append(ElemsDeforS6, DeforsS6)
        # MP_ElemsDeforS1 = np.vstack((MP_ElemsDeforS1, ElemsDeforS1))
        # MP_ElemsDeforS6 = np.vstack((MP_ElemsDeforS6, ElemsDeforS6))

    # DriftTecho_v = np.array(DriftTecho_v)
    maxDriftTecho = np.abs(DriftTecho_v).max()
    maxDriftPisoBdg = np.abs(maxDriftPiso_v).max()
    maxVBasal = np.abs(VBasal_v).max()
    maxRABdg = np.abs(maxRA_v).max()
    maxVuVnBdg = np.abs(maxVu_Vn_v).max()

    md = maxDriftPisoBdg
    print('maxDriftTecho={0:.3f} maxDriftPisoBdg={1:.3f} maxVBasal={2:.3f}'.format(maxDriftTecho, maxDriftPisoBdg,
                                                                                   maxVBasal))

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

    return cIndex, maxDriftTecho, maxDriftPisoBdg, maxVBasal, maxRABdg, maxVuVnBdg, Test, controlTime, Tmax
