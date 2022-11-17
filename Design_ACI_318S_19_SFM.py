global Loc_span, Loc_heigth, ListNodes, Elements, DataBeamDesing, DataColDesing, WDL, WLL, WDLS, Wtotal, \
    cover, num_elems, Beta1B, Beta1C, fcB, fcC, ListNodesDrift, ListNodesBasal, Ta, num_beams, num_cols
# Function: Reads Beams design data from table that allows the user to modify the default design from TAB2 of GU
def data_beams_table(self):
    self.registros_beams = []

    for DB in DataBeamDesing:
        b = DB.b / cm
        h = DB.h / cm
        L_As_top = DB.Ast1 / cm ** 2
        L_As_bot = DB.Asb1 / cm ** 2
        R_As_top = DB.Ast2 / cm ** 2
        R_As_bot = DB.Asb2 / cm ** 2
        L_Leg_n = DB.ns1
        R_Leg_n = DB.ns2
        L_Sstirrup = DB.ss1 / cm
        R_Sstirrup = DB.ss2 / cm
        registro = RegistroBeams(self.ui.tbl_data_design_beams, DB.EleTag, b, h, L_As_top, L_As_bot, L_Leg_n,
                                 L_Sstirrup, R_As_top, R_As_bot, R_Leg_n, R_Sstirrup)
        self.registros_beams.append(registro)

# Function: Reads Columns design data from table that allows the user to modify the default design from TAB2 of GUI.
def data_columns_table(self):
    self.registros_cols = []
    for DC in DataColDesing:
        b = DC.b / cm
        h = DC.h / cm
        roCol = DC.ro
        db = DC.db / mm
        de = DC.de / mm
        nbH = DC.nbH
        nbB = DC.nbB
        nsH = DC.nsH
        nsB = DC.nsB
        sst = DC.sst / cm
        registro = RegistroColumns(self.ui.tbl_data_design_columns, DC.EleTag, b, h, roCol, db, de, nbH, nbB, nsH, nsB,
                                   sst)
        self.registros_cols.append(registro)

def beta1(fc):
    if fc <= 28 * MPa:
        Beta1 = 0.85
    else:
        Beta1 = max([0.85 - 0.05 * (fc - 28.) / 7., 0.65])
    return Beta1

# Design load combinations
def Combo_ACI(DL, LL, E):
    U1 = 1.2 * DL + 1.6 * LL
    U2 = 1.2 * DL + 1.0 * LL + 1.0 * E
    U3 = 1.2 * DL + 1.0 * LL - 1.0 * E
    U4 = 0.9 * DL + 1.0 * E
    U5 = 0.9 * DL - 1.0 * E
    return U1, U2, U3, U4, U5

# Flexural beams design
def AsBeam(Mu, EleTag, cover, ro_min_b, ro_max_b, dst, fy, BBeam, HBeam):
    b, h = BBeam, HBeam
    Mu = abs(Mu)
    db_v = np.array([4, 5, 6, 7, 8, 10])
    for ndb in db_v:
        db = ndb / 8. * inch
        d = h - cover - dst - 0.5 * db
        if Mu == 0.0:
            ro_req = ro_min_b
        else:
            ro_req = 0.85 * fcB / fy * (1. - sqrt(1. - 2. * (Mu / 0.9 / b / d ** 2) / 0.85 / fcB))
        if ro_req < ro_min_b:
            ro_req = ro_min_b
        As_req = ro_req * b * d
        Ab = pi * db ** 2 / 4.
        nb = max(2., ceil(As_req / Ab))
        As_con = nb * Ab
        slb = (b - 2 * cover - 2 * dst - nb * db) / (nb - 1.)  # free clear bars
        if slb >= max(1. * inch, db):
            break
        if ro_req > ro_max_b:
            print("Steel percentage greater than the maximum in Beam " + str(EleTag))
    if slb < min(1. * inch, db):
        print("Bar separation is not ok in Beam " + str(EleTag))
    a = fy * As_con / 0.85 / fcB / b
    Mn = fy * As_con * (d - a / 2.)
    apr = 1.25 * fy * As_con / 0.85 / fcB / b
    Mpr = 1.25 * fy * As_con * (d - apr / 2.)
    return As_con, d, Mn, db, Mpr, nb

# Shear beams design
def AvBeam(Vu, db, d, EleTag, fys, dst, Ast, BBeam, nb, Vc):
    Vs = (Vu - 0.75 * Vc) / 0.75
    if Vs > 4. * Vc:
        print("reshape by shear in Beam " + str(EleTag))
    se_1 = min(d / 4., 6. * db, 150. * mm)
    nr_v = np.array([2, 3, 4, 5, 6])  # vector de numero de ramas
    if Vs <= 0.:
        se = se_1
        nra = 2.
    else:
        for nra in nr_v:
            Ave = Ast * nra  # area transversal del estribo
            se_2 = Ave * fys * d / Vs
            se = min(se_1, se_2)
            if se >= 60. * mm:
                break
    se = floor(se / cm) * cm
    if se < 60. * mm:
        print("Stirrup spacing is less than 6 cm in beam " + str(EleTag))
    return nra, se

# Colmuns P-M design
def AsColumn(b, h, EleTag, cover, dst, fy, Beta1C, Pu_v, Mu_v, Sum_Mn_B, FactorColBeamStr, ncolsn, Nu_min):
    ro_min = 0.01
    ro_max = 0.06
    npts = 20
    ncom = 10
    ecu = 0.003
    Es = 200. * GPa
    verif = False
    nbB = ceil(b * 10)  # bars numbers along B
    nbH = ceil(h * 10)  # bars numbers along H
    D_c = 1.1 * h / npts
    nbH_v = np.array([nbH - 1, nbH, nbH + 1])
    nbB_v = np.array([nbB - 1, nbB, nbB + 1])
    db_v = np.array([4, 5, 6, 7, 8, 9, 10, 11, 14, 18])  # vector bar diameters
    while verif == False:
        for ndb in db_v:
            db = ndb / 8. * inch
            Ab = pi * db ** 2. / 4.
            dp = cover + dst + 0.5 * db
            d = h - dp
            for nbH in nbH_v:
                for nbB in nbB_v:
                    nbT = 2. * (nbB + nbH - 2.)  # numero total de barras
                    Ast = nbT * Ab
                    ro = Ast / b / h
                    As = np.hstack([nbB * Ab, np.ones(nbH - 2) * 2 * Ab, nbB * Ab])
                    dist = np.linspace(dp, h - dp, nbH)
                    if ro >= ro_min:
                        Pn_max = 0.80 * (0.85 * fcC * (b * h - Ast) + fy * Ast)
                        Tn_max = -fy * Ast
                        c = np.linspace(1.1 * h / npts, 1.1 * h, npts)
                        a = Beta1C * c
                        Pconc = 0.85 * fcC * a * b
                        Mconc = Pconc * (h - a) / 2.
                        et = ecu * (d - c) / c
                        fiv = np.copy(et)
                        fiv = np.where(fiv >= 0.005, 0.9, fiv)
                        fiv = np.where(fiv <= 0.002, 0.65, fiv)
                        fiv = np.where((fiv > 0.002) & (fiv < 0.005), (0.65 + 0.25 * (fiv - 0.002) / 0.003),
                                       fiv)
                        c = c[:, np.newaxis]
                        es = ecu * (c - dist) / c
                        fs = Es * es
                        fs = np.where(fs > fy, fy, fs)
                        fs = np.where(fs < -fy, -fy, fs)
                        Pacer = np.sum(fs * As, axis=1)
                        Macer = np.sum(fs * As * (h / 2. - dist), axis=1)
                        Pn = np.hstack(
                            [Tn_max, np.where(Pconc + Pacer > Pn_max, Pn_max, Pconc + Pacer), Pn_max])
                        Mn = np.hstack([0, Mconc + Macer, 0])
                        fiv = np.hstack([0.9, fiv, 0.65])
                        fiPn = fiv * Pn
                        fiMn = fiv * Mn

                        Pn_max_pr = 0.80 * (0.85 * fcC * (b * h - Ast) + 1.25*fy * Ast)
                        Tn_max_pr = -1.25*fy * Ast
                        fspr = Es * es
                        fspr = np.where(fspr > 1.25*fy, 1.25*fy, fspr)
                        fspr = np.where(fspr < -1.25*fy, -1.25*fy, fspr)
                        Pacer_pr = np.sum(fspr * As, axis=1)
                        Macer_pr = np.sum(fspr * As * (h / 2. - dist), axis=1)
                        Ppr = np.hstack(
                            [Tn_max_pr, np.where(Pconc + Pacer_pr > Pn_max_pr, Pn_max_pr, Pconc + Pacer_pr), Pn_max_pr])
                        Mpr = np.hstack([0, Mconc + Macer_pr, 0])
                        if np.all((Pu_v >= min(fiPn)) & (Pu_v <= max(fiPn))):
                            Mu_i = np.interp(Pu_v, fiPn, fiMn)
                            Mn_i = np.interp(Pu_v, Pn, Mn)
                            Mpr_i = np.interp(Pu_v, Ppr, Mpr)
                            Mns = np.interp(Nu_min, Pn, Mn)
                            Col_to_beam_str_ratio = ncolsn * Mns / Sum_Mn_B
                            if ncolsn == 1:
                                if np.all(Mu_i >= Mu_v) == True:
                                    verif = True
                                    break
                            else:
                                if np.all(Mu_i >= Mu_v) == True and Col_to_beam_str_ratio >= FactorColBeamStr:
                                    verif = True
                                    break

                if verif == True:
                    break
            if verif == True:
                break
        if ndb == db_v[-1] and ro > ro_max:
            print('column ' + str(EleTag) + 'needs to be resized by reinforcement ratio')
            break
    return nbH, nbB, db, As, fiPn, fiMn, Mn_i, Mpr_i, d, dist, ro, Mu_i, Col_to_beam_str_ratio

# Shear columns design
def AvColumn(EleTag, Vu, b, h, nbH, nbB, dst, Vc, db, fys, Nu_min):
    fiv = 0.75
    Ag = b * h
    dp = cover + dst + db / 2
    d = h - dp
    neH_v = np.arange(floor(nbH / 2) + 1,nbH+1)
    neB_v = np.arange(floor(nbB / 2) + 1,nbB+1)
    bc1, bc2 = h - 2*cover, b - 2*cover
    de_v = np.array([3, 4])
    for nde in de_v:
        de = nde / 8. * inch
        Ast = pi * de ** 2. / 4.
        for neB in neB_v:
            for neH in neH_v:
                Ash_H = neH * Ast
                Ash_B = neB * Ast
                hx = (h - 2 * dp) / neH
                so = max(100 * mm, min(100 * mm + (350 * mm - hx) / 3, 150 * mm))
                se_1 = min(6. * db, b / 4., h / 4., so)  # minimum spacing c.18.7.5.3 ACI-19
                if Nu_min <= 0.3*b*h*fcC:
                    se_2 = min(Ash_H / bc1 / (0.3 * (Ag / (bc1 * bc2) - 1) * fcC / fys),
                               Ash_B / bc2 / (0.3 * (Ag / (bc1 * bc2) - 1) * fcC / fys),
                               Ash_H / bc1 / (0.09 * fcC / fys),
                               Ash_B / bc1 / (0.09 * fcC / fys))
                elif Nu_min > 0.3*b*h*fcC:
                    kf = max(fcC/(175*MPa)+0.6, 1.0)
                    nl = neH*2+(neB-2)*2
                    kn = nl/(nl-2)
                    se_2 = min(Ash_H / bc1 / (0.3 * (Ag / (bc1 * bc2) - 1) * fcC / fys),
                               Ash_B / bc2 / (0.3 * (Ag / (bc1 * bc2) - 1) * fcC / fys),
                               Ash_H / bc1 / (0.09 * fcC / fys),
                               Ash_B / bc1 / (0.09 * fcC / fys),
                               Ash_H/bc1/((0.2*kf*kn*Nu_min)/(fys*bc1*bc2)),
                               Ash_B/bc2/((0.2*kf*kn*Nu_min)/(fys*bc1*bc2)))
                se_1 = min(se_1, se_2)
                Vs = (Vu - fiv * Vc) / fiv
                if Vs <= 1 / 3 * sqrt(fcC * MPa) * b * d:
                    se_1 = se_1
                elif Vs >= 1 / 3 * sqrt(fcC * MPa) * b * d:
                    se_1 = min(se_1, h / 4)
                if Vs > 0.66 * sqrt(fcC * MPa) * b * d:
                    print('Resize the column' + str(EleTag) + ' by shear ')
                Ave = Ash_B  # area transversal del estribo
                if Vs <= 0.:
                    se = se_1
                else:
                    se_2 = Ave * fys * d / Vs
                    se = min([se_1, se_2])
                if se >= 80.*mm:
                    break
                if se < 60. * mm:
                    print('Minimum spacing of stirrups is not met in column ' + str(EleTag))
            if se >= 80.*mm:
                break
        if se >= 80.*mm:
            break
    Vn = Vc + Ave*fys*d/se
    return se, neB, neH, Vn, de

# Compression block parameters beta as function f'c

# Input geometric, materials and seismic design parameters from TAB1 of GUI
Lafg = float(self.ui.Lafg.text())
Lafs = float(self.ui.Lafs.text())
DL = float(self.ui.DL.text())
LL = float(self.ui.LL.text())
HColi = float(self.ui.HColi.text())  # Column inside Depth
BColi = float(self.ui.BColi.text())  # Column inside Width
HCole = float(self.ui.HCole.text())  # Column outside Depth
BCole = float(self.ui.BCole.text())  # Column outside Width
HBeam = float(self.ui.HBeam.text())
BBeam = float(self.ui.BBeam.text())
IFC = float(self.ui.InertiaColumnsFactor.text())
IFB = float(self.ui.InertiaBeamsFactor.text())
heigth_v = self.ui.heigth_v.text()
heigth_v = heigth_v.split(',')
heigth_v = np.array(heigth_v, dtype=float)
span_v = self.ui.span_v.text()
span_v = span_v.split(',')
span_v = np.array(span_v, dtype=float)
fy = float(self.ui.fy.text()) * MPa
fys = float(self.ui.fys.text()) * MPa
fcB = float(self.ui.fcB.text()) * MPa
fcC = float(self.ui.fcC.text()) * MPa
WDL = Lafg * DL
WDLS = Lafs * DL
WLL = Lafg * LL

plt.close('all')
op.wipe()
op.model('Basic', '-ndm', 2, '-ndf', 3)

# Nodes Creations
Loc_span = np.append(0, np.cumsum(span_v))
Loc_heigth = np.append(0, np.cumsum(heigth_v))
n_col_axes = len(Loc_span)
xn_v, yn_v = np.meshgrid(Loc_span, Loc_heigth)
xn_vf = np.ravel(xn_v)
yn_vf = np.ravel(yn_v)
num_nodes = len(Loc_span) * len(Loc_heigth)
ListNodes = np.empty([num_nodes, 3])
nodeTag = 0
for (xn, yn) in zip(xn_vf, yn_vf):
    ListNodes[nodeTag, :] = [nodeTag, xn, yn]
    op.node(nodeTag, xn, yn)
    if yn == 0.:
        op.fix(nodeTag, 1, 1, 1)
    nodeTag += 1
for node in ListNodes:
    if node[2] > 0. and node[1] == 0.:
        MasterNode = node[0]
    if node[2] > 0. and node[1] != 0.:
        op.equalDOF(int(MasterNode), int(node[0]), 1)

ListNodesDrift = ListNodes[np.where(ListNodes[:, 1] == 0.)]
ListNodesBasal = ListNodes[np.where(ListNodes[:, 2] == 0.)]
MassType = "-lMass"  # -lMass, -cMass

# Columns creation for elastic analysis
op.geomTransf('Linear', 1, '-jntOffset', 0, 0, 0, -HBeam / 2)
op.geomTransf('Linear', 2, '-jntOffset', 0, HBeam / 2, 0, -HBeam / 2)
AColi = BColi * HColi  # cross-sectional area
ACole = BCole * HCole  # cross-sectional area
EcC = 4700 * sqrt(fcC * MPa)
IzColi = 1. / 12. * BColi * HColi ** 3  # Column moment of inertia
IzCole = 1. / 12. * BCole * HCole ** 3  # Column moment of inertia
EleTag = 1
Elements = []
for Nod_ini in range(num_nodes):
    if ListNodes[Nod_ini, 2] != Loc_heigth[-1]:
        Nod_end = Nod_ini + n_col_axes
        if ListNodes[Nod_ini, 2] == 0.:
            gTr = 1
            RZi = 0
            RZe = HBeam / 2
            LCol = ListNodes[Nod_end, 2] - ListNodes[Nod_ini, 2] - RZi - RZe
        else:
            gTr = 2
            RZi = HBeam / 2
            RZe = HBeam / 2
            LCol = ListNodes[Nod_end, 2] - ListNodes[Nod_ini, 2] - RZi - RZe
        if ListNodes[Nod_ini, 1] == 0. or ListNodes[Nod_ini, 1] == Loc_span[-1]:
            BCol, HCol = BCole, HCole
            ACol = ACole
            IzCol = IFC * IzCole
        else:
            BCol, HCol = BColi, HColi
            ACol = AColi
            IzCol = IFC * IzColi
        MassDens = ACol * GConc / g
        Elements.append(BeamElasticElement(EleTag, Nod_ini, Nod_end, ACol, EcC, IzCol, LCol, BCol, HCol, gTr,
                                           RZi, RZe))
        op.element('elasticBeamColumn', EleTag, Nod_ini, Nod_end, ACol, EcC, IzCol, gTr, '-mass', MassDens,
                   MassType)
        EleTag += 1
num_cols = EleTag

# Beams creation for elastic analysis
op.geomTransf('Linear', 3, '-jntOffset', HColi / 2., 0, -HColi / 2., 0)
op.geomTransf('Linear', 4, '-jntOffset', HCole / 2., 0, -HColi / 2., 0)
op.geomTransf('Linear', 5, '-jntOffset', HColi / 2., 0, -HCole / 2., 0)
ABeam = BBeam * HBeam
EcB = 4700 * sqrt(fcB * MPa)
IzBeam = IFB * BBeam * HBeam ** 3 / 12
MassDens = ABeam * GConc / g + WDLS / g
for Nod_ini in range(num_nodes):
    if ListNodes[Nod_ini, 1] != Loc_span[-1] and ListNodes[Nod_ini, 2] != 0.:
        Nod_end = Nod_ini + 1
        if ListNodes[Nod_ini, 1] == 0.:
            gTr = 4
            RZi = HCole / 2.
            RZe = HColi / 2.
            LBeam = ListNodes[Nod_end, 1] - ListNodes[Nod_ini, 1] - RZi - RZe
        elif ListNodes[Nod_ini, 1] == Loc_span[-2]:
            gTr = 5
            RZi = HColi / 2.
            RZe = HCole / 2.
            LBeam = ListNodes[Nod_end, 1] - ListNodes[Nod_ini, 1] - RZi - RZe
        else:
            gTr = 3
            RZi = HColi / 2.
            RZe = HColi / 2.
            LBeam = ListNodes[Nod_end, 1] - ListNodes[Nod_ini, 1] - RZi - RZe
        Elements.append(BeamElasticElement(EleTag, Nod_ini, Nod_end, ABeam, EcB, IzBeam, LBeam, BBeam, HBeam,
                                           gTr, RZi, RZe))
        op.element('elasticBeamColumn', EleTag, Nod_ini, Nod_end, ABeam, EcB, IzBeam, gTr,
                   '-mass', MassDens, MassType)
        EleTag += 1
num_elems = EleTag
num_beams = num_elems - num_cols

# Create a Plain load pattern for gravity loading with a Linear TimeSeries
Pvig = ABeam * GConc
PColi = AColi * GConc
PCole = ACole * GConc
PColSlabD, PColSlabL = 0, 0
if self.ui.radioButtonSpatial.isChecked() and Lafs > Lafg:
    PColSlabD = (Lafs-Lafg)*WDL+Pvig*Lafs
    PColSlabL = (Lafs-Lafg)*WLL

op.timeSeries('Linear', 1)
op.pattern('Plain', 1, 1)
for Element in Elements:
    if ListNodes[Element.Nod_ini, 1] == ListNodes[Element.Nod_end, 1]:
        if ListNodes[Element.Nod_ini, 1] == 0. or ListNodes[Element.Nod_ini, 1] == Loc_span[-1]:
            PCol = PCole
        else:
            PCol = PColi
        op.eleLoad('-ele', Element.EleTag, '-type', '-beamUniform', 0, -PCol)
        op.eleLoad('-ele', Element.EleTag, '-type', '-beamPoint', 0, 1, -PColSlabD)
    if ListNodes[Element.Nod_ini, 2] == ListNodes[Element.Nod_end, 2]:
        op.eleLoad('-ele', Element.EleTag, '-type', '-beamUniform', -Pvig - WDL)

op.system('UmfPack')
op.numberer('Plain')
op.constraints('Plain')
op.integrator('LoadControl', 1.0)
op.algorithm('Linear')
op.analysis('Static')
op.analyze(1)
ElemnsForceD = []
for Element in Elements:
    Forces = op.eleForce(Element.EleTag)
    Forces.insert(0, Element.EleTag)
    ElemnsForceD.append(Forces)
ElemnsForceD = np.array(ElemnsForceD)
Wtotal = np.sum(ElemnsForceD[:len(Loc_span), 2]) * Lafs/Lafg #debo mirar esto

op.loadConst('-time', 0.0)
op.timeSeries('Linear', 2)
op.pattern('Plain', 2, 1)
for Element in Elements:
    if ListNodes[Element.Nod_ini, 1] == ListNodes[Element.Nod_end, 1]:
        op.eleLoad('-ele', Element.EleTag, '-type', '-beamPoint', 0, 1, -PColSlabL)
    if ListNodes[Element.Nod_ini, 2] == ListNodes[Element.Nod_end, 2]:
        op.eleLoad('-ele', Element.EleTag, '-type', '-beamUniform', -WLL)
op.analyze(1)

ElemnsForceDL = []
for Element in Elements:
    Forces = op.eleForce(Element.EleTag)
    Forces.insert(0, Element.EleTag)
    ElemnsForceDL.append(Forces)
ElemnsForceDL = np.array(ElemnsForceDL)

# Create a Plain load pattern for seismic loading with a Linear TimeSeries (LLEF)
op.loadConst('-time', 0.0)
Htotal = Loc_heigth[-1]
SeismicLoadCode = self.ui.comboBoxSeismicLoadCode.currentText()
Ie = float(self.ui.Ie.text())
R = float(self.ui.R.text())
Cd = float(self.ui.Cd.text())
Omo = float(self.ui.Omo.text())
if SeismicLoadCode == 'ASCE 7-16':
    Sds = float(self.ui.Sds.text())
    Sd1 = float(self.ui.Sd1.text())
    Tl = float(self.ui.Tl.text())
    Ct = 0.0466
    x = 0.9
    Ta = Ct * Htotal ** x
    print('Ta =', Ta)
    Cu = np.interp(Sd1, [0, 0.1, 0.15, 0.2, 0.3, 5], [1.7, 1.7, 1.6, 1.5, 1.4, 1.4])
    print('Cu =', Cu)
    T = Cu * Ta
    Ts = Sd1 / Sds
    if T <= Ts:
        Sa = max(Sds * Ie / R, 0.044 * Sds * Ie, 0.01)
    elif T <= Tl:
        Sa = max(Sd1 * Ie / T / R, 0.044 * Sds * Ie, 0.01)
    else:
        Sa = max(Sd1 * Tl * Ie / (T ** 2) / R, 0.044 * Sds * Ie, 0.01)
elif SeismicLoadCode == 'NSR-10':
    Cd = R
    Aa = float(self.ui.Aa_10.text())
    Av = float(self.ui.Av_10.text())
    Fa = float(self.ui.Fa_10.text())
    Fv = float(self.ui.Fv_10.text())
    Sds = 2.5*Aa
    Ct = 0.047
    x = 0.9
    Ta = Ct * Htotal ** x
    print('Ta =', Ta)
    Cu = max(1.75-1.2*Av*Fv, 1.2)
    print('Cu =', Cu)
    T = Cu * Ta
    Tc = 0.48*(Av*Fv)/(Aa*Fa)
    Tl = 2.4*Fv
    if T <= Tc:
        Sa = 2.5*Aa*Fa*Ie/R
    elif T <= Tl:
        Sa = 1.2*Av*Fv*Ie/T/R
    else:
        Sa = 1.2*Av*Fv*Tl*Ie/T**2/R
elif SeismicLoadCode == 'NSR-98':
    Cd = R
    Aa = float(self.ui.Aa_98.text())
    S = float(self.ui.S_98.text())
    Sds = 2.5*Aa
    Ct = 0.08
    x = 3/4
    Ta = Ct * Htotal ** x
    print('Ta =', Ta)
    Cu = 1.2
    print('Cu =', Cu)
    T = Cu * Ta
    Tc = 0.48*S
    Tl = 2.4*S
    if T <= Tc:
        Sa = 2.5*Aa*Ie/R
    elif T <= Tl:
        Sa = 1.2*Aa*S*Ie/T/R
    else:
        Sa = Aa*Ie/R
elif SeismicLoadCode == 'CCCSR-84':
    Aa = float(self.ui.Aa_84.text())
    Av = float(self.ui.Av_84.text())
    S = float(self.ui.S_84.text())
    Sds = 2.5*Aa
    Ct = 0.08
    x = 3/4
    Ta = Ct * Htotal ** x
    print('Ta =', Ta)
    Cu = 1.2
    print('Cu =', Cu)
    T = Cu * Ta
    Tc = 0.48*S
    Tl = 2.4*S
    if S >= 1.5 and Aa >= 0.30:
        Sa = min(2.0*Aa*Ie/R, 1.2*Av*S*Ie/T**(2/3)/R)
    else:
        Sa = min(2.5*Aa*Ie/R, 1.2*Av*S*Ie/T**(2/3)/R)
elif SeismicLoadCode == 'Weight Percentage':
    Cd = 1
    Sa = float(self.ui.WP.text())/100
    Sds = Sa
    T = 0.5  # an arbitrary period is taken provided that the distribution is essentially triangular k = 1
    Omo = float(self.ui.Omo.text())
elif SeismicLoadCode == 'Spectra file':
    SpecFile = self.ui.SpectraFile.text()
    Spec = np.loadtxt('Spectra/' + SpecFile + '.txt')
    Ct = 0.08
    x = 3/4
    Ta = Ct * Htotal ** x
    print('Ta =', Ta)
    Cu = 1.2
    print('Cu =', Cu)
    T = Cu * Ta
    Sds = 2.5*Spec[0, 1]
    Sa = np.interp(T, Spec[:, 0], Spec[:, 1])/R
if T <= 0.5:
    k = 1.
elif T <= 2.5:
    k = 0.75 + 0.5 * T
else:
    k = 2.
sumH = np.sum(np.power(Loc_heigth, k))

op.timeSeries('Linear', 3)
op.pattern('Plain', 3, 1)
print('Wtotal =', Wtotal)
print('Sa =', Sa)
Fp = Sa * Wtotal * np.power(Loc_heigth, k) / sumH
print('FSis =', Fp)
for (fp, ind) in zip(Fp, range(len(Loc_heigth))):
    op.load(int(ListNodesDrift[ind, 0]), fp, 0.0, 0.0)
Vbasal = Sa * Wtotal

op.analyze(1)
ElemnsForceDLE = []
for Element in Elements:
    Forces = op.eleForce(Element.EleTag)
    Forces.insert(0, Element.EleTag)
    ElemnsForceDLE.append(Forces)
ElemnsForceDLE = np.array(ElemnsForceDLE)
np.set_printoptions(precision=6)
np.set_printoptions(suppress=True)

# Story drift caculations
DriftMax = 0.02
nodesDisp = []
Desp_x = []
Desp_y = []
Id_Node_Drift = ListNodesDrift[:, 0]
Id_Node_Drift = np.int64(Id_Node_Drift)
Id_Node_Drift = Id_Node_Drift.tolist()
Id_Node = ListNodes[:, 0]
Id_Node = np.int64(Id_Node)
Id_Node = Id_Node.tolist()
for nodo in Id_Node_Drift:
    nodesDisp.append([nodo, op.nodeDisp(nodo, 1)])
nodesDisp = np.array(nodesDisp) * Cd
for nodo in Id_Node:
    Desp_x.append([nodo, op.nodeDisp(nodo, 1)])
    Desp_y.append([nodo, op.nodeDisp(nodo, 2)])
Desp_x = np.array(Desp_x) * Cd
Desp_y = np.array(Desp_y) * Cd
drift = nodesDisp[1:, 1] - nodesDisp[:-1, 1]
drift_p = np.divide(drift, np.array(heigth_v))
ver_drift = np.where(drift_p < DriftMax, 'ok', 'no ok')
Id_Floor = np.arange(1, len(Loc_heigth))
drift_table = pd.DataFrame({"1.Floor": Id_Floor, "2.Drift": drift_p * 100, "3.": ver_drift})
print(drift_table)


# self.ui.progressBarBeamDesign.setFormat('designing elements . . .')
# self.ui.progressBarBeamDesign.setStyleSheet('text-align: center')

# subprocess.call('ACI_318S_19_IFM.py', shell=True)
# os.system("ACI_318S_19_IFM.py")

# runpy.run_path(path_name='ACI_318S_19_IFM.py')
# exec(open("ACI_318S_19_IFM.py").read())
# subprocess.Popen(['python', 'ACI_318S_19_IFM.py'])
# import ACI_318S_19_IFM

self.ui.progressBarBeamDesign.show()
# Beams and columns design procedures
Beta1B = beta1(fcB)
cover = 4 * cm
dst = 3 / 8 * inch
Ast = pi * dst ** 2 / 4.  # area de la barra del estribo
ro_max_b = min(0.85 * Beta1B * fcB * 3. / fy / 8., 0.025)  # maximun steel percentage
ro_min_b = max(0.25 * sqrt(fcB / MPa) * MPa / fy, 1.4 * MPa / fy)  # minimun steel percentage
DataBeamDesing = []
nprog = 0
nelems = num_elems - 1
for (Ele, EleForceD, EleForceDL, EleForceDLE) in zip(Elements, ElemnsForceD, ElemnsForceDL, ElemnsForceDLE):
    self.ui.progressBarBeamDesign.setValue(int(100 * nprog / nelems))
    nprog = nprog + 1
    if ListNodes[Ele.Nod_ini, 2] == ListNodes[Ele.Nod_end, 2]:
        VID = EleForceD[2]
        VIL = EleForceDL[2] - VID
        VIE = EleForceDLE[2] - VID - VIL
        VED = abs(EleForceD[5])
        VEL = abs(EleForceDL[5]) - VED
        VEE = abs(EleForceDLE[5]) - VED - VEL

        MID = EleForceD[3] - EleForceD[2] * Ele.RZi
        MIL = EleForceDL[3] - EleForceDL[2] * Ele.RZi - MID
        MIE = EleForceDLE[3] - EleForceDLE[2] * Ele.RZi - MID - MIL
        MED = EleForceD[6] + EleForceD[5] * Ele.RZe
        MEL = EleForceDL[6] + EleForceDL[5] * Ele.RZe - MED
        MEE = EleForceDLE[6] + EleForceDLE[5] * Ele.RZe - MED - MEL
        MED, MEL, MEE = -MED, -MEL, -MEE
        # print('MID ', MID, 'MED', MED, 'MIL ', MIL, 'MEL', MEL, 'MIE ', MIE, 'MEE', MEE)
        MI1, MI2, MI3, MI4, MI5 = Combo_ACI(MID, MIL, MIE)
        MNU1 = max([MI1, MI2, MI3, MI4, MI5, 0.])  # Negative initial design node moment
        MPU1 = min([MI1, MI2, MI3, MI4, MI5, abs(MNU1) / 2])  # Positive initial design node momentum
        ME1, ME2, ME3, ME4, ME5 = Combo_ACI(MED, MEL, MEE)
        MNU2 = max([ME1, ME2, ME3, ME4, ME5, 0.])  # Negative moment final design
        MPU2 = min([ME1, ME2, ME3, ME4, ME5, abs(MNU2) / 2])  # Positive moment final design
        Mmax = max([MNU1, -MPU1, MNU2, -MPU2])
        MNU1 = max([MNU1, Mmax / 4])
        MPU1 = min([MPU1, -Mmax / 4])
        MNU2 = max([MNU2, Mmax / 4])
        MPU2 = min([MPU2, -Mmax / 4])
        # print('MNU1 ', MNU1, 'MPU1', MPU1, 'MNU2 ', MNU2, 'MPU2', MPU2)
        Ast1, dt1, Mn_N1, db_t1, Mpr_N1, nbt1 = AsBeam(MNU1, Ele.EleTag, cover, ro_min_b, ro_max_b, dst, fy, BBeam, HBeam)
        Asb1, db1, Mn_P1, db_b1, Mpr_P1, nbb1 = AsBeam(MPU1, Ele.EleTag, cover, ro_min_b, ro_max_b, dst, fy, BBeam, HBeam)
        Ast2, dt2, Mn_N2, db_t2, Mpr_N2, nbt2 = AsBeam(MNU2, Ele.EleTag, cover, ro_min_b, ro_max_b, dst, fy, BBeam, HBeam)
        Asb2, db2, Mn_P2, db_b2, Mpr_P2, nbb2 = AsBeam(MPU2, Ele.EleTag, cover, ro_min_b, ro_max_b, dst, fy, BBeam, HBeam)

        VI1 = 1.2 * VID + 1.6 * VIL
        VI2 = 1.2 * VID + 1.0 * VIL - 1.0 * VIE
        VI3 = 0.9 * VID - 1.0 * VIE
        VI4 = abs(-(Mpr_P1 + Mpr_N2) / Ele.LEle + ((1.2 + 0.2*Sds)*WDL + 1.0*WLL) * Ele.LEle / 2.)
        VI5 = (Mpr_N1 + Mpr_P2) / Ele.LEle + ((1.2 + 0.2*Sds)*WDL + 1.0*WLL) * Ele.LEle / 2.

        VU1 = max(VI1, VI2, VI3, VI4, VI5)  # Negative shear node final design

        VE1 = 1.2 * VED + 1.6 * VEL
        VE2 = 1.2 * VED + 1.0 * VEL + 1.0 * VEE
        VE3 = 0.9 * VED + 1.0 * VEE
        VE4 = (Mpr_P1 + Mpr_N2) / Ele.LEle + ((1.2 + 0.2*Sds)*WDL + 1.0*WLL) * Ele.LEle / 2.
        VE5 = abs(-(Mpr_N1 + Mpr_P2) / Ele.LEle + ((1.2 + 0.2*Sds)*WDL + 1.0*WLL) * Ele.LEle / 2.)

        VU2 = max(VE1, VE2, VE3, VE4, VE5)  # Negative shear node final design

        Vpr_1 = (Mpr_P1 + Mpr_N2) / Ele.LEle + (1.2 * WDL + WLL) * Ele.LEle / 2.
        Vpr_2 = (Mpr_N1 + Mpr_P2) / Ele.LEle + (1.2 * WDL + WLL) * Ele.LEle / 2.
        Vpr = max(Vpr_1, Vpr_2)

        nb1 = max(nbt1, nbb1)
        nb2 = max(nbt2, nbb2)

        PID = EleForceD[2]
        PIL = EleForceDL[2] - PID
        PIE = EleForceDLE[2] - PID - PIL
        PED = -EleForceD[5]
        PEL = -EleForceDL[5] - PED
        PEE = -EleForceDLE[5] - PED - PEL
        PI1, PI2, PI3, PI4, PI5 = Combo_ACI(PID, PIL, PIE)
        PE1, PE2, PE3, PE4, PE5 = Combo_ACI(PED, PEL, PEE)
        PU1 = min(PI2, PI3, PI4, PI5)
        PU2 = min(PE2, PE3, PE4, PE5)
        if max(VI4, VI5) >= 0.5*VU1 and PU1 < BBeam*HBeam*fcB/20:
            Vc1 = 0
        else:
            Vc1 = 0.17 * sqrt(fcB / 1000.) * MPa * BBeam * d
        if max(VE4, VE5) >= 0.5*VU2 and PU2 < BBeam*HBeam*fcB/20:
            Vc2 = 0
        else:
            Vc2 = 0.17 * sqrt(fcB / 1000.) * MPa * BBeam * d

        nst1, sst1 = AvBeam(VU1, db_t1, dt1, Ele.EleTag, fys, dst, Ast, BBeam, nb1, Vc1)
        nst2, sst2 = AvBeam(VU2, db_t2, dt2, Ele.EleTag, fys, dst, Ast, BBeam, nb2, Vc2)

        DataBeamDesing.append(BeamDesing(Ele.EleTag, BBeam, HBeam, Ast1, dt1, Mn_N1, Asb1, db1, Mn_P1, nst1,
                                         sst1, Ast2, dt2, Mn_N2, Asb2, db2, Mn_P2, nst2, sst2, Ele.Nod_ini,
                                         Ele.Nod_end, db_t1, db_b1, db_t2, db_b2, Vpr, VU1, VU2))
        self.ui.tbl_data_design_beams.setRowCount(0)
        data_beams_table(self)

self.ui.progressBarBeamDesign.hide()
# self.QProgressBar.reset()

# Column design procedure

self.ui.progressBarColumnDesign.show()
# self.ui.progressBarBeamDesign.setValue(50)
# self.ui.progressBarColumnDesign.setFormat('designing columns . . .')
# self.ui.progressBarColumnDesign.setStyleSheet('text-align: center')
Beta1C = beta1(fcC)
DataColDesing = []
nprog = 0
for (Ele, EleForceD, EleForceDL, EleForceDLE) in zip(Elements, ElemnsForceD, ElemnsForceDL, ElemnsForceDLE):
    self.ui.progressBarColumnDesign.setValue(int(100 * nprog / nelems))
    nprog = nprog + 1
    if ListNodes[Ele.Nod_ini, 1] == ListNodes[Ele.Nod_end, 1]:
        if ListNodes[Ele.Nod_end, 2] == Loc_heigth[-1]:
            ncolsn = 1
        else:
            ncolsn = 2
        Mn_N_R, Mn_P_R, Mn_N_L, Mn_P_L = 0, 0, 0, 0
        for DB in DataBeamDesing:
            if Ele.Nod_end == DB.Nod_ini:
                Mn_N_R, Mn_P_R = DB.Mn_n1, DB.Mn_p1
            if Ele.Nod_end == DB.Nod_end:
                Mn_N_L, Mn_P_L = DB.Mn_n2, DB.Mn_p2
        Sum_Mn_B = max(Mn_P_R + Mn_N_L, Mn_N_R + Mn_P_L)
        # print('Node =', Ele.Nod_end, 'Sum_Mn_Beams =', Sum_Mn_B)
        b, h = Ele.BEle, Ele.HEle

        MID = EleForceD[3]
        MIL = EleForceDL[3] - MID
        MIE = EleForceDLE[3] - MID - MIL

        PID = EleForceD[2]
        PIL = EleForceDL[2] - PID
        PIE = EleForceDLE[2] - PID - PIL

        MI1, MI2, MI3, MI4, MI5 = Combo_ACI(MID, MIL, MIE)
        PI1, PI2, PI3, PI4, PI5 = Combo_ACI(PID, PIL, PIE)

        MED = -EleForceD[6]
        MEL = -EleForceDL[6] - MED
        MEE = -EleForceDLE[6] - MED - MEL
        # print('MID ', MID, 'MED', MED, 'MIL ', MIL, 'MEL', MEL, 'MIE ', MIE, 'MEE', MEE)

        PED = -EleForceD[5]
        PEL = -EleForceDL[5] - PED
        PEE = -EleForceDLE[5] - PED - PEL

        ME1, ME2, ME3, ME4, ME5 = Combo_ACI(MED, MEL, MEE)
        PE1, PE2, PE3, PE4, PE5 = Combo_ACI(PED, PEL, PEE)

        Nu_min = min([PI2, PI3, PI4, PI5, PE2, PE3, PE4, PE5])

        Pu_v = np.array([PI1, PI2, PI3, PI4, PI5, PE1, PE2, PE3, PE4, PE5])
        Mu_v = np.array([MI1, MI2, MI3, MI4, MI5, ME1, ME2, ME3, ME4, ME5])
        Mu_v = np.absolute(Mu_v)
        FactorColBeamStr = self.ui.Col_to_beam_str_ratio.value()

        nbH, nbB, db, As, fiPn, fiMn, Mn_i, Mpr_i, d, dist, ro, Mu_i, ColBeamStr = AsColumn(b, h, EleTag, cover, dst,
                                                                                            fy, Beta1C, Pu_v, Mu_v,
                                                                                            Sum_Mn_B, FactorColBeamStr,
                                                                                            ncolsn, Nu_min)

        VID = EleForceD[1]
        VIL = EleForceDL[1] - VID
        VIE = EleForceDLE[1] - VID - VIL
        VID, VIL, VIE = abs(VID), abs(VIL), abs(VIE)

        Mn_is = Mn_i[[1, 2, 3, 4, 6, 7, 8, 9]]
        Mn_max = np.max(Mn_is)  # Maximum moment nominal of all seismic combo
        Mpr_is = Mpr_i[[1, 2, 3, 4, 6, 7, 8, 9]]
        Mpr_max = np.max(Mpr_is)  # Maximum moment probable of all seismic combo

        VI1, VI2, VI3, VI4, VI5 = Combo_ACI(VID, VIL, VIE)
        VI6 = 2.0 * Mpr_max / Ele.LEle
        Vu = max([VI1, VI2, VI3, VI4, VI5, VI6])
        if VI6 >= 0.5*Vu and Nu_min < BCol*HCol*fcB/20:
            Vc = 0
        else:
            dp = cover + dst + db / 2
            d = h - dp
            Vc = (0.17 * sqrt(fcC * MPa) + Nu_min / 6 / b / h) * b * d
        sst, nsB, nsH, Vn, de = AvColumn(EleTag, Vu, b, h, nbH, nbB, dst, Vc, db, fys, Nu_min)
        NUG1 = abs(PID + 0.25 * PIL)
        NUG2 = abs(PED + 0.25 * PEL)
        NUD1 = abs(PID + 0.25 * PIL + PIE)
        NUD2 = abs(PED + 0.25 * PEL + PEE)
        MUD1 = abs(MID + 0.25 * MIL + MIE)
        MUD2 = abs(MED + 0.25 * MEL + MEE)
        VUD1 = abs(VID + 0.25 * VIL + VIE)
        VUD2 = abs(VED + 0.25 * VEL + VEE)
        DataColDesing.append(ColDesing(Ele.EleTag, b, h, nbH, nbB, db, de, As, Pu_v, Mu_v, fiPn, fiMn, Mn_i, d,
                                       dist, ro, Mu_i, sst, nsB, nsH, Ele.Nod_ini, Ele.Nod_end, NUD1, NUD2,
                                       NUG1, NUG2, MUD1, MUD2, VUD1, VUD2, ColBeamStr, Vn))

    self.ui.tbl_data_design_columns.setRowCount(0)
    data_columns_table(self)
    self.ui.tabWidget.setCurrentIndex(1)

self.ui.progressBarColumnDesign.hide()

# Frame Geometry plot
fig = self.ui.DataFrame.canvas.axes
fig.clear()
ax = fig.add_axes([0, 0, 0.9, 1])

ax.plot(ListNodes[:, 1], ListNodes[:, 2], 'ks')
# print('ListNodes[:, 1]', ListNodes[:, 1])
# print('Desp_x[-1, :]', Desp_x[:, -1])
Nodes_desp_x = ListNodes[:, 1] + 20 * Desp_x[:, -1]
Nodes_desp_y = ListNodes[:, 2] + 5 * Desp_y[:, -1]
ax.plot(Nodes_desp_x, Nodes_desp_y, 's', color='red', alpha=0.5)
ax.axis('off')
ind = 0

for Ele in Elements:
    xi = ListNodes[Ele.Nod_ini, 1]
    yi = ListNodes[Ele.Nod_ini, 2]
    xe = ListNodes[Ele.Nod_end, 1]
    ye = ListNodes[Ele.Nod_end, 2]
    ax.plot([xi, xe], [yi, ye], 'k-', alpha=.8)

    xid = Nodes_desp_x[Ele.Nod_ini]
    yid = Nodes_desp_y[Ele.Nod_ini]
    xed = Nodes_desp_x[Ele.Nod_end]
    yed = Nodes_desp_y[Ele.Nod_end]
    xd = np.array([xid, xed])
    yd = np.array([yid, yed])
    ax.plot(xd, yd, 'r-', alpha=.2)

    if xi == xe:
        ax.text(xi, (ye + yi) / 2, r'C{}'.format(Ele.EleTag), style='italic', fontsize=8, rotation='vertical',
                verticalalignment='center')
    if yi == ye:
        ax.text((xe + xi) / 2, yi, r'B{}'.format(Ele.EleTag), style='italic', fontsize=8,
                horizontalalignment='center')
        if xe == Loc_span[-1]:
            Delta_x = xed - xid
            ax.text(xed + .05 * Delta_x, yed, r'$\Delta = {:.2f} %$'.format(drift_p[ind] * 100),
                    style='italic', fontsize=8)
            ind += 1

for DC in DataColDesing:
    xi = ListNodes[DC.Nod_ini, 1]
    yi = ListNodes[DC.Nod_ini, 2]
    xe = ListNodes[DC.Nod_end, 1]
    ye = ListNodes[DC.Nod_end, 2]
    Delta_x = xe - xi
    Delta_y = ye - yi
    ax.text(xe + 0.05*Delta_y, ye + 0.10*Delta_y, r'{:.1f} '.format(DC.ColBeamStr), style='italic', fontsize=8,
                va='top', ha='right', multialignment="right", bbox=dict(boxstyle='round', fc="w", ec="k"))
ax.axis('equal')
fig.set_tight_layout(False)
self.ui.DataFrame.canvas.draw()
self.ui.DataFrame.canvas.show()
