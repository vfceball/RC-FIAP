global Loc_span, Loc_heigth, ListNodes, Elements, DataBeamDesign, DataColDesign, WDL, WLL, WDLS, Wtotal, \
    cover, num_elems, Beta1B, Beta1C, fcB, fcC, ListNodesDrift, ListNodesBasal, Ta, num_beams, num_cols,\
    DataWallDesign, ListNodesW1, ListNodesW2, nnt, ListEleTagCols, ListEleTagBeams, ListEleTagW1, ListEleTagW2, \
    ZeroLengthElement, ZeroLengths, net, num_walls, PWall1SlabD, PWall1SlabL, PWall2SlabD, PWall2SlabL, PColSlabD,\
    PColSlabL, FN, infill_m, infill_lw


# Function: Reads Beams design data from table that allows the user to modify the default design from TAB2 of GU
def data_beams_table(self):
    self.registros_beams = []
    for DB in DataBeamDesign:
        b = DB.b / cm
        h = DB.h / cm
        L_As_top = DB.Ast1 / (DB.b * DB.dt1) * 100
        L_As_bot = DB.Asb1 / (DB.b * DB.db1) * 100
        R_As_top = DB.Ast2 / (DB.b * DB.dt2) * 100
        R_As_bot = DB.Asb2 / (DB.b * DB.db2) * 100
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
    for DC in DataColDesign:
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
        Vu_Vn = DC.Vu_Vn
        registro = RegistroColumns(self.ui.tbl_data_design_columns, DC.EleTag, b, h, roCol, db, de, nbH, nbB, nsH, nsB,
                                   sst, Vu_Vn)
        self.registros_cols.append(registro)

def data_walls_table(self):
    self.registros_walls = []
    for DW in DataWallDesign:
        tw = DW.b / cm
        lw = DW.h / cm
        ro_l = DW.ro_l
        ro_t = DW.ro_t
        db = DW.db / mm
        dst = DW.dst / mm
        sst = DW.sst / cm
        cMaxS = DW.cMaxS / DW.h
        sigma_c = DW.sigma_c / fcC
        BE = DW.BE
        registro = RegistroWalls(self.ui.tbl_data_design_walls, DW.EleTag, tw, lw, ro_l, ro_t, db, dst, sst, cMaxS,
                                 sigma_c, BE)
        self.registros_walls.append(registro)

def beta1(fc):
    if fc <= 28 * MPa:
        Beta1 = 0.85
    else:
        Beta1 = max([0.85 - 0.05 * (fc - 28.) / 7., 0.65])
    return Beta1

# Design load combinations
def Combo_ACI(DL, LL, E):
    U1 = 1.4 * DL + 1.7 * LL
    U2 = 1.05 * DL + 1.28 * LL + 1.0 * E
    U3 = 1.05 * DL + 1.28 * LL - 1.0 * E
    U4 = 0.9 * DL + 1.0 * E
    U5 = 0.9 * DL - 1.0 * E
    return U1, U2, U3, U4, U5

# Flexural beams design
def AsBeam(Mu, EleTag, cover, ro_min_b, ro_max_b, dst, fy, BBeam, HBeam):
    b, h = BBeam, HBeam
    Mu = abs(Mu)
    db_v = np.array([4, 5, 6, 7, 8, 9, 10])
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
    Vc1 = 0.53 * sqrt(fcB * kgf/cm**2) * BBeam * d
    Vs = (Vu - 0.85 * Vc) / 0.85
    if Vs > 4. * Vc1:
        print("reshape by shear in Beam " + str(EleTag))
    elif 2. * Vc < Vs <= 4. * Vc:
        se_1 = d / 4
    elif Vs <= 2. * Vc:
        se_1 = d / 2
    nr_v = np.array([2, 3, 4])  # vector de numero de ramas
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
    # nra1 = ceil((BBeam - 2 * cover)/(15*cm))
    nra1 = 2
    nra = max(nra, nra1)
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
                        fiv = np.where(fiv >= 0.005, 0.90, fiv)
                        fiv = np.where(fiv <= 0.002, 0.70, fiv)
                        fiv = np.where((fiv > 0.002) & (fiv < 0.005), (0.70 + 0.25 * (fiv - 0.002) / 0.003), fiv)
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
                        fiv = np.hstack([0.9, fiv, 0.70])
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

# Walls P-M design
def AsWall(b, h, EleTag, cover, fy, Beta1C, Pu_v, Mu_v, Nu_min, PS, ro_min):
    # ro_min = 0.0015
    ro_max = 0.06
    dst = 3 / 8 * inch
    sl = np.arange(0.45, 0.05, -0.05)
    npts = 20
    # ncom = 10
    ecu = 0.003
    Es = 200. * GPa
    verif = False
    nbB = 2  # bars numbers along B
    # nbH = ceil(h * 10)  # bars numbers along H
    D_c = 1.1 * h / npts
    # nbB_v = np.array([nbB - 1, nbB, nbB + 1])
    db_v = np.array([2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 14, 18])  # vector bar diameters
    ro_b = ro_min
    while verif == False:
        for ndb in db_v:
            db = ndb / 8. * inch
            Ab = pi * db ** 2. / 4.
            dp = cover + dst + 0.5 * db
            d = h - dp
            nbH_v = np.ceil((h - 2 * dp) / sl)
            nbH_v = np.unique(nbH_v)
            nbH_v = np.int_(nbH_v)
            for nbH in nbH_v:
                nbT = 2. * (nbB + nbH - 2.)  # numero total de barras
                Ast = nbT * Ab
                ro = Ast / b / h
                As = np.hstack([nbB * Ab, np.ones(nbH - 2) * 2 * Ab, nbB * Ab])
                dist = np.linspace(dp, h - dp, nbH)
                if ro >= ro_min and ro >= ro_b:
                    # print('ro', ro, '#', ndb, 'nbH', nbH)
                    Pn_max = 0.80 * (0.85 * fcC * (b * h - Ast) + fy * Ast)
                    Tn_max = -fy * Ast
                    c = np.linspace(1.1 * h / npts, 1.1 * h, npts)
                    cv = np.copy(c)
                    a = Beta1C * c
                    Pconc = 0.85 * fcC * a * b
                    Mconc = Pconc * (h - a) / 2.
                    et = ecu * (d - c) / c
                    fiv = np.copy(et)
                    fiv = np.where(fiv >= 0.005, 0.9, fiv)
                    fiv = np.where(fiv <= 0.002, 0.7, fiv)
                    fiv = np.where((fiv > 0.002) & (fiv < 0.005), (0.7 + 0.25 * (fiv - 0.002) / 0.003), fiv)
                    c = c[:, np.newaxis]
                    es = ecu * (c - dist) / c
                    fs = Es * es
                    fs = np.where(fs > fy, fy, fs)
                    fs = np.where(fs < -fy, -fy, fs)
                    Pacer = np.sum(fs * As, axis=1)
                    Macer = np.sum(fs * As * (h / 2. - dist), axis=1)
                    Pn = np.hstack([Tn_max, np.where(Pconc + Pacer > Pn_max, Pn_max, Pconc + Pacer), Pn_max])
                    Mn = np.hstack([0, Mconc + Macer, 0])
                    fiv = np.hstack([0.9, fiv, 0.7])
                    cv = np.hstack([0, cv, 1e15])
                    fiPn = fiv * Pn
                    fiMn = fiv * Mn
                    ro_b = ro
                    if np.all((Pu_v >= min(fiPn)) & (Pu_v <= max(fiPn))):
                        Mu_i = np.interp(Pu_v, fiPn, fiMn)
                        Mn_i = np.interp(Pu_v, Pn, Mn)
                        Mns = np.interp(Nu_min, Pn, Mn)
                        if np.all(Mu_i >= Mu_v) == True:
                            # print('ok')
                            cMaxS = np.max(np.interp(PS, Pn, cv))
                            verif = True
                            break
            if verif == True:
                break
        if ndb == db_v[-1] and ro > ro_max:
            print('wall ' + str(EleTag) + 'needs to be resized by reinforcement ratio')
            break
    return nbH, nbB, db, As, fiPn, fiMn, Mn_i, d, dist, ro, Mu_i, cMaxS


# Shear columns design
def AvColumn(EleTag, Vu, b, h, nbH, nbB, dst, Vc, db, fys, Nu_min, Vc1):
    fiv = 0.85
    Ag = b * h
    Ast = pi * dst**2/4
    dp = cover + dst + db / 2
    d, dbc = h - dp, b - dp
    neH = floor(nbH / 2) + 1
    neB = floor(nbB / 2) + 1

    Ash_H = neH * Ast
    Ash_B = neB * Ast

    Vc = (0.17 * sqrt(fcC * MPa) + Nu_min / 6 / Ag) * b * d
    Vs = (Vu - fiv * Vc) / fiv
    if Vs > 4. * Vc:
        print("reshape by shear in Column " + str(EleTag))
    elif 2. * Vc < Vs <= 4. * Vc:
        se_1 = min(h / 4, b / 2)
    elif Vs <= 2. * Vc:
        se_1 = min(h / 2, b / 2)
    Ave = Ash_B  # area transversal del estribo
    if Vs <= 0.:
        se = se_1
    else:
        se_2 = Ave * fys * d / Vs
        se = min([se_1, se_2])
    if se < 60. * mm:
        print('Minimum spacing of stirrups is not met in column ' + str(EleTag))
    Vn = Vc + Ave * fys * d / se
    return se, neB, neH, Vn

def AvWall(EleTag, Vu, b, h, nbH, nbB, dst, Ast, Nu_min, db, fys, ro, cover):
    hw = Loc_heigth[-1]
    fiv = 0.85
    Ag = b * h
    dp = cover + dst + db / 2
    d = h - dp
    if hw / h <= 1.0:
        alfa_c = 0.25
    elif hw / h >= 2.0:
        alfa_c = 0.17
    else:
        alfa_c = np.interp(hw / h, [1.0, 2.0], [0.25, 0.17])
    Vc = (alfa_c * sqrt(fcC * MPa)) * Ag
    Vs = (Vu - fiv * Vc) / fiv
    Ave = 2 * Ast  # area transversal del estribo
    se_1 = Ave * fys * h / Vs
    ro_t = max(h / se_1 * Ave / Ag, 0.0020)
    se = h / (ro_t * Ag / Ave)
    ro_l = ro
    nbH = ceil(ro_l * Ag / (pi * db ** 2 / 4) / 2)
    Vn = Vc + Ave * fys * d / se
    if Vn > 0.66 * sqrt(fcC * MPa) * Ag:
        print('Wall ' + str(EleTag) + 'needs to be resized by shear limit')
    if se < 60. * mm:
        print('Minimum spacing of stirrups is not met in column ' + str(EleTag))
    return se, Vn, nbH, ro_t, ro_l

# Compression block parameters beta as function f'c

# Input geometric, materials and seismic design parameters from TAB1 of GUI
FN = float(self.ui.frames_numbers.text())
Lafg = float(self.ui.Lafg.text()) * FN
Lafs = float(self.ui.Lafs.text()) * FN
DL = float(self.ui.DL.text())
LL = float(self.ui.LL.text())
HColi = float(self.ui.HColi.text())  # Column inside Depth
BColi = float(self.ui.BColi.text()) * FN  # Column inside Width
HCole = float(self.ui.HCole.text())  # Column outside Depth
BCole = float(self.ui.BCole.text()) * FN  # Column outside Width
HBeam = float(self.ui.HBeam.text())
BBeam = float(self.ui.BBeam.text()) * FN
tw1, tw2 = float(self.ui.tw1.text()), float(self.ui.tw2.text())
lw1, lw2 = float(self.ui.lw1.text()), float(self.ui.lw2.text())
Afw1, Afw2 = float(self.ui.Af1.text()), float(self.ui.Af2.text())
IFC = float(self.ui.InertiaColumnsFactor.text())
IFB = float(self.ui.InertiaBeamsFactor.text())
IFW = float(self.ui.InertiaWallsFactor.text())
heigth_v = self.ui.heigth_v.text()
heigth_v = heigth_v.split(',')
# heigth_v = np.array(heigth_v, dtype=float)
span_v = self.ui.span_v.text()
span_v = span_v.split(',')
# span_v = np.array(span_v, dtype=float)
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

# Creation infill matrix
infill_v, infill_s = np.zeros(len(heigth_v)), np.zeros(len(span_v))
for ind in range(len(heigth_v)):
    infill_v[ind] = heigth_v[ind].count('i')
for ind in range(len(span_v)):
    infill_s[ind] = span_v[ind].count('i')
heigth_v = list(map(lambda elem: float(elem.split()[0]), heigth_v))
span_v = list(map(lambda elem: float(elem.split()[0]), span_v))
infill_lw = np.append(span_v, 0)
infill_lw, dummy = np.meshgrid(infill_lw, infill_v)
infill_lw = np.ravel(infill_lw)
infill_s = np.append(infill_s, 0)
infill_s, infill_v = np.meshgrid(infill_s, infill_v)
infill_m = infill_s*infill_v
infill_m = np.ravel(infill_m)
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
ListNodesW1, ListNodesW2 = None, None
if tw1 != 0:
    ListNodesW1 = np.empty([len(Loc_heigth), 3])
if tw2 != 0:
    ListNodesW2 = np.empty([len(Loc_heigth), 3])
nodeTag = int(ListNodes[-1, 0]) + 1
indW = 0
# print('ListNodes', ListNodes)
for node in ListNodes:
    if node[1] == Loc_span[-1] and tw1 != 0:
        xn, yn = node[1], node[2]
        op.node(nodeTag, xn, yn)
        if node[2] == 0.:
            # print(nodeTag, 1, 1, 1)
            op.fix(nodeTag, 1, 1, 1)
        ListNodesW1[indW, :] = [nodeTag, xn, yn]
        indW += 1
        nodeTag += 1
indW = 0
for node in ListNodes:
    if node[1] == Loc_span[-1] and tw2 != 0:
        xn, yn = node[1], node[2]
        op.node(nodeTag, xn, yn)
        if node[2] == 0.:
            # print(nodeTag, 1, 1, 1)
            op.fix(nodeTag, 1, 1, 1)
        ListNodesW2[indW, :] = [nodeTag, xn, yn]
        indW += 1
        nodeTag += 1
# print('ListNodesW1', ListNodesW1)
# print('ListNodesW2', ListNodesW2)
num_nodes_tot = num_nodes
if tw1 != 0 and tw2 != 0:
    num_nodes_tot = num_nodes + len(ListNodesW1[:, 0]) + len(ListNodesW2[:, 0])
    ListNodesTagW1 = np.int_(ListNodesW1[:, 0])
    ListNodesTagW2 = np.int_(ListNodesW2[:, 0])
elif tw1 != 0 and tw2 == 0:
    num_nodes_tot = num_nodes + len(ListNodesW1[:, 0])
    ListNodesTagW1 = np.int_(ListNodesW1[:, 0])
elif tw2 != 0 and tw1 == 0:
    num_nodes_tot = num_nodes + len(ListNodesW2[:, 0])
    ListNodesTagW2 = np.int_(ListNodesW2[:, 0])
ListNodesDrift = ListNodes[np.where(ListNodes[:, 1] == 0.)]
ListNodesBasal = ListNodes[np.where(ListNodes[:, 2] == 0.)]
nnt = num_nodes_tot
# print('nnt', nnt)
if tw1 != 0:
    for node in ListNodesW1[:-1, :]:
        # print(int(node[0] + nnt))
        nodeTag = int(node[0] + nnt)
        op.node(nodeTag, node[1], node[2])
        if node[2] == 0.:
            ListNodesBasal = np.vstack((ListNodesBasal, np.array([node[0], node[1], node[2]])))
if tw2 != 0:
    for node in ListNodesW2[:-1, :]:
        nodeTag = int(node[0] + nnt)
        op.node(nodeTag, node[1], node[2])
        if node[2] == 0.:
            ListNodesBasal = np.vstack((ListNodesBasal, np.array([node[0], node[1], node[2]])))
# print('ListNodesBasal', ListNodesBasal)
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
Elements, ZeroLengths = [], []
ListEleTagCols, ListEleTagBeams = [], []
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
        Elements.append(BeamElasticElement(EleTag, Nod_ini, Nod_end, ACol/FN, EcC, IzCol/FN, LCol, BCol/FN, HCol, gTr,
                                           RZi, RZe))
        # op.element('ElasticTimoshenkoBeam', EleTag, Nod_ini, Nod_end, EcC, 0.4*EcC, ACol, IzCol, ACol, gTr, '-mass',
        #            MassDens, MassType)
        # print('elasticBeamColumn', EleTag, Nod_ini, Nod_end, ACol, EcC, IzCol, gTr)
        op.element('elasticBeamColumn', EleTag, Nod_ini, Nod_end, ACol, EcC, IzCol, gTr, '-mass', MassDens,
                   MassType)
        ListEleTagCols.append(EleTag)
        EleTag += 1
num_cols = EleTag
# Beams creation for elastic analysis
op.geomTransf('Linear', 3, '-jntOffset', HColi / 2., 0, -HColi / 2., 0)
op.geomTransf('Linear', 4, '-jntOffset', HCole / 2., 0, -HColi / 2., 0)
op.geomTransf('Linear', 5, '-jntOffset', HColi / 2., 0, -HCole / 2., 0)
op.geomTransf('Linear', 6)
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
        Elements.append(BeamElasticElement(EleTag, Nod_ini, Nod_end, ABeam/FN, EcB, IzBeam/FN, LBeam, BBeam/FN, HBeam,
                                           gTr, RZi, RZe))
        # op.element('ElasticTimoshenkoBeam', EleTag, Nod_ini, Nod_end, EcB, 0.4*EcB, ABeam, IzBeam, ABeam, gTr, '-mass',
        #            MassDens, MassType)
        # print('elasticBeamColumn', EleTag, Nod_ini, Nod_end, ABeam, EcB, IzBeam, gTr)
        op.element('elasticBeamColumn', EleTag, Nod_ini, Nod_end, ABeam, EcB, IzBeam, gTr,
                   '-mass', MassDens, MassType)
        ListEleTagBeams.append(EleTag)
        EleTag += 1
# print('EleTag', EleTag)
num_elems = EleTag
num_beams = num_elems - num_cols

# print('ListNodesW1', ListNodesW1)
ListEleTagW1, ListEleTagW2 = [], []
RZi, RZe = 0, 0
if tw1 != 0:
    AWall1 = tw1 * lw1  # cross-sectional area
    IzWall1 = IFW * 1. / 12. * tw1 * lw1 ** 3  # Column moment of inertia
    # print('AWall1', AWall1, 'IzWall1', IzWall1)
    for Nod_ini in range(len(ListNodesW1[:, 0]) - 1):
        Nod_end = Nod_ini + 1
        gTr = 6
        LWall = ListNodesW1[Nod_end, 2] - ListNodesW1[Nod_ini, 2] - RZi - RZe
        MassDens = AWall1 * GConc / g
        Nod_iniW, Nod_endW = int(ListNodesW1[Nod_ini, 0] + nnt), int(ListNodesW1[Nod_end, 0])
        # print(Nod_iniW, Nod_endW)
        Elements.append(BeamElasticElement(EleTag, Nod_iniW, Nod_endW, AWall1, EcC, IzWall1, LWall, tw1, lw1, gTr,
                                           RZi, RZe))
        op.element('elasticBeamColumn', EleTag, Nod_iniW, Nod_endW, AWall1, EcC, IzWall1, gTr, '-mass', MassDens,
                   MassType)
        # op.element('elasticTimoshenkoBeam', EleTag, Nod_iniW, Nod_endW, EcC, 0.4*EcC, AWall1, IzWall1, AWall1, gTr,
        #            '-mass', MassDens, MassType)
        ListEleTagW1.append(EleTag)
        EleTag += 1
        # print('EleTag', EleTag)

if tw2 != 0:
    AWall2 = tw2 * lw2  # cross-sectional area
    IzWall2 = IFW * 1. / 12. * tw2 * lw2 ** 3  # Column moment of inertia
    # print('AWall2', AWall2, 'IzWall1', IzWall2)
    for Nod_ini in range(len(ListNodesW2) - 1):
        Nod_end = Nod_ini + 1
        gTr = 6
        LWall = ListNodesW1[Nod_end, 2] - ListNodesW1[Nod_ini, 2] - RZi - RZe
        MassDens = AWall2 * GConc / g
        Nod_iniW, Nod_endW = int(ListNodesW2[Nod_ini, 0] + nnt), int(ListNodesW2[Nod_end, 0])
        # print(Nod_iniW, Nod_endW)
        Elements.append(BeamElasticElement(EleTag, Nod_iniW, Nod_endW, AWall2, EcC, IzWall2, LWall, tw2, lw2, gTr,
                                           RZi, RZe))
        # print('elasticBeamColumn', EleTag, Nod_iniW, Nod_endW, AWall2, EcC, IzWall2, gTr)
        op.element('elasticBeamColumn', EleTag, Nod_iniW, Nod_endW, AWall2, EcC, IzWall2, gTr, '-mass', MassDens,
                   MassType)
        # op.element('elasticTimoshenkoBeam', EleTag, Nod_iniW, Nod_endW, EcC, 0.4*EcC, AWall2, IzWall2, AWall2, gTr,
        #            '-mass', MassDens, MassType)
        ListEleTagW2.append(EleTag)
        EleTag += 1
        # print('EleTag', EleTag)

# print('ListEleTagW1', ListEleTagW1)
# print('ListEleTagW2', ListEleTagW2)
num_elems = EleTag
num_walls = num_elems - num_cols - num_beams
# print('num_walls', num_walls)
op.uniaxialMaterial('Elastic', 1, 1e15)
ind = 0
for Nod_ini in range(num_nodes):
    if ListNodes[Nod_ini, 1] == Loc_span[-1] and ListNodes[Nod_ini, 2] != 0.:
        ind += 1
        if tw1 != 0:
            ZeroLengths.append(ZeroLengthElement(EleTag, Nod_ini, int(ListNodesTagW1[ind])))
            op.element('zeroLength', EleTag, Nod_ini, int(ListNodesTagW1[ind]), '-mat', 1, '-dir', 1)
            # op.element("corotTruss", EleTag, Nod_ini, int(ListNodesTagW1[ind]), 1e3, int(1e4))
            EleTag += 1
        if tw2 != 0:
            # print('EleTag', EleTag, Nod_ini, ListNodesTagW2[ind])
            ZeroLengths.append(ZeroLengthElement(EleTag, Nod_ini, int(ListNodesTagW2[ind])))
            op.element('zeroLength', EleTag, Nod_ini, int(ListNodesTagW2[ind]), '-mat', 1, '-dir', 1)
            # op.element("corotTruss", EleTag, Nod_ini, int(ListNodesTagW2[ind]), 1e3, int(1e4))
            EleTag += 1

num_elems_tot = EleTag - 1
net = num_elems_tot
# print('net', net)
for Ele in Elements:
    if Ele.EleTag in ListEleTagW1:
        EleTag = int(Ele.EleTag + net)
        # print(int(Ele.Nod_ini - nnt), int(Ele.Nod_ini))
        Nod_ini = int(Ele.Nod_ini - nnt)
        Ks = 0.4 * EcC * Ele.AEle / 1.2 / Ele.LEle
        # Ks = 1e15
        op.uniaxialMaterial('Elastic', EleTag, Ks)
        op.element('zeroLength', EleTag, Nod_ini, int(Ele.Nod_ini), '-mat', EleTag, 1, 1, '-dir', 1, 2, 3)
    if Ele.EleTag in ListEleTagW2:
        EleTag = int(Ele.EleTag + net)
        Nod_ini = int(Ele.Nod_ini - nnt)
        Ks = 0.4 * EcC * Ele.AEle / 1.2 / Ele.LEle
        # Ks = 1e15
        op.uniaxialMaterial('Elastic', EleTag, Ks)
        op.element('zeroLength', EleTag, Nod_ini, int(Ele.Nod_ini), '-mat', EleTag, 1, 1, '-dir', 1, 2, 3)

# print('num_elems_tot', num_elems_tot)
# Create a Plain load pattern for gravity loading with a Linear TimeSeries
Pvig = ABeam * GConc
PColi = AColi * GConc
PCole = ACole * GConc
PWall1SlabD, PWall2SlabD, PWall1SlabL, PWall2SlabL = 0, 0, 0, 0
if tw1 != 0:
    PWall1 = AWall1 * GConc
    PWall1SlabD, PWall1SlabL = Afw1 * DL, Afw1 * LL
if tw2 != 0:
    PWall2 = AWall2 * GConc
    PWall2SlabD, PWall2SlabL = Afw2 * DL, Afw2 * LL
PColSlabD, PColSlabL = 0, 0  #Pvig / FN * Lafs, 0 es para tener en cuenta el peso perpendicular de la viga
if self.ui.radioButtonSpatial.isChecked() and Lafs > Lafg:
    PColSlabD = (Lafs - Lafg) * DL + Pvig / FN * Lafs
    PColSlabL = (Lafs - Lafg) * LL

op.timeSeries('Linear', 1)
op.pattern('Plain', 1, 1)
for Element in Elements:
    if Element.EleTag in ListEleTagCols:
        LCol = Element.LEle + Element.RZi + Element.RZe
        if ListNodes[Element.Nod_ini, 1] == 0. or ListNodes[Element.Nod_ini, 1] == Loc_span[-1]:
            PCol = PCole #* LCol / Element.LEle
            PColSlabD1 = PColSlabD #+ HCole * Lafg * DL
        else:
            PCol = PColi #* LCol / Element.LEle
            PColSlabD1 = PColSlabD #+ HColi * Lafg * DL
        # print('-ele', Element.EleTag, '-type', '-beamUniform', 0, -PCol)
        # print('-ele', Element.EleTag, '-type', '-beamPoint', 0, 1, -PColSlabD)
        op.eleLoad('-ele', Element.EleTag, '-type', '-beamUniform', 0, -PCol)
        op.eleLoad('-ele', Element.EleTag, '-type', '-beamPoint', 0, 1, -PColSlabD1)
    if Element.EleTag in ListEleTagBeams:
        # print('-ele', Element.EleTag, '-type', '-beamUniform', -Pvig - WDL)
        op.eleLoad('-ele', Element.EleTag, '-type', '-beamUniform', -Pvig - WDL)
    if Element.EleTag in ListEleTagW1:
        # print('-ele', Element.EleTag, '-type', '-beamUniform', 0, -PWall1)
        # print('-ele', Element.EleTag, '-type', '-beamPoint', 0, 1, -PWall1SlabD)
        op.eleLoad('-ele', Element.EleTag, '-type', '-beamUniform', 0, -PWall1)
        op.eleLoad('-ele', Element.EleTag, '-type', '-beamPoint', 0, 1, -PWall1SlabD)
    if Element.EleTag in ListEleTagW2:
        # print('-ele', Element.EleTag, '-type', '-beamUniform', 0, -PWall2)
        # print('-ele', Element.EleTag, '-type', '-beamPoint', 0, 1, -PWall2SlabD)
        op.eleLoad('-ele', Element.EleTag, '-type', '-beamUniform', 0, -PWall2)
        op.eleLoad('-ele', Element.EleTag, '-type', '-beamPoint', 0, 1, -PWall2SlabD)
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
# print('ElemnsForceD', ElemnsForceD)
Wtotal = np.sum(ElemnsForceD[:len(Loc_span), 2]) # debo mirar esto
if tw1 != 0:
    Wtotal = Wtotal + ElemnsForceD[ListEleTagW1[0], 2]
if tw2 != 0:
    Wtotal = Wtotal + ElemnsForceD[ListEleTagW2[0], 2]
# print('Wtotal', Wtotal)
op.loadConst('-time', 0.0)
op.timeSeries('Linear', 2)
op.pattern('Plain', 2, 1)
for Element in Elements:
    if Element.EleTag in ListEleTagCols:
        if ListNodes[Element.Nod_ini, 1] == 0. or ListNodes[Element.Nod_ini, 1] == Loc_span[-1]:
            PColSlabL1 = PColSlabL #+ HCole / 2 * Lafg * LL
        else:
            PColSlabL1 = PColSlabL #+ HColi * Lafg * LL
        op.eleLoad('-ele', Element.EleTag, '-type', '-beamPoint', 0, 1, -PColSlabL1)
    if Element.EleTag in ListEleTagBeams:
        op.eleLoad('-ele', Element.EleTag, '-type', '-beamUniform', -WLL)
    if Element.EleTag in ListEleTagW1:
        # print('EleTag', Element.EleTag, '-beamPoint', -PWall1SlabL)
        op.eleLoad('-ele', Element.EleTag, '-type', '-beamPoint', 0, 1, -PWall1SlabL)
    if Element.EleTag in ListEleTagW2:
        # print('EleTag', Element.EleTag, '-beamPoint', -PWall2SlabL)
        op.eleLoad('-ele', Element.EleTag, '-type', '-beamPoint', 0, 1, -PWall2SlabL)
op.analyze(1)

ElemnsForceDL = []
for Element in Elements:
    Forces = op.eleForce(Element.EleTag)
    Forces.insert(0, Element.EleTag)
    ElemnsForceDL.append(Forces)
ElemnsForceDL = np.array(ElemnsForceDL)
# print('ElemnsForceDL', ElemnsForceDL)
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
    if tw1 != 0 or tw2 != 0:
        Ct, x = 0.0488, 0.75
    Ta = Ct * Htotal ** x
    print('Ta =', Ta)
    Cu = np.interp(Sd1, [0, 0.1, 0.15, 0.2, 0.3, 5], [1.7, 1.7, 1.6, 1.5, 1.4, 1.4])
    print('Cu =', Cu)
    T = Cu * Ta
    print('T =', T)
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
    Sds = 2.5 * Aa
    Ct = 0.047
    x = 0.9
    if tw1 != 0 or tw2 != 0:
        Ct, x = 0.049, 0.75
    Ta = Ct * Htotal ** x
    print('Ta =', Ta)
    Cu = max(1.75 - 1.2 * Av * Fv, 1.2)
    print('Cu =', Cu)
    T = Cu * Ta
    Tc = 0.48 * (Av * Fv) / (Aa * Fa)
    Tl = 2.4 * Fv
    if T <= Tc:
        Sa = 2.5 * Aa * Fa * Ie / R
    elif T <= Tl:
        Sa = 1.2 * Av * Fv * Ie / T / R
    else:
        Sa = 1.2 * Av * Fv * Tl * Ie / T ** 2 / R
elif SeismicLoadCode == 'NSR-98':
    Cd = R
    Aa = float(self.ui.Aa_98.text())
    S = float(self.ui.S_98.text())
    Sds = 2.5 * Aa
    Ct = 0.08
    x = 3 / 4
    if tw1 != 0 or tw2 != 0:
        Ct, x = 0.05, 3 / 4
    Ta = Ct * Htotal ** x
    print('Ta =', Ta)
    Cu = 1.2
    print('Cu =', Cu)
    T = Cu * Ta
    Tc = 0.48 * S
    Tl = 2.4 * S
    if T <= Tc:
        Sa = 2.5 * Aa * Ie / R
    elif T <= Tl:
        Sa = 1.2 * Aa * S * Ie / T / R
    else:
        Sa = Aa * Ie / R
elif SeismicLoadCode == 'CCCSR-84':
    Aa = float(self.ui.Aa_84.text())
    Av = float(self.ui.Av_84.text())
    S = float(self.ui.S_84.text())
    Sds = 2.5 * Aa
    Ct = 0.08
    x = 3 / 4
    if tw1 != 0 or tw2 != 0:
        Ct, x = 0.05, 3 / 4
    Ta = Ct * Htotal ** x
    print('Ta =', Ta)
    Cu = 1.2
    print('Cu =', Cu)
    T = Cu * Ta
    Tc = 0.48 * S
    Tl = 2.4 * S
    if S >= 1.5 and Aa >= 0.30:
        Sa = min(2.0 * Aa * Ie / R, 1.2 * Av * S * Ie / T ** (2 / 3) / R)
    else:
        Sa = min(2.5 * Aa * Ie / R, 1.2 * Av * S * Ie / T ** (2 / 3) / R)
elif SeismicLoadCode == 'Weight Percentage':
    Cd = 1
    Sa = float(self.ui.WP.text()) / 100
    Sds = Sa
    T = 0.5  # an arbitrary period is taken provided that the distribution is essentially triangular k = 1
elif SeismicLoadCode == 'Spectra file':
    SpecFile = self.ui.SpectraFile.text()
    Spec = np.loadtxt('Spectra/' + SpecFile + '.txt')
    Ct = 0.08
    x = 3 / 4
    if tw1 != 0 or tw2 != 0:
        Ct, x = 0.05, 3 / 4
    Ta = Ct * Htotal ** x
    print('Ta =', Ta)
    Cu = 1.2
    print('Cu =', Cu)
    T = Cu * Ta
    Sds = 2.5 * Spec[0, 1]
    Sa = np.interp(T, Spec[:, 0], Spec[:, 1]) / R
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
# print('ElemnsForceDLE', ElemnsForceDLE)
np.set_printoptions(precision=6)
np.set_printoptions(suppress=True)

# Story drift caculations
DriftMax = 0.02
nodesDisp = []
Desp_x = []
Desp_y = []
if tw1 != 0:
    Desp_xW1 = []
    Desp_yW1 = []
if tw2 != 0:
    Desp_xW2 = []
    Desp_yW2 = []
Id_Node_Drift = ListNodesDrift[:, 0]
Id_Node_Drift = np.int64(Id_Node_Drift)
Id_Node_Drift = Id_Node_Drift.tolist()
Id_Node = ListNodes[:, 0]
Id_Node = np.int64(Id_Node)
Id_Node = Id_Node.tolist()
if tw1 != 0:
    Id_NodeW1 = ListNodesW1[:, 0]
    Id_NodeW1 = np.int64(Id_NodeW1)
    Id_NodeW1 = Id_NodeW1.tolist()
if tw2 != 0:
    Id_NodeW2 = ListNodesW2[:, 0]
    Id_NodeW2 = np.int64(Id_NodeW2)
    Id_NodeW2 = Id_NodeW2.tolist()
for nodo in Id_Node_Drift:
    nodesDisp.append([nodo, op.nodeDisp(nodo, 1)])
nodesDisp = np.array(nodesDisp) * Cd
for nodo in Id_Node:
    Desp_x.append([nodo, op.nodeDisp(nodo, 1)])
    Desp_y.append([nodo, op.nodeDisp(nodo, 2)])
Desp_x = np.array(Desp_x) * Cd
Desp_y = np.array(Desp_y) * Cd
if tw1 != 0:
    for nodo in Id_NodeW1:
        Desp_xW1.append([nodo, op.nodeDisp(nodo, 1)])
        Desp_yW1.append([nodo, op.nodeDisp(nodo, 2)])
    Desp_xW1 = np.array(Desp_xW1) * Cd
    Desp_yW1 = np.array(Desp_yW1) * Cd
if tw2 != 0:
    for nodo in Id_NodeW2:
        Desp_xW2.append([nodo, op.nodeDisp(nodo, 1)])
        Desp_yW2.append([nodo, op.nodeDisp(nodo, 2)])
    Desp_xW2 = np.array(Desp_xW2) * Cd
    Desp_yW2 = np.array(Desp_yW2) * Cd
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
self.ui.progressBarBeamDesign.setFormat('designing beams . . .')
self.ui.progressBarBeamDesign.setStyleSheet('text-align: center')

# Beams and columns design procedures
Beta1B = beta1(fcB)
cover = 4 * cm
coverW = 2 * cm
dst = 3 / 8 * inch
Ast = pi * dst ** 2 / 4.  # area de la barra del estribo
ro_max_b = min(0.75 * 0.85 * Beta1B * fcB * 3. / fy / 5., 0.025)  # maximun steel percentage
ro_min_b = max(0.25 * sqrt(fcB / MPa) * MPa / fy, 1.4 * MPa / fy)  # minimun steel percentage
DataBeamDesign = []
nprog = 0
nelems = num_elems - 1
# print('num_beams', num_beams)
# QApplication.processEvents()
for (Ele, EleForceD, EleForceDL, EleForceDLE) in zip(Elements, ElemnsForceD, ElemnsForceDL, ElemnsForceDLE):
    # self.ui.progressBarBeamDesign.setValue(int(100 * nprog / num_elems))
    self.ui.progressBarBeamDesign.setValue(int(100 * nprog / num_beams))
    # QApplication.processEvents()
    nprog = nprog + 1
    if Ele.EleTag in ListEleTagBeams:
        b, h = Ele.BEle, Ele.HEle
        VID = EleForceD[2]/FN
        VIL = EleForceDL[2]/FN - VID
        VIE = EleForceDLE[2]/FN - VID - VIL
        VED = abs(EleForceD[5]/FN)
        VEL = abs(EleForceDL[5]/FN) - VED
        VEE = abs(EleForceDLE[5]/FN) - VED - VEL

        MID = EleForceD[3]/FN - EleForceD[2]/FN * Ele.RZi
        MIL = EleForceDL[3]/FN - EleForceDL[2]/FN * Ele.RZi - MID
        MIE = EleForceDLE[3]/FN - EleForceDLE[2]/FN * Ele.RZi - MID - MIL
        MED = EleForceD[6]/FN + EleForceD[5]/FN * Ele.RZe
        MEL = EleForceDL[6]/FN + EleForceDL[5]/FN * Ele.RZe - MED
        MEE = EleForceDLE[6]/FN + EleForceDLE[5]/FN * Ele.RZe - MED - MEL
        MED, MEL, MEE = -MED, -MEL, -MEE
        # print('MID ', MID, 'MED', MED, 'MIL ', MIL, 'MEL', MEL, 'MIE ', MIE, 'MEE', MEE)
        MI1, MI2, MI3, MI4, MI5 = Combo_ACI(MID, MIL, MIE)
        MNU1 = max([MI1, MI2, MI3, MI4, MI5])  # Negative initial design node moment
        MPU1 = min([MI1, MI2, MI3, MI4, MI5])  # Positive initial design node momentum
        ME1, ME2, ME3, ME4, ME5 = Combo_ACI(MED, MEL, MEE)
        MNU2 = max([ME1, ME2, ME3, ME4, ME5])  # Negative moment final design
        MPU2 = min([ME1, ME2, ME3, ME4, ME5])  # Positive moment final design
        # Mmax = max([MNU1, -MPU1, MNU2, -MPU2])
        # MNU1 = max([MNU1, Mmax / 4])
        # MPU1 = min([MPU1, -Mmax / 4])
        # MNU2 = max([MNU2, Mmax / 4])
        # MPU2 = min([MPU2, -Mmax / 4])

        # print('MNU1 ', MNU1, 'MPU1', MPU1, 'MNU2 ', MNU2, 'MPU2', MPU2)
        Ast1, dt1, Mn_N1, db_t1, Mpr_N1, nbt1 = AsBeam(MNU1, Ele.EleTag, cover, ro_min_b, ro_max_b, dst, fy, b, h)
        Asb1, db1, Mn_P1, db_b1, Mpr_P1, nbb1 = AsBeam(-MPU1, Ele.EleTag, cover, ro_min_b, ro_max_b, dst, fy, b, h)
        Ast2, dt2, Mn_N2, db_t2, Mpr_N2, nbt2 = AsBeam(MNU2, Ele.EleTag, cover, ro_min_b, ro_max_b, dst, fy, b, h)
        Asb2, db2, Mn_P2, db_b2, Mpr_P2, nbb2 = AsBeam(-MPU2, Ele.EleTag, cover, ro_min_b, ro_max_b, dst, fy, b, h)

        VI1 = 1.4 * VID + 1.7 * VIL
        VI2 = 1.05 * VID + 1.28 * VIL - 1.0 * VIE
        VI3 = 0.9 * VID - 1.0 * VIE
        VI4 = abs(-(Mn_P1 + Mn_N2) / Ele.LEle + (1.05 * WDL + 1.28 * WLL)/FN * Ele.LEle / 2.)
        VI5 = (Mn_N1 + Mn_P2) / Ele.LEle + (1.05 * WDL + 1.28 * WLL)/FN * Ele.LEle / 2.

        VI6 = 1.05 * VID + 1.28 * VIL - 2 * VIE
        VI7 = 0.9 * VID - 2 * VIE

        VU1a = np.max(np.abs([VI1, VI2, VI3]))
        VU1b = max(VI4, VI5)
        VU1c = max(VI6, VI7)

        VU1 = VU1a # Negative shear node final design

        VE1 = 1.4 * VED + 1.7 * VEL
        VE2 = 1.05 * VED + 1.28 * VEL + 1.0 * VEE
        VE3 = 0.9 * VED + 1.0 * VEE
        VE4 = (Mn_P1 + Mn_N2) / Ele.LEle + (1.05 * WDL + 1.28 * WLL)/FN * Ele.LEle / 2.
        VE5 = abs(-(Mn_N1 + Mn_P2) / Ele.LEle + (1.05 * WDL + 1.28 * WLL)/FN * Ele.LEle / 2.)

        VE6 = 1.05 * VED + 1.28 * VEL + 2 * VEE
        VE7 = 0.9 * VED + 2 * VEE

        VU2a = np.max(np.abs([VE1, VE2, VE3]))
        VU2b = max(VE4, VE5)
        VU2c = max(VE6, VE7)

        VU2 = VU2a  # Negative shear node final design

        Vpr_1 = (Mpr_P1 + Mpr_N2) / Ele.LEle + (1.2 * WDL + WLL)/FN * Ele.LEle / 2.
        Vpr_2 = (Mpr_N1 + Mpr_P2) / Ele.LEle + (1.2 * WDL + WLL)/FN * Ele.LEle / 2.
        Vpr = max(Vpr_1, Vpr_2)
        Vun = VU2b

        nb1 = max(nbt1, nbb1)
        nb2 = max(nbt2, nbb2)

        PID = EleForceD[2]/FN
        PIL = EleForceDL[2]/FN - PID
        PIE = EleForceDLE[2]/FN - PID - PIL
        PED = -EleForceD[5]/FN
        PEL = -EleForceDL[5]/FN - PED
        PEE = -EleForceDLE[5]/FN - PED - PEL
        PI1, PI2, PI3, PI4, PI5 = Combo_ACI(PID, PIL, PIE)
        PE1, PE2, PE3, PE4, PE5 = Combo_ACI(PED, PEL, PEE)
        PU1 = min(PI2, PI3, PI4, PI5)
        PU2 = min(PE2, PE3, PE4, PE5)
        Vc1 = 0.17 * sqrt(fcB * MPa) * b * dt1
        Vc2 = 0.17 * sqrt(fcB * MPa) * b * dt2

        nst1, sst1 = AvBeam(VU1, db_t1, dt1, Ele.EleTag, fys, dst, Ast, b, nb1, Vc1)
        nst2, sst2 = AvBeam(VU2, db_t2, dt2, Ele.EleTag, fys, dst, Ast, b, nb2, Vc2)

        DataBeamDesign.append(BeamDesign(Ele.EleTag, b, h, Ast1, dt1, Mn_N1, Asb1, db1, Mn_P1, nst1,
                                         sst1, Ast2, dt2, Mn_N2, Asb2, db2, Mn_P2, nst2, sst2, Ele.Nod_ini,
                                         Ele.Nod_end, db_t1, db_b1, db_t2, db_b2, Vpr, VU1, VU2, Mpr_N1, Mpr_P1,
                                         Mpr_N2, Mpr_P2, Vun))
        self.ui.tbl_data_design_beams.setRowCount(0)
        data_beams_table(self)
self.ui.progressBarBeamDesign.reset()
self.ui.progressBarBeamDesign.hide()
# self.QProgressBar.reset()

# Column design procedure

self.ui.progressBarBeamDesign.show()
# self.ui.progressBarBeamDesign.setValue(50)
self.ui.progressBarBeamDesign.setFormat('designing columns . . .')
self.ui.progressBarBeamDesign.setStyleSheet('text-align: center')
Beta1C = beta1(fcC)
DataColDesign = []
DataWallDesign = []
nprog = 0
# print('num_cols', num_cols)
# QApplication.processEvents()
for (Ele, EleForceD, EleForceDL, EleForceDLE) in zip(Elements, ElemnsForceD, ElemnsForceDL, ElemnsForceDLE):
    # self.ui.progressBarColumnDesign.setValue(float(100 * nprog / num_elems))
    self.ui.progressBarBeamDesign.setValue(int(100 * nprog / num_cols))
    nprog = nprog + 1
    if Ele.EleTag in ListEleTagCols:
        # QApplication.processEvents()
        if ListNodes[Ele.Nod_end, 2] == Loc_heigth[-1]:
            ncolsn = 1
        else:
            ncolsn = 2
        Mn_N_R, Mn_P_R, Mn_N_L, Mn_P_L = 0, 0, 0, 0
        Mpr_N_R, Mpr_P_R, Mpr_N_L, Mpr_P_L = 0, 0, 0, 0
        for DB in DataBeamDesign:
            if Ele.Nod_end == DB.Nod_ini:
                Mn_N_R, Mn_P_R = DB.Mn_n1, DB.Mn_p1
                Mpr_N_R, Mpr_P_R = DB.Mpr_n1, DB.Mpr_p1
            if Ele.Nod_end == DB.Nod_end:
                Mn_N_L, Mn_P_L = DB.Mn_n2, DB.Mn_p2
                Mpr_N_L, Mpr_P_L = DB.Mpr_n2, DB.Mpr_p2
        Sum_Mn_B = max(Mn_P_R + Mn_N_L, Mn_N_R + Mn_P_L)
        Sum_Mpr_B = max(Mpr_P_R + Mpr_N_L, Mpr_N_R + Mpr_P_L)
        # print('Node =', Ele.Nod_end, 'Sum_Mn_Beams =', Sum_Mn_B)
        b, h = Ele.BEle, Ele.HEle

        MID = EleForceD[3]/FN
        MIL = EleForceDL[3]/FN - MID
        MIE = EleForceDLE[3]/FN - MID - MIL

        PID = EleForceD[2]/FN
        PIL = EleForceDL[2]/FN - PID
        PIE = EleForceDLE[2]/FN - PID - PIL

        MI1, MI2, MI3, MI4, MI5 = Combo_ACI(MID, MIL, MIE)
        PI1, PI2, PI3, PI4, PI5 = Combo_ACI(PID, PIL, PIE)

        MED = -EleForceD[6]/FN
        MEL = -EleForceDL[6]/FN - MED
        MEE = -EleForceDLE[6]/FN - MED - MEL
        # print('MID ', MID, 'MED', MED, 'MIL ', MIL, 'MEL', MEL, 'MIE ', MIE, 'MEE', MEE)

        PED = -EleForceD[5]/FN
        PEL = -EleForceDL[5]/FN - PED
        PEE = -EleForceDLE[5]/FN - PED - PEL

        ME1, ME2, ME3, ME4, ME5 = Combo_ACI(MED, MEL, MEE)
        PE1, PE2, PE3, PE4, PE5 = Combo_ACI(PED, PEL, PEE)

        Nu_min = min([PI2, PI3, PI4, PI5, PE2, PE3, PE4, PE5])

        Pu_v = np.array([PI1, PI2, PI3, PI4, PI5, PE1, PE2, PE3, PE4, PE5])
        Mu_v = np.array([MI1, MI2, MI3, MI4, MI5, ME1, ME2, ME3, ME4, ME5])
        Mu_v = np.absolute(Mu_v)
        FactorColBeamStr = self.ui.Col_to_beam_str_ratio.value()

        nbH, nbB, db, As, fiPn, fiMn, Mn_i, Mpr_i, d, dist, ro,\
            Mu_i, ColBeamStr = AsColumn(b, h, Ele.EleTag, cover, dst, fy, Beta1C, Pu_v, Mu_v, Sum_Mn_B,
                                        FactorColBeamStr, ncolsn, Nu_min)

        # print('Mn_i =', Mn_i)
        VID = EleForceD[1]/FN
        VIL = EleForceDL[1]/FN - VID
        VIE = EleForceDLE[1]/FN - VID - VIL
        VID, VIL, VIE = abs(VID), abs(VIL), abs(VIE)
        # print('MID ', MID, 'MED', MED, 'MIL ', MIL, 'MEL', MEL, 'MIE ', MIE, 'MEE', MEE)

        Mn_is = Mn_i[[1, 2, 3, 4, 6, 7, 8, 9]]
        Mn_max = np.max(Mn_is)  # Maximum moment nominal of all seismic combo
        Mpr_is = Mpr_i[[1, 2, 3, 4, 6, 7, 8, 9]]
        Mpr_max = np.max(Mpr_is)  # Maximum moment probable of all seismic combo

        VI1, VI2, VI3, VI4, VI5 = Combo_ACI(VID, VIL, VIE)
        VI6 = 2.0 * Mn_max / Ele.LEle
        VI7 = 1.05 * abs(VID) + 1.28 * abs(VIL) + 2 * abs(VIE)
        VI8 = 0.9 * abs(VID) + 2 * abs(VIE)

        VU1a = np.max(np.abs([VI1, VI2, VI3]))
        VU1b = VI6
        VU1c = max(VI7, VI8)
        Vu = VU1a  # Negative shear node final design
        Vc = (0.17 * sqrt(fcC * MPa) + Nu_min / 6 / b / h) * b * d
        # print('Vu =', Vu, 'Vc =', Vc, 'Nu_min', Nu_min, 'Ele.EleTag', Ele.EleTag)
        Vc1 = (0.17 * sqrt(fcC * MPa) + Nu_min / 6 / b / h) * b * d
        sst, nsB, nsH, Vn = AvColumn(Ele.EleTag, Vu, b, h, nbH, nbB, dst, Vc, db, fys, Nu_min, Vc1)
        de = dst
        NUG1 = abs(PID + 0.25 * PIL)
        NUG2 = abs(PED + 0.25 * PEL)
        NUD1 = abs(PID + 0.25 * PIL + PIE)
        NUD2 = abs(PED + 0.25 * PEL + PEE)
        MUD1 = abs(MID + 0.25 * MIL + MIE)
        MUD2 = abs(MED + 0.25 * MEL + MEE)
        VUD1 = abs(VID + 0.25 * VIL + VIE)
        VUD2 = abs(VED + 0.25 * VEL + VEE)
        if Cd <= 2.0:
            knl = 1.0
        elif Cd >= 6.0:
            knl = 0.7
        else:
            knl = np.interp(Cd, [2.0, 6.0], [1.0, 0.7])
        sem = sst
        if sem / d <= 0.75:
            alfa_col = 1.0
        elif sem / d >= 1.0:
            alfa_col = 0.0
        else:
            alfa_col = np.interp(sem / d, [0.75, 1.0], [1.0, 0.0])
        ro_t = min(max(0.075, Ast * nsB / (b * d)), 0.0175)
        yi = ListNodes[Ele.Nod_ini, 2]
        ye = ListNodes[Ele.Nod_end, 2]
        LCol = ye - yi
        VnCol = knl * (alfa_col * Ast * nsB * fy * d / sst +
                        0.5 * sqrt(fcC * MPa) / (min(max(2, MUD1 / VUD1), 4) * d) *
                        sqrt(1 + NUD1 / (0.5 * b * h * sqrt(fcC * MPa))) * 0.8 * b * h)
        Vupr = min(2.0 * Mpr_max / Ele.LEle, Sum_Mpr_B)
        Vu_Vn = Vupr / VnCol
        sigma_IF, wIF = 0, 0
        # print("Vu %5.2f Vn %5.2f VnCol %5.2f Mpr_max %5.2f Vupr %5.2f MUD1 %5.2f VUD1 %5.2f NUD1 %5.2f"
        #       % (Vu, Vn, VnCol, Mpr_max, Vupr, MUD1, VUD1, NUD1))
        # print('knl', knl, 'alfa_col', alfa_col, 'Ast', Ast, 'nsB', nsB, 'sem', sem)

        DataColDesign.append(ColDesign(Ele.EleTag, b, h, nbH, nbB, db, dst, As, Pu_v, Mu_v, fiPn, fiMn, Mn_i, d,
                                       dist, ro, Mu_i, sst, nsB, nsH, Ele.Nod_ini, Ele.Nod_end, NUD1, NUD2,
                                       NUG1, NUG2, MUD1, MUD2, VUD1, VUD2, ColBeamStr, Vn, VI6, Vu_Vn, VnCol, sem,
                                       sigma_IF, wIF))
    self.ui.tbl_data_design_columns.setRowCount(0)
    data_columns_table(self)
    self.ui.tabWidget.setCurrentIndex(1)
self.ui.progressBarBeamDesign.reset()
self.ui.progressBarBeamDesign.hide()
# self.ui.progressBarBeamDesign.setValue(50)
nprog = 0
if tw1 != 0 or tw2 != 0:
    self.ui.progressBarBeamDesign.show()
    self.ui.progressBarBeamDesign.setFormat('designing walls . . .')
    self.ui.progressBarBeamDesign.setStyleSheet('text-align: center')
RFD = self.ui.Factor_red_steel_Wall.value()
for (Ele, EleForceD, EleForceDL, EleForceDLE) in zip(Elements, ElemnsForceD, ElemnsForceDL, ElemnsForceDLE):
    # self.ui.progressBarBeamDesign.setValue(int(100 * nprog / num_elems))
    if Ele.EleTag in ListEleTagW1 or Ele.EleTag in ListEleTagW2:
        if tw1 != 0 and Ele.EleTag == ListEleTagW1[0]:
            ro_min = 0.0015
        if tw2 != 0 and Ele.EleTag == ListEleTagW2[0]:
            ro_min = 0.0015
        self.ui.progressBarBeamDesign.setValue(int(100 * nprog / num_walls))
        # QApplication.processEvents()
        nprog = nprog + 1
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
        PS = np.array([PI2, PI3, PI4, PI5, PE2, PE3, PE4, PE5])
        MS = np.array([MI2, MI3, MI4, MI5, ME2, ME3, ME4, ME5])
        Pu_v = np.array([PI1, PI2, PI3, PI4, PI5, PE1, PE2, PE3, PE4, PE5])
        Mu_v = np.array([MI1, MI2, MI3, MI4, MI5, ME1, ME2, ME3, ME4, ME5])
        Mu_v = np.absolute(Mu_v)
        PS, MS = np.absolute(PS), np.absolute(MS)
        nbH, nbB, db, As, fiPn, fiMn, Mn_i, d, dist, ro, Mu_i, cMaxS = AsWall(b, h, Ele.EleTag, coverW, fy, Beta1C,
                                                                              Pu_v, Mu_v, Nu_min, PS, ro_min)
        ro_min = max(RFD*ro, 0.0020)
        # print('Mn_i =', Mn_i)
        VID = EleForceD[1]
        VIL = EleForceDL[1] - VID
        VIE = EleForceDLE[1] - VID - VIL
        VID, VIL, VIE = abs(VID), abs(VIL), abs(VIE)

        Mn_is = Mn_i[[1, 2, 3, 4, 6, 7, 8, 9]]
        Mn_max = np.max(Mn_is)  # Maximum moment of all seismic combo
        VI1, VI2, VI3, VI4, VI5 = Combo_ACI(VID, VIL, VIE)
        VI6 = 2.0 * Mn_max / Ele.LEle
        VI7 = 1.2 * VID + 1.0 * VIL + Omo * VIE
        VI8 = 1.2 * VID + 1.0 * VIL - Omo * VIE
        VI9 = 0.9 * VID + Omo * VIE
        VI10 = 0.9 * VID - Omo * VIE
        VUa = max([VI1, VI2, VI3, VI4, VI5])
        VUb = VI6
        VUc = max([VI7, VI8, VI9, VI10])
        Vu = VUa  # Cortante maximo de diseño
        sst, Vn, nbH, ro_t, ro_l = AvWall(Ele.EleTag, Vu, b, h, nbH, nbB, dst, Ast, Nu_min, db, fys, ro, coverW)
        NUG1 = abs(PID + 0.25 * PIL)
        NUG2 = abs(PED + 0.25 * PEL)
        NUD1 = abs(PID + 0.25 * PIL + PIE)
        NUD2 = abs(PED + 0.25 * PEL + PEE)
        MUD1 = abs(MID + 0.25 * MIL + MIE)
        MUD2 = abs(MED + 0.25 * MEL + MEE)
        VUD1 = abs(VID + 0.25 * VIL + VIE)
        VUD2 = abs(VED + 0.25 * VEL + VEE)
        sigma_c = np.max(PS/(b*h) + MS*(h/2)/(b*h**3/12))
        lbe = max(cMaxS/2, cMaxS-0.1*h)
        Ag_BE = b * lbe
        bc1, bc2 = lbe - coverW, b - 2 * coverW
        dp_BE = coverW + dst + db / 2
        nsH_BE = ceil(bc1/(20*cm))+1
        nsB_BE = ceil(bc2/(20*cm))+1
        Ast = pi * dst ** 2. / 4.
        Ash_H, Ash_B = nsH_BE * Ast, nsB_BE * Ast
        hx = (lbe - 2 * dp_BE) / nsH_BE

        if WallDetailing == 'OMF':
            BE = 'No'
            se_BE = sst
        elif WallDetailing == 'IMF':
            if sigma_c >= 0.3*fcC:
                BE = 'Yes'
                se_1 = min(8. * db, 16 * dst, b / 2., 150*mm)  # minimum spacing c.18.7.5.3 ACI-19
                se_2 = min(Ash_H / bc1 / (0.06 * fcC / fys), Ash_B / bc2 / (0.06 * fcC / fys))
                se_BE = min([se_1, se_2])
            else:
                BE = 'NO'
                se_BE = sst
        elif WallDetailing == 'SMF':
            if sigma_c >= 0.2*fcC:
                BE = 'Yes'
                se_1 = min(6. * db, b / 3.)  # minimum spacing c.18.7.5.3 ACI-19
                se_2 = min(Ash_H / bc1 / (0.09 * fcC / fys), Ash_B / bc2 / (0.09 * fcC / fys))
                se_BE = min([se_1, se_2])
            else:
                BE = 'NO'
                se_BE = sst
        DataWallDesign.append(WallDesign(Ele.EleTag, b, h, nbH, nbB, db, dst, As, Pu_v, Mu_v, fiPn, fiMn, Mn_i, d, dist,
                                         ro_t, ro_l, Mu_i, sst, Ele.Nod_ini, Ele.Nod_end, NUD1, NUD2, NUG1, NUG2, MUD1,
                                         MUD2, VUD1, VUD2, Vn, cMaxS, sigma_c, BE, lbe, nsH_BE, nsB_BE, se_BE))

    self.ui.tbl_data_design_walls.setRowCount(0)
    data_walls_table(self)
    self.ui.tabWidget.setCurrentIndex(1)
self.ui.progressBarBeamDesign.reset()
self.ui.progressBarBeamDesign.hide()
# self.ui.progressBarBeamDesign.reset()


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
if tw1 != 0 and tw2 != 0:
    sepW = 2 + lw1
    numW1 = len(ListNodesW1[:, 0])
elif tw1 == 0 and tw2 != 0:
    numW1 = 0
    sepW = 0.
if tw1 != 0:
    Nodes_desp_xW1 = 2 + ListNodesW1[:, 1] + 20 * Desp_xW1[:, -1]
    Nodes_desp_yW1 = ListNodesW1[:, 2] + 5 * Desp_yW1[:, -1]
if tw2 != 0:
    Nodes_desp_xW2 = 2 + ListNodesW2[:, 1] + sepW + 20 * Desp_xW2[:, -1]
    Nodes_desp_yW2 = ListNodesW2[:, 2] + 5 * Desp_yW2[:, -1]

# if tw1 != 0:
#     ax.plot(Nodes_desp_xW1, Nodes_desp_yW1, 's', color='red', alpha=0.5)
#     ax.plot(Nodes_desp_xW1 + lw1, Nodes_desp_yW1, 's', color='red', alpha=0.5)
# if tw2 != 0:
#     ax.plot(Nodes_desp_xW2, Nodes_desp_yW2, 's', color='red', alpha=0.5)
#     ax.plot(Nodes_desp_xW2 + lw2, Nodes_desp_yW2, 's', color='red', alpha=0.5)
ax.axis('off')
ind = 0
from matplotlib.patches import Rectangle, Polygon

xdF, ydF = np.empty(shape=0), np.empty(shape=0)
for Ele in Elements:
    if Ele.EleTag in ListEleTagCols or Ele.EleTag in ListEleTagBeams:
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
            if xe == Loc_span[-1] and tw1 == 0 and tw2 == 0:
                Delta_x = xed - xid
                ax.text(xed + .05 * Delta_x, yed, r'$\Delta = {:.2f} %$'.format(drift_p[ind] * 100),
                        style='italic', fontsize=8)
                ind += 1
            if xe == Loc_span[-1]:
                # print('xed', xed)
                xdF, ydF = np.append(xed, xdF), np.append(yed, ydF)

    if Ele.EleTag in ListEleTagW1:
        ind_ini, ind_end = Ele.Nod_ini - num_nodes - nnt, Ele.Nod_end - num_nodes
        xiW1, yiW1 = ListNodesW1[ind_ini, 1] + 2, ListNodesW1[ind_ini, 2]
        xeW1, yeW1 = ListNodesW1[ind_end, 1] + 2, ListNodesW1[ind_end, 2]
        hW1 = ListNodesW1[ind_end, 2] - ListNodesW1[ind_ini, 2]
        ax.add_patch(Rectangle((xiW1, yiW1), lw1, hW1, ec="black", fc="silver"))
        ax.text(xiW1 + lw1 / 2, (yeW1 + yiW1) / 2, r'W{}'.format(Ele.EleTag), style='italic', fontsize=8,
                rotation='vertical', verticalalignment='center')

    if Ele.EleTag in ListEleTagW2:
        ind_ini, ind_end = Ele.Nod_ini - num_nodes - numW1 - nnt, Ele.Nod_end - num_nodes - numW1
        xiW2, yiW2 = sepW + ListNodesW2[ind_ini, 1] + 2, ListNodesW2[ind_ini, 2]
        xeW2, yeW2 = sepW + ListNodesW2[ind_end, 1] + 2, ListNodesW2[ind_end, 2]
        hW2 = ListNodesW2[ind_end, 2] - ListNodesW2[ind_ini, 2]
        ax.add_patch(Rectangle((xiW2, yiW2), lw2, hW2, ec="black", fc="silver"))
        ax.text(xiW2 + lw2 / 2, (yeW2 + yiW2) / 2, r'W{}'.format(Ele.EleTag), style='italic', fontsize=8,
                rotation='vertical', verticalalignment='center')
ind = 0
for Ele in Elements:
    if Ele.EleTag in ListEleTagW1:
        ind_ini, ind_end = Ele.Nod_ini - num_nodes - nnt, Ele.Nod_end - num_nodes
        xid, yid = Nodes_desp_xW1[ind_ini], Nodes_desp_yW1[ind_ini]
        xed, yed = Nodes_desp_xW1[ind_end], Nodes_desp_yW1[ind_end]
        lat = np.array([xid, xid + lw1, xed + lw1, xed, xid])
        long = np.array([yid, yid, yed, yed, yid])
        ax.add_patch(Polygon(np.c_[lat, long], facecolor='red', alpha=0.5))
        if tw2 == 0:
            Delta_x, ind = xed - xid, ind_ini
            ax.text(xed + lw1 + .1 * lw1, yed, r'$\Delta = {:.2f} %$'.format(drift_p[ind] * 100), style='italic',
                    fontsize=8)
        fpos = 0.1
        # print('xdF', xdF)
        # print('ind_ini', ind_ini)
        # xi, xe, yi, ye = xdF[ind_ini], xed, ydF[ind_ini], yed
        # Delta_x = xe - xi
        # Delta_y = ye - yi
        # xi_R = xi + fpos * Delta_x
        # yi_R = yi + fpos * Delta_y
        # xe_R = xe - fpos * Delta_x
        # ye_R = ye - fpos * Delta_y
        # ax.plot([xi, xe], [yi, ye], 'r-', alpha=.2)
        # ax.scatter([xi_R, xe_R], [yi_R, ye_R], s=5, c='red', alpha=0.2)
        ind += 1
    if Ele.EleTag in ListEleTagW2:
        ind_ini, ind_end = Ele.Nod_ini - num_nodes - numW1 - nnt, Ele.Nod_end - num_nodes - numW1
        xid, yid = Nodes_desp_xW2[ind_ini], Nodes_desp_yW2[ind_ini]
        xed, yed = Nodes_desp_xW2[ind_end], Nodes_desp_yW2[ind_end]
        lat = np.array([xid, xid + lw2, xed + lw2, xed, xid])
        long = np.array([yid, yid, yed, yed, yid])
        ax.add_patch(Polygon(np.c_[lat, long], facecolor='red', alpha=0.5))
        Delta_x, ind = xed - xid, ind_ini
        ax.text(xed + lw2 + .1 * lw2, yed, r'$\Delta = {:.2f} %$'.format(drift_p[ind] * 100), style='italic',
                fontsize=8)
        ind += 1

        # lat = [20.2, 30.2, 40.3, 50.3]
        # long = [21.3, 22.3, 22.7, 22.9]
        # lat = np.append(lat, lat[0])
        # long = np.append(long, long[0])
        #
        # plt.gca().add_patch(plt.Polygon(np.c_[lat, long], facecolor='red', edgecolor="blue"))
        # plt.gca().autoscale()
        # plt.show()

for DC in DataColDesign:
    xi = ListNodes[DC.Nod_ini, 1]
    yi = ListNodes[DC.Nod_ini, 2]
    xe = ListNodes[DC.Nod_end, 1]
    ye = ListNodes[DC.Nod_end, 2]
    Delta_x = xe - xi
    Delta_y = ye - yi
    ax.text(xe + 0.05 * Delta_y, ye + 0.10 * Delta_y, r'{:.1f} '.format(DC.ColBeamStr), style='italic', fontsize=8,
            va='top', ha='right', multialignment="right", bbox=dict(boxstyle='round', fc="w", ec="k"))
ax.axis('equal')
fig.set_tight_layout(False)
self.ui.DataFrame.canvas.draw()
self.ui.DataFrame.canvas.show()
