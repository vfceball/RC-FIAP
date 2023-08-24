## ############################################################### ##
## RC_FIAP (Reinforced Concrete Frame Inelastic Analysis Platform) ##
##                                                                 ##
## Developed by:                                                   ##
##       Victor F. Ceballos (vceballos@uninorte.edu.co)            ##
##       Carlos A. Arteta (carteta@uninorte.edu.co)                ##
## RCFIAPMain.py : this is the main script that calls              ##
## GUIFrameNonLinearACI.py : graphical environment                 ##
## mplwidget.py : cript to help plot the plastic hinge projector   ##
## ############################################################### ##

global Loc_span, Loc_heigth, ListNodes, Elements, DataBeamDesign, DataColDesign, Wtotal, num_elems, ListNodesDrift, \
    cIndex, ListNodesBasal, T1m, IM, Sa_max, RDR_max, SDR_max, nrecs, DataWallDesign, net

import sys
from math import pi, sqrt, ceil, floor, atan, degrees
from scipy import interpolate
from scipy.stats import norm

import openseespy.opensees as op
from PyQt5.QtWidgets import *
# from PyQt5.uic import loadUi
from PyQt5.QtGui import QDoubleValidator, QIntValidator
from PyQt5.QtCore import Qt

from PyQt5.QtWidgets import QDialog, QApplication
from GUIFrameNonLinearACI import *
import numpy as np  # load the numpy module, calling it np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import colors, colorbar
from matplotlib.ticker import MaxNLocator
from matplotlib.patches import Rectangle, Polygon

from threading import Thread

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure

matplotlib.use("Qt5Agg")
import pandas as pd
import os
import subprocess
import runpy

from mpl_toolkits.axes_grid1 import make_axes_locatable

# Definition of units
m = 1.  # define basic units -- output units
kN = 1.  # define basic units -- output units
sec = 1.  # define basic units -- output units
mm = m / 1000.  # define engineering units
cm = m / 100.
N = kN / 1000.
MPa = N / mm ** 2
GPa = MPa * 1000
m2 = m ** 2  # m^2
m3 = m ** 3  # m^3
m4 = m ** 4  # m^4
inch = cm * 2.54
ft = 12. * inch
g = 9.81 * m / sec ** 2  # gravitational acceleration
kip = 4.448 * kN
ksi = kip / inch ** 2
psi = ksi / 1000.
lbf = psi * inch ** 2  # pounds force
pcf = lbf / ft ** 3  # pounds per cubic foot
psf = lbf / ft ** 3  # pounds per square foot
in2 = inch ** 2  # inch^2
in4 = inch ** 4  # inch^4
GConc = 24. * kN / m ** 3  # Specific gravity of concrete
kgf = g * N
tonf = 1000 * kgf
cbar = False
np.set_printoptions(precision=6)


class RegistroBeams:
    def __init__(self, tbl_data_design_beams, id_, b, h, L_As_top, L_As_bot, L_Leg_n, L_Sstirrup, R_As_top, R_As_bot,
                 R_Leg_n, R_Sstirrup):
        fila = tbl_data_design_beams.rowCount()
        tbl_data_design_beams.insertRow(fila)

        self.spx_id = QLineEdit(tbl_data_design_beams)  # setattr(self, 'spx_id', QLineEdit(tbl_data_design_beams))
        self.spx_id.setValidator(QIntValidator(0, 100))
        self.spx_id.setText(f'B{id_}')
        # self.spx_id.setStyleSheet('border-top: none; border-right: none; border-bottom: none')
        # self.spx_id.setFont(('Times', 10))

        self.spx_b = QLineEdit(tbl_data_design_beams)
        self.spx_b.setValidator(QIntValidator(20, 1000))
        self.spx_b.setText('{:d}'.format(int(b)))
        # self.spx_b.setStyleSheet('border-top: none; border-right: none; border-bottom: none')

        self.spx_h = QLineEdit(tbl_data_design_beams)
        self.spx_h.setValidator(QIntValidator(20, 1000))
        self.spx_h.setText('{:d}'.format(int(h)))
        # self.spx_h.setStyleSheet('border-top: none; border-right: none; border-bottom: none')

        self.spx_L_As_top = QLineEdit(tbl_data_design_beams)
        self.spx_L_As_top.setValidator(QDoubleValidator(2., 400., 2))
        self.spx_L_As_top.setText('{:.3f}'.format(L_As_top))
        # self.spx_L_As_top.setStyleSheet('border-top: none; border-right: none; border-bottom: none')

        self.spx_L_As_bot = QLineEdit(tbl_data_design_beams)
        self.spx_L_As_bot.setValidator(QDoubleValidator(2., 400., 2))
        self.spx_L_As_bot.setText('{:.3f}'.format(L_As_bot))

        self.spx_L_Leg_n = QLineEdit(tbl_data_design_beams)
        self.spx_L_Leg_n.setValidator(QIntValidator(2, 10))
        self.spx_L_Leg_n.setText('{:d}'.format(int(L_Leg_n)))

        self.spx_L_Sstirrup = QLineEdit(tbl_data_design_beams)
        self.spx_L_Sstirrup.setValidator(QIntValidator(4, 30))
        self.spx_L_Sstirrup.setText('{:d}'.format(int(L_Sstirrup)))

        self.spx_R_As_top = QLineEdit(tbl_data_design_beams)
        self.spx_R_As_top.setValidator(QDoubleValidator(2., 400., 2))
        self.spx_R_As_top.setText('{:.3f}'.format(R_As_top))

        self.spx_R_As_bot = QLineEdit(tbl_data_design_beams)
        self.spx_R_As_bot.setValidator(QDoubleValidator(2., 400., 2))
        self.spx_R_As_bot.setText('{:.3f}'.format(R_As_bot))

        self.spx_R_Leg_n = QLineEdit(tbl_data_design_beams)
        self.spx_R_Leg_n.setValidator(QIntValidator(2, 10))
        self.spx_R_Leg_n.setText('{:d}'.format(int(R_Leg_n)))

        self.spx_R_Sstirrup = QLineEdit(tbl_data_design_beams)
        self.spx_R_Sstirrup.setValidator(QIntValidator(4, 30))
        self.spx_R_Sstirrup.setText('{:d}'.format(int(R_Sstirrup)))

        tbl_data_design_beams.setCellWidget(fila, 0, self.spx_id)
        tbl_data_design_beams.setCellWidget(fila, 1, self.spx_b)
        tbl_data_design_beams.setCellWidget(fila, 2, self.spx_h)

        tbl_data_design_beams.setCellWidget(fila, 3, self.spx_L_As_top)
        tbl_data_design_beams.setCellWidget(fila, 4, self.spx_L_As_bot)
        tbl_data_design_beams.setCellWidget(fila, 5, self.spx_L_Leg_n)
        tbl_data_design_beams.setCellWidget(fila, 6, self.spx_L_Sstirrup)

        tbl_data_design_beams.setCellWidget(fila, 7, self.spx_R_As_top)
        tbl_data_design_beams.setCellWidget(fila, 8, self.spx_R_As_bot)
        tbl_data_design_beams.setCellWidget(fila, 9, self.spx_R_Leg_n)
        tbl_data_design_beams.setCellWidget(fila, 10, self.spx_R_Sstirrup)

        tbl_data_design_beams.setColumnWidth(0, 40)
        tbl_data_design_beams.setColumnWidth(1, 40)
        tbl_data_design_beams.setColumnWidth(2, 40)
        tbl_data_design_beams.setColumnWidth(3, 50)
        tbl_data_design_beams.setColumnWidth(4, 50)
        tbl_data_design_beams.setColumnWidth(5, 40)
        tbl_data_design_beams.setColumnWidth(6, 50)
        tbl_data_design_beams.setColumnWidth(7, 50)
        tbl_data_design_beams.setColumnWidth(8, 50)
        tbl_data_design_beams.setColumnWidth(9, 40)
        tbl_data_design_beams.setColumnWidth(10, 50)

        stylesheet = "::section{border-style: solid;" \
                     "border-width: 1px;}"
        tbl_data_design_beams.horizontalHeader().setStyleSheet(stylesheet)


class RegistroColumns:
    def __init__(self, tbl_data_design_columns, id_, b, h, ro, db, de, nbH, nbB, Leg_n_H, Leg_n_B, Sstirrup, Vu_Vn):
        fila = tbl_data_design_columns.rowCount()
        tbl_data_design_columns.insertRow(fila)

        self.spx_id = QLineEdit(tbl_data_design_columns)
        self.spx_id.setValidator(QIntValidator(0, 1000))
        self.spx_id.setText(f'C{id_}')

        self.spx_b = QLineEdit(tbl_data_design_columns)
        self.spx_b.setValidator(QIntValidator(20, 1000))
        self.spx_b.setText('{:d}'.format(int(b)))

        self.spx_h = QLineEdit(tbl_data_design_columns)
        self.spx_h.setValidator(QIntValidator(20, 1000))
        self.spx_h.setText('{:d}'.format(int(h)))

        self.spx_ro = QLineEdit(tbl_data_design_columns)
        self.spx_ro.setValidator(QDoubleValidator(2., 400., 2))
        self.spx_ro.setText('{:.2f}'.format(ro * 100))

        self.spx_db = QLineEdit(tbl_data_design_columns)
        self.spx_db.setValidator(QDoubleValidator(1., 10., 2))
        self.spx_db.setText('{:.2f}'.format(db))

        self.spx_de = QLineEdit(tbl_data_design_columns)
        self.spx_de.setValidator(QDoubleValidator(1., 10., 2))
        self.spx_de.setText('{:.2f}'.format(de))

        self.spx_nbH = QLineEdit(tbl_data_design_columns)
        self.spx_nbH.setValidator(QIntValidator(2, 100))
        self.spx_nbH.setText('{:d}'.format(int(nbH)))

        self.spx_nbB = QLineEdit(tbl_data_design_columns)
        self.spx_nbB.setValidator(QIntValidator(2, 100))
        self.spx_nbB.setText('{:d}'.format(int(nbB)))

        self.spx_Leg_n_H = QLineEdit(tbl_data_design_columns)
        self.spx_Leg_n_H.setValidator(QIntValidator(2, 100))
        self.spx_Leg_n_H.setText('{:d}'.format(int(Leg_n_H)))

        self.spx_Leg_n_B = QLineEdit(tbl_data_design_columns)
        self.spx_Leg_n_B.setValidator(QIntValidator(2, 100))
        self.spx_Leg_n_B.setText('{:d}'.format(int(Leg_n_B)))

        self.spx_Sstirrup = QLineEdit(tbl_data_design_columns)
        self.spx_Sstirrup.setValidator(QIntValidator(2, 100))
        self.spx_Sstirrup.setText('{:d}'.format(int(Sstirrup)))

        self.spx_Vu_Vn = QLineEdit(tbl_data_design_columns)
        self.spx_Vu_Vn.setText('{:.2f}'.format(Vu_Vn))


        tbl_data_design_columns.setCellWidget(fila, 0, self.spx_id)
        tbl_data_design_columns.setCellWidget(fila, 1, self.spx_b)
        tbl_data_design_columns.setCellWidget(fila, 2, self.spx_h)
        tbl_data_design_columns.setCellWidget(fila, 3, self.spx_ro)
        tbl_data_design_columns.setCellWidget(fila, 4, self.spx_db)
        tbl_data_design_columns.setCellWidget(fila, 5, self.spx_de)
        tbl_data_design_columns.setCellWidget(fila, 6, self.spx_nbH)
        tbl_data_design_columns.setCellWidget(fila, 7, self.spx_nbB)
        tbl_data_design_columns.setCellWidget(fila, 8, self.spx_Leg_n_H)
        tbl_data_design_columns.setCellWidget(fila, 9, self.spx_Leg_n_B)
        tbl_data_design_columns.setCellWidget(fila, 10, self.spx_Sstirrup)
        tbl_data_design_columns.setCellWidget(fila, 11, self.spx_Vu_Vn)

        tbl_data_design_columns.setColumnWidth(0, 40)
        tbl_data_design_columns.setColumnWidth(1, 40)
        tbl_data_design_columns.setColumnWidth(2, 40)
        tbl_data_design_columns.setColumnWidth(3, 40)
        tbl_data_design_columns.setColumnWidth(4, 60)
        tbl_data_design_columns.setColumnWidth(5, 60)
        tbl_data_design_columns.setColumnWidth(6, 40)
        tbl_data_design_columns.setColumnWidth(7, 40)
        tbl_data_design_columns.setColumnWidth(8, 60)
        tbl_data_design_columns.setColumnWidth(9, 60)
        tbl_data_design_columns.setColumnWidth(10, 60)
        tbl_data_design_columns.setColumnWidth(11, 60)


        stylesheet = "::section{border-style: solid;" \
                     "border-width: 1px;}"
        tbl_data_design_columns.horizontalHeader().setStyleSheet(stylesheet)


class RegistroWalls:
    def __init__(self, tbl_data_design_walls, id_, tw, lw, ro_l, ro_t, db, de, Sst, cMaxS, sigma_c, BE):
        fila = tbl_data_design_walls.rowCount()
        tbl_data_design_walls.insertRow(fila)

        self.spx_id = QLineEdit(tbl_data_design_walls)
        self.spx_id.setValidator(QIntValidator(0, 1000))
        self.spx_id.setText(f'W{id_}')

        self.spx_b = QLineEdit(tbl_data_design_walls)
        self.spx_b.setValidator(QIntValidator(15, 1000))
        self.spx_b.setText('{:d}'.format(int(tw)))

        self.spx_h = QLineEdit(tbl_data_design_walls)
        self.spx_h.setValidator(QIntValidator(20, 2000))
        self.spx_h.setText('{:d}'.format(int(lw)))

        self.spx_ro_l = QLineEdit(tbl_data_design_walls)
        self.spx_ro_l.setValidator(QDoubleValidator(0.1, 8., 3))
        self.spx_ro_l.setText('{:.2f}'.format(ro_l * 100))

        self.spx_ro_t = QLineEdit(tbl_data_design_walls)
        self.spx_ro_t.setValidator(QDoubleValidator(0.1, 8., 3))
        self.spx_ro_t.setText('{:.2f}'.format(ro_t * 100))

        self.spx_db = QLineEdit(tbl_data_design_walls)
        self.spx_db.setValidator(QDoubleValidator(1., 100., 2))
        self.spx_db.setText('{:.2f}'.format(db))

        self.spx_de = QLineEdit(tbl_data_design_walls)
        self.spx_de.setValidator(QDoubleValidator(1., 100., 2))
        self.spx_de.setText('{:.2f}'.format(de))

        self.spx_Sst = QLineEdit(tbl_data_design_walls)
        self.spx_Sst.setValidator(QIntValidator(2, 100))
        self.spx_Sst.setText('{:d}'.format(int(Sst)))

        self.spx_cMaxS = QLineEdit(tbl_data_design_walls)
        # self.spx_cMaxS.setValidator(QDoubleValidator(0.1, 100))
        self.spx_cMaxS.setAlignment(QtCore.Qt.AlignCenter)
        self.spx_cMaxS.setText('{:1.3f}'.format(cMaxS))

        self.spx_sigma_c = QLineEdit(tbl_data_design_walls)
        # self.spx_sigma_c.setValidator(QDoubleValidator(0.01, 100))
        # self.spx_sigma_c.setAlignment(Qt.AlignHCenter)
        self.spx_sigma_c.setAlignment(QtCore.Qt.AlignCenter)
        self.spx_sigma_c.setText('{:1.3f}'.format(sigma_c))

        self.spx_BE = QLineEdit(tbl_data_design_walls)
        # self.spx_id.setValidator(QIntValidator(0, 1000))
        self.spx_BE.setText(BE)

        tbl_data_design_walls.setCellWidget(fila, 0, self.spx_id)
        tbl_data_design_walls.setCellWidget(fila, 1, self.spx_b)
        tbl_data_design_walls.setCellWidget(fila, 2, self.spx_h)
        tbl_data_design_walls.setCellWidget(fila, 3, self.spx_ro_l)
        tbl_data_design_walls.setCellWidget(fila, 4, self.spx_ro_t)
        tbl_data_design_walls.setCellWidget(fila, 5, self.spx_db)
        tbl_data_design_walls.setCellWidget(fila, 6, self.spx_de)
        tbl_data_design_walls.setCellWidget(fila, 7, self.spx_cMaxS)
        tbl_data_design_walls.setCellWidget(fila, 8, self.spx_sigma_c)
        tbl_data_design_walls.setCellWidget(fila, 9, self.spx_BE)

        tbl_data_design_walls.setColumnWidth(0, 40)
        tbl_data_design_walls.setColumnWidth(1, 40)
        tbl_data_design_walls.setColumnWidth(2, 40)
        tbl_data_design_walls.setColumnWidth(3, 40)
        tbl_data_design_walls.setColumnWidth(4, 60)
        tbl_data_design_walls.setColumnWidth(5, 60)
        tbl_data_design_walls.setColumnWidth(6, 40)
        tbl_data_design_walls.setColumnWidth(7, 60)
        tbl_data_design_walls.setColumnWidth(8, 60)
        tbl_data_design_walls.setColumnWidth(9, 60)

        stylesheet = "::section{border-style: solid;" \
                     "border-width: 1px;}"
        tbl_data_design_walls.horizontalHeader().setStyleSheet(stylesheet)


class BeamElasticElement:
    def __init__(self, EleTag, Nod_ini, Nod_end, AEle, EcEle, IzEle, LEle, BEle, HEle, ElegTr, RZi, RZe):
        self.EleTag = EleTag
        self.Nod_ini = Nod_ini
        self.Nod_end = Nod_end
        self.AEle = AEle
        self.EcEle = EcEle
        self.IzEle = IzEle
        self.LEle = LEle
        self.BEle = BEle
        self.HEle = HEle
        self.ElegTr = ElegTr
        self.RZi = RZi
        self.RZe = RZe


class ZeroLengthElement:
    def __init__(self, EleTag, Nod_ini, Nod_end):
        self.EleTag = EleTag
        self.Nod_ini = Nod_ini
        self.Nod_end = Nod_end


class BeamDesign:
    def __init__(self, EleTag, b, h, Ast1, dt1, Mn_n1, Asb1, db1, Mn_p1, ns1, ss1, Ast2, dt2, Mn_n2, Asb2, db2, Mn_p2,
                 ns2, ss2, Nod_ini, Nod_end, db_t1, db_b1, db_t2, db_b2, Vpr, VU1, VU2, Mpr_n1, Mpr_p1, Mpr_n2, Mpr_p2,
                 Vun):
        self.EleTag = EleTag
        self.b = b
        self.h = h
        self.Ast1 = Ast1
        self.dt1 = dt1
        self.Mn_n1 = Mn_n1
        self.Asb1 = Asb1
        self.db1 = db1
        self.Mn_p1 = Mn_p1
        self.ns1 = ns1
        self.ss1 = ss1
        self.Ast2 = Ast2
        self.dt2 = dt2
        self.Mn_n2 = Mn_n2
        self.Asb2 = Asb2
        self.db2 = db2
        self.Mn_p2 = Mn_p2
        self.ns2 = ns2
        self.ss2 = ss2
        self.Nod_ini = Nod_ini
        self.Nod_end = Nod_end
        self.db_t1 = db_t1
        self.db_b1 = db_b1
        self.db_t2 = db_t2
        self.db_b2 = db_b2
        self.Vpr = Vpr
        self.VU1 = VU1
        self.VU2 = VU2
        self.Mpr_n1 = Mpr_n1
        self.Mpr_p1 = Mpr_p1
        self.Mpr_n2 = Mpr_n2
        self.Mpr_p2 = Mpr_p2
        self.Vun = Vun
        
        
class ColDesign:
    def __init__(self, EleTag, b, h, nbH, nbB, db, de, As, Pu_v, Mu_v, fiPn, fiMn, Mn_i, d, dist, ro, Mu_i,
                 sst, nsB, nsH, Nod_ini, Nod_end, NUD1, NUD2, NUG1, NUG2, MUD1, MUD2, VUD1, VUD2, ColBeamStr, Vn, VI6,
                 Vu_Vn, VnCol, sem, sigma_IF, wIF):
        self.EleTag = EleTag
        self.b = b
        self.h = h
        self.nbH = nbH
        self.nbB = nbB
        self.db = db
        self.de = de
        self.As = As
        self.Pu_v = Pu_v
        self.Mu_v = Mu_v
        self.fiPn = fiPn
        self.fiMn = fiMn
        self.Mn_i = Mn_i
        self.d = d
        self.dist = dist
        self.ro = ro
        self.Mu_i = Mu_i
        self.sst = sst
        self.nsB = nsB
        self.nsH = nsH
        self.Nod_ini = Nod_ini
        self.Nod_end = Nod_end
        self.NUD1 = NUD1
        self.NUD2 = NUD2
        self.NUG1 = NUG1
        self.NUG2 = NUG2
        self.MUD1 = MUD1
        self.MUD2 = MUD2
        self.VUD1 = VUD1
        self.VUD2 = VUD2
        self.ColBeamStr = ColBeamStr
        self.Vn = Vn
        self.Vy = VI6
        self.Vu_Vn = Vu_Vn
        self.VnCol = VnCol
        self.sem = sem
        self.sigma_IF = sigma_IF
        self.wIF = wIF


class WallDesign:
    def __init__(self, EleTag, b, h, nbH, nbB, db, dst, As, Pu_v, Mu_v, fiPn, fiMn, Mn_i, d, dist, ro_t, ro_l, Mu_i,
                 sst, Nod_ini, Nod_end, NUD1, NUD2, NUG1, NUG2, MUD1, MUD2, VUD1, VUD2, Vn, cMaxS, sigma_c, BE, lbe,
                 nsH_BE, nsB_BE, se_BE):
        self.EleTag = EleTag
        self.b = b
        self.h = h
        self.nbH = nbH
        self.nbB = nbB
        self.db = db
        self.dst = dst
        self.As = As
        self.Pu_v = Pu_v
        self.Mu_v = Mu_v
        self.fiPn = fiPn
        self.fiMn = fiMn
        self.Mn_i = Mn_i
        self.d = d
        self.dist = dist
        self.ro_t = ro_t
        self.ro_l = ro_l
        self.Mu_i = Mu_i
        self.sst = sst
        self.Nod_ini = Nod_ini
        self.Nod_end = Nod_end
        self.NUD1 = NUD1
        self.NUD2 = NUD2
        self.NUG1 = NUG1
        self.NUG2 = NUG2
        self.MUD1 = MUD1
        self.MUD2 = MUD2
        self.VUD1 = VUD1
        self.VUD2 = VUD2
        self.Vn = Vn
        self.cMaxS = cMaxS
        self.sigma_c = sigma_c
        self.BE = BE
        self.lbe = lbe
        self.nsH_BE = nsH_BE
        self.nsB_BE = nsB_BE
        self.se_BE = se_BE


# class DuctilityCurve:
#     def __init__(self, xi, xe, yi, ye, CD_i, CD_e):
#         self.xi = xi
#         self.xe = xe
#         self.yi = yi
#         self.ye = ye
#         self.CD_i = CD_i
#         self.CD_e = CD_e


class PlasticRotationAngle:
    def __init__(self, xi, xe, yi, ye, PRA_i, PRA_e):
        self.xi = xi
        self.xe = xe
        self.yi = yi
        self.ye = ye
        self.PRA_i = PRA_i
        self.PRA_e = PRA_e


class AcceptanceCriteria:
    def __init__(self, IO_1, LS_1, CP_1, IO_2, LS_2, CP_2):
        self.IO_1 = IO_1
        self.LS_1 = LS_1
        self.CP_1 = CP_1
        self.IO_2 = IO_2
        self.LS_2 = LS_2
        self.CP_2 = CP_2


class PlasticHingeLength:
    def __init__(self, phl1, phl2):
        self.phl1 = phl1
        self.phl2 = phl2


class MyForm(QDialog):
    def __init__(self):
        super().__init__()
        self.ui = Ui_NonLinearFrameDialog()
        # self.setStyleSheet("QLineEdit {border: none}")
        self.ui.setupUi(self)
        self.ui.Design.clicked.connect(self.Design)
        self.ui.AcceptDesign.clicked.connect(self.AcceptDesign)
        self.ui.CreateNLM.clicked.connect(self.CreateNLM)
        self.ui.Pushover.clicked.connect(self.Pushover)
        self.ui.IDA.clicked.connect(self.IDA)
        self.ui.CSS.clicked.connect(self.CSS)
        self.ui.PlotIDA.clicked.connect(self.PlotIDA)
        self.ui.PlotCSS.clicked.connect(self.PlotCSS)
        self.ui.PlotCSSMaxFloors.clicked.connect(self.PlotCSSMaxFloors)
        self.ui.progressBarPushover.hide()
        self.ui.progressBarBeamDesign.hide()
        self.ui.progressBarColumnDesign.hide()
        self.ui.Exit.clicked.connect(self.Exit)
        self.show()

    def Exit(self):
        self.close()

    def Design(self):
        global Loc_span, Loc_heigth, ListNodes, Elements, DataBeamDesign, DataColDesign, WDL, WLL, WDLS, Wtotal, \
            cover, num_elems, Beta1B, Beta1C, fcB, fcC, CodeDesign, DataWallDesign, PWall1SlabD, PWall1SlabL, \
            PWall2SlabD, PWall2SlabL, PColSlabD, PColSlabL, ListNodesW1, ListNodesW2, FrameDetailing, WallDetailing,\
            FN, infill_m, infill_lw

        CodeDesign = self.ui.comboBoxDesignCode.currentText()
        FrameDetailing = self.ui.comboBoxSeismicFrameDesign.currentText()
        WallDetailing = self.ui.comboBoxSeismicWallDesign.currentText()
        if CodeDesign == 'ACI 318-19':
            if FrameDetailing == 'IMF':
                exec(open("Design_ACI_318S_19_IFM.py").read())
            elif FrameDetailing == 'SMF':
                exec(open("Design_ACI_318S_19_SFM.py").read())
            elif FrameDetailing == 'OMF':
                exec(open("Design_ACI_318S_19_OFM.py").read())
        elif CodeDesign == 'CCCSR-84':
            if FrameDetailing == 'IMF':
                exec(open("Design_CCCSR_84_IFM.py").read())
            elif FrameDetailing == 'SMF':
                exec(open("Design_CCCSR_84_SFM.py").read())
            elif FrameDetailing == 'OMF':
                exec(open("Design_CCCSR_84_OFM.py").read())
        elif CodeDesign == 'Gravity Load Design':
            exec(open("Design_GLD.py").read())

    def AcceptDesign(self):

        global T1m, T2m, EleCol, EleBeam, DataBeamPhl, DataColPhl, ListNodesLC, DataColPhl, list_beams, list_cols, \
            ListNodes, num_nodes, PDG_Beams, PDG_Cols, net, EleWall, ListNodesW1, ListNodesW2, PDG_Walls
        # Validation of beam and column design table data

        def validate_data(self):
            cover = 4 * cm
            dst = 3 / 8 * inch
            for (r, DB) in zip(self.registros_beams, DataBeamDesign):
                DB.b = float(r.spx_b.text()) * cm
                DB.h = float(r.spx_h.text()) * cm
                DB.Ast1 = float(r.spx_L_As_top.text()) * (DB.b * DB.dt1) / 100
                DB.Asb1 = float(r.spx_L_As_bot.text()) * (DB.b * DB.db1) / 100
                DB.Ast2 = float(r.spx_R_As_top.text()) * (DB.b * DB.dt2) / 100
                DB.Asb2 = float(r.spx_R_As_bot.text()) * (DB.b * DB.db2) / 100
                DB.ns1 = int(r.spx_L_Leg_n.text())
                DB.ns2 = int(r.spx_R_Leg_n.text())
                DB.ss1 = float(r.spx_L_Sstirrup.text()) * cm
                DB.ss2 = float(r.spx_R_Sstirrup.text()) * cm

            for (r, DC) in zip(self.registros_cols, DataColDesign):
                DC.b = float(r.spx_b.text()) * cm
                DC.h = float(r.spx_h.text()) * cm
                DC.db = float(r.spx_db.text()) * mm
                DC.de = float(r.spx_de.text()) * mm
                DC.nbH = int(r.spx_nbH.text())
                DC.nbB = int(r.spx_nbB.text())
                DC.nsH = int(r.spx_Leg_n_H.text())
                DC.nsB = int(r.spx_Leg_n_B.text())
                DC.sst = float(r.spx_Sstirrup.text()) * cm
                dp = cover + dst + 0.5 * DC.db
                DC.dist = np.linspace(dp, DC.h - dp, DC.nbH)
                Ab = pi * DC.db ** 2. / 4.
                DC.As = np.hstack([DC.nbB * Ab, np.ones(DC.nbH - 2) * 2 * Ab, DC.nbB * Ab])
                DC.sem = float(r.spx_Sstirrup.text()) * cm

            for (r, DW) in zip(self.registros_walls, DataWallDesign):
                DW.b = float(r.spx_b.text()) * cm
                DW.h = float(r.spx_h.text()) * cm
                DW.ro_l = float(r.spx_ro_l.text()) / 100
                # print('DW.ro_l', DW.ro_l)
                DW.ro_t = float(r.spx_ro_t.text()) / 100
                DW.db = float(r.spx_db.text()) * mm
                DW.dst = float(r.spx_de.text()) * mm

                Ag = DW.b * DW.h
                Ab = pi * DW.db ** 2 / 4
                DW.nbH = ceil(DW.ro_l * Ag / Ab / 2) * 2
                DW.nbB = 2
                Ast = pi * DW.dst ** 2 / 4.  # area de la barra del estribo
                Ave = 2 * Ast  # area transversal del estribo
                DW.sst = DW.h / DW.ro_t * Ave / Ag
                dp = cover + DW.dst + 0.5 * DW.db
                DW.dist = np.linspace(dp, DW.h - dp, DW.nbH)
                DW.As = np.hstack([DW.nbB * Ab, np.ones(DW.nbH - 2) * 2 * Ab, DW.nbB * Ab])
        validate_data(self)
        self.ui.tabWidget.setCurrentIndex(2)


    # Creation of the nonlinear model
    def CreateNLM(self):
        global T1m, T2m, EleCol, EleBeam, DataBeamPhl, DataColPhl, ListNodesLC, DataColPhl, list_beams, list_cols, \
            ListNodes, num_nodes, PDG_Beams, PDG_Cols, net, EleWall, ListNodesW1, ListNodesW2, PDG_Walls

        # Function: Parameters of regularized unconfined concrete
        def con_inconf_regu():
            fpc = -fc
            epsc0 = 2 * fpc / Ec
            Gfc = max(2.0 * (-fpc / MPa) * N / mm, 25.0 * N / mm)
            epscu = min(Gfc / (0.6 * fpc * phl) - 0.8 * fpc / Ec + epsc0, 1.5 * epsc0)
            fcu = 0.2 * fpc
            lambdaU = 0.10
            ft = 0.33 * sqrt(-fpc * MPa)
            Ets = ft / 0.002
            # print('fpc, epsc0, fcu, epscu, lambdaU, ft, Ets', fpc, epsc0, fcu, epscu, lambdaU, ft, Ets)
            return fpc, epsc0, fcu, epscu, lambdaU, ft, Ets

        # Function: Parameters of regularized confined concrete
        def con_conf_regu(b, h, nsB, nsH, sst):
            fpc = -fc
            bcx = h - 2. * cover - dst
            bcy = b - 2. * cover - dst
            Asx = nsB * Ast
            Asy = nsH * Ast
            Asvt = Asx + Asy
            flx = Asvt * fys / sst / bcx
            fly = Asvt * fys / sst / bcy
            slx = bcx / (nsB - 1)
            sly = bcy / (nsH - 1)
            k2x = min(0.26 * sqrt((bcx / sst) * (bcx / slx) * (1000. / flx)), 1)
            k2y = min(0.26 * sqrt((bcy / sst) * (bcy / sly) * (1000. / fly)), 1)
            flex = k2x * flx
            fley = k2y * fly
            fle = (flex * bcx + fley * bcy) / (bcx + bcy)
            k1 = 6.7 * (fle / 1000.) ** (-0.17)
            fcc = fc + k1 * fle
            fpcc = -fcc
            Ecc = Ec
            Gfc = max(2.0 * (-fpc / MPa) * N / mm, 25.0 * N / mm)
            # Kp = fcc / fc
            # Gfcc = max(5 * (Kp - 0.85) * Gfc, 2.5 * Gfc)
            K = k1 * fle / fc
            epscc0 = -eo1 * (1. + 5. * K)
            Gfcc = 1.7 * Gfc
            epsccu = min(Gfcc / (0.6 * fpcc * phl) - 0.8 * fpcc / Ecc + epscc0, 1.5 * epscc0)
            fccu = 0.2 * fpcc
            lambdaC = 0.10
            ft = 0.33 * sqrt(-fpc * MPa)
            Ets = ft / 0.002
            # print('fpcc, epscc0, fccu, epsccu, lambdaC, ft, Ets', fpcc, epscc0, fccu, epsccu, lambdaC, ft, Ets)
            return fpcc, epscc0, fccu, epsccu, lambdaC, ft, Ets

        def con_conf_reguW(b, h, nsB, nsH, sst):
            fpc = -fc
            bcx = h - 2. * cover - dst
            bcy = b - 2. * cover - dst
            Asx = nsB * Ast
            Asy = nsH * Ast
            Asvt = Asx + Asy
            flx = Asvt * fys / sst / bcx
            fly = Asvt * fys / sst / bcy
            slx = bcx / (nsB - 1)
            sly = bcy / (nsH - 1)
            k2x = min(0.26 * sqrt((bcx / sst) * (bcx / slx) * (1000. / flx)), 1)
            k2y = min(0.26 * sqrt((bcy / sst) * (bcy / sly) * (1000. / fly)), 1)
            flex = k2x * flx
            fley = k2y * fly
            fle = (flex * bcx + fley * bcy) / (bcx + bcy)
            k1 = 6.7 * (fle / 1000.) ** (-0.17)
            fcc = fc + k1 * fle
            fpcc = -fcc
            Ecc = Ec
            Gfc = max(2.0 * (-fpc / MPa) * N / mm, 25.0 * N / mm)
            Kp = fcc / fc
            Gfcc = max(5 * (Kp - 0.85) * Gfc, 2.5 * Gfc)
            K = k1 * fle / fc
            epscc0 = -eo1 * (1. + 5. * K)
            # Gfcc = 1.7 * Gfc
            epsccu = min(Gfcc / (0.6 * fpcc * phl) - 0.8 * fpcc / Ecc + epscc0, 1.5 * epscc0)
            fccu = 0.2 * fpcc
            lambdaC = 0.10
            ft = 0.33 * sqrt(-fpc * MPa)
            Ets = ft / 0.002
            # print('fpcc, epscc0, fccu, epsccu, lambdaC, ft, Ets', fpcc, epscc0, fccu, epsccu, lambdaC, ft, Ets)
            return fpcc, epscc0, fccu, epsccu, lambdaC, ft, Ets

        # Function: Regularized steel parameters
        def steel_mat_regu():
            FyTestN4 = fy
            FsuTestN4 = 1.25 * FyTestN4
            epsuTestN4 = 0.05
            LgageTestN4 = 200.0 * mm
            Es = 200.0 * GPa
            FyPosN4 = FyTestN4
            epsyPosN4 = FyPosN4 / Es
            FyNegN4 = FyTestN4
            epsyNegN4 = FyNegN4 / Es
            FsuPosN4 = FsuTestN4
            epsuPosN4 = max(epsyPosN4 + LgageTestN4 / phl * (epsuTestN4 - epsyPosN4), epsuTestN4)
            bPosN4 = (FsuPosN4 - FyPosN4) / (Es * (epsuPosN4 - epsyPosN4))
            if CodeDesign == 'ACI 318-19' or CodeDesign == 'CCCSR-84':
                if FrameDetailing == 'IMF' or FrameDetailing == 'SMF':
                    epsuNegN4 = max(-epsccu, 0.005)
                else:
                    epsuNegN4 = 0.005
            if CodeDesign == 'Gravity Load Design':
                epsuNegN4 = 0.005
            if tw1 != 0 or tw2 != 0:
                epsuNegN4 = 0.005
            #     if DW.BE == 'Yes':
            #         epsuNegN4 = max(-epsccu, 0.005)
            #     else:
            #         epsuNegN4 = 0.005
            bNegN4 = bPosN4
            # FsuNegN4 = FsuTestN4
            FsuNegN4 = FyNegN4 + bNegN4 * (Es * (epsuNegN4 - epsyNegN4))
            FsrPosN4 = 0.15 * FyPosN4
            epsrPosN4 = 0.15
            FsrNegN4 = 0.15 * FyNegN4
            if CodeDesign == 'ACI 318-19' or CodeDesign == 'CCCSR-84':
                if FrameDetailing == 'IMF' or FrameDetailing == 'SMF':
                    epsrNegN4 = max(-epsccu,  1.5 * epsuNegN4)
                else:
                    epsrNegN4 = max(-epscu, 1.5 * epsuNegN4)
            if CodeDesign == 'Gravity Load Design':
                epsrNegN4 = max(-epscu, 1.5 * epsuNegN4)
            if tw1 != 0 or tw2 != 0:
                epsrNegN4 = max(-epscu, 1.5 * epsuNegN4)
            #     if DW.BE == 'Yes':
            #         epsrNegN4 = min(-epsccu, 0.02)
            #     else:
            #         epsrNegN4 = -epscu
            pinchX = 0.2
            pinchY = 0.8
            damage1 = 0.0
            damage2 = 0.0
            beta = 0.0
            if self.ui.radioButtonHysteretic.isChecked() == True:
                # print('Hysteretic', Ele.EleTag * 20 + 4 + pos, FyPosN4, epsyPosN4, FsuPosN4, epsuPosN4
                #                     , FsrPosN4, epsrPosN4, -FyNegN4, -epsyNegN4, -FsuNegN4, -epsuNegN4, -FsrNegN4
                #                     , -epsrNegN4, pinchX, pinchY, damage1, damage2, beta)
                op.uniaxialMaterial('Hysteretic', Ele.EleTag * 20 + 4 + pos, FyPosN4, epsyPosN4, FsuPosN4, epsuPosN4
                                    , FsrPosN4, epsrPosN4, -FyNegN4, -epsyNegN4, -FsuNegN4, -epsuNegN4, -FsrNegN4
                                    , -epsrNegN4, pinchX, pinchY, damage1, damage2, beta)
            elif self.ui.radioButtonSteelMPF.isChecked() == True:
                SteelN4Mat = Ele.EleTag * 20 + 4 + pos
                SteelMPFTag = int(1e6 * SteelN4Mat)
                R0 = 20.0
                cR1 = 0.925
                cR2 = 0.15
                a1 = 0.0
                a2 = 1.0
                a3 = 0.0
                a4 = 0.0
                # print(Ele.EleTag, 'SteelMPF', int(SteelMPFTag), FyPosN4/MPa, FyNegN4/MPa, Es/GPa, bPosN4, bNegN4, R0, cR1, cR2, a1, a2, a3, a4)
                op.uniaxialMaterial('SteelMPF', SteelMPFTag, FyPosN4, FyNegN4, Es, bPosN4, bNegN4, R0, cR1, cR2, a1, a2, a3, a4)
                # print('MinMax', int(SteelN4Mat), int(SteelMPFTag), '-min', -epsuNegN4, '-max', epsuPosN4)
                op.uniaxialMaterial('MinMax', SteelN4Mat, SteelMPFTag, '-min', -epsuNegN4, '-max', epsuPosN4)

        # Function: Parameters of non-regularized confined concrete
        def con_conf(b, h, nsB, nsH, sst):
            fpc = -fc
            bcx = h - 2. * cover - dst
            bcy = b - 2. * cover - dst
            Asx = nsB * Ast
            Asy = nsH * Ast
            Asvt = Asx + Asy
            flx = Asvt * fys / sst / bcx
            fly = Asvt * fys / sst / bcy
            slx = bcx / (nsB - 1)
            sly = bcy / (nsH - 1)
            k2x = min(0.26 * sqrt((bcx / sst) * (bcx / slx) * (1000. / flx)), 1)
            k2y = min(0.26 * sqrt((bcy / sst) * (bcy / sly) * (1000. / fly)), 1)
            flex = k2x * flx
            fley = k2y * fly
            fle = (flex * bcx + fley * bcy) / (bcx + bcy)
            k1 = 6.7 * (fle / 1000.) ** (-0.17)
            fcc = -fc + k1 * fle
            fpcc = -fcc
            # print('fpcc/fpc =', fpcc/fpc)
            K = k1 * fle / fc
            epscc0 = -eo1 * (1. + 5. * K)
            rov = Asvt / sst / (bcx + bcy)
            e85 = 260 * rov * epscc0 - eo85
            # print('e85', e85)
            epsccu = (e85 - epscc0) * (0.2 * fcc - fcc) / (0.85 * fcc - fcc) + epscc0
            fccu = 0.2 * fpcc
            lambdaC = 0.10
            ft = 0.33 * sqrt(-fpc * MPa)
            Ets = ft / 0.002
            # print('fpcc, epscc0, fccu, epsccu, lambdaC, ft, Ets', fpcc, epscc0, fccu, epsccu, lambdaC, ft, Ets)
            return fpcc, epscc0, fccu, epsccu, lambdaC, ft, Ets

        # Function: Parameters of non-regularized steel
        def steel_mat():
            FyTestN4 = fy
            FsuTestN4 = 1.25 * FyTestN4
            epsuTestN4 = 0.10
            LgageTestN4 = phl
            Es = 200.0 * GPa
            FyPosN4 = FyTestN4
            epsyPosN4 = FyPosN4 / Es
            FyNegN4 = FyTestN4
            epsyNegN4 = FyNegN4 / Es
            FsuPosN4 = FsuTestN4
            epsuPosN4 = epsyPosN4 + LgageTestN4 / phl * (epsuTestN4 - epsyPosN4)
            bPosN4 = (FsuPosN4 - FyPosN4) / (Es * (epsuPosN4 - epsyPosN4))
            if CodeDesign == 'ACI 318-19' or CodeDesign == 'CCCSR-84':
                if FrameDetailing == 'IMF' or FrameDetailing == 'SMF':
                    epsuNegN4 = max(-epscc0, 0.003)
                else:
                    epsuNegN4 = 0.003
            if CodeDesign == 'Gravity Load Design':
                epsuNegN4 = 0.003
            bNegN4 = bPosN4
            FsuNegN4 = FyNegN4 + bNegN4 * (Es * (epsuNegN4 - epsyNegN4))
            FsrPosN4 = 0.15 * FyPosN4
            epsrPosN4 = 0.15
            FsrNegN4 = 0.15 * FsuNegN4
            if CodeDesign == 'ACI 318-19' or CodeDesign == 'CCCSR-84':
                if FrameDetailing == 'IMF' or FrameDetailing == 'SMF':
                    epsrNegN4 = min(-epsccu, 0.02)
                else:
                    epsrNegN4 = eo20
            if CodeDesign == 'Gravity Load Design':
                epsrNegN4 = eo20
            pinchX = 0.2
            pinchY = 0.8
            damage1 = 0.0
            damage2 = 0.0
            beta = 0.0
            if self.ui.radioButtonHysteretic.isChecked() == True:
                # print('Hysteretic', Ele.EleTag * 20 + 4 + pos, FyPosN4, epsyPosN4, FsuPosN4, epsuPosN4
                #                     , FsrPosN4, epsrPosN4, -FyNegN4, -epsyNegN4, -FsuNegN4, -epsuNegN4, -FsrNegN4
                #                     , -epsrNegN4, pinchX, pinchY, damage1, damage2, beta)
                op.uniaxialMaterial('Hysteretic', Ele.EleTag * 20 + 4 + pos, FyPosN4, epsyPosN4, FsuPosN4, epsuPosN4
                                    , FsrPosN4, epsrPosN4, -FyNegN4, -epsyNegN4, -FsuNegN4, -epsuNegN4, -FsrNegN4
                                    , -epsrNegN4, pinchX, pinchY, damage1, damage2, beta)
            elif self.ui.radioButtonSteelMPF.isChecked() == True:
                SteelN4Mat = Ele.EleTag * 20 + 4 + pos
                SteelMPFTag = 1e6 * SteelN4Mat
                R0 = 20.0
                cR1 = 0.925
                cR2 = 0.15
                a1 = 0.0
                a2 = 1.0
                a3 = 0.0
                a4 = 0.0
                # print('SteelMPF', int(SteelMPFTag), FyPosN4, FyNegN4, Es, bPosN4, bNegN4, R0, cR1, cR2, a1, a2, a3, a4)
                op.uniaxialMaterial('SteelMPF', SteelMPFTag, FyPosN4, FyNegN4, Es, bPosN4, bNegN4, R0, cR1, cR2, a1, a2,
                                    a3,
                                    a4)
                # print('MinMax', int(SteelN4Mat), int(SteelMPFTag), '-min', -epsuNegN4, '-max', epsuPosN4)
                op.uniaxialMaterial('MinMax', SteelN4Mat, SteelMPFTag, '-min', -epsuNegN4, '-max', epsuPosN4)

        # Function: Creation of fibers in beams
        def fiber_beam(Ast, Asb, pos):
            op.section('Fiber', Ele.EleTag * 2 + pos)
            op.patch('rect', Ele.EleTag * 20 + 2 + pos, NFBeams, 1, -y2 + dp, -z2 + dp1, y2 - dp, z2 - dp1)
            op.patch('rect', Ele.EleTag * 20 + pos, NFBeams, 1, -y2 + dp, z2 - dp1, y2 - dp, z2)
            op.patch('rect', Ele.EleTag * 20 + pos, NFBeams, 1, -y2 + dp, -z2, y2 - dp, -z2 + dp1)
            op.patch('rect', Ele.EleTag * 20 + pos, 2, 1, -y2, -z2, -y2 + dp, z2)
            op.patch('rect', Ele.EleTag * 20 + pos, 2, 1, y2 - dp, -z2, y2, z2)
            # print(Ele.EleTag * 20 + 4 + pos, 1, Ast/cm**2, y2 - dp, z2 - dp, y2 - dp, -z2 + dp)
            op.layer('straight', Ele.EleTag * 20 + 4 + pos, 1, Ast, y2 - dp, z2 - dp1, y2 - dp, -z2 + dp1)
            # print(Ele.EleTag * 20 + 4 + pos, 1, Asb/cm**2, -y2 + dp, z2 - dp, -y2 + dp, -z2 + dp)
            op.layer('straight', Ele.EleTag * 20 + 4 + pos, 1, Asb, -y2 + dp, z2 - dp1, -y2 + dp, -z2 + dp1)

        op.wipe()  # The models is restarted in opensees
        op.model('Basic', '-ndm', 2, '-ndf', 3)
        ListNodesLC = np.empty([len(Loc_heigth), 3])
        nodeTag = nnt + 1
        indLC = 0
        for node in ListNodes:
            # print(int(node[0]), node[1], node[2])
            op.node(int(node[0]), node[1], node[2])
            if node[2] == 0.:
                op.fix(int(node[0]), 1, 1, 1)
            if node[1] == Loc_span[-1] and self.ui.radioButtonYesLC.isChecked() == True:
                xn, yn = node[1] + 1 * m, node[2]
                op.node(nodeTag, xn, yn)
                if node[2] == 0.:
                    op.fix(nodeTag, 1, 1, 0)
                    # print(nodeTag, 1, 1, 0)
                else:
                    # op.fix(nodeTag, 0, 0, 1)
                    op.fix(nodeTag, 0, 0, 0)
                ListNodesLC[indLC, :] = [nodeTag, xn, yn]
                indLC += 1
                nodeTag += 1
        nntot = nodeTag - 1
        if np.any(ListNodesW1):
            for node in ListNodesW1:
                # print(int(node[0]), node[1], node[2])
                op.node(int(node[0]), node[1], node[2])
                if node[2] == 0.:
                    op.fix(int(node[0]), 1, 1, 1)
        if np.any(ListNodesW2):
            for node in ListNodesW2:
                op.node(int(node[0]), node[1], node[2])
                if node[2] == 0.:
                    op.fix(int(node[0]), 1, 1, 1)
        if np.any(ListNodesW1):
            for node in ListNodesW1[:-1, :]:
                nodeTag = int(node[0] + nntot)
                # print(nodeTag, node[1], node[2])
                op.node(nodeTag, node[1], node[2])
        if np.any(ListNodesW2):
            for node in ListNodesW2[:-1, :]:
                nodeTag = int(node[0] + nntot)
                op.node(nodeTag, node[1], node[2])
        # Non-linear spring nodes
        for node in ListNodes:
            if node[2] == 0.:
                op.node(int(node[0]) + nntot * 2, node[1], node[2])
            if node[1] == 0. and node[2] != Loc_heigth[-1] and node[2] != 0.:
                # op.node(int(node[0]) + nntot * 1, node[1], node[2])
                op.node(int(node[0]) + nntot * 2, node[1], node[2])
                op.node(int(node[0]) + nntot * 3, node[1], node[2])
                # op.node(int(node[0]) + nntot * 4, node[1], node[2])
            elif node[1] == Loc_span[-1] and node[2] != Loc_heigth[-1] and node[2] != 0.:
                # op.node(int(node[0]) + nntot * 1, node[1], node[2])
                op.node(int(node[0]) + nntot * 2, node[1], node[2])
                op.node(int(node[0]) + nntot * 3, node[1], node[2])
                # op.node(int(node[0]) + nntot * 5, node[1], node[2])
            elif node[1] == 0. and node[2] == Loc_heigth[-1] and node[2] != 0.:
                # op.node(int(node[0]) + nntot * 1, node[1], node[2])
                op.node(int(node[0]) + nntot * 3, node[1], node[2])
                # op.node(int(node[0]) + nntot * 4, node[1], node[2])
            elif node[1] == Loc_span[-1] and node[2] == Loc_heigth[-1] and node[2] != 0.:
                # op.node(int(node[0]) + nntot * 1, node[1], node[2])
                op.node(int(node[0]) + nntot * 3, node[1], node[2])
                # op.node(int(node[0]) + nntot * 5, node[1], node[2])
            elif node[2] == Loc_heigth[-1] and node[2] != 0.:
                # op.node(int(node[0]) + nntot * 1, node[1], node[2])
                op.node(int(node[0]) + nntot * 3, node[1], node[2])
                # op.node(int(node[0]) + nntot * 4, node[1], node[2])
                # op.node(int(node[0]) + nntot * 5, node[1], node[2])
            elif node[2] != 0.:
                # op.node(int(node[0]) + nntot * 1, node[1], node[2])
                op.node(int(node[0]) + nntot * 2, node[1], node[2])
                op.node(int(node[0]) + nntot * 3, node[1], node[2])
                # op.node(int(node[0]) + nntot * 4, node[1], node[2])
                # op.node(int(node[0]) + nntot * 5, node[1], node[2])

        # print('ListNodes', ListNodes)
        # if self.ui.radioButtonYesLC.isChecked() == True:
        # print('ListNodesLC', ListNodesLC)
        cover = 4 * cm
        dst = 3 / 8 * inch
        Ast = pi * dst ** 2 / 4.  # area de la barra del estribo
        FN = float(self.ui.frames_numbers.text())
        # creacion de columnas
        HBeam = float(self.ui.HBeam.text())
        HColi = float(self.ui.HColi.text())  # Column inside Depth            # if node[2] > 0 and node[1] == 0:
        #     MasterNode = node[0]
        # if node[2] > 0 and node[1] != 0:
        #     op.equalDOF(int(MasterNode), int(node[0]), 1)

        HCole = float(self.ui.HCole.text())  # Column outside Depth
        fy = float(self.ui.fy.text()) * MPa
        fys = float(self.ui.fys.text()) * MPa
        fy, fys = 1.1 * fy, 1.1 * fys
        Es = 200.0 * GPa
        Cd = float(self.ui.Cd.text())
        fcB = float(self.ui.fcB.text()) * MPa
        fcC = float(self.ui.fcC.text()) * MPa
        tw1, tw2 = float(self.ui.tw1.text()), float(self.ui.tw2.text())
        lw1, lw2 = float(self.ui.lw1.text()), float(self.ui.lw2.text())
        NFBeams = int(self.ui.nFBeamsLineEdit.text())
        NFCols = int(self.ui.nFColsLineEdit.text())
        NFWalls = int(self.ui.nFWallsLineEdit.text())
        NFMVLEM = int(self.ui.nFMVLEM.value())
        Lafg = float(self.ui.Lafg.text()) * FN
        Lafs = float(self.ui.Lafs.text()) * FN

        twIF = float(self.ui.twIF.text()) * mm
        CASM = float(self.ui.CASM.text()) / 100  #Central strut area percentage
        fmIF = float(self.ui.fmIF.text()) * MPa
        EmIF = float(self.ui.EmIF.text()) * MPa
        CIF = float(self.ui.CohesionIF.text()) * MPa
        fmtetaIF = float(self.ui.fmtetaIF.text()) * MPa

        if self.ui.radioButtonJS1.isChecked():
            op.geomTransf('PDelta', 1, '-jntOffset', 0, 0, 0, -HBeam / 2)
            op.geomTransf('PDelta', 2, '-jntOffset', 0, HBeam / 2, 0, -HBeam / 2)
            op.geomTransf('Corotational', 3, '-jntOffset', HColi / 2., 0, -HColi / 2., 0)
            op.geomTransf('Corotational', 4, '-jntOffset', HCole / 2., 0, -HColi / 2., 0)
            op.geomTransf('Corotational', 5, '-jntOffset', HColi / 2., 0, -HCole / 2., 0)
        elif self.ui.radioButtonJS2.isChecked():
            op.geomTransf('PDelta', 1, '-jntOffset', 0, 0, 0, -HBeam / 2)
            op.geomTransf('PDelta', 2, '-jntOffset', 0, HBeam / 2, 0, -HBeam / 2)
            op.geomTransf('Corotational', 3, '-jntOffset', 0, 0, 0, 0)
            op.geomTransf('Corotational', 4, '-jntOffset', 0, 0, 0, 0)
            op.geomTransf('Corotational', 5, '-jntOffset', 0, 0, 0, 0)
        elif self.ui.radioButtonJS3.isChecked():
            op.geomTransf('PDelta', 1, '-jntOffset', 0, 0, 0, 0)
            op.geomTransf('PDelta', 2, '-jntOffset', 0, 0, 0, 0)
            op.geomTransf('Corotational', 3, '-jntOffset', HColi / 2., 0, -HColi / 2., 0)
            op.geomTransf('Corotational', 4, '-jntOffset', HCole / 2., 0, -HColi / 2., 0)
            op.geomTransf('Corotational', 5, '-jntOffset', HColi / 2., 0, -HCole / 2., 0)
        elif self.ui.radioButtonJS4.isChecked():
            op.geomTransf('PDelta', 1, '-jntOffset', 0, 0, 0, -HBeam / 4)
            op.geomTransf('PDelta', 2, '-jntOffset', 0, HBeam / 4, 0, -HBeam / 4)
            op.geomTransf('Corotational', 3, '-jntOffset', HColi / 4., 0, -HColi / 4., 0)
            op.geomTransf('Corotational', 4, '-jntOffset', HCole / 4., 0, -HColi / 4., 0)
            op.geomTransf('Corotational', 5, '-jntOffset', HColi / 4., 0, -HCole / 4., 0)
        op.geomTransf('Corotational', 6)

        num_nodes = len(Loc_span) * len(Loc_heigth)
        if self.ui.radioButtonYesLC.isChecked():
            EleTag = net + 1
            Iner_rat = float(self.ui.InertiaRatio.text())  # Column outside Depth
            fcC = float(self.ui.fcC.text()) * MPa
            HColi = float(self.ui.HColi.text())  # Column inside Depth
            BColi = float(self.ui.BColi.text())  # Column inside Width
            IFC = float(self.ui.InertiaColumnsFactor.text())
            EcC = 4700 * sqrt(fcC * MPa)
            IzColi = 1. / 12. * BColi * HColi ** 3  # Column moment of inertia
            IzCol = IFC * IzColi
            for Nod_ini in range(int(ListNodesLC[0, 0]), int(ListNodesLC[-1, 0])):
                Nod_end = Nod_ini + 1
                # op.element("corotTruss", EleTag, Nod_ini, Nod_end, 1e3, int(1e4))
                # print("corotTruss", EleTag, Nod_ini, Nod_end, 1e3, int(1e4))
                op.element('elasticBeamColumn', EleTag, Nod_ini, Nod_end, 1e3, EcC, Iner_rat * IzCol, 6)
                # print('elasticBeamColumn', EleTag, Nod_ini, Nod_end, 1e3, 1e3, 1e+1/2, 6)
                EleTag += 1
            for Nod_ini in range(num_nodes):
                if ListNodes[Nod_ini, 1] == Loc_span[-1] and ListNodes[Nod_ini, 2] != 0.:
                    for indLC in range(ListNodesLC[:, 0].size):
                        if ListNodesLC[indLC, 2] == ListNodes[Nod_ini, 2]:
                            Nod_end = int(ListNodesLC[indLC, 0])
                            op.element("corotTruss", EleTag, Nod_ini, Nod_end, 1e3, int(1e4))
                            # print("corotTruss", EleTag, Nod_ini, Nod_end, 1e3, int(1e4))
                            EleTag += 1
            net = EleTag - 1

        EleCol, EleBeam, EleWall = [], [], []
        for Ele in Elements:
            if Ele.EleTag in ListEleTagCols:
                EleCol.append(Ele)
            if Ele.EleTag in ListEleTagBeams:
                EleBeam.append(Ele)
            if Ele.EleTag in ListEleTagW1 or Ele.EleTag in ListEleTagW2:
                EleWall.append(Ele)

        # Creation of non-linear elements (beams and columns)
        conc_fact = 1.2
        SeismicLoadCode = self.ui.comboBoxSeismicLoadCode.currentText()
        if SeismicLoadCode == 'CCCSR-84' or SeismicLoadCode == 'Weight Percentage':
            conc_fact = 1.0
        fcC, fcB = conc_fact * fcC, conc_fact * fcB
        eo1, eo85, eo20, lambdaU = 0.002, 0.0038, 0.006, 0.1
        # Creation of non-linear columns
        Krigid = 1e15
        DataColPhl = []
        for (Ele, DC, IFE, lwIF) in zip(EleCol, DataColDesign, infill_m, infill_lw):
            fc, Ec = fcC, Ele.EcEle
            if self.ui.radioButton05H.isChecked():
                phl = 0.5 * DC.h
            elif self.ui.radioButtonPark.isChecked():
                phl = 0.08 * Ele.LEle + 0.022 * fy / MPa * DC.db
            elif self.ui.radioButtonBerry.isChecked():
                phl = 0.05 * Ele.LEle + 0.1 * fy / MPa * DC.db / sqrt(fc * MPa)
            DataColPhl.append(PlasticHingeLength(phl, phl))
            if self.ui.checkBox_ReguConcrete.isChecked() == True:
                fpc, epsc0, fcu, epscu, lambdaU, ft, Ets = con_inconf_regu()
                # print('fpc, epsc0, fcu, epscu, lambdaU, ft, Ets', fpc/MPa, epsc0, fcu/MPa, epscu, lambdaU, ft/MPa, Ets/MPa)
                op.uniaxialMaterial('Concrete02', Ele.EleTag * 20, fpc, epsc0, fcu, epscu, lambdaU, ft, Ets)
                op.uniaxialMaterial('Concrete02', Ele.EleTag * 20 + 1, fpc, epsc0, fcu, epscu, lambdaU, ft, Ets)
                if CodeDesign == 'ACI 318-19' or CodeDesign == 'CCCSR-84' or CodeDesign == 'Gravity Load Design':
                    if FrameDetailing == 'IMF' or FrameDetailing == 'SMF':
                        fpcc, epscc0, fccu, epsccu, lambdaC, ft, Ets = con_conf_regu(DC.b, DC.h, DC.nsB, DC.nsH, DC.sst)
                        # print('fpcc, epscc0, fccu, epsccu, lambdaC, ft, Ets', fpcc/MPa, epscc0, fccu/MPa, epsccu, lambdaC, ft/MPa, Ets/MPa)
                        op.uniaxialMaterial('Concrete02', Ele.EleTag * 20 + 2, fpcc, epscc0, fccu, epsccu, lambdaC, ft,
                                            Ets)
                        op.uniaxialMaterial('Concrete02', Ele.EleTag * 20 + 3, fpcc, epscc0, fccu, epsccu, lambdaC, ft,
                                            Ets)
                    else:
                        op.uniaxialMaterial('Concrete02', Ele.EleTag * 20 + 2, fpc, epsc0, fcu, epscu, lambdaU, ft, Ets)
                        op.uniaxialMaterial('Concrete02', Ele.EleTag * 20 + 3, fpc, epsc0, fcu, epscu, lambdaU, ft, Ets)
                        epsccu = epscu
            if self.ui.checkBox_ReguSteel.isChecked() == True:
                pos = 0
                steel_mat_regu()
                pos = 1
                steel_mat_regu()
            if self.ui.checkBox_ReguConcrete.isChecked() == False:
                ft = 0.33 * sqrt(fcC * MPa)
                Ets = ft / 0.002
                # print('Concrete02', Ele.EleTag * 20, -fcC, eo1, -0.2 * fcC, eo20, lambdaU, ft, Ets)
                op.uniaxialMaterial('Concrete02', Ele.EleTag * 20, -fcC, -eo1, -0.2 * fcC, -eo20, lambdaU, ft, Ets)
                op.uniaxialMaterial('Concrete02', Ele.EleTag * 20 + 1, -fcC, -eo1, -0.2 * fcC, -eo20, lambdaU, ft, Ets)
                if CodeDesign == 'ACI 318-19' or CodeDesign == 'CCCSR-84' or CodeDesign == 'Gravity Load Design':
                    if FrameDetailing == 'IMF' or FrameDetailing == 'SMF':
                        fpcc, epscc0, fccu, epsccu, lambdaC, ft, Ets = con_conf(DC.b, DC.h, DC.nsB, DC.nsH, DC.sst)
                        op.uniaxialMaterial('Concrete02', Ele.EleTag * 20 + 2, fpcc, epscc0, fccu, epsccu, lambdaC, ft,
                                            Ets)
                        op.uniaxialMaterial('Concrete02', Ele.EleTag * 20 + 3, fpcc, epscc0, fccu, epsccu, lambdaC, ft,
                                            Ets)
                    else:
                        op.uniaxialMaterial('Concrete02', Ele.EleTag * 20 + 2, -fcC, -eo1, -0.2 * fcC, -eo20, lambdaU,
                                            ft, Ets)
                        op.uniaxialMaterial('Concrete02', Ele.EleTag * 20 + 3, -fcC, -eo1, -0.2 * fcC, -eo20, lambdaU,
                                            ft, Ets)
                        epsccu = eo20
            if self.ui.checkBox_ReguSteel.isChecked() == False:
                pos = 0
                steel_mat()
                pos = 1
                steel_mat()
            dp = DC.dist[0]
            dp1 = dp * FN
            y1 = DC.h / 2.0
            z1 = DC.b / 2.0 * FN
            op.section('Fiber', Ele.EleTag * 2)
            op.patch('rect', Ele.EleTag * 20 + 2, NFCols, 1, -y1 + dp, -z1 + dp1, y1 - dp, z1 - dp1)
            op.patch('rect', Ele.EleTag * 20, NFCols, 1, -y1 + dp, z1 - dp1, y1 - dp, z1)
            op.patch('rect', Ele.EleTag * 20, NFCols, 1, -y1 + dp, -z1, y1 - dp, -z1 + dp1)
            op.patch('rect', Ele.EleTag * 20, 2, 1, -y1, -z1, -y1 + dp, z1)
            op.patch('rect', Ele.EleTag * 20, 2, 1, y1 - dp, -z1, y1, z1)
            for dist, As in zip(DC.dist, DC.As):
                op.layer('straight', Ele.EleTag * 20 + 4, 1, As * FN, -y1 + dist, z1 - dp1, -y1 + dist, -z1 + dp1)
            MassDens = Ele.AEle * GConc / g * FN
            Nod_ini = Ele.Nod_ini + nntot * 2
            Nod_end = Ele.Nod_end + nntot * 3
            op.beamIntegration('HingeRadau', Ele.EleTag, Ele.EleTag * 2, phl, Ele.EleTag * 2, phl, Ele.EleTag * 2)
            op.element('forceBeamColumn', Ele.EleTag * 10, Nod_ini, Nod_end, Ele.ElegTr, Ele.EleTag, '-mass', MassDens)
            EleTag = int(Ele.EleTag * 10 + 2)
            Nod_ini = Ele.Nod_end + nntot * 3
            Nod_end = Ele.Nod_end
            MatTagX, MatTagY, MatTagZ = EleTag * 20 + 7, EleTag * 20 + 9, EleTag * 20 + 11
            op.uniaxialMaterial('Elastic', MatTagX, Krigid)
            if Cd <= 2.0:
                knl = 1.0
            elif Cd >= 6.0:
                knl = 0.7
            else:
                knl = np.interp(Cd, [2.0, 6.0], [1.0, 0.7])
            if DC.sem / DC.d <= 0.75:
                alfa_col = 1.0
            elif DC.sem / DC.d >= 1.0:
                alfa_col = 0.0
            else:
                alfa_col = np.interp(DC.sem / DC.d, [0.75, 1.0], [1.0, 0.0])
            Ast = pi * DC.de ** 2 / 4
            Linf = Ele.LEle / 2
            VnCol = knl * (alfa_col * Ast * DC.nsB * fys * DC.d / DC.sem +
                           0.5 * sqrt(fcC * MPa) / (min(max(2, Linf / DC.d), 4)) *
                           sqrt(1 + DC.NUG1 / (0.5 * DC.b * DC.h * sqrt(fcC * MPa))) * 0.8 * DC.b * DC.h)

            # VnCol = knl * (alfa_col * Ast * DC.nsB * fys * DC.d / DC.sem +
            #                0.5 * sqrt(fcC * MPa) / (min(max(2, DC.MUD1 / DC.VUD1), 4) * DC.d) *
            #                sqrt(1 + DC.NUD1 / (0.5 * DC.b * DC.h * sqrt(fcC * MPa))) * 0.8 * DC.b * DC.h)
            DC.VnCol = VnCol
            if self.ui.radioButtonSMFNone.isChecked():
                Ks = Krigid
                op.uniaxialMaterial('Elastic', MatTagY, Ks)
            elif self.ui.radioButtonSMFLinear.isChecked():
                Ks = 0.4 * Ele.EcEle * DC.b * DC.d / 1.2 / Ele.LEle
                op.uniaxialMaterial('Elastic', MatTagY, Ks)
            elif self.ui.radioButtonSMFNonlinear.isChecked():
                fcr = 1 / 3 * sqrt(fcC * MPa)
                Vcr = fcr / (Ele.LEle / DC.h) * sqrt(1 + DC.NUD1 / fcr / Ele.AEle) * 0.8 * Ele.AEle
                # print('Vcr', Vcr)
                Gc = 0.42 * Ele.EcEle
                GA0 = 0.8 * Ele.AEle * Gc
                Gamma_cr = Vcr / GA0
                d_cr = Gamma_cr * Ele.LEle
                Gamma_vn = VnCol / GA0
                # VnCol = 0.17 * sqrt(fcC * MPa) * 0.8 * DC.b * DC.h
                # print('VnCol', VnCol)
                # Vultm = VnCol + (0.0015 - Gamma_vn) * 0.10 * GA0
                # Vcolap = Vultm + (0.15-0.0015) * 0.01 * GA0
                # Vcolap = 1.01 * Vultm
                f1 = VnCol
                f2 = 1.03 * f1
                f3 = 0.1 * f2
                # f1 = 0.1*f1
                VcACI = min((0.17 * sqrt(fcC * MPa) + min(DC.NUG1 / 6 / Ele.AEle, 0.05 * fcC)) * DC.b * DC.d,
                            0.42 * sqrt(fcC * MPa) * DC.b * DC.d)
                VnACI = VcACI + Ast * DC.nsB * fys * DC.d / DC.sem
                d1 = Gamma_vn * Ele.LEle / 2
                d2 = (4 - 12 * (VnACI/(DC.b * DC.d))/fcC) * d1
                d3 = 2 * d2
                op.uniaxialMaterial('Hysteretic', MatTagY, f1, d1, f2, d2, f3, d3, -f1, -d1, -f2, -d2, -f3, -d3, 0.8, 0.2, 0.0, 0.02, 0.0)
                # op.uniaxialMaterial('Hysteretic', MatTagY, f1, d1, f2, d2, -f1, -d1, -f2, -d2, 0.8, 0.2, 0.0, 0.02, 0.0)
                # op.uniaxialMaterial('Hysteretic', MatTagY, f1, d1, 1.01*f1, d3, -f1, -d1, -1.01*f1, -d3, 0.8, 0.2, 0.0, 0.02, 0.0)
            op.uniaxialMaterial('Elastic', MatTagZ, Krigid)
            # print('zeroLength', EleTag, Nod_ini, int(Ele.Nod_ini), '-mat', MatTagX, MatTagY, MatTagZ, '-dir', 1, 2, 3)
            op.element('zeroLength', EleTag, Nod_ini, Nod_end, '-mat', MatTagY, MatTagX, MatTagZ, '-dir', 1, 2, 3)
            EleTag = int(Ele.EleTag * 10 + 1)
            Nod_ini = Ele.Nod_ini
            Nod_end = Ele.Nod_ini + nntot * 2
            op.element('zeroLength', EleTag, Nod_ini, Nod_end, '-mat', MatTagY, MatTagX, MatTagZ, '-dir', 1, 2, 3)

            if IFE != 0:
                EleTag = int(Ele.EleTag * 10 + 4)
                MatTag = EleTag * 20 + 13
                hwIF = Ele.LEle
                if ListNodes[Ele.Nod_ini, 1] == 0. or ListNodes[Ele.Nod_ini, 1] == Loc_span[-1]:
                    lwIF = lwIF - HColi / 2 - HCole / 2
                else:
                    lwIF = lwIF - HColi
                ldIF = sqrt(lwIF**2 + hwIF**2)
                wIF = 0.25 * ldIF
                ADIF = wIF * twIF
                if self.ui.radioButtonIMWASCENonlinear.isChecked():
                    VmaxIF = CIF * lwIF * twIF
                    PmaxIF = VmaxIF / (lwIF / ldIF)
                    sigma_max = PmaxIF / ADIF
                elif self.ui.radioButtonIMWManualNonlinear.isChecked():
                    sigma_max = fmtetaIF
                DC.sigma_IF = sigma_max
                DC.wIF = wIF
                fpc = -sigma_max
                fpcu = 0.2 * fpc
                epsc0, epsU = -0.002, -0.0032
                op.uniaxialMaterial('Concrete01', MatTag, fpc, epsc0, fpcu, epsU)
                print('CASM', CASM)
                if CASM < 1:
                    ADIFS = (1-CASM) * ADIF
                    Nod_ini = Ele.Nod_end + nntot * 3
                    Nod_end = (Ele.Nod_ini + 1) + nntot * 2
                    op.element("corotTruss", EleTag + 1, Nod_ini, Nod_end, ADIFS, MatTag)
                    Nod_ini = Ele.Nod_ini + nntot * 2
                    Nod_end = (Ele.Nod_end + 1) + nntot * 3
                    op.element("corotTruss", EleTag + 3, Nod_ini, Nod_end, ADIFS, MatTag)
                ADIFC = CASM * ADIF
                Nod_ini = Ele.Nod_end
                Nod_end = Ele.Nod_ini + 1
                op.element("corotTruss", EleTag, Nod_ini, Nod_end, ADIFC, MatTag)
                Nod_ini = Ele.Nod_ini
                Nod_end = Ele.Nod_end + 1
                op.element("corotTruss", EleTag + 2, Nod_ini, Nod_end, ADIFC, MatTag)

        op.uniaxialMaterial("Elastic", int(1e4), 200 * GPa)

        DataWallPhl = []
        if tw1 != 0 and tw2 != 0:
            ListNodesW = np.vstack((ListNodesW1, ListNodesW2))
        elif tw1 == 0 and tw2 != 0:
            ListNodesW = ListNodesW2
        else:
            ListNodesW = ListNodesW1
        indW = 0
        for (Ele, DW) in zip(EleWall, DataWallDesign):
            fc, Ec = fcC, Ele.EcEle
            phl1 = 0.2 * DW.h + 0.044 * (DW.MUD1 / DW.VUD1)
            phl2 = min((0.2 * DW.h + 0.05 * (DW.MUD1 / DW.VUD1)) * (1 - 1.5 * DW.NUD1 / fc / DW.b / DW.h), 0.8 * DW.h)
            phl3 = 0.27 * DW.h * (1 - DW.NUD1 / fc / DW.b / DW.h) * (1 - fy * DW.ro_t / fc) * (
                    (DW.MUD1 / DW.VUD1) / DW.h) ** 0.45
            if self.ui.radioButtonPriestley.isChecked():
                phl = phl1
            elif self.ui.radioButtonBohl.isChecked():
                phl = phl2
            elif self.ui.radioButtonKazaz.isChecked():
                phl = phl3
            elif self.ui.radioButtonMean.isChecked():
                phl = np.mean([phl1, phl2, phl3])
            if ListNodesW[indW, 2] == 0:
                phl = 0.111 * Ele.LEle
            else:
                phl = 0.35 * Ele.LEle
            if self.ui.radioButtonMVLEM.isChecked():
                phl = Ele.LEle
            DataWallPhl.append(PlasticHingeLength(phl, phl))
            if self.ui.checkBox_ReguConcrete.isChecked() == True:
                fpc, epsc0, fcu, epscu, lambdaU, ft, Ets = con_inconf_regu()
                # print('fpc, epsc0, fcu, epscu, lambdaU, ft, Ets', fpc/MPa, epsc0, fcu/MPa, epscu, lambdaU, ft/MPa, Ets/MPa)
                op.uniaxialMaterial('Concrete02', Ele.EleTag * 20, fpc, epsc0, fcu, epscu, lambdaU, ft, Ets)
                op.uniaxialMaterial('Concrete02', Ele.EleTag * 20 + 1, fpc, epsc0, fcu, epscu, lambdaU, ft, Ets)
                if CodeDesign == 'ACI 318-19' or CodeDesign == 'CCCSR-84' or CodeDesign == 'Gravity Load Design':
                    if WallDetailing == 'IMF' or WallDetailing == 'SMF':
                        fpcc, epscc0, fccu, epsccu, lambdaC, ft, Ets = con_conf_reguW(DW.b, DW.lbe, DW.nsB_BE,
                                                                                      DW.nsH_BE, DW.se_BE)
                        if DW.BE == 'Yes':
                            op.uniaxialMaterial('Concrete02', Ele.EleTag * 20 + 2, fpcc, epscc0, fccu, epsccu, lambdaC,
                                                ft, Ets)
                            op.uniaxialMaterial('Concrete02', Ele.EleTag * 20 + 3, fpcc, epscc0, fccu, epsccu, lambdaC,
                                                ft, Ets)
                            # print('Concrete02', Ele.EleTag * 20 + 3, fpcc, epscc0, fccu, epsccu, lambdaC, ft, Ets)
                        else:
                            op.uniaxialMaterial('Concrete02', Ele.EleTag * 20 + 2, fpc, epsc0, fcu, epscu, lambdaU, ft,
                                                Ets)
                            op.uniaxialMaterial('Concrete02', Ele.EleTag * 20 + 3, fpc, epsc0, fcu, epscu, lambdaU, ft,
                                                Ets)
                            epsccu = epscu
                            # print('Concrete02', Ele.EleTag * 20 + 3, fpc, epsc0, fcu, epscu, lambdaU, ft, Ets)
                    else:
                        op.uniaxialMaterial('Concrete02', Ele.EleTag * 20 + 2, fpc, epsc0, fcu, epscu, lambdaU, ft, Ets)
                        op.uniaxialMaterial('Concrete02', Ele.EleTag * 20 + 3, fpc, epsc0, fcu, epscu, lambdaU, ft, Ets)
                        epsccu = epscu
                        # print('Concrete02', Ele.EleTag * 20 + 3, fpc, epsc0, fcu, epscu, lambdaU, ft, Ets)
            if self.ui.checkBox_ReguSteel.isChecked() == True:
                pos = 0
                steel_mat_regu()
                pos = 1
                steel_mat_regu()
            if self.ui.checkBox_ReguConcrete.isChecked() == False:
                ft = 0.33 * sqrt(fcC * MPa)
                Ets = ft / 0.002
                # print('Concrete02', Ele.EleTag * 20, -fcC, eo1, -0.2 * fcC, eo20, lambdaU, ft, Ets)
                op.uniaxialMaterial('Concrete02', Ele.EleTag * 20, -fcC, -eo1, -0.2 * fcC, -eo20, lambdaU, ft, Ets)
                op.uniaxialMaterial('Concrete02', Ele.EleTag * 20 + 1, -fcC, -eo1, -0.2 * fcC, -eo20, lambdaU, ft, Ets)
                if CodeDesign == 'ACI 318-19' or CodeDesign == 'CCCSR-84' or CodeDesign == 'Gravity Load Design':
                    if WallDetailing == 'IMF' or WallDetailing == 'SMF':
                        fpcc, epscc0, fccu, epsccu, lambdaC, ft, Ets = con_conf(DW.b, DW.lbe, DW.nsB_BE, DW.nsH_BE,
                                                                                DW.se_BE)
                        if DW.BE == 'Yes':
                            op.uniaxialMaterial('Concrete02', Ele.EleTag * 20 + 2, fpcc, epscc0, fccu, epsccu, lambdaC,
                                                ft, Ets)
                            op.uniaxialMaterial('Concrete02', Ele.EleTag * 20 + 3, fpcc, epscc0, fccu, epsccu, lambdaC,
                                                ft, Ets)
                            # print('Concrete02', Ele.EleTag * 20 + 3, fpcc, epscc0, fccu, epsccu, lambdaC,ft, Ets)

                        else:
                            op.uniaxialMaterial('Concrete02', Ele.EleTag * 20 + 2, -fcC, -eo1, -0.2 * fcC, -eo20,
                                                lambdaU, ft, Ets)
                            op.uniaxialMaterial('Concrete02', Ele.EleTag * 20 + 3, -fcC, -eo1, -0.2 * fcC, -eo20,
                                                lambdaU, ft, Ets)
                            epsccu = eo20
                            # print('Concrete02', Ele.EleTag * 20 + 3, -fcC, -eo1, -0.2 * fcC, -eo20, lambdaU, ft, Ets)
                    else:
                        op.uniaxialMaterial('Concrete02', Ele.EleTag * 20 + 2, -fcC, -eo1, -0.2 * fcC, -eo20, lambdaU,
                                            ft, Ets)
                        op.uniaxialMaterial('Concrete02', Ele.EleTag * 20 + 3, -fcC, -eo1, -0.2 * fcC, -eo20, lambdaU,
                                            ft, Ets)
                        epsccu = eo20
                        # print('Concrete02', Ele.EleTag * 20 + 3, -fcC, -eo1, -0.2 * fcC, -eo20, lambdaU, ft, Ets)
            if self.ui.checkBox_ReguSteel.isChecked() == False:
                pos = 0
                steel_mat()
                pos = 1
                steel_mat()
            dp = DW.dist[0]
            y1 = DW.h / 2.0
            z1 = DW.b / 2.0
            if DW.BE == 'Yes':
                op.section('Fiber', Ele.EleTag * 2)
                op.patch('rect', Ele.EleTag * 20 + 2, NFCols, 1, -y1 + dp, -z1 + dp, -y1 + DW.lbe, z1 - dp)
                op.patch('rect', Ele.EleTag * 20 + 2, NFCols, 1, y1 - DW.lbe, -z1 + dp, y1 - dp, z1 - dp)
                op.patch('rect', Ele.EleTag * 20, NFCols, 1, -y1 + dp, z1 - dp, -y1 + DW.lbe, z1)
                op.patch('rect', Ele.EleTag * 20, NFCols, 1, y1 - DW.lbe, z1 - dp, y1 - dp, z1)
                op.patch('rect', Ele.EleTag * 20, NFCols, 1, -y1 + dp, -z1, -y1 + DW.lbe, -z1 + dp)
                op.patch('rect', Ele.EleTag * 20, NFCols, 1, y1 - DW.lbe, -z1, y1 - dp, -z1 + dp)
                op.patch('rect', Ele.EleTag * 20, 2, 1, -y1, -z1, -y1 + dp, z1)
                op.patch('rect', Ele.EleTag * 20, 2, 1, y1 - dp, -z1, y1, z1)
                op.patch('rect', Ele.EleTag * 20, NFWalls, 1, -y1 + DW.lbe, -z1, y1 - DW.lbe, z1)
                widthsMVLEM = np.hstack([DW.lbe, np.ones(NFMVLEM-2)*((DW.h-2*DW.lbe)/(NFMVLEM-2)), DW.lbe])
                ConcMVLEM = np.hstack([Ele.EleTag * 20 + 2, np.ones(NFMVLEM-2) * Ele.EleTag * 20,
                                          Ele.EleTag * 20 + 2])
            else:
                op.section('Fiber', Ele.EleTag * 2)
                op.patch('rect', Ele.EleTag * 20 + 2, NFWalls, 1, -y1 + dp, -z1 + dp, y1 - dp, z1 - dp)
                op.patch('rect', Ele.EleTag * 20, NFWalls, 1, -y1 + dp, z1 - dp, y1 - dp, z1)
                op.patch('rect', Ele.EleTag * 20, NFWalls, 1, -y1 + dp, -z1, y1 - dp, -z1 + dp)
                op.patch('rect', Ele.EleTag * 20, 2, 1, -y1, -z1, -y1 + dp, z1)
                op.patch('rect', Ele.EleTag * 20, 2, 1, y1 - dp, -z1, y1, z1)
                widthsMVLEM = np.ones(NFMVLEM) * DW.h / NFMVLEM
                ConcMVLEM = np.ones(NFMVLEM) * Ele.EleTag * 20
            for dist, As in zip(DW.dist, DW.As):
                # print('Wall ', Ele.EleTag * 20 + 4, 1, As, -y1 + dist, z1 - dp, -y1 + dist, -z1 + dp)
                op.layer('straight', Ele.EleTag * 20 + 4, 1, As, -y1 + dist, z1 - dp, -y1 + dist, -z1 + dp)
            SteelMVLEM = np.ones(NFMVLEM) * (Ele.EleTag * 20 + 4)
            rhoMVLEM = np.ones(NFMVLEM) * DW.ro_l
            thickMVLEM = np.ones(NFMVLEM) * DW.b
            MassDens = Ele.AEle * GConc / g
            # if ListNodesW[indW, 2] == 0:
            #     op.beamIntegration('Radau', Ele.EleTag, Ele.EleTag * 2, 3)
            # else:
            #     op.beamIntegration('Legendre', Ele.EleTag, Ele.EleTag * 2, 3)
            indW += 1
            Nod_ini = Ele.Nod_ini
            Nod_end = Ele.Nod_end
            eleNodes = [Nod_ini, Nod_end]
            if self.ui.radioButtonFBC.isChecked():
                # op.beamIntegration('HingeRadau', Ele.EleTag, Ele.EleTag*2, phl_hr, Ele.EleTag*2, phl_hr, Ele.EleTag*2)
                if ListNodesW[indW, 2] == 0:
                    op.beamIntegration('Radau', Ele.EleTag, Ele.EleTag * 2, 3)
                else:
                    op.beamIntegration('Legendre', Ele.EleTag, Ele.EleTag * 2, 3)
                op.element('forceBeamColumn', Ele.EleTag * 10, Nod_ini, Nod_end, Ele.ElegTr, Ele.EleTag, '-mass', MassDens)
            MatTagX, MatTagY, MatTagZ = Ele.EleTag * 20 + 6, Ele.EleTag * 20 + 8, Ele.EleTag * 20 + 10
            op.uniaxialMaterial('Elastic', MatTagX, Krigid)
            if self.ui.radioButtonFBC.isChecked():
                if self.ui.radioButtonSMWLinear.isChecked():
                    Ks = 0.4 * Ele.EcEle * DW.b * DW.d / 1.2 / Ele.LEle
                    op.uniaxialMaterial('Elastic', MatTagY, Ks)
                elif self.ui.radioButtonSMWNonlinear.isChecked():
                    Vcr = min(0.5 * DW.Vn, 0.25 * sqrt(fcC * MPa) * DW.b * DW.d)
                    Gamma_cr = Vcr / (0.4 * Ele.EcEle * DW.b * DW.d)
                    f1, f2, f3 = Vcr, 1.5 * DW.Vn, 1.03 * 1.5 * DW.Vn
                    d1, d2, d3 = Gamma_cr * Ele.LEle, 0.0015 * Ele.LEle, 0.015 * Ele.LEle
                    # f1 = 0.1*f1
                    # op.uniaxialMaterial('Hysteretic', MatTagY, f1, d1, 1.1*f1, d2, 1.5*f1, d3, -f1, -d1, -1.1*f1, -d2, -1.5*f1, -d3, 0.8, 0.2, 0.0, 0.02, 0.0)
                    op.uniaxialMaterial('Hysteretic', MatTagY, f1, d1, f2, d2, f3, d3, -f1, -d1, -f2, -d2, -f3, -d3,
                                        0.8, 0.2, 0.0, 0.02, 0.0)
            elif self.ui.radioButtonMVLEM.isChecked():
                MatTagSMVLEM = Ele.EleTag * 20 + 12
                Ks = 0.4 * Ele.EcEle * DW.b * DW.d / 1.2  # / Ele.LEle
                op.uniaxialMaterial('Elastic', MatTagSMVLEM, Ks)
                ShearMVLEM = np.ones(NFMVLEM) * MatTagSMVLEM
                thickMVLEM, widthsMVLEM, rhoMVLEM, ConcMVLEM, SteelMVLEM, ShearMVLEM = np.float64(thickMVLEM),\
                    np.float64(widthsMVLEM), np.float64(rhoMVLEM), np.int64(ConcMVLEM), np.int64(SteelMVLEM),\
                    np.int64(ShearMVLEM)
                thickMVLEM, widthsMVLEM, rhoMVLEM, ConcMVLEM, SteelMVLEM, ShearMVLEM = thickMVLEM.tolist(),\
                    widthsMVLEM.tolist(), rhoMVLEM.tolist(), ConcMVLEM.tolist(), SteelMVLEM.tolist(),\
                    ShearMVLEM.tolist()
                op.element('MVLEM', Ele.EleTag * 10, MassDens, *eleNodes, NFMVLEM, 0.4, '-thick', *thickMVLEM, '-width',
                           *widthsMVLEM, '-rho', *rhoMVLEM, '-matConcrete', *ConcMVLEM, '-matSteel', *SteelMVLEM,
                           '-matShear', *ShearMVLEM)
                # print('MVLEM', Ele.EleTag * 10, MassDens, *eleNodes, NFMVLEM, 0.4, '-thick', *thickMVLEM, '-width',
                #            *widthsMVLEM, '-rho', *rhoMVLEM, '-matConcrete', *ConcMVLEM, '-matSteel', *SteelMVLEM,
                #            '-matShear', *ShearMVLEM)

                op.uniaxialMaterial('Elastic', MatTagY, Krigid)
            op.uniaxialMaterial('Elastic', MatTagZ, Krigid)
            # print('zeroLength', EleTag, Nod_ini, int(Ele.Nod_ini), '-mat', MatTagX, MatTagY, MatTagZ, '-dir', 1, 2, 3)
            EleTag = int(Ele.EleTag * 10 + 1)
            Nod_ini = int(Ele.Nod_ini - nnt)
            Nod_end = int(Ele.Nod_ini)
            op.element('zeroLength', EleTag, Nod_ini, Nod_end, '-mat', MatTagY, MatTagX, MatTagZ, '-dir', 1, 2, 3)
        # print('DataColPhl =', DataColPhl)
        # Creation of non-linear beams
        DataBeamPhl = []
        for (Ele, DB) in zip(EleBeam, DataBeamDesign):
            fc, Ec, nsH = fcB, Ele.EcEle, 2
            if self.ui.radioButton05H.isChecked():
                phl1 = 0.5 * DB.h
                phl2 = 0.5 * DB.h
            elif self.ui.radioButtonPark.isChecked():
                phl1 = 0.08 * Ele.LEle + 0.022 * fy / MPa * DB.db_t1
                phl2 = 0.08 * Ele.LEle + 0.022 * fy / MPa * DB.db_t2
            elif self.ui.radioButtonBerry.isChecked():
                phl1 = 0.05 * Ele.LEle + 0.1 * fy / MPa * DB.db_t1 / sqrt(fc * MPa)
                phl2 = 0.05 * Ele.LEle + 0.1 * fy / MPa * DB.db_t2 / sqrt(fc * MPa)
            DataBeamPhl.append(PlasticHingeLength(phl1, phl2))
            if self.ui.checkBox_ReguConcrete.isChecked() == True and self.ui.checkBox_ReguSteel.isChecked() == True:
                phl = phl1
                fpc, epsc0, fcu, epscu, lambdaU, ft, Ets = con_inconf_regu()
                op.uniaxialMaterial('Concrete02', Ele.EleTag * 20, fpc, epsc0, fcu, epscu, lambdaU, ft, Ets)
                phl = phl2
                fpc, epsc0, fcu, epscu, lambdaU, ft, Ets = con_inconf_regu()
                op.uniaxialMaterial('Concrete02', Ele.EleTag * 20 + 1, fpc, epsc0, fcu, epscu, lambdaU, ft, Ets)
                if CodeDesign == 'ACI 318-19' or CodeDesign == 'CCCSR-84' or CodeDesign == 'Gravity Load Design':
                    if FrameDetailing == 'IMF' or FrameDetailing == 'SMF':
                        phl, pos = phl1, 0
                        fpcc, epscc0, fccu, epsccu, lambdaC, ft, Ets = con_conf_regu(DB.b, DB.h, DB.ns1, nsH, DB.ss1)
                        op.uniaxialMaterial('Concrete02', Ele.EleTag * 20 + 2, fpcc, epscc0, fccu, epsccu, lambdaC, ft,
                                            Ets)
                        steel_mat_regu()
                        phl, pos = phl2, 1
                        fpcc, epscc0, fccu, epsccu, lambdaC, ft, Ets = con_conf_regu(DB.b, DB.h, DB.ns2, nsH, DB.ss2)
                        op.uniaxialMaterial('Concrete02', Ele.EleTag * 20 + 3, fpcc, epscc0, fccu, epsccu, lambdaC, ft,
                                            Ets)
                        steel_mat_regu()
                    else:
                        phl, pos = phl1, 0
                        fpc, epsc0, fcu, epscu, lambdaU, ft, Ets = con_inconf_regu()
                        op.uniaxialMaterial('Concrete02', Ele.EleTag * 20 + 2, fpc, epsc0, fcu, epscu, lambdaU, ft, Ets)
                        steel_mat_regu()
                        phl, pos = phl2, 1
                        fpc, epsc0, fcu, epscu, lambdaU, ft, Ets = con_inconf_regu()
                        op.uniaxialMaterial('Concrete02', Ele.EleTag * 20 + 3, fpc, epsc0, fcu, epscu, lambdaU, ft, Ets)
                        steel_mat_regu()
            if self.ui.checkBox_ReguConcrete.isChecked() == False and self.ui.checkBox_ReguSteel.isChecked() == False:
                ft = 0.33 * sqrt(fcB * MPa)
                Ets = ft / 0.002
                op.uniaxialMaterial('Concrete02', Ele.EleTag * 20, -fcB, -eo1, -0.2 * fcB, -eo20, lambdaU, ft, Ets)
                op.uniaxialMaterial('Concrete02', Ele.EleTag * 20 + 1, -fcB, -eo1, -0.2 * fcB, -eo20, lambdaU, ft, Ets)
                if CodeDesign == 'ACI 318-19' or CodeDesign == 'CCCSR-84' or CodeDesign == 'Gravity Load Design':
                    if FrameDetailing == 'IMF' or FrameDetailing == 'SMF':
                        fpcc, epscc0, fccu, epsccu, lambdaC, ft, Ets = con_conf(DB.b, DB.h, DB.ns1, nsH, DB.ss1)
                        op.uniaxialMaterial('Concrete02', Ele.EleTag * 20 + 2, fpcc, epscc0, fccu, epsccu, lambdaC, ft,
                                            Ets)
                        fpcc, epscc0, fccu, epsccu, lambdaC, ft, Ets = con_conf(DB.b, DB.h, DB.ns2, nsH, DB.ss2)
                        op.uniaxialMaterial('Concrete02', Ele.EleTag * 20 + 3, fpcc, epscc0, fccu, epsccu, lambdaC, ft,
                                            Ets)
                    else:
                        op.uniaxialMaterial('Concrete02', Ele.EleTag * 20 + 2, -fcB, -eo1, -0.2 * fcB, -eo20, lambdaU,
                                            ft, Ets)
                        op.uniaxialMaterial('Concrete02', Ele.EleTag * 20 + 3, -fcB, -eo1, -0.2 * fcB, -eo20, lambdaU,
                                            ft, Ets)
                pos = 0
                steel_mat()
                pos = 1
                steel_mat()
            y2 = DB.h / 2.0
            z2 = DB.b / 2.0 * FN
            dp = DB.h - min(DB.db1, DB.dt1)
            dp1 = dp * FN
            # print('dp1 =', dp)
            pos = 0
            fiber_beam(DB.Ast1 * FN, DB.Asb1 * FN, pos)
            dp = DB.h - min(DB.db2, DB.dt2)
            dp1 = dp * FN
            # print('dp2 =', dp)
            pos = 1
            fiber_beam(DB.Ast2 * FN, DB.Asb2 * FN, pos)
            LColT = Ele.LEle + Ele.RZi + Ele.RZe
            MassDens = ((Ele.AEle * GConc / g) * FN + WDLS / g) #* (LColT / Ele.LEle)
            Nod_ini = Ele.Nod_ini  # + nntot * 4
            Nod_end = Ele.Nod_end  #+ nntot * 5
            op.beamIntegration('HingeRadau', Ele.EleTag, Ele.EleTag * 2, phl1, Ele.EleTag * 2 + 1, phl2, Ele.EleTag * 2)
            op.element('forceBeamColumn', Ele.EleTag * 10, Nod_ini, Nod_end, Ele.ElegTr, Ele.EleTag, '-mass', MassDens)
            # print('DataBeamPhl =', DataBeamPhl)

        for Ele in ZeroLengths:
            # print('Elastic',  Ele.EleTag * 20, 1e15)
            op.uniaxialMaterial('Elastic', Ele.EleTag * 20, Krigid)
            # print('zeroLength', Ele.EleTag, Ele.Nod_ini, Ele.Nod_end, '-mat',  Ele.EleTag * 20, '-dir', 1)
            op.element('zeroLength', Ele.EleTag * 10, Ele.Nod_ini, Ele.Nod_end, '-mat', Ele.EleTag * 20, '-dir', 1)
        list_beams = [Ele.EleTag for Ele in EleBeam]
        list_cols = [Ele.EleTag for Ele in EleCol]
        list_walls = [Ele.EleTag * 10 for Ele in EleWall]
        list_walls = np.int64(list_walls)
        list_walls = list_walls.tolist()
        # print('list_beams =', list_beams)
        # print('list_cols =', list_cols)
        if not os.path.exists("Pushover"):
            os.mkdir("Pushover")
        if self.ui.radioButtonMVLEM.isChecked():
            op.recorder('Element', '-file', 'Pushover/walls_curvature.out', '-time', '-ele', *list_walls, 'Curvature')
            # print('Element', '-file', 'Pushover/walls_curvature.out', '-time', '-ele', *list_walls, 'Curvature')

        # Recording of forces and deformations from nonlinear analysis
        # op.recorder('Element', '-file', 'Pushover/beams_force_1.out',
        #             '-time', '-ele', *list_beams, 'section', 1, 'force')
        # op.recorder('Element', '-file', 'Pushover/beams_def_1.out',
        #             '-time', '-ele', *list_beams, 'section', 1, 'deformation')
        # op.recorder('Element', '-file', 'Pushover/beams_force_6.out',
        #             '-time', '-ele', *list_beams, 'section', 6, 'force')
        # op.recorder('Element', '-file', 'Pushover/beams_def_6.out',
        #             '-time', '-ele', *list_beams, 'section', 6, 'deformation')
        # op.recorder('Element', '-file', 'Pushover/cols_force_1.out',
        #             '-time', '-ele', *list_cols, 'section', 1, 'force')
        # op.recorder('Element', '-file', 'Pushover/cols_def_1.out',
        #             '-time', '-ele', *list_cols, 'section', 1, 'deformation')
        # op.recorder('Element', '-file', 'Pushover/cols_force_6.out',
        #             '-time', '-ele', *list_cols, 'section', 6, 'force')
        # op.recorder('Element', '-file', 'Pushover/cols_def_6.out',
        #             '-time', '-ele', *list_cols, 'section', 6, 'deformation')
        op.recorder('Node', '-file', 'Pushover/HoriNodes.out', '-time', '-node', *ListNodes, '-dof', 1, 'disp')
        op.recorder('Node', '-file', 'Pushover/VertNodes.out', '-time', '-node', *ListNodes, '-dof', 2, 'disp')
        if tw1 != 0:
            Id_NodeW1 = ListNodesW1[:, 0]
            Id_NodeW1 = np.int64(Id_NodeW1)
            Id_NodeW1 = Id_NodeW1.tolist()

        if tw2 != 0:
            Id_NodeW2 = ListNodesW2[:, 0]
            Id_NodeW2 = np.int64(Id_NodeW2)
            Id_NodeW2 = Id_NodeW2.tolist()
        if np.any(ListNodesW1):
            op.recorder('Node', '-file', 'Pushover/HoriNodesW1.out', '-time', '-node', *Id_NodeW1, '-dof', 1, 'disp')
            op.recorder('Node', '-file', 'Pushover/VertNodesW1.out', '-time', '-node', *Id_NodeW1, '-dof', 2, 'disp')
        if np.any(ListNodesW2):
            op.recorder('Node', '-file', 'Pushover/HoriNodesW2.out', '-time', '-node', *Id_NodeW2, '-dof', 1, 'disp')
            op.recorder('Node', '-file', 'Pushover/VertNodesW2.out', '-time', '-node', *Id_NodeW2, '-dof', 2, 'disp')
        # Create a Plain load pattern for gravity loading with a Linear TimeSeries
        op.timeSeries('Linear', 1)
        op.pattern('Plain', 1, 1)
        for Ele in EleCol:
            LCol = Ele.LEle + Ele.RZi + Ele.RZe
            PCol = Ele.AEle * GConc #* LCol / Ele.LEle
            PColSlabD1 = PColSlabD + HCole * WDL
            op.eleLoad('-ele', Ele.EleTag * 10, '-type', '-beamUniform', 0, -PCol * FN)
            op.eleLoad('-ele', Ele.EleTag * 10, '-type', '-beamPoint', 0, 1, -PColSlabD1 - 0.25 * PColSlabL)
        for Ele in EleWall:
            if self.ui.radioButtonFBC.isChecked():
                op.eleLoad('-ele', Ele.EleTag * 10, '-type', '-beamUniform', 0, -Ele.AEle * GConc)
                if np.any(ListNodesW1):
                    op.eleLoad('-ele', Ele.EleTag * 10, '-type', '-beamPoint', 0, 1, -PWall1SlabD - 0.25 * PWall1SlabL)
                if np.any(ListNodesW2):
                    op.eleLoad('-ele', Ele.EleTag * 10, '-type', '-beamPoint', 0, 1, -PWall2SlabD - 0.25 * PWall2SlabL)
            elif self.ui.radioButtonMVLEM.isChecked():
                op.load(Ele.Nod_end, 0., -Ele.AEle * Ele.LEle * GConc, 0)
                if np.any(ListNodesW1):
                    op.load(Ele.Nod_end, 0., -PWall1SlabD - 0.25 * PWall1SlabL, 0)
                if np.any(ListNodesW2):
                    op.load(Ele.Nod_end, 0., -PWall2SlabD - 0.25 * PWall2SlabL, 0)
        for Ele in EleBeam:
            op.eleLoad('-ele', Ele.EleTag * 10, '-type', '-beamUniform', (-Ele.AEle * GConc) * FN - WDL - 0.25 * WLL)
        if self.ui.radioButtonYesLC.isChecked():
            PD_LC = (WDLS - WDL) * Loc_span[-1]
            PL_LC = (WDLS - WDL) / WDL * WLL * Loc_span[-1]
            for Nod_load in ListNodesLC[1:, 0]:
                op.load(int(Nod_load), 0., -PD_LC - 0.25 * PL_LC, 0)
                # print('Nod_load', Nod_load, 0., -PD_LC-0.25*PL_LC, 0)
        # op.printModel()
        Tol = 1.0e-6  # convergence tolerance for test
        op.constraints('Plain')  # how it handles boundary conditions
        op.numberer('Plain')  # renumber dof to minimize band-width (optimization), if you want to
        op.system('BandGeneral')  # how to store and solve the system of equations in the analysis
        op.test('NormDispIncr', Tol, 100)  # determine if convergence has been achieved at the end of an iteration step
        op.algorithm('KrylovNewton')  # use Newton solution algorithm: updates tangent stiffness at every iteration
        NstepGravity = 10  # apply gravity in 10 steps
        DGravity = 1. / NstepGravity  # first load increment;
        op.integrator('LoadControl', DGravity)  # determine the next time step for an analysis
        op.analysis('Static')  # define type of analysis static or transient
        nele = num_elems - 1
        # MG_ElemsForceS1, MG_ElemsDeforS1, MG_ElemsForceS6, MG_ElemsDeforS6 = np.zeros(2 * nele), np.zeros(2 * nele), \
        #     np.zeros(2 * nele), np.zeros(2 * nele)
        step = 1
        loadf = 1.0
        nodeTag = int(ListNodes[-1, 0]) + 1
        while step <= NstepGravity and loadf > 0:
            # ElemsForceS1, ElemsDeforS1, ElemsForceS6, ElemsDeforS6 = [], [], [], []
            op.analyze(1)
            # for Element in Elements:
            #     ForcesS1 = np.array(op.eleResponse(Element.EleTag, 'section', 1, 'force'))
            #     ForcesS6 = np.array(op.eleResponse(Element.EleTag, 'section', 6, 'force'))
            #     DeforsS1 = np.array(op.eleResponse(Element.EleTag, 'section', 1, 'deformation'))
            #     DeforsS6 = np.array(op.eleResponse(Element.EleTag, 'section', 6, 'deformation'))
            #     ElemsForceS1 = np.append(ElemsForceS1, ForcesS1)
            #     ElemsDeforS1 = np.append(ElemsDeforS1, DeforsS1)
            #     ElemsForceS6 = np.append(ElemsForceS6, ForcesS6)
            #     ElemsDeforS6 = np.append(ElemsDeforS6, DeforsS6)
            # MG_ElemsForceS1 = np.vstack((MG_ElemsForceS1, ElemsForceS1))
            # MG_ElemsDeforS1 = np.vstack((MG_ElemsDeforS1, ElemsDeforS1))
            # MG_ElemsForceS6 = np.vstack((MG_ElemsForceS6, ElemsForceS6))
            # MG_ElemsDeforS6 = np.vstack((MG_ElemsDeforS6, ElemsDeforS6))
            if self.ui.radioButtonYesLC.isChecked() == True:
                op.reactions()
                ReaLC = op.nodeReaction(nodeTag, 1)
                # print('nodeTag', nodeTag, 'ReaLC', ReaLC)
            loadf = op.getTime()
            step += 1
        num_beams, num_cols, num_walls = len(EleBeam), len(EleCol), len(EleWall)
        PDG_Beams, PDG_Cols, PDG_Walls = np.zeros([num_beams, 2]), np.zeros([num_cols, 2]), np.zeros([num_walls, 2])
        for (Ele, ind) in zip(EleBeam, range(num_beams)):
            PDG = op.eleResponse(Ele.EleTag * 10, 'plasticDeformation')
            PDG_Beams[ind, :] = [PDG[1], PDG[2]]
        for (Ele, ind) in zip(EleCol, range(num_cols)):
            PDG = op.eleResponse(Ele.EleTag * 10, 'plasticDeformation')
            PDG_Cols[ind, :] = [PDG[1], PDG[2]]
        if self.ui.radioButtonFBC.isChecked():
            for (Ele, ind) in zip(EleWall, range(num_walls)):
                PDG = op.eleResponse(Ele.EleTag * 10, 'plasticDeformation')
                # print('PDG', PDG)
                PDG_Walls[ind, :] = [PDG[1], PDG[2]]

        # print('PDG_Cols', PDG_Cols)
        # print('PDG_Beams', PDG_Beams)
        op.loadConst('-time', 0.0)
        print("Model Nonlinear Built")

        # xi = 0.025  # damping ratio
        xi = float(self.ui.Damping.text()) / 100
        floors_num = len(Loc_heigth) - 1
        # print('floors_num =', floors_num)
        if floors_num > 2 and self.ui.radioButtonFBC.isChecked():
            from ModalAnalysis import ModalAnalysis
            T, Mratios, Mfactors, Mtots = ModalAnalysis(floors_num, pflag=0, outname='Prueba2')
            MPart = Mratios[1]
            # print('MPart', MPart)
            CumMPart = np.cumsum(MPart)
            # print('CumMPart', CumMPart)
            indN90 = np.where(CumMPart >= 90)[0][0]
            # print('indN90', indN90)
            # print('Mratios', Mratios[1])
            self.ui.T1.setText(str(round(T[0], 2)) + ' sec')
            self.ui.T2.setText(str(round(T[1], 2)) + ' sec')
            self.ui.T3.setText(str(round(T[2], 2)) + ' sec')

            wn = 2 * pi / T  # eigenvalue analysis for nEigenJ modes
            if self.ui.radioButtonT1_T3.isChecked():
                title = r'Damping Model:$\omega_1, \omega_3$'
                wi, wj = wn[0], wn[2]
                A = np.array([[1 / wi, wi], [1 / wj, wj]])
                b = np.array([xi, xi])
                x = np.linalg.solve(A, 2 * b)
                xiv = x[0] / (2 * wn) + x[1] * wn / 2

            elif self.ui.radioButtonASCE_41_17.isChecked():
                # ASCE Damping for NLD
                title = 'Damping Model: ASCE 41-17'
                wi = wn[0] / 2
                for ind in range(1, floors_num):
                    wj = wn[ind]
                    A = np.array([[1 / wi, wi], [1 / wj, wj]])
                    b = np.array([xi, xi])
                    x = np.linalg.solve(A, 2 * b)
                    xiv = x[0] / (2 * wn) + x[1] * wn / 2
                    # print('xiv', xiv)
                    if 8 * xiv[0] > xiv[indN90] and xi >= np.sum(MPart[:indN90] * xiv[:indN90] / 100):
                        break

            if self.ui.radioButtonCommitted.isChecked():
                op.rayleigh(x[0], 0, 0, x[1])  # RAYLEIGH damping
            elif self.ui.radioButtonCurrent.isChecked():
                op.rayleigh(x[0], x[1], 0, 0)  # RAYLEIGH damping
            elif self.ui.radioButtonInitial.isChecked():
                op.rayleigh(x[0], 0, x[1], 0)  # RAYLEIGH damping
            T1m = T[0]
            cmap = plt.cm.get_cmap('jet')
            fig = self.ui.DampingCurve.canvas.axes
            fig.clear()
            ax = fig.add_axes([0.15, 0.2, 0.70, 0.70])
            # ax = fig.add_axes([0, 0.1, 1, 0.9])

            ax.scatter(wn, xiv, s=5 * MPart, c=MPart, cmap=cmap, alpha=1)

            wn_v = np.linspace(wn[0] / 2, wn[-1], num=200)
            xi_v = x[0] / (2 * wn_v) + x[1] * wn_v / 2
            ax.plot(wn_v, xi_v, '-')
            ax.axhline(y=xi, c='red', linestyle='dashed')

            ax.set_ylabel(r'$\zeta_n$')
            # r'$j = %1.1f \/\/\theta = %1.2f \/\/\beta = %1.2f$'
            ax.set_xlabel(r'$\omega_n(rad/sec)$')
            ax.set_title(title)
            ax.grid(True)
            ax2 = fig.add_axes([0.9, 0.2, 0.02, 0.7])
            # ax2 = fig.add_axes([0.05, 0.05, 0.7, 0.03])
            norm = plt.Normalize(vmin=np.min(MPart), vmax=np.max(MPart))
            cbar = colorbar.ColorbarBase(ax2, cmap=cmap, norm=norm, orientation='vertical')  # ,
            # label='Mass(%)', labelpad=1)
            cbar.set_label('M*(%)', labelpad=-20, y=1.10, rotation=0, fontstyle='italic')
            fig.set_tight_layout(False)
            self.ui.DampingCurve.canvas.draw()
            self.ui.DampingCurve.canvas.show()

            if not os.path.exists("NonLinearModel"):
                os.mkdir("NonLinearModel")
            DampingModel = np.vstack((wn, xiv, MPart))
            np.savetxt('NonLinearModel/' + 'DampingModel.txt', DampingModel.T)

        elif floors_num == 1:
            nEigenI = 1  # mode 1
            # nEigenI2 = 2  # mode 2
            lambdaN = op.eigen(nEigenI)  # eigenvalue analysis for nEigenJ modes
            lambdaI = lambdaN[nEigenI - 1]  # eigenvalue mode i
            # lambdaI2 = lambdaN[nEigenI2 - 1]  # eigenvalue mode i2
            # print('lambdaN ', lambdaN)
            omegaI = pow(lambdaI, 0.5)
            # omegaI2 = pow(lambdaI2, 0.5)
            T1m = 2. * pi / omegaI
            self.ui.T1.setText(str(round(T1m, 2)) + ' sec')
            self.ui.T2.setText(str(' '))
            self.ui.T3.setText(str(' '))

            # T2m = 2. * pi / omegaI2
            # print('Ta1=', T1m, 'seg')
        elif self.ui.radioButtonMVLEM.isChecked():
            Lambda = np.asarray(op.eigen(3))
            Omega = Lambda ** 0.5
            T = 2 * np.pi / Omega
            self.ui.T1.setText(str(round(T[0], 2)) + ' sec')
            self.ui.T2.setText(str(round(T[1], 2)) + ' sec')
            self.ui.T3.setText(str(round(T[2], 2)) + ' sec')
            wn = 2 * pi / T  # eigenvalue analysis for nEigenJ modes
            wi, wj = wn[0], wn[2]
            A = np.array([[1 / wi, wi], [1 / wj, wj]])
            b = np.array([xi, xi])
            x = np.linalg.solve(A, 2 * b)
            if self.ui.radioButtonCommitted.isChecked():
                op.rayleigh(x[0], 0, 0, x[1])  # RAYLEIGH damping
            elif self.ui.radioButtonCurrent.isChecked():
                op.rayleigh(x[0], x[1], 0, 0)  # RAYLEIGH damping
            elif self.ui.radioButtonInitial.isChecked():
                op.rayleigh(x[0], 0, x[1], 0)  # RAYLEIGH damping
            T1m = T[0]


        elif floors_num == 2:
            nEigenI = 1  # mode 1
            nEigenI2 = 2  # mode 2
            lambdaN = op.eigen(nEigenI2)  # eigenvalue analysis for nEigenJ modes
            lambdaI = lambdaN[nEigenI - 1]  # eigenvalue mode i
            lambdaI2 = lambdaN[nEigenI2 - 1]  # eigenvalue mode i2
            # print('lambdaN ', lambdaN)
            omegaI = pow(lambdaI, 0.5)
            omegaI2 = pow(lambdaI2, 0.5)
            T1m = 2. * pi / omegaI
            T2m = 2. * pi / omegaI2
            self.ui.T1.setText(str(round(T1m, 2)) + ' sec')
            self.ui.T2.setText(str(round(T2m, 2)) + ' sec')
            self.ui.T3.setText(str(' '))

            # print('Ta1=', T1m, 'seg')

        self.ui.tabWidget.setCurrentIndex(3)

    # Pushover function
    def Pushover(self):
        global cbar, DataColPhl, num_nodes

        def singlePush1(dref, mu, ctrlNode, dispDir, nSteps):
            IOflag = 2
            testType = 'NormDispIncr'
            # set testType	EnergyIncr;					# Dont use with Penalty constraints
            # set testType	RelativeNormUnbalance;		# Dont use with Penalty constraints
            # set testType	RelativeNormDispIncr;		# Dont use with Lagrange constraints
            # set testType	RelativeTotalNormDispIncr;	# Dont use with Lagrange constraints
            # set testType	RelativeEnergyIncr;			# Dont use with Penalty constraints
            tolInit = 1.0e-7  # the initial Tolerance, so it can be referred back to
            iterInit = 50  # the initial Max Number of Iterations
            algorithmType = 'KrylovNewton'  # the algorithm type

            op.test(testType, tolInit
                    , iterInit)  # determine if convergence has been achieved at the end of an iteration step
            op.algorithm(algorithmType)  # use Newton solution algorithm: updates tangent stiffness at every iteration
            disp = dref * mu
            dU = disp / (1.0 * nSteps)
            # print('dref ', dref, 'mu ', mu, 'dU ', dU, 'disp ', disp)
            op.integrator('DisplacementControl', ctrlNode, dispDir, dU)  # determine the next time step for an analysis
            op.analysis('Static')  # define type of analysis static or transient

            # Print values
            if IOflag >= 1:
                print('singlePush: Push ', ctrlNode, ' to ', mu)

            #      the initial values to start the while loop
            nele = num_elems - 1
            # MP_ElemsForceS1, MP_ElemsDeforS1, MP_ElemsForceS6, MP_ElemsDeforS6 = np.zeros(2 * nele), \
            #     np.zeros(2 * nele), \
            #     np.zeros(2 * nele), np.zeros(2 * nele)
            ok = 0
            step = 1
            loadf = 1.0
            # This feature of disabling the possibility of having a negative loading has been included.
            # This has been adapted from a similar script by Prof. Garbaggio
            htot = op.nodeCoord(ctrlNode, 2)
            Der_obj = dref * nSteps / htot
            maxDriftPiso = 0.0
            VBasal_v = [0]
            DriftTecho_v = [0]
            MaxDriftStory_v = [0]
            DriftStory_m = np.zeros(floors_num)
            num_beams, num_cols, num_walls = len(EleBeam), len(EleCol), len(EleWall)
            PD_Beams, PD_Cols, PD_Walls = np.zeros([num_beams, 2]), np.zeros([num_cols, 2]), np.zeros([num_walls, 2])
            Vu_Vn_Cols = np.zeros(num_cols)
            ratio_infill = np.zeros(num_cols)
            T1_v, T2_v = [], []
            while step <= nSteps and ok == 0 and loadf > 0:
                self.ui.progressBarPushover.setValue(100 * step / nSteps)
                ElemsForceS1, ElemsDeforS1, ElemsForceS6, ElemsDeforS6 = [], [], [], []
                ok = op.analyze(1)
                loadf = op.getTime()
                temp = op.nodeDisp(ctrlNode, dispDir)
                # Print the current displacement
                if IOflag >= 2:
                    print('Pushed ', ctrlNode, ' in ', dispDir, ' to ', temp, ' with ', loadf, 'step', step)

                # If the analysis fails, try the following changes to achieve convergence
                # Analysis will be slower in here though...
                if ok != 0:
                    print('Trying relaxed convergence..')
                    op.test(testType, tolInit * 0.01,
                            iterInit * 50)  # determine if convergence has been achieved at the end of an iteration step
                    ok = op.analyze(1)
                    op.test(testType, tolInit,
                            iterInit)  # determine if convergence has been achieved at the end of an iteration step
                if ok != 0:
                    print('Trying Newton with initial then current .')
                    op.test(testType, tolInit * 0.01,
                            iterInit * 50)  # determine if convergence has been achieved at the end of an iteration step
                    op.algorithm('Newton', '-initialThenCurrent')
                    ok = op.analyze(1)
                    op.algorithm(algorithmType)
                    op.test(testType, tolInit,
                            iterInit)  # determine if convergence has been achieved at the end of an iteration step
                if ok != 0:
                    print('Trying ModifiedNewton with initial ..')
                    op.test(testType, tolInit * 0.01,
                            iterInit * 50)  # determine if convergence has been achieved at the end of an iteration step
                    op.algorithm('ModifiedNewton', '-initial')
                    ok = op.analyze(1)
                    op.algorithm(algorithmType)
                    op.test(testType, tolInit,
                            iterInit)  # determine if convergence has been achieved at the end of an iteration step
                if ok != 0:
                    print('Trying KrylovNewton ..')
                    op.test(testType, tolInit * 0.01,
                            iterInit * 50)  # determine if convergence has been achieved at the end of an iteration step
                    op.algorithm('KrylovNewton')
                    ok = op.analyze(1)
                    op.algorithm(algorithmType)
                    op.test(testType, tolInit,
                            iterInit)  # determine if convergence has been achieved at the end of an iteration step
                if ok != 0:
                    print('Trying FixedNumIter .. ....')
                    op.test('FixedNumIter',
                            iterInit)  # determine if convergence has been achieved at the end of an iteration step
                    ok = op.analyze(1)
                drift_piso_v = np.zeros(floors_num)
                for (ii, nod_ini, nod_end) in zip(range(floors_num), ListNodesDrift[:-1, 0], ListNodesDrift[1:, 0]):
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
                    drift_piso_v[ii] = drift_piso
                    if drift_piso >= maxDriftPiso:
                        maxDriftPiso = drift_piso

                VBasal = 0.
                op.reactions()
                for node in ListNodesBasal:
                    # print('ind Basal ', node[0], 'VBasal', VBasal)
                    VBasal = VBasal + op.nodeReaction(int(node[0]), 1)
                if self.ui.radioButtonYesLC.isChecked():
                    VBasal = VBasal + op.nodeReaction(int(ListNodesLC[0, 0]), 1)
                # w2 = op.eigen('-fullGenLapack', 3)
                # T1 = 2. * pi / sqrt(abs(w2[0]))
                # T2 = 2. * pi / sqrt(abs(w2[1]))
                # T1, T2 = w2[0], w2[1]
                # T1_v = np.append(T1_v, T1)
                # T2_v = np.append(T2_v, T2)
                VBasal_v = np.append(VBasal_v, VBasal)
                DriftTecho = op.nodeDisp(ctrlNode, dispDir) / htot
                DriftTecho_v = np.append(DriftTecho_v, DriftTecho)
                MaxDriftStory_v = np.append(MaxDriftStory_v, maxDriftPiso)
                DriftStory_m = np.vstack((DriftStory_m, drift_piso_v))

                for (Ele, DC, ind, IFE) in zip(EleCol, DataColDesign, range(num_cols), infill_m):
                    ForcesS1 = np.array(op.eleResponse(Ele.EleTag * 10, 'force'))
                    Vu = abs(ForcesS1[0]) / FN
                    Vu_Vn = Vu / DC.VnCol
                    Vu_Vn_Cols[ind] = max(Vu_Vn_Cols[ind], Vu_Vn)
                    if IFE != 0:
                        EleTag = int(Ele.EleTag * 10 + 4)
                        ForcesS1 = np.array(op.eleResponse(EleTag, 'force'))
                        PInfillC = sqrt(ForcesS1[0]**2 + ForcesS1[1]**2)
                        EleTag = int(Ele.EleTag * 10 + 5)
                        ForcesS1 = np.array(op.eleResponse(EleTag, 'force'))
                        PInfillS = sqrt(ForcesS1[0]**2 + ForcesS1[1]**2)
                        PInfill = PInfillC + PInfillS
                        sigma_infill = PInfill / DC.wIF / twIF
                        ratio_sigma = sigma_infill / DC.sigma_IF
                        ratio_infill[ind] = max(ratio_infill[ind], ratio_sigma)

                # for Element in Elements:
                #     ForcesS1 = np.array(op.eleResponse(Element.EleTag, 'section', 1, 'force'))
                #     ForcesS6 = np.array(op.eleResponse(Element.EleTag, 'section', 6, 'force'))
                #     DeforsS1 = np.array(op.eleResponse(Element.EleTag, 'section', 1, 'deformation'))
                #     DeforsS6 = np.array(op.eleResponse(Element.EleTag, 'section', 6, 'deformation'))
                #     ElemsForceS1 = np.append(ElemsForceS1, ForcesS1)
                #     ElemsDeforS1 = np.append(ElemsDeforS1, DeforsS1)
                #     ElemsForceS6 = np.append(ElemsForceS6, ForcesS6)
                #     ElemsDeforS6 = np.append(ElemsDeforS6, DeforsS6)
                # MP_ElemsForceS1 = np.vstack((MP_ElemsForceS1, ElemsForceS1))
                # MP_ElemsDeforS1 = np.vstack((MP_ElemsDeforS1, ElemsDeforS1))
                # MP_ElemsForceS6 = np.vstack((MP_ElemsForceS6, ElemsForceS6))
                # MP_ElemsDeforS6 = np.vstack((MP_ElemsDeforS6, ElemsDeforS6))

                if maxDriftPiso >= Der_obj:
                    break

                loadf = op.getTime()
                step += 1
            # print('MP_ElemsForceS1 =', MP_ElemsForceS1)
            # print('Tamñano de MP_ElemsForceS1 =', np.shape(MP_ElemsForceS1))
            maxDriftTecho = dU * step / htot
            maxDriftTecho2 = op.nodeDisp(ctrlNode, dispDir) / htot
            for (Ele, ind) in zip(EleBeam, range(num_beams)):
                PD = op.eleResponse(Ele.EleTag * 10, 'plasticDeformation')
                PD_Beams[ind, :] = [PD[1], PD[2]]
            for (Ele, ind) in zip(EleCol, range(num_cols)):
                PD = op.eleResponse(Ele.EleTag * 10, 'plasticDeformation')
                PD_Cols[ind, :] = [PD[1], PD[2]]
                # ForcesS1 = np.array(op.eleResponse(Ele.EleTag * 10, 'force'))
            if self.ui.radioButtonFBC.isChecked():
                for (Ele, ind) in zip(EleWall, range(num_walls)):
                    PD = op.eleResponse(Ele.EleTag * 10, 'plasticDeformation')
                    # PD_Walls[ind, :] = [PD[1], PD[2]]
                    PD_Walls[ind, :] = [PD[1], 0]


            # print('PD_Cols', PD_Cols)
            # print('PD_Beams', PD_Beams)

            PD_Beams = abs(np.array(PDG_Beams + PD_Beams))
            PD_Cols = abs(np.array(PDG_Cols + PD_Cols))
            if self.ui.radioButtonFBC.isChecked():
                PD_Walls = abs(np.array(PDG_Walls + PD_Walls))

            # T_v = np.vstack((T1_v, T2_v))
            # print('PD_Cols', PD_Cols)
            # print('PD_Beams', PD_Beams)

            if ok != 0:
                print('DispControl Analysis FAILED')
            else:
                print('DispControl Analysis SUCCESSFUL')
            if loadf <= 0:
                print('Stopped because of Load factor below zero: ', loadf)
            #    if PrintFlag == 0:
            #        os.remove("singlePush.txt")
            #        print singlePush.txt
            return maxDriftPiso, maxDriftTecho, maxDriftTecho2, VBasal_v, MaxDriftStory_v, DriftTecho_v, \
                DriftStory_m, PD_Beams, PD_Cols, PD_Walls, Vu_Vn_Cols, ratio_infill

            # return maxDriftPiso, maxDriftTecho, maxDriftTecho2, VBasal_v, MaxDriftStory_v, DriftTecho_v, \
            #     MP_ElemsForceS1, MP_ElemsDeforS1, MP_ElemsForceS6, MP_ElemsDeforS6, DriftStory_m, PD_Beams, PD_Cols

        # Pushover function varying tests and algorithms
        def singlePush(dref, mu, ctrlNode, dispDir, nSteps):
            # --------------------------------------------------
            # Description of Parameters
            # --------------------------------------------------
            # dref:			Reference displacement to which cycles are run. Corresponds to yield or equivalent other, such as 1mm
            # mu:			Multiple of dref to which the push is run. So pushover can be run to a specifived ductility or displacement
            # ctrlNode:		Node to control with the displacement integrator.
            # dispDir:		DOF the loading is applied.
            # nSteps:		Number of steps.
            # IOflag:		Option to print details on screen. 2 for print of each step, 1 for basic info (default), 0 for off
            # ---------------------------------------------------
            test = {1: 'NormDispIncr', 2: 'RelativeEnergyIncr', 3: 'EnergyIncr',
                    4: 'RelativeNormUnbalance', 5: 'RelativeNormDispIncr',
                    6: 'NormUnbalance', 7: 'FixedNumIter'}
            alg = {1: 'KrylovNewton', 2: 'SecantNewton', 3: 'ModifiedNewton',
                   4: 'RaphsonNewton', 5: 'PeriodicNewton', 6: 'BFGS',
                   7: 'Broyden', 8: 'NewtonLineSearch'}

            # test = {1:'NormDispIncr', 2: 'RelativeEnergyIncr', 3:'EnergyIncr'}
            # alg = {1:'KrylovNewton', 2:'ModifiedNewton'}

            IOflag = 2
            PrintFlag = 0
            testType = 'RelativeNormDispIncr'  # Dont use with Penalty constraints

            tolInit = 1.0e-7  # the initial Tolerance, so it can be referred back to
            iterInit = 50  # the initial Max Number of Iterations
            algorithmType = 'KrylovNewton'  # the algorithm type
            #      	algorithmType Newton;		#      the algorithm type
            #      	algorithmType Newton;		#      the algorithm type

            # op.constraints('Transformation') # how it handles boundary conditions
            # op.numberer('RCM')    # renumber dof to minimize band-width (optimization), if you want to
            # op.system('BandGeneral') # how to store and solve the system of equations in the analysis

            op.test(testType, tolInit,
                    iterInit)  # determine if convergence has been achieved at the end of an iteration step
            op.algorithm(algorithmType)  # use Newton solution algorithm: updates tangent stiffness at every iteration
            disp = dref * mu
            dU = disp / (1.0 * nSteps)
            # print('dref ', dref, 'mu ', mu, 'dU ', dU, 'disp ', disp, 'nSteps ', nSteps)
            op.integrator('DisplacementControl', ctrlNode, dispDir, dU)  # determine the next time step for an analysis
            op.analysis('Static')  # defivne type of analysis static or transient

            # Print values
            if IOflag >= 1:
                print('singlePush: Push ', ctrlNode, ' to ', mu)
            nele = num_elems - 1
            # MP_ElemsForceS1, MP_ElemsDeforS1, MP_ElemsForceS6, MP_ElemsDeforS6 = np.zeros(2 * nele), \
            #     np.zeros(2 * nele), \
            #     np.zeros(2 * nele), np.zeros(2 * nele)
            #      the initial values to start the while loop
            ok = 0
            step = 1
            loadf = 1.0
            # This feature of disabling the possibility of having a negative loading has been included.
            # This has been adapted from a similar script by Prof. Garbaggio
            maxDriftPiso = 0.0
            htot = op.nodeCoord(ctrlNode, 2)
            maxDriftPiso = 0.0
            VBasal_v = [0]
            DriftTecho_v = [0]
            MaxDriftStory_v = [0]
            DriftStory_m = np.zeros(floors_num)
            PD_Beams, PD_Cols, PD_Walls = np.zeros([num_beams, 2]), np.zeros([num_cols, 2]), np.zeros([num_walls, 2])
            Vu_Vn_Cols = np.zeros(num_cols)
            ratio_infill = np.zeros(num_cols)
            T1_v, T2_v = [], []
            # factor_v = np.array([1,0.75,0.5,0.25,0.1,2,3,5,10])
            # fact_v = np.array([50,100,500])
            # factor = 100
            # fact = 1.
            while step <= nSteps and ok == 0 and loadf > 0:
                self.ui.progressBarPushover.setValue(100 * step / nSteps)
                # ElemsForceS1, ElemsDeforS1, ElemsForceS6, ElemsDeforS6 = [], [], [], []
                ok = op.analyze(1)
                loadf = op.getTime()
                temp = op.nodeDisp(ctrlNode, dispDir)
                if IOflag >= 2:
                    print('Pushed ', ctrlNode, ' in ', dispDir, ' to ', temp, ' with ', loadf, 'step ', step)
                # for factor in factor_v:
                # op.integrator('DisplacementControl',ctrlNode,dispDir,factor*dU)  # determine the next time step for an analysis
                # for fact in fact_v:
                for j in alg:
                    for i in test:
                        for fact in [1, 20, 50]:
                            if ok != 0 and j >= 4 and i != 7:
                                # print('Trying ',str(alg[j]))
                                op.test(test[i], tolInit * .01, iterInit * fact)
                                op.algorithm(alg[j])
                                ok = op.analyze(1)
                                op.algorithm(algorithmType)
                                op.test(testType, tolInit, iterInit)
                            elif ok != 0 and j < 4 and i != 7:
                                # print('Trying ',str(alg[j]))
                                op.test(test[i], tolInit, iterInit * fact)
                                op.algorithm(alg[j], '-initial')
                                ok = op.analyze(1)
                                op.algorithm(algorithmType)
                                op.test(testType, tolInit, iterInit)
                            if ok == 0:
                                break
                        if ok != 0 and i == 7:
                            op.test(test[i], iterInit)
                            op.algorithm(alg[j])
                            ok = op.analyze(1)
                        if ok == 0:
                            break
                    if ok == 0:
                        break
                    # if ok == 0:
                    #     break
                    # if ok == 0:
                    #     break
                # op.integrator('DisplacementControl',ctrlNode,dispDir,dU)  # determine the next time step for an analysis
                # Calculation of maximum Drift between floors
                drift_piso_v = np.zeros(floors_num)
                for (ii, nod_ini, nod_end) in zip(range(floors_num), ListNodesDrift[:-1, 0], ListNodesDrift[1:, 0]):
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
                    drift_piso_v[ii] = drift_piso
                    if drift_piso >= maxDriftPiso:
                        maxDriftPiso = drift_piso

                VBasal = 0.
                op.reactions()
                for node in ListNodesBasal:
                    # print('ind Basal ', node[0], 'VBasal', VBasal)
                    VBasal = VBasal + op.nodeReaction(int(node[0]), 1)
                if self.ui.radioButtonYesLC.isChecked():
                    VBasal = VBasal + op.nodeReaction(int(ListNodesLC[0, 0]), 1)
                # w2 = op.eigen('-fullGenLapack', 3)
                # T1 = 2. * pi / sqrt(abs(w2[0]))
                # T2 = 2. * pi / sqrt(abs(w2[1]))
                # T1, T2 = w2[0], w2[1]
                # T1_v = np.append(T1_v, T1)
                # T2_v = np.append(T2_v, T2)
                VBasal_v = np.append(VBasal_v, VBasal)
                DriftTecho = op.nodeDisp(ctrlNode, dispDir) / htot
                DriftTecho_v = np.append(DriftTecho_v, DriftTecho)
                MaxDriftStory_v = np.append(MaxDriftStory_v, maxDriftPiso)
                # print('MaxDriftStory_v', MaxDriftStory_v)
                DriftStory_m = np.vstack((DriftStory_m, drift_piso_v))
                for (Ele, DC, ind, IFE) in zip(EleCol, DataColDesign, range(num_cols), infill_m):

                    ForcesS1 = np.array(op.eleResponse(Ele.EleTag * 10, 'force'))
                    Vu = abs(ForcesS1[0]) / FN
                    Vu_Vn = Vu / DC.VnCol
                    Vu_Vn_Cols[ind] = max(Vu_Vn_Cols[ind], Vu_Vn)
                    if IFE != 0:
                        EleTag = int(Ele.EleTag * 10 + 4)
                        ForcesS1 = np.array(op.eleResponse(EleTag, 'force'))
                        PInfillC = sqrt(ForcesS1[0]**2 + ForcesS1[1]**2)
                        EleTag = int(Ele.EleTag * 10 + 6)
                        ForcesS1 = np.array(op.eleResponse(EleTag, 'force'))
                        PInfillS = sqrt(ForcesS1[0]**2 + ForcesS1[1]**2)
                        PInfill = PInfillC + PInfillS
                        sigma_infill = PInfill / DC.wIF / twIF
                        ratio_sigma = sigma_infill / DC.sigma_IF
                        ratio_infill[ind] = max(ratio_infill[ind], ratio_sigma)

                # for Element in Elements:
                #     ForcesS1 = np.array(op.eleResponse(Element.EleTag, 'section', 1, 'force'))
                #     ForcesS6 = np.array(op.eleResponse(Element.EleTag, 'section', 6, 'force'))
                #     DeforsS1 = np.array(op.eleResponse(Element.EleTag, 'section', 1, 'deformation'))
                #     DeforsS6 = np.array(op.eleResponse(Element.EleTag, 'section', 6, 'deformation'))
                #     ElemsForceS1 = np.append(ElemsForceS1, ForcesS1)
                #     ElemsDeforS1 = np.append(ElemsDeforS1, DeforsS1)
                #     ElemsForceS6 = np.append(ElemsForceS6, ForcesS6)
                #     ElemsDeforS6 = np.append(ElemsDeforS6, DeforsS6)
                # MP_ElemsForceS1 = np.vstack((MP_ElemsForceS1, ElemsForceS1))
                # MP_ElemsDeforS1 = np.vstack((MP_ElemsDeforS1, ElemsDeforS1))
                # MP_ElemsForceS6 = np.vstack((MP_ElemsForceS6, ElemsForceS6))
                # MP_ElemsDeforS6 = np.vstack((MP_ElemsDeforS6, ElemsDeforS6))
                if maxDriftPiso >= Der_obj:
                    break
                loadf = op.getTime()
                step += 1
            maxDriftTecho = dU * step / htot
            maxDriftTecho2 = op.nodeDisp(ctrlNode, dispDir) / htot
            for (Ele, ind) in zip(EleBeam, range(num_beams)):
                PD = op.eleResponse(Ele.EleTag * 10, 'plasticDeformation')
                PD_Beams[ind, :] = [PD[1], PD[2]]
            for (Ele, ind) in zip(EleCol, range(num_cols)):
                PD = op.eleResponse(Ele.EleTag * 10, 'plasticDeformation')
                PD_Cols[ind, :] = [PD[1], PD[2]]
            if self.ui.radioButtonFBC.isChecked():
                for (Ele, ind) in zip(EleWall, range(num_walls)):
                    PD = op.eleResponse(Ele.EleTag * 10, 'plasticDeformation')
                    # PD_Walls[ind, :] = [PD[1], PD[2]]
                    PD_Walls[ind, :] = [PD[1], 0]
            PD_Beams = abs(np.array(PDG_Beams + PD_Beams))
            PD_Cols = abs(np.array(PDG_Cols + PD_Cols))
            if self.ui.radioButtonFBC.isChecked():
                PD_Walls = abs(np.array(PDG_Walls + PD_Walls))

            if ok != 0:
                print('DispControl Analysis FAILED')
            else:
                print('DispControl Analysis SUCCESSFUL')
            if loadf <= 0:
                print('Stopped because of Load factor below zero: ', loadf)
            #    if PrintFlag == 0:
            #        os.remove("singlePush.txt")
            #        print singlePush.txt
            return maxDriftPiso, maxDriftTecho, maxDriftTecho2, VBasal_v, MaxDriftStory_v, DriftTecho_v, \
                DriftStory_m, PD_Beams, PD_Cols, PD_Walls, Vu_Vn_Cols, ratio_infill

            # return maxDriftPiso, maxDriftTecho, maxDriftTecho2, VBasal_v, MaxDriftStory_v, DriftTecho_v, \
            #     MP_ElemsForceS1, MP_ElemsDeforS1, MP_ElemsForceS6, MP_ElemsDeforS6, DriftStory_m, PD_Beams, PD_Cols,\
            #     PD_Walls

        ListNodesDrift = ListNodes[np.where(ListNodes[:, 1] == 0.)]
        # ListNodesBasal = ListNodes[np.where(ListNodes[:, 2] == 0.)]
        if T1m <= 0.5:
            k = 1.
        elif T1m <= 2.5:
            k = 0.75 + 0.5 * T1m
        else:
            k = 2.

        sumH = np.sum(np.power(Loc_heigth, k))
        floors_num = len(Loc_heigth)
        num_beams, num_cols, num_walls = len(EleBeam), len(EleCol), len(EleWall)
        # Defining the pushover lateral distribution type
        if self.ui.radioButtonTriangular.isChecked():
            Fp = np.power(Loc_heigth, k) / sumH
        if self.ui.radioButtonUniform.isChecked():
            Fp = 1. / floors_num * np.ones(floors_num + 1)
        print('Fp =', Fp)
        op.loadConst('-time', 0.0)
        op.timeSeries('Linear', 2)
        op.pattern('Plain', 2, 1)
        # for (node, fp, ind) in zip(ListNodesDrift, Fp, range(floors_num)):
        #     op.load(int(node[0]), fp, 0.0, 0.0)
        num_nodes = len(Loc_span) * len(Loc_heigth)
        for (h_floor, fp) in zip(Loc_heigth[1:], Fp):
            for ind_node in range(num_nodes):
                if ListNodes[ind_node, 2] == h_floor:
                    op.load(ind_node, fp, 0.0, 0.0)
        Htotal = Loc_heigth[-1]
        Der_obj = float(self.ui.Der_obj.text())
        twIF = float(self.ui.twIF.text()) * mm
        Des_obj = Der_obj * Htotal  # Desplazamiento objetivo
        nSteps = int(self.ui.nSteps.text())
        dref = Des_obj / nSteps
        mu = nSteps
        IDctrlNode = int(ListNodesDrift[-1, 0])  # Node where displacement is read
        # print('IDctrlNode =', IDctrlNode)
        IDctrlDOF = 1  # DOF x=1, y=2
        Tol = 1.0e-4  # Tolerance
        self.ui.progressBarPushover.show()
        if self.ui.radioButtonFast.isChecked():
            # maxDriftPiso, maxDriftTecho, maxDriftTecho2, VBasal_v, MaxDriftStory_v, DriftTecho_v, MP_ElemsForceS1, \
            #     MP_ElemsDeforS1, MP_ElemsForceS6, MP_ElemsDeforS6, DriftStory_m, PD_Beams, PD_Cols, PD_Walls\
            #     = singlePush(dref, mu, IDctrlNode, IDctrlDOF, nSteps)
            maxDriftPiso, maxDriftTecho, maxDriftTecho2, VBasal_v, MaxDriftStory_v, DriftTecho_v, DriftStory_m, \
                PD_Beams, PD_Cols, PD_Walls, Vu_Vn_Cols, ratio_infill = singlePush1(dref, mu, IDctrlNode, IDctrlDOF,
                                                                                    nSteps)

        if self.ui.radioButtonForced.isChecked():
            # maxDriftPiso, maxDriftTecho, maxDriftTecho2, VBasal_v, MaxDriftStory_v, DriftTecho_v, MP_ElemsForceS1, \
            #     MP_ElemsDeforS1, MP_ElemsForceS6, MP_ElemsDeforS6, DriftStory_m, PD_Beams, PD_Cols, PD_Walls\
            #     = singlePush(dref, mu, IDctrlNode, IDctrlDOF, nSteps)
            maxDriftPiso, maxDriftTecho, maxDriftTecho2, VBasal_v, MaxDriftStory_v, DriftTecho_v, DriftStory_m, \
                PD_Beams, PD_Cols, PD_Walls, Vu_Vn_Cols, ratio_infill = singlePush(dref, mu, IDctrlNode, IDctrlDOF,
                                                                                   nSteps)
        self.ui.progressBarPushover.hide()

        op.wipe()
        # Plot pushover curve
        # print('VBasal_v', VBasal_v)
        # print('Vu_Vn_Cols', Vu_Vn_Cols)
        fig = self.ui.PushCurve.canvas.axes
        fig.clear()
        ax = fig.add_axes([0.2, 0.2, 0.75, 0.75])
        if self.ui.checkBox_SDR.isChecked() == True:
            ax.plot(MaxDriftStory_v * 100, -VBasal_v / Wtotal, '-', c='red', linewidth=3, label=r'$SDR$', marker='None')
        if self.ui.checkBox_RDR.isChecked() == True:
            ax.plot(DriftTecho_v * 100, -VBasal_v / Wtotal, '-', c='blue', linewidth=3, label=r'$RDR$', marker='None')
        ax.set_xlabel('Drift %')
        ax.set_ylabel('Vb/Ws')
        ax.set_title('Pushover Curve')
        ax.legend(loc="lower right", fancybox=True, shadow=True, ncol=1, fontsize=8)
        ax.grid(True)
        self.ui.PushCurve.canvas.draw()
        self.ui.PushCurve.canvas.show()
        OutputPushFile = self.ui.OutputPushFile.text()
        if not os.path.exists("Pushover"):
            os.mkdir("Pushover")

        # PushData = np.concatenate((maxDriftTecho2*100, MaxDriftStory_v*100, -VBasal_v / Wtotal), axis=1)
        # print('PushData =', DriftTecho_v * 100)
        # print('PushData =', MaxDriftStory_v * 100)
        # print('PushData =', -VBasal_v / Wtotal)
        PushData = np.vstack((DriftTecho_v * 100, MaxDriftStory_v * 100, -VBasal_v / Wtotal))
        # print('PushData =', PushData.T)
        np.savetxt('Pushover/' + OutputPushFile + 'PushData.txt', PushData.T)
        # print('num_cols =', num_cols)
        # Reading of forces and deflections of beams and columns from recorders
        # print('M_BeamsForceS1', M_BeamsForceS1)
        # print('Tamñano de M_BeamsForceS1 =', np.shape(M_BeamsForceS1))
        # print('M_ColsForceS1', M_ColsForceS1)
        # print('Tamñano de M_ColsForceS1 =', np.shape(M_ColsForceS1))
        # beams_force_1 = np.loadtxt('Pushover/beams_force_1.out')
        # beams_def_1 = np.loadtxt('Pushover/beams_def_1.out')
        # beams_force_6 = np.loadtxt('Pushover/beams_force_6.out')
        # beams_def_6 = np.loadtxt('Pushover/beams_def_6.out')
        # cols_force_1 = np.loadtxt('Pushover/cols_force_1.out')
        # cols_def_1 = np.loadtxt('Pushover/cols_def_1.out')
        # cols_force_6 = np.loadtxt('Pushover/cols_force_6.out')
        # cols_def_6 = np.loadtxt('Pushover/cols_def_6.out')
        # print('cols_def_1', cols_def_1)
        if self.ui.radioButtonMVLEM.isChecked():
            walls_curvature = np.loadtxt('Pushover/walls_curvature.out')
            # print('walls_curvature', walls_curvature)
            HElev = np.diff(Loc_heigth)
            PD_Walls[:, 0] = np.absolute(walls_curvature[-1, 1:]) * HElev
            # print('PD_Walls', PD_Walls)
            # print('num_walls', num_walls)
            # PD_Walls = np.concatenate((PD_Walls.T,  np.zeros((num_walls, 1))), axis=1)
            # print('PD_Walls', PD_Walls)

        fy = float(self.ui.fy.text()) * MPa
        tw1, tw2 = float(self.ui.tw1.text()), float(self.ui.tw2.text())
        lw1, lw2 = float(self.ui.lw1.text()), float(self.ui.lw2.text())
        twIF = float(self.ui.twIF.text()) * mm
        Es = 200.0 * GPa
        ey = fy / Es
        # CD_Beams = np.zeros([num_beams, 2])  # Curvature Ductility - Beams
        # PRA_Beams = np.zeros([num_beams, 2])  # Plastic Rotation Angle - Beams
        # My_Beams = np.zeros([num_beams, 2])

        # Calculation of curvature ductility of beams and columns
        # for (ind, DB, DBPhl) in zip(range(num_beams), DataBeamDesign, DataBeamPhl):
        #     ets_beam_1 = M_BeamsDeforS1[:-2, 2 * ind] + M_BeamsDeforS1[:-2, 2 * ind + 1] * (DB.dt1 - DB.h / 2)
        #     ebs_beam_1 = M_BeamsDeforS1[:-2, 2 * ind] + M_BeamsDeforS1[:-2, 2 * ind + 1] * (DB.h / 2 - DB.db1)
        #     ets_beam_6 = M_BeamsDeforS6[:-2, 2 * ind] + M_BeamsDeforS6[:-2, 2 * ind + 1] * (DB.dt2 - DB.h / 2)
        #     ebs_beam_6 = M_BeamsDeforS6[:-2, 2 * ind] + M_BeamsDeforS6[:-2, 2 * ind + 1] * (DB.h / 2 - DB.db1)
        #     fi_1 = np.absolute(M_BeamsDeforS1[:-2, 2 * ind + 1])
        #     M_beam_1 = np.absolute(M_BeamsForceS1[:-2, 2 * ind + 1])
        #     fi_6 = np.absolute(M_BeamsDeforS6[:-2, 2 * ind + 1])
        #     M_beam_6 = np.absolute(M_BeamsForceS6[:-2, 2 * ind + 1])
        #
        #     # es_beam_1 = np.maximum(np.absolute(ets_beam_1), np.absolute(ebs_beam_1))
        #     # es_beam_6 = np.maximum(np.absolute(ets_beam_6), np.absolute(ebs_beam_6))
        #     # print('es_beam_1', es_beam_1, 'es_beam_6', es_beam_6)
        #     if np.max(ets_beam_1) <= ey and np.max(ets_beam_1) <= ey:
        #         CD_1 = 0
        #         My_1 = 0
        #         PRA1 = DBPhl.phl1 * np.max(fi_1)
        #     else:
        #         if np.max(ets_beam_1) >= ey:
        #             ft = interpolate.interp1d(ets_beam_1, M_beam_1, kind='nearest')
        #             My_1t = ft(ey)
        #         else:
        #             My_1t = float('inf')
        #         if np.max(ebs_beam_1) >= ey:
        #             fb = interpolate.interp1d(ebs_beam_1, M_beam_1, kind='nearest')
        #             My_1b = fb(ey)
        #         else:
        #             My_1b = float('inf')
        #         # print('ind', ind, 'My_1t', My_1t, 'My_1b', My_1b)
        #         My_1 = min(My_1t, My_1b)
        #         f = interpolate.interp1d(M_beam_1, fi_1, kind='nearest')
        #         fiy_1 = f(My_1)
        #         CD_1 = np.max(fi_1) / fiy_1
        #         PRA1 = DBPhl.phl1 * np.max(fi_1)
        #     if np.max(ets_beam_6) <= ey and np.max(ebs_beam_6) <= ey:
        #         CD_6 = 0
        #         My_6 = 0
        #         PRA6 = DBPhl.phl2 * np.max(fi_6)
        #     else:
        #         if np.max(ets_beam_6) >= ey:
        #             ft = interpolate.interp1d(ets_beam_6, M_beam_6, kind='nearest')
        #             My_6t = ft(ey)
        #         else:
        #             My_6t = float('inf')
        #         if np.max(ebs_beam_6) >= ey:
        #             fb = interpolate.interp1d(ebs_beam_6, M_beam_6, kind='nearest')
        #             My_6b = fb(ey)
        #         else:
        #             My_6b = float('inf')
        #         My_6 = min(My_6t, My_6b)
        #         f = interpolate.interp1d(M_beam_6, fi_6, kind='nearest')
        #         fiy_6 = f(My_6)
        #         CD_6 = np.max(fi_6) / fiy_6
        #         PRA6 = DBPhl.phl2 * np.max(fi_6)
        #     CD_Beams[ind, :] = [CD_1, CD_6]
        #     PRA_Beams[ind, :] = [PRA1, PRA6]
        #     My_Beams[ind, :] = [My_1, My_6]
        #
        #     # print('CD_Beams =', CD_Beams)
        #
        # CD_Cols = np.zeros([num_cols, 2])
        # PRA_Cols = np.zeros([num_cols, 2])
        # My_Cols = np.zeros([num_cols, 2])
        # for (ind, DC, DCPhl) in zip(range(num_cols), DataColDesign, DataColPhl):
        #     ets_col_1 = np.absolute(M_ColsDeforS1[:-2, 2 * ind] + M_ColsDeforS1[:-2, 2 * ind + 1] * (DC.d - DC.h / 2))
        #     ebs_col_1 = np.absolute(M_ColsDeforS1[:-2, 2 * ind] + M_ColsDeforS1[:-2, 2 * ind + 1] * (DC.h / 2 - DC.d))
        #     ets_col_6 = np.absolute(M_ColsDeforS6[:-2, 2 * ind] + M_ColsDeforS6[:-2, 2 * ind + 1] * (DC.d - DC.h / 2))
        #     ebs_col_6 = np.absolute(M_ColsDeforS6[:-2, 2 * ind] + M_ColsDeforS6[:-2, 2 * ind + 1] * (DC.h / 2 - DC.d))
        #     fi_1 = np.absolute(M_ColsDeforS1[:-2, 2 * ind + 1])
        #     M_col_1 = np.absolute(M_ColsForceS1[:-2, 2 * ind + 1])
        #     fi_6 = np.absolute(M_ColsDeforS6[:-2, 2 * ind + 1])
        #     M_col_6 = np.absolute(M_ColsForceS6[:-2, 2 * ind + 1])
        #     # es_col_1 = np.maximum(np.absolute(ets_col_1), np.absolute(ebs_col_1))
        #     # es_col_6 = np.maximum(np.absolute(ets_col_6), np.absolute(ebs_col_6))
        #     # print('es_col_1', es_col_1, 'es_col_6', es_col_6)
        #     if np.max(ets_col_1) <= ey and np.max(ebs_col_1) <= ey:
        #         CD_1 = 0
        #         Mfy_1 = 0
        #         PRA1 = DCPhl.phl1 * np.max(fi_1)
        #     else:
        #         if np.max(ets_col_1) >= ey:
        #             ft = interpolate.interp1d(ets_col_1, M_col_1, kind='nearest')
        #             Mfy_1t = ft(ey)
        #         else:
        #             Mfy_1t = float('inf')
        #         if np.max(ebs_col_1) >= ey:
        #             fb = interpolate.interp1d(ebs_col_1, M_col_1, kind='nearest')
        #             Mfy_1b = fb(ey)
        #         else:
        #             Mfy_1b = float('inf')
        #         Mfy_1 = min(Mfy_1t, Mfy_1b)
        #         f = interpolate.interp1d(M_col_1, fi_1, kind='nearest')
        #         fify_1 = f(Mfy_1)
        #         My_1 = np.max(M_col_1)
        # fiy_1 = My_1 / Mfy_1 * fify_1
        #         CD_1 = np.max(fi_1) / fiy_1
        #         PRA1 = DCPhl.phl1 * np.max(fi_1)
        #     if np.max(ets_col_6) <= ey and np.max(ebs_col_6) <= ey:
        #         CD_6 = 0
        #         Mfy_6 = 0
        #         PRA6 = DCPhl.phl2 * np.max(fi_6)
        #     else:
        #         if np.max(ets_col_6) >= ey:
        #             ft = interpolate.interp1d(ets_col_6, M_col_6, kind='nearest')
        #             Mfy_6t = ft(ey)
        #         else:
        #             Mfy_6t = float('inf')
        #         if np.max(ebs_col_6) >= ey:
        #             fb = interpolate.interp1d(ebs_col_6, M_col_6, kind='nearest')
        #             Mfy_6b = fb(ey)
        #         else:
        #             Mfy_6b = float('inf')
        #         Mfy_6 = min(Mfy_6t, Mfy_6b)
        #         f = interpolate.interp1d(M_col_6, fi_6, kind='nearest')
        #         fify_6 = f(Mfy_6)
        #         My_6 = np.max(M_col_6)
        #         fiy_6 = My_6 / Mfy_6 * fify_6
        #         CD_6 = np.max(fi_6) / fiy_6
        #         PRA6 = DCPhl.phl2 * np.max(fi_6)
        #     CD_Cols[ind, :] = [CD_1, CD_6]
        #     PRA_Cols[ind, :] = [PRA1, PRA6]
        #     My_Cols[ind, :] = [Mfy_1, Mfy_6]

        # print('CD_Cols =', CD_Cols)
        # CD_Ele = np.concatenate((CD_Cols, CD_Beams), axis=0)
        # PRA_Ele = np.concatenate((PRA_Cols, PRA_Beams), axis=0)
        # print('PD_Cols =', PD_Cols)
        # print('PD_Beams =', PD_Beams)
        PD_Ele = np.concatenate((PD_Cols, PD_Beams, PD_Walls), axis=0)
        PRA_Ele = PD_Ele  # the plastic rotation calculated with the opensees is taken
        # print('CD_Ele =', CD_Ele)
        # print('PRA_Ele =', PRA_Ele)

        # Drawing of curvature ductility in the plastic hinge projector
        # self.ui.PHP.canvas.figure.clear()
        fig = self.ui.PHP.canvas.axes

        if cbar:
            cbar.remove()
            fig.clear()
        # cmap = plt.cm.get_cmap("jet")
        Desp_x = np.loadtxt('Pushover/HoriNodes.out')
        Desp_y = np.loadtxt('Pushover/VertNodes.out')
        # print('Desp_x', Desp_x)
        num_nodes = np.shape(ListNodes)[0]
        Nodes_desp_x = ListNodes[:, 1] + 3 * Desp_x[-1, 1:num_nodes + 1]
        Nodes_desp_y = ListNodes[:, 2] + 3 * Desp_y[-1, 1:num_nodes + 1]
        # self.figure = plt.figure()
        ax1 = fig.add_axes([0, 0.1, 1, 0.9])
        ax1.plot(Nodes_desp_x, Nodes_desp_y, 'ks')
        ax1.axis('off')
        if tw1 != 0 and tw2 != 0:
            sepW = 2 + lw1
        elif tw1 == 0 and tw2 != 0:
            numW1 = 0
            sepW = 0.
        if tw1 != 0:
            Desp_xW1 = np.loadtxt('Pushover/HoriNodesW1.out')
            Desp_yW1 = np.loadtxt('Pushover/VertNodesW1.out')
            # print('Desp_xW1', Desp_xW1)
            numW1 = len(ListNodesW1[:, 0])
            Nodes_desp_xW1 = 2 + ListNodesW1[:, 1] + 3 * Desp_xW1[-1, 1:numW1 + 1]
            Nodes_desp_yW1 = ListNodesW1[:, 2] + 3 * Desp_yW1[-1, 1:numW1 + 1]
        if tw2 != 0:
            Desp_xW2 = np.loadtxt('Pushover/HoriNodesW2.out')
            Desp_yW2 = np.loadtxt('Pushover/VertNodesW2.out')
            numW2 = len(ListNodesW2[:, 0])
            Nodes_desp_xW2 = 2 + ListNodesW2[:, 1] + sepW + 3 * Desp_xW2[-1, 1:numW2 + 1]
            Nodes_desp_yW2 = ListNodesW2[:, 2] + 3 * Desp_yW2[-1, 1:numW2 + 1]
        fpos = 0.1
        fsize = 1
        DataDC = []
        DataPRA = []
        for (ind, Ele, IFE, DC) in zip(range(num_cols), EleCol, infill_m, DataColDesign):
            xi = Nodes_desp_x[Ele.Nod_ini]
            yi = Nodes_desp_y[Ele.Nod_ini]
            xe = Nodes_desp_x[Ele.Nod_end]
            ye = Nodes_desp_y[Ele.Nod_end]
            x = np.array([xi, xe])
            y = np.array([yi, ye])
            # print(xi, yi, xe, ye)
            ax1.plot(x, y, 'k-', alpha=.3)
            Delta_x = xe - xi
            Delta_y = ye - yi
            xi_CD = xi + fpos * Delta_x
            yi_CD = yi + fpos * Delta_y
            xe_CD = xe - fpos * Delta_x
            ye_CD = ye - fpos * Delta_y
            # CD_i = CD_Ele[Ele.EleTag - 1, 0]
            # CD_e = CD_Ele[Ele.EleTag - 1, 1]
            PRA_i = PRA_Ele[Ele.EleTag - 1, 0]
            PRA_e = PRA_Ele[Ele.EleTag - 1, 1]
            # DataDC.append(DuctilityCurve(xi_CD, xe_CD, yi_CD, ye_CD, fsize * CD_i, fsize * CD_e))
            DataPRA.append(PlasticRotationAngle(xi_CD, xe_CD, yi_CD, ye_CD, PRA_i, PRA_e))
            #Plot Shear Arrow
            xi_SA = xi + 0.5 * Delta_x + Vu_Vn_Cols[ind] / 2
            yi_SA = yi + 0.5 * Delta_y
            if Vu_Vn_Cols[ind] <= 0.7:
                colorArrow = 'green'
            elif 0.7 < Vu_Vn_Cols[ind] <= 1.0:
                colorArrow = 'orange'
            else:
                colorArrow = 'red'
            ax1.arrow(xi_SA, yi_SA, -Vu_Vn_Cols[ind], 0, head_width=0.2, width=0.05, ec=colorArrow)
            ax1.text(xi_SA, yi_SA, r'{:.2f}'.format(Vu_Vn_Cols[ind]), style='italic', fontsize=8)

            if IFE != 0:
                if ratio_infill[ind] <= 0.4:
                    colorInfill = 'green'
                elif 0.4 < ratio_infill[ind] <= 0.7:
                    colorInfill = 'orange'
                else:
                    colorInfill = 'red'
                ind_ini, ind_end = Ele.Nod_end, Ele.Nod_ini + 1
                xid, yid = Nodes_desp_x[ind_ini], Nodes_desp_y[ind_ini]
                xed, yed = Nodes_desp_x[ind_end], Nodes_desp_y[ind_end]
                ind_ini2, ind_end2 = Ele.Nod_ini, Ele.Nod_end + 1
                xid2, yid2 = Nodes_desp_x[ind_ini2], Nodes_desp_y[ind_ini2]
                xed2, yed2 = Nodes_desp_x[ind_end2], Nodes_desp_y[ind_end2]
                fposInfill = DC.wIF/Ele.LEle
                Delta_x1 = xid - xid2
                Delta_y1 = yid - yid2
                Delta_x2 = xed2 - xed
                Delta_y2 = yed2 - yed
                lat = np.array([xid, xid - fposInfill * Delta_x1, xed, xed + fposInfill * Delta_x2])
                long = np.array([yid, yid - fposInfill * Delta_y1, yed, yed + fposInfill * Delta_y2])
                ax1.add_patch(Polygon(np.c_[lat, long], facecolor=colorInfill, alpha=0.8))
                theta = degrees(atan((xed + fposInfill * Delta_x2 - xid)/(yed + fposInfill * Delta_y2 - yid))) + 90
                ax1.text((xid+xed + fposInfill * Delta_x2)/2, (yid+yed + fposInfill * Delta_y2)/2 -
                         (fposInfill * Delta_y2) / 2, r'{:.2f}'.format(ratio_infill[ind]), rotation=-theta,
                         style='italic', rotation_mode='anchor', transform_rotates_text=True, fontsize=8,
                         horizontalalignment='center', verticalalignment='center')

        for Ele in EleBeam:
            xi = Nodes_desp_x[Ele.Nod_ini]
            yi = Nodes_desp_y[Ele.Nod_ini]
            xe = Nodes_desp_x[Ele.Nod_end]
            ye = Nodes_desp_y[Ele.Nod_end]
            x = np.array([xi, xe])
            y = np.array([yi, ye])
            # print(xi, yi, xe, ye)
            ax1.plot(x, y, 'k-', alpha=.3)
            Delta_x = xe - xi
            Delta_y = ye - yi
            xi_CD = xi + fpos * Delta_x
            yi_CD = yi + fpos * Delta_y
            xe_CD = xe - fpos * Delta_x
            ye_CD = ye - fpos * Delta_y
            # CD_i = CD_Ele[Ele.EleTag - 1, 0]
            # CD_e = CD_Ele[Ele.EleTag - 1, 1]
            PRA_i = PRA_Ele[Ele.EleTag - 1, 0]
            PRA_e = PRA_Ele[Ele.EleTag - 1, 1]
            # DataDC.append(DuctilityCurve(xi_CD, xe_CD, yi_CD, ye_CD, fsize * CD_i, fsize * CD_e))
            DataPRA.append(PlasticRotationAngle(xi_CD, xe_CD, yi_CD, ye_CD, PRA_i, PRA_e))

        for Ele in EleWall:
            if Ele.EleTag in ListEleTagW1:
                # print('Ele.EleTag', Ele.EleTag)
                ind_ini, ind_end = Ele.Nod_ini - num_nodes - nnt, Ele.Nod_end - num_nodes
                xid, yid = Nodes_desp_xW1[ind_ini], Nodes_desp_yW1[ind_ini]
                xed, yed = Nodes_desp_xW1[ind_end], Nodes_desp_yW1[ind_end]
                lat = np.array([xid, xid + lw1, xed + lw1, xed, xid])
                long = np.array([yid, yid, yed, yed, yid])
                ax1.add_patch(Polygon(np.c_[lat, long], facecolor='black', alpha=0.1))
                Delta_x = xed - xid
                Delta_y = yed - yid
                xi_CD = xid + lw1 / 2 + fpos * Delta_x
                yi_CD = yid + fpos * Delta_y
                xe_CD = xed + lw1 / 2 - fpos * Delta_x
                ye_CD = yed - fpos * Delta_y
                PRA_i = PRA_Ele[Ele.EleTag - 1, 0]
                PRA_e = PRA_Ele[Ele.EleTag - 1, 1]
                DataPRA.append(PlasticRotationAngle(xi_CD, xe_CD, yi_CD, ye_CD, PRA_i, PRA_e))

            if Ele.EleTag in ListEleTagW2:
                ind_ini, ind_end = Ele.Nod_ini - num_nodes - numW1 - nnt, Ele.Nod_end - num_nodes - numW1
                xid, yid = Nodes_desp_xW2[ind_ini], Nodes_desp_yW2[ind_ini]
                xed, yed = Nodes_desp_xW2[ind_end], Nodes_desp_yW2[ind_end]
                lat = np.array([xid, xid + lw2, xed + lw2, xed, xid])
                long = np.array([yid, yid, yed, yed, yid])
                ax1.add_patch(Polygon(np.c_[lat, long], facecolor='black', alpha=0.1))
                Delta_x = xed - xid
                Delta_y = yed - yid
                xi_CD = xid + lw2 / 2 + fpos * Delta_x
                yi_CD = yid + fpos * Delta_y
                xe_CD = xed + lw2 / 2 - fpos * Delta_x
                ye_CD = yed - fpos * Delta_y
                PRA_i = PRA_Ele[Ele.EleTag - 1, 0]
                PRA_e = PRA_Ele[Ele.EleTag - 1, 1]
                DataPRA.append(PlasticRotationAngle(xi_CD, xe_CD, yi_CD, ye_CD, PRA_i, PRA_e))
        ax1.axis('equal')
        PRA_x, PRA_y, PRA_size = [], [], []
        for PRA in DataPRA:
            PRA_x.append([PRA.xi, PRA.xe])
            PRA_y.append([PRA.yi, PRA.ye])
            PRA_size.append([PRA.PRA_i, PRA.PRA_e])
        PRA_x = np.array(PRA_x)
        PRA_x = PRA_x.flatten()
        PRA_y = np.array(PRA_y)
        PRA_y = PRA_y.flatten()
        PRA_size = np.array(PRA_size)
        PRA_size = PRA_size.flatten()

        if self.ui.radioButtonPHP_RA.isChecked():
            cmap = plt.cm.get_cmap('jet')

            sc = ax1.scatter(PRA_x, PRA_y, s=3000 * PRA_size, c=PRA_size, cmap=cmap, alpha=1)

            ax2 = fig.add_axes([0.05, 0.1, 0.9, 0.03])
            norm = plt.Normalize(vmin=np.min(PRA_size), vmax=np.max(PRA_size))
            cbar = colorbar.ColorbarBase(ax2, cmap=cmap, norm=norm, orientation='horizontal',
                                         label='Rotation angle (radians)')
            self.ui.PHP.canvas.draw()
            self.ui.PHP.canvas.show()

        PRA_size = []
        for PRA in DataPRA:
            PRA_size.append([PRA.PRA_i, PRA.PRA_e])
        PRA_size = np.array(PRA_size)
        if self.ui.radioButtonAC.isChecked():
            # Calculation of plastic rotation angles (radians) - Performance levels
            ecu = 0.003
            ro_bal = 0.85 * fcB / fy * ecu / (ecu + fy / Es) * Beta1B
            dst = 3 / 8 * inch
            Ast = pi * dst ** 2 / 4.  # area de la barra del estribo
            DataAC = []

            Cd = float(self.ui.Cd.text())
            for (ind, DC, Ele) in zip(range(num_cols), DataColDesign, EleCol):
                if Cd <= 2.0:
                    knl = 1.0
                elif Cd >= 6.0:
                    knl = 0.7
                else:
                    knl = np.interp(Cd, [2.0, 6.0], [1.0, 0.7])
                if DC.sst / DC.d <= 0.75:
                    alfa_col = 1.0
                elif DC.sst / DC.d >= 1.0:
                    alfa_col = 0.0
                else:
                    alfa_col = np.interp(DC.sst / DC.d, [0.75, 1.0], [1.0, 0.0])
                ro_t = min(max(0.075, Ast * DC.nsB / (DC.b * DC.d)), 0.0175)
                yi = ListNodes[Ele.Nod_ini, 2]
                ye = ListNodes[Ele.Nod_end, 2]
                LCol = ye - yi
                VColOE = knl * (alfa_col * Ast * DC.nsB * fy * DC.d / DC.sst +
                                0.5 * sqrt(fcC * MPa) / (min(max(2, DC.MUD1 / DC.VUD1), 4) * DC.d) *
                                sqrt(1 + DC.NUD1 / (0.5 * DC.b * DC.h * sqrt(fcC * MPa))) * 0.8 * DC.b * DC.h)
                Vy = DC.Vy
                a = max(0.042 - 0.043 * DC.NUD1 / (DC.b * DC.h * fcC) + 0.63 * ro_t - 0.023 * max(Vy / VColOE, 0.2), 0)
                bc = max(0.5 / (5 + DC.NUD1 / (0.8 * DC.b * DC.h * fcC) * 1 / ro_t * fcC / fy) - 0.01, a)
                if 0.1 <= DC.NUD1 / (DC.b * DC.h * fcC) <= 0.5:
                    b = bc
                elif 0.5 < DC.NUD1 / (DC.b * DC.h * fcC) <= 0.7:
                    b = np.interp(DC.NUD1 / (DC.b * DC.h * fcC), [0.5, 0.7], [bc, max(0, a)])
                elif DC.NUD1 / (DC.b * DC.h * fcC) > 0.7:
                    b = a
                elif DC.NUD1 / (DC.b * DC.h * fcC) < 0.1:
                    b = max(0.5 / (5 + 0.1 / 0.8 * 1 / ro_t * fcC / fy) - 0.01, a)
                IO_1 = min(0.15 * a, 0.005)
                LS_1 = 0.5 * b
                CP_1 = 0.7 * b

                VColOE = knl * (alfa_col * Ast * DC.nsB * fy * DC.d / DC.sst +
                                0.5 * sqrt(fcC * MPa) / (min(max(2, DC.MUD2 / DC.VUD2), 4) / DC.d) *
                                sqrt(1 + DC.NUD2 / (0.5 * DC.b * DC.h * sqrt(fcC * MPa))) * 0.8 * DC.b * DC.h)
                # Vy = 2 * My_Cols[ind, 1] / LCol
                a = max(0.042 - 0.043 * DC.NUD2 / (DC.b * DC.h * fcC) + 0.63 * ro_t - 0.023 * max(Vy / VColOE, 0.2), 0)
                bc = max(0.5 / (5 + DC.NUD2 / (0.8 * DC.b * DC.h * fcC) * 1 / ro_t * fcC / fy) - 0.01, a)
                if 0.1 <= DC.NUD2 / (DC.b * DC.h * fcC) <= 0.5:
                    b = bc
                elif 0.5 < DC.NUD2 / (DC.b * DC.h * fcC) <= 0.7:
                    b = np.interp(DC.NUD1 / (DC.b * DC.h * fcC), [0.5, 0.7], [bc, max(0, a)])
                elif DC.NUD2 / (DC.b * DC.h * fcC) > 0.7:
                    b = a
                elif DC.NUD2 / (DC.b * DC.h * fcC) < 0.1:
                    b = max(0.5 / (5 + 0.1 / 0.8 * 1 / ro_t * fcC / fy) - 0.01, a)
                IO_2 = min(0.15 * a, 0.005)
                LS_2 = 0.5 * b
                CP_2 = 0.7 * b
                # print('IO_1, LS_1, CP_1, IO_2, LS_2, CP_2', IO_1, LS_1, CP_1, IO_2, LS_2, CP_2)
                DataAC.append(AcceptanceCriteria(IO_1, LS_1, CP_1, IO_2, LS_2, CP_2))

            for DB in DataBeamDesign:
                ro1 = DB.Asb1 / DB.b / DB.db1
                ro1p = DB.Ast1 / DB.b / DB.db1
                Vs = DB.ns1 * Ast * fy * DB.db1 / DB.ss1
                if DB.VU1 / (DB.b * DB.db1 * sqrt(fcB * MPa)) <= 0.25:
                    IO1 = 0.010
                    LS1 = 0.025
                    CP1 = 0.050
                elif DB.VU1 / (DB.b * DB.db1 * sqrt(fcB * MPa)) >= 0.50:
                    IO1 = 0.005
                    LS1 = 0.020
                    CP1 = 0.040
                else:
                    IO1 = np.interp(DB.VU1 / (DB.b * DB.db1 * sqrt(fcB * MPa)), [0.25, 0.5], [0.010, 0.005])
                    LS1 = np.interp(DB.VU1 / (DB.b * DB.db1 * sqrt(fcB * MPa)), [0.25, 0.5], [0.025, 0.020])
                    CP1 = np.interp(DB.VU1 / (DB.b * DB.db1 * sqrt(fcB * MPa)), [0.25, 0.5], [0.050, 0.040])
                if DB.VU1 / (DB.b * DB.db1 * sqrt(fcB * MPa)) <= 0.25:
                    IO2 = 0.005
                    LS2 = 0.020
                    CP2 = 0.030
                elif DB.VU1 / (DB.b * DB.db1 * sqrt(fcB * MPa)) >= 0.50:
                    IO2 = 0.005
                    LS2 = 0.015
                    CP2 = 0.020
                else:
                    IO2 = np.interp(DB.VU1 / (DB.b * DB.db1 * sqrt(fcB * MPa)), [0.25, 0.5], [0.005, 0.005])
                    LS2 = np.interp(DB.VU1 / (DB.b * DB.db1 * sqrt(fcB * MPa)), [0.25, 0.5], [0.020, 0.015])
                    CP2 = np.interp(DB.VU1 / (DB.b * DB.db1 * sqrt(fcB * MPa)), [0.25, 0.5], [0.030, 0.020])
                if DB.VU1 / (DB.b * DB.db1 * sqrt(fcB * MPa)) <= 0.25:
                    IO3 = 0.005
                    LS3 = 0.020
                    CP3 = 0.030
                elif DB.VU1 / (DB.b * DB.db1 * sqrt(fcB * MPa)) >= 0.50:
                    IO3 = 0.0015
                    LS3 = 0.010
                    CP3 = 0.015
                else:
                    IO3 = np.interp(DB.VU1 / (DB.b * DB.db1 * sqrt(fcB * MPa)), [0.25, 0.5], [0.005, 0.0015])
                    LS3 = np.interp(DB.VU1 / (DB.b * DB.db1 * sqrt(fcB * MPa)), [0.25, 0.5], [0.020, 0.010])
                    CP3 = np.interp(DB.VU1 / (DB.b * DB.db1 * sqrt(fcB * MPa)), [0.25, 0.5], [0.030, 0.015])
                if DB.VU1 / (DB.b * DB.db1 * sqrt(fcB * MPa)) <= 0.25:
                    IO4 = 0.005
                    LS4 = 0.010
                    CP4 = 0.015
                elif DB.VU1 / (DB.b * DB.db1 * sqrt(fcB * MPa)) >= 0.50:
                    IO4 = 0.0015
                    LS4 = 0.005
                    CP4 = 0.010
                else:
                    IO4 = np.interp(DB.VU1 / (DB.b * DB.db1 * sqrt(fcB * MPa)), [0.25, 0.5], [0.005, 0.0015])
                    LS4 = np.interp(DB.VU1 / (DB.b * DB.db1 * sqrt(fcB * MPa)), [0.25, 0.5], [0.010, 0.005])
                    CP4 = np.interp(DB.VU1 / (DB.b * DB.db1 * sqrt(fcB * MPa)), [0.25, 0.5], [0.015, 0.010])
                if DB.ss1 <= DB.db1 / 3 and Vs >= 3 / 4 * DB.VU1:
                    if (ro1 - ro1p) / ro_bal <= 0:
                        IO_1 = IO1
                        LS_1 = LS1
                        CP_1 = CP1
                    elif (ro1 - ro1p) / ro_bal >= 0.5:
                        IO_1 = IO2
                        LS_1 = LS2
                        CP_1 = CP2
                    else:
                        IO_1 = np.interp((ro1 - ro1p) / ro_bal, [0.0, 0.5], [IO1, IO2])
                        LS_1 = np.interp((ro1 - ro1p) / ro_bal, [0.0, 0.5], [LS1, LS2])
                        CP_1 = np.interp((ro1 - ro1p) / ro_bal, [0.0, 0.5], [CP1, CP2])
                if DB.ss1 > DB.db1 / 3 or Vs < 3 / 4 * DB.VU1:
                    if (ro1 - ro1p) / ro_bal <= 0:
                        IO_1 = IO3
                        LS_1 = LS3
                        CP_1 = CP3
                    elif (ro1 - ro1p) / ro_bal >= 0.5:
                        IO_1 = IO4
                        LS_1 = LS4
                        CP_1 = CP4
                    else:
                        IO_1 = np.interp((ro1 - ro1p) / ro_bal, [0.0, 0.5], [IO3, IO4])
                        LS_1 = np.interp((ro1 - ro1p) / ro_bal, [0.0, 0.5], [LS3, LS4])
                        CP_1 = np.interp((ro1 - ro1p) / ro_bal, [0.0, 0.5], [CP3, CP4])
                ro2 = DB.Ast2 / DB.b / DB.db2
                ro2p = DB.Asb2 / DB.b / DB.db2
                Vs = DB.ns2 * Ast * fy * DB.db2 / DB.ss2
                if DB.VU2 / (DB.b * DB.db2 * sqrt(fcB * MPa)) <= 0.25:
                    IO1 = 0.010
                    LS1 = 0.025
                    CP1 = 0.050
                elif DB.VU1 / (DB.b * DB.db2 * sqrt(fcB * MPa)) >= 0.50:
                    IO1 = 0.005
                    LS1 = 0.020
                    CP1 = 0.040
                else:
                    IO1 = np.interp(DB.VU1 / (DB.b * DB.db2 * sqrt(fcB * MPa)), [0.25, 0.5], [0.010, 0.005])
                    LS1 = np.interp(DB.VU1 / (DB.b * DB.db2 * sqrt(fcB * MPa)), [0.25, 0.5], [0.025, 0.020])
                    CP1 = np.interp(DB.VU1 / (DB.b * DB.db2 * sqrt(fcB * MPa)), [0.25, 0.5], [0.050, 0.040])
                if DB.VU1 / (DB.b * DB.db2 * sqrt(fcB * MPa)) <= 0.25:
                    IO2 = 0.005
                    LS2 = 0.020
                    CP2 = 0.030
                elif DB.VU1 / (DB.b * DB.db2 * sqrt(fcB * MPa)) >= 0.50:
                    IO2 = 0.005
                    LS2 = 0.015
                    CP2 = 0.020
                else:
                    IO2 = np.interp(DB.VU1 / (DB.b * DB.db2 * sqrt(fcB * MPa)), [0.25, 0.5], [0.005, 0.005])
                    LS2 = np.interp(DB.VU1 / (DB.b * DB.db2 * sqrt(fcB * MPa)), [0.25, 0.5], [0.020, 0.015])
                    CP2 = np.interp(DB.VU1 / (DB.b * DB.db2 * sqrt(fcB * MPa)), [0.25, 0.5], [0.030, 0.020])
                if DB.VU1 / (DB.b * DB.db2 * sqrt(fcB * MPa)) <= 0.25:
                    IO3 = 0.005
                    LS3 = 0.020
                    CP3 = 0.030
                elif DB.VU1 / (DB.b * DB.db2 * sqrt(fcB * MPa)) >= 0.50:
                    IO3 = 0.0015
                    LS3 = 0.010
                    CP3 = 0.015
                else:
                    IO3 = np.interp(DB.VU1 / (DB.b * DB.db2 * sqrt(fcB * MPa)), [0.25, 0.5], [0.005, 0.0015])
                    LS3 = np.interp(DB.VU1 / (DB.b * DB.db2 * sqrt(fcB * MPa)), [0.25, 0.5], [0.020, 0.010])
                    CP3 = np.interp(DB.VU1 / (DB.b * DB.db2 * sqrt(fcB * MPa)), [0.25, 0.5], [0.030, 0.015])
                if DB.VU1 / (DB.b * DB.db2 * sqrt(fcB * MPa)) <= 0.25:
                    IO4 = 0.005
                    LS4 = 0.010
                    CP4 = 0.015
                elif DB.VU1 / (DB.b * DB.db2 * sqrt(fcB * MPa)) >= 0.50:
                    IO4 = 0.0015
                    LS4 = 0.005
                    CP4 = 0.010
                else:
                    IO4 = np.interp(DB.VU1 / (DB.b * DB.db2 * sqrt(fcB * MPa)), [0.25, 0.5], [0.005, 0.0015])
                    LS4 = np.interp(DB.VU1 / (DB.b * DB.db2 * sqrt(fcB * MPa)), [0.25, 0.5], [0.010, 0.005])
                    CP4 = np.interp(DB.VU1 / (DB.b * DB.db2 * sqrt(fcB * MPa)), [0.25, 0.5], [0.015, 0.010])
                if DB.ss2 <= DB.db2 / 3 and Vs >= 3 / 4 * DB.VU2:
                    if (ro2 - ro2p) / ro_bal <= 0:
                        IO_2 = IO1
                        LS_2 = LS1
                        CP_2 = CP1
                    elif (ro2 - ro2p) / ro_bal >= 0.5:
                        IO_2 = IO2
                        LS_2 = LS2
                        CP_2 = CP2
                    else:
                        IO_2 = np.interp((ro2 - ro2p) / ro_bal, [0.0, 0.5], [IO1, IO2])
                        LS_2 = np.interp((ro2 - ro2p) / ro_bal, [0.0, 0.5], [LS1, LS2])
                        CP_2 = np.interp((ro2 - ro2p) / ro_bal, [0.0, 0.5], [CP1, CP2])
                if DB.ss1 > DB.db2 / 3 or Vs < 3 / 4 * DB.VU1:
                    if (ro2 - ro2p) / ro_bal <= 0:
                        IO_2 = IO3
                        LS_2 = LS3
                        CP_2 = CP3
                    elif (ro2 - ro2p) / ro_bal >= 0.5:
                        IO_2 = IO4
                        LS_2 = LS4
                        CP_2 = CP4
                    else:
                        IO_2 = np.interp((ro2 - ro2p) / ro_bal, [0.0, 0.5], [IO3, IO4])
                        LS_2 = np.interp((ro2 - ro2p) / ro_bal, [0.0, 0.5], [LS3, LS4])
                        CP_2 = np.interp((ro2 - ro2p) / ro_bal, [0.0, 0.5], [CP3, CP4])
                # print('IO_1, LS_1, CP_1, IO_2, LS_2, CP_2', IO_1, LS_1, CP_1, IO_2, LS_2, CP_2)
                DataAC.append(AcceptanceCriteria(IO_1, LS_1, CP_1, IO_2, LS_2, CP_2))

            for (ind, DW, Ele) in zip(range(num_walls), DataWallDesign, EleWall):
                P_Rel = DW.NUD1 / (DW.b * DW.h * fcC)
                V_Rel = DW.VUD1 / (DW.b * DW.h * sqrt(fcC * MPa))
                # if P_Rel <= 0.1:
                if V_Rel <= 1 / 3:
                    IO1, LS1, CP1 = 0.002, 0.008, 0.015
                elif V_Rel >= 1 / 2:
                    IO1, LS1, CP1 = 0.002, 0.006, 0.010
                else:
                    IO1 = np.interp(V_Rel, [1 / 3, 1 / 2], [0.002, 0.002])
                    LS1 = np.interp(V_Rel, [1 / 3, 1 / 2], [0.008, 0.006])
                    CP1 = np.interp(V_Rel, [1 / 3, 1 / 2], [0.015, 0.010])
                if V_Rel <= 1 / 3:
                    IO2, LS2, CP2 = 0.001, 0.003, 0.005
                elif V_Rel >= 1 / 2:
                    IO2, LS2, CP2 = 0.001, 0.002, 0.004
                else:
                    IO2 = np.interp(V_Rel, [1 / 3, 1 / 2], [0.001, 0.001])
                    LS2 = np.interp(V_Rel, [1 / 3, 1 / 2], [0.003, 0.002])
                    CP2 = np.interp(V_Rel, [1 / 3, 1 / 2], [0.005, 0.004])
                if P_Rel <= 0.1:
                    IO_1 = IO1
                    LS_1 = LS1
                    CP_1 = CP1
                elif P_Rel >= 0.25:
                    IO_1 = IO2
                    LS_1 = LS2
                    CP_1 = CP2
                else:
                    IO_1 = np.interp(P_Rel, [0.10, 0.25], [IO1, IO2])
                    LS_1 = np.interp(P_Rel, [0.10, 0.25], [LS1, LS2])
                    CP_1 = np.interp(P_Rel, [0.10, 0.25], [CP1, CP2])
                # print('IO_1, LS_1, CP_1, IO_2, LS_2, CP_2', IO_1, LS_1, CP_1, IO_2, LS_2, CP_2)
                DataAC.append(AcceptanceCriteria(IO_1, LS_1, CP_1, IO_2, LS_2, CP_2))

            # print('DataAC =\n', DataAC)

            AC_size = []
            for AC in DataAC:
                AC_size.append([AC.IO_1, AC.LS_1, AC.CP_1, AC.IO_2, AC.LS_2, AC.CP_2])
            AC_size = np.array(AC_size)
            np.savetxt('Pushover/' + OutputPushFile + 'AC_Data.txt', AC_size, fmt='%0.5f')

            # print('AC_size', AC_size)
            Color_PRA = []
            for ind in range(num_elems - 1):
                # print('ind =', ind)
                if PRA_size[ind, 0] <= AC_size[ind, 0]:
                    Color_PRA.append('green')
                elif AC_size[ind, 0] < PRA_size[ind, 0] <= AC_size[ind, 1]:
                    Color_PRA.append('yellow')
                elif AC_size[ind, 1] < PRA_size[ind, 0] <= AC_size[ind, 2]:
                    Color_PRA.append('orange')
                elif PRA_size[ind, 0] > AC_size[ind, 2]:
                    Color_PRA.append('red')
                if PRA_size[ind, 1] <= AC_size[ind, 3]:
                    Color_PRA.append('green')
                elif AC_size[ind, 3] < PRA_size[ind, 1] <= AC_size[ind, 4]:
                    Color_PRA.append('yellow')
                elif AC_size[ind, 4] < PRA_size[ind, 1] <= AC_size[ind, 5]:
                    Color_PRA.append('orange')
                elif PRA_size[ind, 1] > AC_size[ind, 5]:
                    Color_PRA.append('red')
            # cm = plt.cm.get_cmap('jet')

            from matplotlib.colors import ListedColormap
            # cMap = ListedColormap(['green', 'yellow', 'orange', 'red'])

            cmap = colors.ListedColormap(['green', 'yellow', 'orange', 'red'])
            norm = colors.BoundaryNorm([0, 1, 2, 3, 4], cmap.N)

            PRA_size = PRA_size.flatten()
            ax1.scatter(PRA_x, PRA_y, s=30, c=Color_PRA, cmap=cmap, norm=norm, alpha=1)
            # divider = make_axes_locatable(self.ui.PHP.canvas.axes)
            # cax = divider.new_vertical(size="5%", pad=0.5, pack_start=True)
            ax2 = fig.add_axes([0.05, 0.1, 0.9, 0.03])
            # cmap = colors.ListedColormap(['green', 'yellow', 'orange', 'red'])
            # norm = colors.BoundaryNorm([0, 1, 2, 3, 4], cmap.N)
            cmap.set_over('0.25')
            cmap.set_under('0.75')
            bounds = [0, 1, 2, 3, 4]
            boundsTick = ['0', 'IO', 'LS', 'CP']
            norm = colors.BoundaryNorm(bounds, cmap.N)
            cbar = colorbar.ColorbarBase(ax2, cmap=cmap, norm=norm,
                                         orientation='horizontal')  # , ticks=boundsTick, spacing='proportional',
            # )
            cbar.set_label('Acceptance Criteria ASCE 41-17')
            #
            # cbar = self.ui.PHP.canvas.figure.colorbar(mappable=sc, orientation="horizontal") #, aspect=80)
            # cbar.set_label('Curvature ductility $\mu_{\phi}=\phi_{u}/\phi_{y}$')
            # # cbar.set_ticks([0, 1, 2, 3, 4])
            cbar.set_ticklabels(['', 'IO', 'LS', 'CP', ''])

            # cbar.ax2.text((ii - vmin) / (vmax - vmin), 1.5, str(ii), transform=cbar.ax.transAxes, va='bottom',
            #              ha='center')
            self.ui.PHP.canvas.draw()
            self.ui.PHP.canvas.show()

            # self.ui.PHP.canvas.axes.set_position([0.2, 0.3, 0.9, 0.9])  # left,bottom,width,height
            # self.ui.tabWidget.setCurrentIndex(3)

        # Plot of the lower moment-curvature diagram of the first column
        indb = 1  # + 2 * (naxes - 1)
        indc = 2

        # fig = self.ui.MomFi.canvas.axes
        # fig.clear()
        # ax = fig.add_axes([0.1, 0.2, 0.85, 0.7])
        # ax.plot(np.absolute(M_ColsDeforS6[:-2, 2*(indc-1)+1]), np.absolute(M_ColsForceS6[:-2, 2*(indc-1)+1]), 'b.-')
        # ax.plot(np.absolute(M_BeamsDeforS6[:-2, 2*(indb-1)+1]), np.absolute(M_BeamsForceS6[:-2, 2*(indb-1)+1]), 'r.-')
        # ax.set_xlabel('Curvature (1/m)')
        # ax.set_ylabel('Moment (kN*m)')
        # ax.set_title('Column')
        # ax.grid(True)
        # self.ui.MomFi.canvas.draw()
        # self.ui.MomFi.canvas.show()
        # self.CreateNLM()
        # self.ui.tabWidget.setCurrentIndex(3)

        # fig = self.ui.MomFi.canvas.axes
        # fig.clear()
        # ax = fig.add_axes([0.1, 0.2, 0.85, 0.7])
        # ax.plot(MaxDriftStory_v*100, T_v[0, :], 'b.-', label='T1')
        # ax.plot(MaxDriftStory_v*100, T_v[1, :], 'r.-', label='T2')
        # ax.set_xlabel('Maximum Story Drift %')
        # ax.set_ylabel('Period (sec)')
        # ax.set_title('Period variation curve')
        # ax.grid(True)
        # ax.legend(loc='upper left')
        # self.ui.MomFi.canvas.draw()
        # self.ui.MomFi.canvas.show()

        AlphaCurves = self.ui.AlphaCurves.text()
        AlphaCurves = AlphaCurves.split(',')
        AlphaCurves = np.array(AlphaCurves, dtype=float) / 100
        if AlphaCurves[-1] > MaxDriftStory_v[-1]:
            AlphaCurves[-1] = np.round(MaxDriftStory_v[-1] * 100, 1) / 100

        nSDR = AlphaCurves.size
        Alpha = np.zeros((nSDR, floors_num - 1))

        for jj in range(nSDR):
            RDR = np.interp(AlphaCurves[jj], MaxDriftStory_v[1:], DriftTecho_v[1:])
            for ii in range(floors_num - 1):
                SDR = np.interp(AlphaCurves[jj], MaxDriftStory_v[1:], DriftStory_m[1:, ii])
                Alpha[jj, ii] = SDR / RDR

        Floors = np.arange(floors_num - 1) + 1

        fig = self.ui.AlphaCurve.canvas.axes
        fig.clear()
        ax = fig.add_axes([0.25, 0.2, 0.7, 0.7])

        markers = ['^', 's', 'p', 'h', '8']

        for jj in range(nSDR):
            ax.plot(Alpha[jj, :], Floors, marker=markers[jj], label=r'$SDR_{max} = %1.1f$' % (AlphaCurves[jj] * 100))

        ax.set_xlabel(r'$Alpha_j = SDR_j/RDR$', fontsize=8)
        ax.set_ylabel(r'$Story_j$', rotation='vertical', fontsize=8)
        # ax.yaxis.set_label_coords(0.1, 1.02)
        ax.set_xlim(0, None)
        ax.set_ylim(1 - .1, floors_num - 1 + .1)
        ax.xaxis.set_tick_params(labelsize=8)
        ax.yaxis.set_tick_params(labelsize=8)
        ax.yaxis.set_major_locator(MaxNLocator(integer=True))
        ax.set_title(r'Organized by $RDR_{max}$', fontsize=10)
        ax.grid(True)
        ax.legend(loc='upper right', bbox_to_anchor=(1, -0.1), fontsize=8,
                  fancybox=True, shadow=True, ncol=1)

        # labels = ["{0:.0f}".format(y) for y in Floors]
        # ax.set_yticklabels(labels)

        self.ui.AlphaCurve.canvas.draw()
        self.ui.AlphaCurve.canvas.show()

        self.CreateNLM()
        self.ui.tabWidget.setCurrentIndex(4)

    def IDA(self):
        global Loc_span, Loc_heigth, ListNodes, Elements, DataBeamDesign, DataColDesign, Wtotal, num_elems, \
            ListNodesDrift, cIndex, ListNodesBasal, T1m, Wtotal, IMv, Sa_maxv, RDR_maxv, SDR_maxv, nrecs, dCap, \
            list_beams, list_cols
        exec(open("IDA.py").read())

        OutputIDAFile = self.ui.OutputIDAFile.text()
        IMv = np.asmatrix(IMv)
        Sa_maxv = np.asmatrix(Sa_maxv)
        RDR_maxv = np.asmatrix(RDR_maxv)
        SDR_maxv = np.asmatrix(SDR_maxv)

        np.savetxt('IDA/' + OutputIDAFile + '_dCap.txt', np.array([dCap]), fmt='%0.3f')
        with open('IDA/' + OutputIDAFile + '_indv.txt', 'w') as f:
            for line in IMv:
                np.savetxt(f, line, fmt='%.6f')

        with open('IDA/' + OutputIDAFile + '_IMv.txt', 'w') as f:
            for line in IMv:
                np.savetxt(f, line, fmt='%.6f')

        with open('IDA/' + OutputIDAFile + '_SDR_maxv.txt', 'w') as f:
            for line in SDR_maxv:
                np.savetxt(f, line, fmt='%.6f')

        with open('IDA/' + OutputIDAFile + '_Sa_maxv.txt', 'w') as f:
            for line in Sa_maxv:
                np.savetxt(f, line, fmt='%.6f')

        with open('IDA/' + OutputIDAFile + '_RDR_maxv.txt', 'w') as f:
            for line in RDR_maxv:
                np.savetxt(f, line, fmt='%.6f')

        IMv = np.asarray(IMv)
        Sa_maxv = np.asarray(Sa_maxv)
        SDR_maxv = np.asarray(SDR_maxv)
        RDR_maxv = np.asarray(RDR_maxv)
        self.ui.tabWidget.setCurrentIndex(5)

    def CSS(self):
        global Loc_span, Loc_heigth, ListNodes, Elements, DataBeamDesign, DataColDesign, Wtotal, num_elems, \
            ListNodesDrift, cIndex, ListNodesBasal, T1m, Wtotal, IM, Sa_max, RDR_max, SDR_max, nrecs, RA_max, EleCol, \
            EleBeam, DataColPhl, VnVu_max, list_beams, list_cols, DataBeamPhl, \
            maxSDRBdg, maxSDRBdgv, maxAccelBdgv, maxPhRot_Colcv, maxPhRot_Beamcv, MedPhRot_Colmv_v, MedPhRot_Beammv_v, \
            ResiBdg, ResiBdgv, maxCountColapv, SDRBdgVColapv
        if not os.path.exists("CSS"):
            os.mkdir("CSS")
        exec(open("CSS.py").read())
        op.wipeAnalysis()
        OutputCSSFile = self.ui.OutputCSSFile.text()
        np.savetxt('CSS/' + OutputCSSFile + '_IM.txt', IM, fmt='%.6f')
        np.savetxt('CSS/' + OutputCSSFile + '_Sa_max.txt', Sa_max, fmt='%.6f')
        np.savetxt('CSS/' + OutputCSSFile + '_RDR_max.txt', RDR_max, fmt='%.6f')
        np.savetxt('CSS/' + OutputCSSFile + '_SDR_max.txt', SDR_max, fmt='%.6f')
        np.savetxt('CSS/' + OutputCSSFile + '_RA_max.txt', RA_max, fmt='%.6f')
        np.savetxt('CSS/' + OutputCSSFile + '_VuVn_max.txt', VnVu_max, fmt='%.6f')
        # np.savetxt('CSS/' + OutputCSSFile + '_PhRot_Col_max.txt', maxPhRot_Colv, fmt='%.6f')
        # np.savetxt('CSS/' + OutputCSSFile + '_PhRot_Beam_max.txt', maxPhRot_Beamv, fmt='%.6f')
        np.savetxt('CSS/' + OutputCSSFile + '_PhRot1_Col_max.txt', maxPhRot_Colcv, fmt='%.6f')
        np.savetxt('CSS/' + OutputCSSFile + '_PhRot1_Beam_max.txt', maxPhRot_Beamcv, fmt='%.6f')
        np.savetxt('CSS/' + OutputCSSFile + '_SDR.txt', maxSDRBdgv, fmt='%.6f')
        np.savetxt('CSS/' + OutputCSSFile + '_Desp_Floor_resi.txt', ResiBdgv, fmt='%.6f')
        np.savetxt('CSS/' + OutputCSSFile + '_PFA.txt', maxAccelBdgv, fmt='%.6f')
        np.savetxt('CSS/' + OutputCSSFile + '_PhRot_Col_Med.txt', MedPhRot_Colmv_v, fmt='%.6f')
        np.savetxt('CSS/' + OutputCSSFile + '_PhRot_Beam_Med.txt', MedPhRot_Beammv_v, fmt='%.6f')
        np.savetxt('CSS/' + OutputCSSFile + '_Count_Colap_V.txt', maxCountColapv, fmt='%.6f')
        np.savetxt('CSS/' + OutputCSSFile + '_SDR_Colap_V.txt', SDRBdgVColapv, fmt='%.6f')


    def PlotIDA(self):
        exec(open("FragilityFunction.py").read())

    def PlotCSS(self):
        exec(open("FragilityFunctionCSS.py").read())

    def PlotCSSMaxFloors(self):
        exec(open("CSSMaxFloors.py").read())


# from ReadRecord import ReadRecord
# from getSaT import getSaT
# from runNRHA_CSS import runNRHA_CSS
# from Exci_pattern import Exci_pattern
#
#
# def compute_conccurent(EQname, sf, objeto):
#     print('Run: ' + EQname)
#     inFile = EQname + '.AT2'
#     outFile = EQname + '.G4'
#     dt, npts = ReadRecord(inFile, outFile)
#     dur = npts * dt
#     Sa, Sv, Sd, pga, amax = getSaT(outFile, dt, T1m, xi, npts)  # Get the PGA of the record in the X direction
#     IMgeomean = Sa
#     IM = np.append(IM, sf * IMgeomean)
#     objeto.CreateNLM()
#     Exci_pattern(dt, outFile, sf * g)
#     cIndex, mDT, mDPB, maxVB = runNRHA_CSS(dt, dur, Loc_heigth, ListNodesDrift, ListNodesBasal)
#     Sa_max = np.append(Sa_max, maxVB / Wtotal)
#     RDR_max = np.append(RDR_max, mDT)
#     SDR_max = np.append(SDR_max, mDPB)
#     op.wipe


if __name__ == "__main__":
    app = QApplication(sys.argv)
    w = MyForm()
    w.show()
    sys.exit(app.exec_())
