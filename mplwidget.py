# -*- coding: utf-8 -*-
"""
Created on Fri Oct  2 06:46:48 2020

@author: Victor Ceballos
"""

# ------------------------------------------------- -----
# -------------------- mplwidget.py --------------------
# -------------------------------------------------- ----
from PyQt5.QtWidgets import *

from PyQt5 import QtCore

from matplotlib.backends.backend_qt5agg import FigureCanvas

from matplotlib.figure import Figure

from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar

class MplWidget(QWidget):

    def __init__(self, parent=None):
        QWidget.__init__(self, parent)
        self.fig = Figure(tight_layout=True)
        self.canvas = FigureCanvas(self.fig)
        vertical_layout = QVBoxLayout()
        vertical_layout.addWidget(self.canvas)
        # self.canvas.axes = self.canvas.figure.add_subplot(111)
        self.canvas.axes = self.canvas.figure #.add_axes([0, 0, 1, 1])


        self.canvas.fig = FigureCanvas(Figure())
        # self.canvas.fig = FigureCanvas(self.fig)
        self.setLayout(vertical_layout)
        self.toolbar = NavigationToolbar(self.canvas, self, coordinates=False)
        self.toolbar.setMinimumWidth(300)
        self.toolbar.setStyleSheet("QToolBar { border: 0px }")
        self.toolbar.setIconSize(QtCore.QSize(15, 15))


