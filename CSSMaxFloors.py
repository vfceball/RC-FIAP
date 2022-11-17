import numpy as np  # load the numpy module, calling it np
import matplotlib.pyplot as plt

InputCSSFile = self.ui.InputCSSFile_2.text()
SDR_Floor_max = np.loadtxt('CSS/' + InputCSSFile + '_SDR_Floor_max.txt')
Accel_Floor_max = np.loadtxt('CSS/' + InputCSSFile + '_Accel_Floor_max.txt')
PhRot_Beam_max = np.loadtxt('CSS/' + InputCSSFile + '_PhRot_Beam_max.txt')
Floors = np.arange(np.size(SDR_Floor_max, 1)+1)
fig1 = self.ui.Hazard_CSS_SDR_max.canvas.axes
fig1.clear()
ax1 = fig1.add_subplot(111)
fig2 = self.ui.Hazard_CSS_Accel_max.canvas.axes
fig2.clear()
ax2 = fig2.add_subplot(111)
fig3 = self.ui.Hazard_CSS_Rot_max.canvas.axes
fig3.clear()
ax3 = fig3.add_subplot(111)

ax1.plot(SDR_Floor_max.T*100, Floors[1:], color='red', linewidth=0.5, marker='None')

# ax2.scatter(RA_max*100, IM, facecolors='none', edgecolors='b')
# ax3.scatter(VuVn_max, IM,  facecolors='none', edgecolors='b')
#
# ax1.set_ylabel(r'$S_a(T_1)$ [g]')
# ax1.set_xlabel(r'$SDR_{max}$ [%]')
# ax1.set_title(r'Maximum story drift ratio vs spectral acceleration $S_a$')
# ax1.set_xlim([0, 10])
# ax1.set_ylim([0, 5])
# ax1.grid(True)

# ax1.legend(loc='upper right', label=r'$i^{th} run$', fontsize=8)

self.ui.Hazard_CSS_SDR_max.canvas.draw()
self.ui.Hazard_CSS_SDR_max.canvas.show()
# ax2.set_ylabel(r'$S_a(T_1)$ [g]')
# ax2.set_xlabel(r'$Rotation_{max}$ [rad]')
# ax2.set_title(r'Maximun rotation angle vs spectral acceleration $S_a$')
# ax2.set_xlim([0, 10])
# ax2.set_ylim([0, 5])
# ax2.grid(True)
# self.ui.MaxRA_IM_CSS.canvas.draw()
# self.ui.MaxRA_IM_CSS.canvas.show()
# ax3.set_ylabel(r'$S_a(T_1)$ [g]')
# ax3.set_xlabel(r'$(Vu/Vn)_{max}$')
# ax3.set_title(r'Maximum normalized base shear vs spectral acceleration $S_a$')
# ax3.grid(True)
# ax3.set_xlim([0, 1.5])
# ax3.set_ylim(bottom=0)
# self.ui.Vu_Vn_IM_CSS.canvas.draw()
# self.ui.Vu_Vn_IM_CSS.canvas.show()
#
