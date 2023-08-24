import numpy as np  # load the numpy module, calling it np
import matplotlib.pyplot as plt
# colors = np.array(["green", "green", "orange", "red", "magenta", "blue", "gold", "purple", "cyan", "brown",
#                    "lime", "gray", "olive"])
colors = plt.cm.jet(np.linspace(0, 1, 10))
InputCSSFile = self.ui.InputCSSFile_2.text()
t_min_lim = float(self.ui.t_min_2.text())
SDR_lim = float(self.ui.SDR_limit_2.text())
HL_Plot = self.ui.HL_Plot_2.text()
HL_Plot = HL_Plot.split(',')
HL_Plot = np.array(HL_Plot, dtype=int)-1
Test = np.loadtxt('CSS/' + InputCSSFile + '_Test_Run_Earthquake.txt', dtype=np.str, delimiter=',', skiprows=0)
SDR_Floor_max = np.loadtxt('CSS/' + InputCSSFile + '_SDR.txt')
Desp_Floor_resi = np.loadtxt('CSS/' + InputCSSFile + '_Desp_Floor_resi.txt')
Accel_Floor_max = np.loadtxt('CSS/' + InputCSSFile + '_PFA.txt')
PhRot_Beam_max = np.loadtxt('CSS/' + InputCSSFile + '_PhRot1_Beam_max.txt')
RDR_max = np.loadtxt('CSS/' + InputCSSFile + '_RDR_max.txt')

ncolors = Test[:, 0]
t_run, t_earth, maxDriftPisoBdg = Test[:, 2], Test[:, 3], Test[:, 4]
t_run, t_earth, maxDriftPisoBdg = t_run.astype('float'), t_earth.astype('float'), maxDriftPisoBdg.astype('float')
t_min = t_run / t_earth
ncolors = ncolors.astype('int')-1

ind = np.where(ncolors == np.expand_dims(HL_Plot, axis=1))[1]
ncolors = ncolors[ind]
Floors = np.arange(np.size(SDR_Floor_max, 1)+1)
t_min, maxDriftPisoBdg = t_min[ind], maxDriftPisoBdg[ind]
SDR_Floor_max = SDR_Floor_max[ind, :]
Desp_Floor_resi = Desp_Floor_resi[ind, :]
Accel_Floor_max = Accel_Floor_max[ind, :]
PhRot_Beam_max = PhRot_Beam_max[ind, :]
RDR_max = RDR_max[ind]
Alpha_Floor_max = (SDR_Floor_max.T/RDR_max).T

ind = np.where(((t_min < t_min_lim / 100) & (maxDriftPisoBdg >= SDR_lim)) | (t_min >= t_min_lim / 100))
# ind = ind.tolist()
# ind = list(ind)
ind = np.asarray(ind)
ind = ind.flatten()
ncolors = ncolors[ind]
t_min, maxDriftPisoBdg = t_min[ind], maxDriftPisoBdg[ind]
SDR_Floor_max = SDR_Floor_max[ind, :]
Desp_Floor_resi = Desp_Floor_resi[ind, :]
Accel_Floor_max = Accel_Floor_max[ind, :]
PhRot_Beam_max = PhRot_Beam_max[ind, :]
Alpha_Floor_max = Alpha_Floor_max[ind, :]

fig1 = self.ui.Hazard_CSS_SDR_max.canvas.axes
fig1.clear()
ax1 = fig1.add_subplot(111)
fig2 = self.ui.Hazard_CSS_Alpha_max.canvas.axes
fig2.clear()
ax2 = fig2.add_subplot(111)
fig3 = self.ui.Hazard_CSS_Accel_max.canvas.axes
fig3.clear()
ax3 = fig3.add_subplot(111)
fig4 = self.ui.Hazard_CSS_Resid_despl.canvas.axes
fig4.clear()
ax4 = fig4.add_subplot(111)

maxlimit_SDR = 0
maxlimit_Alpha = 0
maxlimit_Accel = 0
maxlimit_Resid = 0

for HL in HL_Plot:
    ind = np.where(ncolors == HL)
    ind = np.asarray(ind)
    ind = ind.flatten()
    Floor_m = Floors[1:]*np.ones((len(ind), len(Floors[1:])))
    ax1.plot(SDR_Floor_max[ind, :].T * 100, Floor_m.T, c=colors[HL], linewidth=0.5, marker='None', alpha=0.25)
    ax1.plot([], [], c=colors[HL], linewidth=0.5, marker='None', alpha=0.25)
    ax1.plot(np.mean(SDR_Floor_max[ind, :], axis=0)*100, Floors[1:], c=colors[HL], linewidth=3,
             label=r'$%1.0f$' % (HL+1), marker='o', ms=5)
    ax2.plot(Alpha_Floor_max[ind, :].T, Floor_m.T, c=colors[HL], linewidth=0.5, marker='None', alpha=0.25)
    ax2.plot([], [], c=colors[HL], linewidth=0.5, marker='None', alpha=0.25)
    ax2.plot(np.mean(Alpha_Floor_max[ind, :], axis=0), Floors[1:], c=colors[HL], linewidth=3, label=r'$%1.0f$' % (HL+1),
             marker='o', ms=5)
    ax3.plot(Accel_Floor_max[ind, :].T, Floor_m.T, c=colors[HL], linewidth=0.5, marker='None', alpha=0.25)
    ax3.plot([], [], c=colors[HL], linewidth=0.5, marker='None', alpha=0.25)
    ax3.plot(np.mean(Accel_Floor_max[ind, :], axis=0), Floors[1:], c=colors[HL], linewidth=3, label=r'$%1.0f$' % (HL+1),
             marker='o', ms=5)
    # print('Desp_Floor_resi', Desp_Floor_resi)
    Desp_Floor_resi_plot = Desp_Floor_resi[ind, :]
    ind2 = np.isfinite(Desp_Floor_resi_plot[:, 0])
    # print('Desp_Floor_resi', Desp_Floor_resi)

    # print(Floor_m.T)
    # print('ind2', ind2)
    # print(Desp_Floor_resi_plot[ind2, :].T)
    ax4.plot(Desp_Floor_resi_plot[ind2, :].T, Floor_m[ind2].T, c=colors[HL], linewidth=0.5, marker='None', alpha=0.25)
    ax4.plot([], [], c=colors[HL], linewidth=0.5, marker='None', alpha=0.25)
    ax4.plot(np.mean(Desp_Floor_resi_plot[ind2, :], axis=0), Floors[1:], c=colors[HL], linewidth=3, label=r'$%1.0f$' % (HL+1),
             marker='o', ms=5)
    xlimit_SDR = np.ceil(np.max(np.mean(SDR_Floor_max[ind, :], axis=0)*100)*10)/10
    xlimit_Alpha = np.ceil(np.max(np.mean(Alpha_Floor_max[ind, :], axis=0))*10)/10
    xlimit_Accel = np.ceil(np.max(np.mean(Accel_Floor_max[ind, :], axis=0))*10)/10
    xlimit_Resid = np.ceil(np.max(np.mean(Desp_Floor_resi_plot[ind2, :], axis=0))*10)/10

    if maxlimit_SDR < xlimit_SDR:
        maxlimit_SDR = xlimit_SDR
    if maxlimit_Alpha < xlimit_Alpha:
        maxlimit_Alpha = xlimit_Alpha
    if maxlimit_Accel < xlimit_Accel:
        maxlimit_Accel = xlimit_Accel
    if maxlimit_Resid < xlimit_Resid:
        maxlimit_Resid = xlimit_Resid

# ax1.plot(np.percentile(SDR_Floor_max, 84, axis=0)*100, Floors[1:], c='k', linestyle='dashed', linewidth=2,
#          label='84th pct')

#     ax1.scatter(SDR_max[ind] * 100, IM[ind], facecolors='none', edgecolors=colors[HL], label=r'$%1.0f$' % HL)
#     ax2.scatter(RA_max[ind] * 100, IM[ind], facecolors='none', edgecolors=colors[HL], label=r'$%1.0f$' % HL)
#     ax3.scatter(VuVn_max[ind], IM[ind], facecolors='none', edgecolors=colors[HL], label=r'$%1.0f$' % HL)
#
#
# ax1.plot(SDR_Floor_max.T*100, Floors[1:], color='red', linewidth=0.5, marker='None')

# ax2.scatter(RA_max*100, IM, facecolors='none', edgecolors='b')
# ax3.scatter(VuVn_max, IM,  facecolors='none', edgecolors='b')
#
ax1.set_ylabel(r'Floors')
ax1.set_xlabel(r'$SDR_{max}$ [%]')
ax1.set_title(r'Maximum story drift')
ax1.set_xlim([0, maxlimit_SDR])
ax1.set_ylim(1 - .1, len(Floors) - 1 + .1)
ax1.yaxis.set_major_locator(MaxNLocator(integer=True))
ax1.legend(loc="upper right", title="HL", title_fontsize='small', fancybox=True, shadow=True, ncol=1,
           fontsize=6)
# ax1.set_ylim([0, 5])
ax1.grid(True)

ax2.set_ylabel(r'Floors')
ax2.set_xlabel(r'$Alpha = SDR/RDR$')
#
# ax2.set_xlabel(r'$SDR_{max}$ [%]')
ax2.set_title(r'Maximum alpha ratio')
ax2.set_xlim([0, maxlimit_Alpha])
ax2.set_ylim(1 - .1, len(Floors) - 1 + .1)
ax2.yaxis.set_major_locator(MaxNLocator(integer=True))
ax2.legend(loc="upper right", title="HL", title_fontsize='small', fancybox=True, shadow=True, ncol=1, fontsize=6)
ax2.grid(True)

ax3.set_ylabel(r'Floors')
ax3.set_xlabel(r'$Accel_{max}$ [g]')
ax3.set_title(r'Maximun story acceleration')
ax3.set_xlim([0, maxlimit_Accel])
ax3.set_ylim(1 - .1, len(Floors) - 1 + .1)
ax3.yaxis.set_major_locator(MaxNLocator(integer=True))
ax3.legend(loc="upper right", title="HL", title_fontsize='small', fancybox=True, shadow=True, ncol=1, fontsize=6)
ax3.grid(True)

ax4.set_ylabel(r'Floors')
ax4.set_xlabel(r'$\Delta_{residual}$ [cm]')
ax4.set_title(r'Residual displacement')
ax4.set_xlim([0, maxlimit_Resid])
ax4.set_ylim(1 - .1, len(Floors) - 1 + .1)
ax4.yaxis.set_major_locator(MaxNLocator(integer=True))
ax4.legend(loc="upper right", title="HL", title_fontsize='small', fancybox=True, shadow=True, ncol=1, fontsize=6)
ax4.grid(True)

# ax1.legend(loc='upper right', label=r'$i^{th} run$', fontsize=8)

self.ui.Hazard_CSS_SDR_max.canvas.draw()
self.ui.Hazard_CSS_SDR_max.canvas.show()
self.ui.Hazard_CSS_Alpha_max.canvas.draw()
self.ui.Hazard_CSS_Alpha_max.canvas.show()
self.ui.Hazard_CSS_Accel_max.canvas.draw()
self.ui.Hazard_CSS_Accel_max.canvas.show()
self.ui.Hazard_CSS_Resid_despl.canvas.draw()
self.ui.Hazard_CSS_Resid_despl.canvas.show()

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
