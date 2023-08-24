global nrecs, Sa_T1, SDR_max, HL_v
import numpy as np  # load the numpy module, calling it np
# import statsmodels.api as sm
import matplotlib.pyplot as plt
from scipy.stats import norm, binom, lognorm
from scipy.optimize import fmin


# def fn_mle_pc_probit(IM, num_gms, num_collapse):
#     import statsmodels.api as sm
#
#     Y = np.array([num_collapse, num_gms])
#     sm_probit_Link = sm.genmod.families.links.probit
#     glm_binom = sm.GLM(sm.add_constant(np.log(IM)), sm.add_constant(Y), family=sm.families.Binomial(link=sm_probit_Link))
#     glm_result = glm_binom.fit()
#     b = glm_result.params
#     theta = np.exp(-b[1] / b[2])
#     beta = 1 / b[2]
#     return theta, beta

# def fn_mle_pc(IMo, num_gms, num_collapse):
#     from scipy.optimize import fmin
#     from scipy.stats import norm, binom
#
#     def mlefit(params, num_gms, num_collapse, IMo):
#         if params[1] < 0:
#             loglik = 1e10
#         else:
#             p = norm.cdf(np.log(IMo), params[0], params[1])
#             likelihood = binom.pmf(num_collapse, num_gms, p)
#             print('likelihood', likelihood)
#             likelihood = np.where(likelihood == 0, np.finfo(float).tiny, likelihood)
#             print('likelihood', likelihood)
#             loglik = -np.sum(np.log(likelihood))
#         return loglik
#     x0 = [np.mean(np.log(IMo)), np.std(np.log(IMo))]
#     print('x0', x0)
#     # options = optimset('MaxFunEvals',1000, 'GradObj', 'off')
#     x = fmin(mlefit, x0, args=(num_gms, num_collapse, IMo))
#     print('x=', x)
#     theta = np.exp(x[0])
#     beta = x[1]
#     return theta, beta


def fn_sse_pc(IMo, num_gms, num_collapse):
    from scipy.optimize import fmin

    def ssefit(x, num_gms, num_collapse, IMo):
        if x[0] < 0:
            x[0] = 0
        p = norm.cdf(np.log(IMo), np.log(x[0]), x[1])
        sse = np.sum((p - num_collapse / num_gms) ** 2)
        return sse

    x0 = np.array([0.8, 0.4])
    x = fmin(ssefit, x0, args=(num_gms, num_collapse, IMo))
    theta = x[0]
    beta = x[1]
    return theta, beta


# colors = np.array(["green", "orange", "red", "pink", "black", "yellow", "purple", "beige", "brown", "gray", "cyan",
#                    "magenta"])
colors = plt.cm.jet(np.linspace(0, 1, 11))

InputCSSFile = self.ui.InputCSSFile.text()
CSS_SDRCurves = self.ui.CSS_LimitStage.text()
CSS_SDRCurves = CSS_SDRCurves.split(',')
CSS_SDRCurves = np.array(CSS_SDRCurves, dtype=float) / 100
t_min_lim = float(self.ui.t_min.text())
SDR_lim = float(self.ui.SDR_limit.text())
HL_Plot = self.ui.HL_Plot.text()
HL_Plot = HL_Plot.split(',')
HL_Plot = np.array(HL_Plot, dtype=int)
dCap = CSS_SDRCurves[-1]
Test = np.loadtxt('CSS/' + InputCSSFile + '_Test_Run_Earthquake.txt', dtype=np.str, delimiter=',', skiprows=0)
IM = np.loadtxt('CSS/' + InputCSSFile + '_IM.txt')
SDR_max = np.loadtxt('CSS/' + InputCSSFile + '_SDR_max.txt')
Sa_max = np.loadtxt('CSS/' + InputCSSFile + '_Sa_max.txt')
RDR_max = np.loadtxt('CSS/' + InputCSSFile + '_RDR_max.txt')
# RA_max = np.loadtxt('CSS/' + InputCSSFile + '_RA_max.txt')
VuVn_max = np.loadtxt('CSS/' + InputCSSFile + '_VuVn_max.txt')
ncolors = Test[:, 0]
t_run, t_earth, maxDriftPisoBdg = Test[:, 2], Test[:, 3], Test[:, 4]
t_run, t_earth, maxDriftPisoBdg = t_run.astype('float'), t_earth.astype('float'), maxDriftPisoBdg.astype('float')
# print('t_earth', t_earth)
# print('t_run', t_run)
t_min = t_run / t_earth
# print('t_min', t_min)
ncolors = ncolors.astype('int')
ind = np.where(ncolors == np.expand_dims(HL_Plot, axis=1))[1]
# print('ind', ind)
ncolors = ncolors[ind]
IM = IM[ind]
SDR_max = SDR_max[ind]
Sa_max = Sa_max[ind]
RDR_max = RDR_max[ind]
# RA_max = RA_max[ind]
VuVn_max = VuVn_max[ind]
t_min, maxDriftPisoBdg = t_min[ind], maxDriftPisoBdg[ind]
# print('SDR_lim', SDR_lim,'t_min_lim', t_min_lim )
ind = np.where(((t_min < t_min_lim / 100) & (maxDriftPisoBdg >= SDR_lim)) | (t_min >= t_min_lim / 100))
# print('ind', ind)
ncolors = ncolors[ind]
IM = IM[ind]
SDR_max = SDR_max[ind]
Sa_max = Sa_max[ind]
RDR_max = RDR_max[ind]
# RA_max = RA_max[ind]
VuVn_max = VuVn_max[ind]
t_min, maxDriftPisoBdg = t_min[ind], maxDriftPisoBdg[ind]
nrecs = IM.size
fig1 = self.ui.SDR_IM_CSS.canvas.axes
fig1.clear()
ax1 = fig1.add_subplot(111)
fig2 = self.ui.MaxRA_IM_CSS.canvas.axes
fig2.clear()
ax2 = fig2.add_subplot(111)
fig3 = self.ui.Vu_Vn_IM_CSS.canvas.axes
fig3.clear()
ax3 = fig3.add_subplot(111)

nSDR = CSS_SDRCurves.size
# if np.max(SDR_max) >= dCap:
#     ind = np.where(SDR_max >= dCap)[0][0]
#     SDR_max = SDR_max[:ind + 1]
#     RDR_max = RDR_max[:ind + 1]
#     Sa_max = Sa_max[:ind + 1]
#     IM = IM[:ind + 1]

for HL in HL_Plot:
    ind = np.where(ncolors == HL)
    ax1.scatter(SDR_max[ind] * 100, IM[ind], facecolors='none', edgecolors=colors[HL], label=r'$%1.0f$' % HL)
    ax2.scatter(Sa_max[ind], IM[ind], facecolors='none', edgecolors=colors[HL], label=r'$%1.0f$' % HL)
    ax3.scatter(VuVn_max[ind], IM[ind], facecolors='none', edgecolors=colors[HL], label=r'$%1.0f$' % HL)

ax1.set_ylabel(r'$S_a(T_1)$ [g]')
ax1.set_xlabel(r'$SDR_{max}$ [%]')
ax1.set_title(r'Maximum story drift ratio vs spectral acceleration $S_a$')
ax1.set_xlim([0, 10])
ax1.set_ylim([0, 5])
ax1.grid(True)
# ax1.legend(loc='lower right', fontsize=8, fancybox=True, shadow=True, title="Hazard levels")

# plt.legend(handles=scatter.legend_elements()[0], labels=classes)
# kw = dict(prop="edgecolors", color=colors[ncolors], fmt="$ {x:.2f}", func=lambda s: np.sqrt(s/.3)/3)
# legend2 = ax.legend(*scatter.legend_elements(**kw),
#                     loc="lower right", title="Price")
# handles, labels = scatter1.legend_elements(prop="colors")
# legend1 = ax1.legend(handles, labels, loc="upper left", title="Hazard levels")

# handles, labels = scatter1.legend_elements()
# legend1 = ax.legend(handles, labels, loc="lower left", title="Colors")
# ax.add_artist(legend1)
ax1.legend(loc="upper left", title="Hazard levels", title_fontsize='small', fancybox=True, shadow=True, ncol=1,
           fontsize=6)
ax2.legend(loc="upper right", title="Hazard levels", title_fontsize='small', fancybox=True, shadow=True, ncol=1,
           fontsize=6)
ax3.legend(loc="upper left", title="Hazard levels", title_fontsize='small', fancybox=True, shadow=True, ncol=1,
           fontsize=6)
self.ui.SDR_IM_CSS.canvas.draw()
self.ui.SDR_IM_CSS.canvas.show()
# ax2.set_ylabel(r'$S_a(T_1)$ [g]')
# ax2.set_xlabel(r'$Rotation_{max}$ [rad]')
# ax2.set_title(r'Maximun rotation angle vs spectral acceleration $S_a$')
# ax2.set_xlim([0, 10])
# ax2.set_ylim([0, 5])
# ax2.grid(True)
ax2.set_ylabel(r'$S_a(T_1)$ [g]')
ax2.set_xlabel(r'$Vbasal_{max}$ [g]')
ax2.set_title(r'Maximun basal shear vs spectral acceleration $S_a$')
ax2.set_xlim([0, 2])
ax2.set_ylim([0, 5])
ax2.grid(True)

self.ui.MaxRA_IM_CSS.canvas.draw()
self.ui.MaxRA_IM_CSS.canvas.show()
ax3.set_ylabel(r'$S_a(T_1)$ [g]')
ax3.set_xlabel(r'$(Vu/Vn)_{max}$')
ax3.set_title(r'Maximum normalized base shear vs spectral acceleration $S_a$')
ax3.grid(True)
ax3.set_xlim([0, 1.5])
ax3.set_ylim(bottom=0)

self.ui.Vu_Vn_IM_CSS.canvas.draw()
self.ui.Vu_Vn_IM_CSS.canvas.show()

fig4 = self.ui.FragilityCurveCSS.canvas.axes
fig4.clear()
ax4 = fig4.add_subplot(111)

# ind = np.where(Sa_max <= 5)[0]
# IM = IM[ind]
hist_tot, bin_edges_tot = np.histogram(IM, bins=np.geomspace(0.0001, 5, 25))
# print('hist_tot=', hist_tot)
# print('bin_edges_tot=', bin_edges_tot)
# print('hist_tot=', len(hist_tot))
# print('bin_edges_tot=', len(bin_edges_tot))
bin_avg_tot = (bin_edges_tot[0:-1] + bin_edges_tot[1:]) / 2
ind_t = np.where(hist_tot > 0)[0]
hist_tot = hist_tot[ind_t]
bin_avg_tot = bin_avg_tot[ind_t]

# print('hist_tot=', hist_tot)
# print('bin_edges_tot=', bin_edges_tot)
markers = ['^', 's', 'p', 'h', '8']

FragCurve = self.ui.comboBoxFragCurve.currentText()
if FragCurve == 'Max SDR':
    FRAG_max = SDR_max
elif FragCurve == 'Max VBasal':
    FRAG_max = Sa_max
elif FragCurve == 'Max Vu/Vn':
    FRAG_max = VuVn_max

colors = np.array(["green", "orange", "red", "brown", "magenta", "cyan", "yellow", "purple", "beige", "black", "gray",
                   "pink"])

for jj in range(nSDR):
    ind = np.where(FRAG_max >= CSS_SDRCurves[jj])[0]
    # IM_bin = np.sort(IM[ind])
    hist, bin_edges = np.histogram(IM[ind], bins=np.geomspace(0.0001, 5, 25))
    hist = hist[ind_t]
    bin_avg = (bin_edges[0:-1] + bin_edges[1:]) / 2
    bin_avg = bin_avg[ind_t]
    # print('bin_avg', bin_avg)
    # theta_hat_probit, beta_hat_probit = fn_mle_pc_probit(bin_avg, num_gms, hist)
    # theta_hat_mle, beta_hat_mle = fn_mle_pc(bin_avg, num_gms, hist)
    theta_hat_sse, beta_hat_sse = fn_sse_pc(bin_avg, hist_tot, hist)
    # x_vals = np.arange(0.01, np.ceil(np.max(IM_bin)), 0.001)
    # p_collapse_mle = norm.cdf((np.log(x_vals/theta_hat_sse))/beta_hat_sse)
    p_collapse_points = hist / hist_tot
    # ax4.plot(bin_avg, p_collapse_points, marker=markers[jj], label=r'$j = %1.1f \theta = %1.3 \beta = %1.3$' %(CSS_SDRCurves[jj] * 100, theta_hat_sse, beta_hat_sse))
    ax4.scatter(bin_avg, p_collapse_points, marker=markers[jj], s=8, c=colors[jj])
    p_collapse_points = np.sort(p_collapse_points)
    # print('hist =', hist)
    # print('hist_tot=', hist_tot)
    # print('p_collapse_points', p_collapse_points)
    # p_collapse_hist = p_collapse_points[1:] - p_collapse_points[0:-1]
    # print('p_collapse_hist', p_collapse_hist)

    # theta_hat_mom = np.exp(np.mean(np.log(p_collapse_hist)))
    # beta_hat_mom = np.std(np.log(p_collapse_hist))
    x_vals = np.geomspace(0.0001, 5, 1000)
    p_collapse = norm.cdf(np.log(x_vals / theta_hat_sse) / beta_hat_sse)
    ax4.plot(x_vals, p_collapse, c=colors[jj], label=r'$j = %1.1f \/\/\theta = %1.2f \/\/\beta = %1.2f$'
                                                     % (CSS_SDRCurves[jj] * 100, theta_hat_sse, beta_hat_sse))
    # shape, loc, scale = lognorm.fit(np.exp(p_collapse_hist))
    # x = np.linspace(0, 5, 100)
    # cdf = lognorm.cdf(x, shape, loc, scale)
    # ax4.plot(x, cdf)
    # ax4.plot(x_vals, p_collapse_mle, 'b-')

if FragCurve == 'Max SDR':
    ax4.set_ylabel(r'$P(SDR_{max}>j|S_a=y)$')
elif FragCurve == 'Max Rotation':
    ax4.set_ylabel(r'$P(Rotation_{max}>j|S_a=y)$')
elif FragCurve == 'Max Vu/Vn':
    ax4.set_ylabel(r'$P((Vu/Vn)_{max}>j|S_a=y)$')
ax4.set_xlabel(r'$S_a(T_1)$ [g]')
ax4.set_xlim([0, 1])
ax4.set_ylim(bottom=0)
ax4.legend(loc='lower right', fontsize=8, fancybox=True, shadow=True, ncol=1)
ax4.grid(True)
self.ui.FragilityCurveCSS.canvas.draw()
self.ui.FragilityCurveCSS.canvas.show()
