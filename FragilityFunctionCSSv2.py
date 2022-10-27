global nrecs, Sa_T1, SDR_max
import numpy as np  # load the numpy module, calling it np
import statsmodels.api as sm
import matplotlib.pyplot as plt
from scipy.stats import norm, binom
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
        sse = np.sum((p - num_collapse/num_gms)**2)
        return sse

    x0 = np.array([0.8, 0.4])
    x = fmin(ssefit, x0, args=(num_gms, num_collapse, IMo))
    theta = x[0]
    beta = x[1]
    return theta, beta



InputCSSFile = self.ui.InputCSSFile.text()
CSS_SDRCurves = self.ui.CSS_SDRCurves.text()
CSS_SDRCurves = CSS_SDRCurves.split(',')
CSS_SDRCurves = np.array(CSS_SDRCurves, dtype=float)/100
dCap = CSS_SDRCurves[-1]
IM = np.loadtxt('CSS/' + InputCSSFile + '_IM.txt')
SDR_max = np.loadtxt('CSS/' + InputCSSFile + '_SDR_max.txt')
Sa_max = np.loadtxt('CSS/' + InputCSSFile + '_Sa_max.txt')
RDR_max = np.loadtxt('CSS/' + InputCSSFile + '_RDR_max.txt')

nrecs = IM.size

print('nrecs =', nrecs)

fig1 = self.ui.SDR_IM_CSS.canvas.axes
fig1.clear()
ax1 = fig1.add_subplot(111)
fig2 = self.ui.RDR_IM_CSS.canvas.axes
fig2.clear()
ax2 = fig2.add_subplot(111)
fig3 = self.ui.Vb_W_IM_CSS.canvas.axes
fig3.clear()
ax3 = fig3.add_subplot(111)

nSDR = CSS_SDRCurves.size
# if np.max(SDR_max) >= dCap:
#     ind = np.where(SDR_max >= dCap)[0][0]
#     SDR_max = SDR_max[:ind + 1]
#     RDR_max = RDR_max[:ind + 1]
#     Sa_max = Sa_max[:ind + 1]
#     IM = IM[:ind + 1]



ax1.plot(SDR_max*100, IM, 'bo')
ax2.plot(RDR_max*100, IM, 'bo')
ax3.plot(IM[:-1], Sa_max[:-1], 'bo')

ax1.set_ylabel(r'$S_a(T_1)$ [g]')
ax1.set_xlabel(r'$SDR_{max}$ [%]')
ax1.set_title(r'Maximum story drift ratio vs spectral acceleration $S_a$')
ax1.set_xlim([0, 10])
ax1.set_ylim([0, 5])
ax1.grid(True)
self.ui.SDR_IM_CSS.canvas.draw()
self.ui.SDR_IM_CSS.canvas.show()
ax2.set_ylabel(r'$S_a(T_1)$ [g]')
ax2.set_xlabel(r'$RDR_{max}$ [%]')
ax2.set_title(r'Maximun roof drift ratio vs spectral acceleration $S_a$')
ax2.set_xlim([0, 10])
ax2.set_ylim([0, 5])
ax2.grid(True)
self.ui.RDR_IM_CSS.canvas.draw()
self.ui.RDR_IM_CSS.canvas.show()
ax3.set_xlabel(r'$S_a(T_1)$ [g]')
ax3.set_ylabel(r'$Vb_{max}/W$')
ax3.set_title(r'Maximum normalized base shear vs spectral acceleration $S_a$')
ax3.grid(True)
ax3.set_xlim([0, 5])
ax3.set_ylim(bottom=0)
self.ui.Vb_W_IM_CSS.canvas.draw()
self.ui.Vb_W_IM_CSS.canvas.show()

# fig4 = self.ui.FragilityCurveCSS.canvas.axes
# fig4.clear()
# ax4 = fig4.add_subplot(111)

ind = np.where(Sa_max <= 5)[0]
IM = IM[ind]
SDR_max = SDR_max[ind]
ind = np.where(SDR_max <= 0.1)[0]
IM = IM[ind]
SDR_max = SDR_max[ind]

# for jj in range(nSDR):
#     ind = np.where((SDR_max >= CSS_SDRCurves[jj]) & (SDR_max <= 0.1))[0]
    # IM_bin = np.sort(IM[ind])
    # hist, bin_edges = np.histogram(IM_bin, bins=np.linspace(0, 5, 26))

H, xedges, yedges = np.histogram2d(IM, SDR_max, bins=(np.linspace(0, 5, 26), np.linspace(0, 0.1, 21)))
print('Histogram2d =', H)
print('SumH =', np.sum(H))
binx_avg = (xedges[0:-1] + xedges[1:]) / 2
biny_avg = (yedges[0:-1] + yedges[1:]) / 2
print('binx_avg', binx_avg)
print('biny_avg', biny_avg)
num_gms = H.sum(axis=0)
num_gms = np.full(H.shape, num_gms)
print('num_gms', num_gms)
print('prob acumulada', np.cumsum(H, axis=0) / num_gms)
    # theta_hat_probit, beta_hat_probit = fn_mle_pc_probit(bin_avg, num_gms, hist)
    # theta_hat_mle, beta_hat_mle = fn_mle_pc(bin_avg, num_gms, hist)
    # theta_hat_sse, beta_hat_sse = fn_sse_pc(bin_avg, num_gms, hist)
    # x_vals = np.arange(0.01, np.ceil(np.max(IM_bin)), 0.001)
    # p_collapse_mle = norm.cdf((np.log(x_vals/theta_hat_sse))/beta_hat_sse)
#     ax4.plot(bin_avg, np.cumsum(hist) / num_gms, 'b.')
#     # ax4.plot(x_vals, p_collapse_mle, 'b-')
# ax4.set_ylabel(r'$P(SDR_{max}>j|S_a=y)$')
# ax4.set_xlabel(r'$S_a(T_1)$ [g]')
# ax4.set_xlim([0, 5])
# ax4.set_ylim(bottom=0)
# ax4.grid(True)
# self.ui.FragilityCurveCSS.canvas.draw()
# self.ui.FragilityCurveCSS.canvas.show()
