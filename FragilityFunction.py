global nrecs, Sa_T1v, SDR_maxv
import numpy as np  # load the numpy module, calling it np
import matplotlib.pyplot as plt
from scipy.stats import norm

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

InputIDAFile = self.ui.InputIDAFile.text()
SDRCurves = self.ui.SDRCurves.text()
SDRCurves = SDRCurves.split(',')
SDRCurves = np.array(SDRCurves, dtype=float)/100

dCap = np.loadtxt('IDA/' + InputIDAFile + '_dCap.txt')
IMv = np.loadtxt('IDA/' + InputIDAFile + '_IMv.txt')
SDR_maxv = np.loadtxt('IDA/' + InputIDAFile + '_SDR_maxv.txt')
Sa_maxv = np.loadtxt('IDA/' + InputIDAFile + '_Sa_maxv.txt')
RDR_maxv = np.loadtxt('IDA/' + InputIDAFile + '_RDR_maxv.txt')
# print('IMv =', IMv)
nrecs = np.array(IMv, ndmin=2).shape[0]

# print('nrecs =', nrecs)

fig1 = self.ui.SDR_Sa.canvas.axes
fig1.clear()
ax1 = fig1.add_subplot(111)
fig2 = self.ui.RDR_Sa.canvas.axes
fig2.clear()
ax2 = fig2.add_subplot(111)
fig3 = self.ui.Vb_W_Sa.canvas.axes
fig3.clear()
ax3 = fig3.add_subplot(111)


nSDR = SDRCurves.size
IM_cap = np.zeros((nSDR, nrecs))
for ii in range(nrecs):
    if nrecs == 1:
        SDR_max = SDR_maxv
        RDR_max = RDR_maxv
        Sa_max = Sa_maxv
        IM = IMv
    else:
        SDR_max = SDR_maxv[ii, :]
        RDR_max = RDR_maxv[ii, :]
        Sa_max = Sa_maxv[ii, :]
        IM = IMv[ii, :]
    if np.max(SDR_max) >= dCap:
        ind = np.where(SDR_max >= dCap)[0][0]
        SDR_max = SDR_max[:ind + 1]
        RDR_max = RDR_max[:ind + 1]
        Sa_max = Sa_max[:ind + 1]
        IM = IM[:ind + 1]

    ax1.plot(SDR_max*100, IM, 'b.-')
    ax2.plot(RDR_max*100, IM, 'b.-')
    ax3.plot(IM[:-1], Sa_max[:-1], 'b.-')

    for jj in range(nSDR):
        if nrecs == 1:
            IM_cap[jj, ii] = np.interp(SDRCurves[jj], SDR_max, IM)
        else:
            IM_cap[jj, ii] = np.interp(SDRCurves[jj], SDR_max, IM)
ax1.set_ylabel(r'$S_a(T_1)$ [g]')
ax1.set_xlabel(r'$SDR_{max}$ [%]')
ax1.set_title(r'Maximum story drift ratio vs spectral acceleration $S_a$')
ax1.set_xlim([0, SDRCurves[-1]*100])
ax1.set_ylim(bottom=0)
ax1.grid(True)
self.ui.SDR_Sa.canvas.draw()
self.ui.SDR_Sa.canvas.show()
ax2.set_ylabel(r'$S_a(T_1)$ [g]')
ax2.set_xlabel(r'$RDR_{max}$ [%]')
ax2.set_title(r'Maximun roof drift ratio vs spectral acceleration $S_a$')
ax2.set_xlim([0, SDRCurves[-1]*100])
ax2.set_ylim(bottom=0)
ax2.grid(True)
self.ui.RDR_Sa.canvas.draw()
self.ui.RDR_Sa.canvas.show()
ax3.set_xlabel(r'$S_a(T_1)$ [g]')
ax3.set_ylabel(r'$Vb_{max}/W$')
ax3.set_title(r'Maximum normalized base shear vs spectral acceleration $S_a$')
ax3.grid(True)
ax3.set_xlim(left=0)
ax3.set_ylim(bottom=0)
self.ui.Vb_W_Sa.canvas.draw()
self.ui.Vb_W_Sa.canvas.show()

IM_cap = np.sort(IM_cap)
fig4 = self.ui.FragilityCurve.canvas.axes
fig4.clear()
ax4 = fig4.add_subplot(111)
markers = ['^', 's', 'p', 'h', '8']
colors = np.array(["green", "yellow", "red", "pink", "black", "orange", "purple", "beige", "brown", "gray", "cyan",
                   "magenta"])
if nrecs >= 5:
    for jj in range(nSDR):
        # theta_hat_mom = np.exp(np.mean(np.log(IM_cap[jj, :])))
        # beta_hat_mom = np.std(np.log(IM_cap[jj, :]))
        theta_hat_sse, beta_hat_sse = fn_sse_pc(IM_cap[jj, :], IM_cap[jj, :].size, np.arange(IM_cap[jj, :].size))
        x_vals = np.arange(0.01, np.ceil(np.max(IM_cap[jj, :])), 0.01)
        p_collapse = norm.cdf(np.log(x_vals/theta_hat_sse)/beta_hat_sse)
        # p_collapse = norm.cdf(np.log(x_vals/theta_hat_mom)/beta_hat_mom)
        ax4.scatter(IM_cap[jj, :], np.arange(IM_cap[jj, :].size)/IM_cap[jj, :].size, marker=markers[jj],
                    s=8, c=colors[jj])
        ax4.plot(x_vals, p_collapse, c=colors[jj],
                 label=r'$j = %1.1f \/\/\theta = %1.2f \/\/\beta = %1.2f$' % (SDRCurves[jj] * 100, theta_hat_sse,
                                                                              beta_hat_sse))
        # Sa_p_05 = np.interp(0.5, p_collapse, x_vals)
        # ax4.annotate('j =' + str(SDRCurves[jj]*100) + '[%]', xy=(Sa_p_05, 0.5), xycoords='data',
        #             xytext=(20, 0), textcoords='offset points',
        #             arrowprops=dict(arrowstyle="->"))
    ax4.set_ylabel(r'$P(SDR_{max}>j|S_a=y)$')
    ax4.set_xlabel(r'$S_a(T_1)$ [g]')
    ax4.set_xlim([0, 5])
    ax4.set_ylim(bottom=0)
    ax4.grid(True)
    ax4.legend(loc='lower right', fontsize=8, fancybox=True, shadow=True, ncol=1)
    self.ui.FragilityCurve.canvas.draw()
    self.ui.FragilityCurve.canvas.show()
