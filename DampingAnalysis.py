# This is a sample Python script.

# Press Mayús+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.

import numpy as np  # load the numpy module, calling it np
import matplotlib.pyplot as plt
from scipy.stats import norm, binom, lognorm


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


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
# 4 PISOS K0
    InputCSSFile = 'ASCE_CSS_4SF_NLC_11_K0_1_0'
    IM_4SF_ASCE_KO_1_0 = np.loadtxt('CSS/' + InputCSSFile + '_IM.txt')
    SDR_max_4SF_ASCE_KO_1_0 = np.loadtxt('CSS/' + InputCSSFile + '_SDR_max.txt')
    Sa_max_4SF_ASCE_KO_1_0 = np.loadtxt('CSS/' + InputCSSFile + '_Sa_max.txt')
    RDR_max_4SF_ASCE_KO_1_0 = np.loadtxt('CSS/' + InputCSSFile + '_RDR_max.txt')
    RA_max_4SF_ASCE_KO_1_0 = np.loadtxt('CSS/' + InputCSSFile + '_RA_max.txt')
    VuVn_max_4SF_ASCE_KO_1_0 = np.loadtxt('CSS/' + InputCSSFile + '_VuVn_max.txt')
    InputCSSFile = 'T1T3_CSS_4SF_NLC_11_K0_1_0'
    IM_4SF_T1T3_KO_1_0 = np.loadtxt('CSS/' + InputCSSFile + '_IM.txt')
    SDR_max_4SF_T1T3_KO_1_0 = np.loadtxt('CSS/' + InputCSSFile + '_SDR_max.txt')
    Sa_max_4SF_T1T3_KO_1_0 = np.loadtxt('CSS/' + InputCSSFile + '_Sa_max.txt')
    RDR_max_4SF_T1T3_KO_1_0 = np.loadtxt('CSS/' + InputCSSFile + '_RDR_max.txt')
    RA_max_4SF_T1T3_KO_1_0 = np.loadtxt('CSS/' + InputCSSFile + '_RA_max.txt')
    VuVn_max_4SF_T1T3_KO_1_0 = np.loadtxt('CSS/' + InputCSSFile + '_VuVn_max.txt')
    # x, y = SDR_max_4SF_T1T3_KO_1_0*100, SDR_max_4SF_ASCE_KO_1_0*100
    # ind = np.where(x <= 6)[0]
    # x, y = x[ind], y[ind]
    # ind = np.where(y <= 6)[0]
    # x, y = x[ind], y[ind]
    # m, b = np.polyfit(x, y, 1)
    # plt.plot(x, b+m*x, '--r')
    # plt.scatter(x, y, alpha=0.5)
    # plt.xlabel(r'$SDR_{max}$ $T_{1}T_{3}$ [%]')
    # plt.ylabel(r'$SDR_{max}$ ASCE [%]')
    # plt.title(r'4 building floors $\zeta =$ 1.0% $K_{initial}$')
    # plt.axis('square')
    # plt.xlim([0, 6])
    # plt.ylim([0, 6])
    # plt.grid(True)
    # plt.show()
    
    InputCSSFile = 'ASCE_CSS_4SF_NLC_11_K0_2_5'
    IM_4SF_ASCE_KO_2_5 = np.loadtxt('CSS/' + InputCSSFile + '_IM.txt')
    SDR_max_4SF_ASCE_KO_2_5 = np.loadtxt('CSS/' + InputCSSFile + '_SDR_max.txt')
    Sa_max_4SF_ASCE_KO_2_5 = np.loadtxt('CSS/' + InputCSSFile + '_Sa_max.txt')
    RDR_max_4SF_ASCE_KO_2_5 = np.loadtxt('CSS/' + InputCSSFile + '_RDR_max.txt')
    RA_max_4SF_ASCE_KO_2_5 = np.loadtxt('CSS/' + InputCSSFile + '_RA_max.txt')
    VuVn_max_4SF_ASCE_KO_2_5 = np.loadtxt('CSS/' + InputCSSFile + '_VuVn_max.txt')
    InputCSSFile = 'T1T3_CSS_4SF_NLC_11_K0_2_5'
    IM_4SF_T1T3_KO_2_5 = np.loadtxt('CSS/' + InputCSSFile + '_IM.txt')
    SDR_max_4SF_T1T3_KO_2_5 = np.loadtxt('CSS/' + InputCSSFile + '_SDR_max.txt')
    Sa_max_4SF_T1T3_KO_2_5 = np.loadtxt('CSS/' + InputCSSFile + '_Sa_max.txt')
    RDR_max_4SF_T1T3_KO_2_5 = np.loadtxt('CSS/' + InputCSSFile + '_RDR_max.txt')
    RA_max_4SF_T1T3_KO_2_5 = np.loadtxt('CSS/' + InputCSSFile + '_RA_max.txt')
    VuVn_max_4SF_T1T3_KO_2_5 = np.loadtxt('CSS/' + InputCSSFile + '_VuVn_max.txt')
    # x, y = SDR_max_4SF_T1T3_KO_2_5*100, SDR_max_4SF_ASCE_KO_2_5*100
    # ind = np.where(x <= 6)[0]
    # x, y = x[ind], y[ind]
    # ind = np.where(y <= 6)[0]
    # x, y = x[ind], y[ind]
    # m, b = np.polyfit(x, y, 1)
    # plt.plot(x, b+m*x, '--r')
    # plt.scatter(x, y, alpha=0.5)
    # plt.xlabel(r'$SDR_{max}$ $T_{1}T_{3}$ [%]')
    # plt.ylabel(r'$SDR_{max}$ ASCE [%]')
    # plt.title(r'4 building floors $\zeta =$ 2.5% $K_{initial}$')
    # plt.axis('square')
    # plt.xlim([0, 6])
    # plt.ylim([0, 6])
    # plt.grid(True)
    # plt.show()
    
    InputCSSFile = 'ASCE_CSS_4SF_NLC_15_K0_5_0'
    IM_4SF_ASCE_KO_5_0 = np.loadtxt('CSS/' + InputCSSFile + '_IM.txt')
    SDR_max_4SF_ASCE_KO_5_0 = np.loadtxt('CSS/' + InputCSSFile + '_SDR_max.txt')
    Sa_max_4SF_ASCE_KO_5_0 = np.loadtxt('CSS/' + InputCSSFile + '_Sa_max.txt')
    RDR_max_4SF_ASCE_KO_5_0 = np.loadtxt('CSS/' + InputCSSFile + '_RDR_max.txt')
    RA_max_4SF_ASCE_KO_5_0 = np.loadtxt('CSS/' + InputCSSFile + '_RA_max.txt')
    VuVn_max_4SF_ASCE_KO_5_0 = np.loadtxt('CSS/' + InputCSSFile + '_VuVn_max.txt')
    InputCSSFile = 'T1T3_CSS_4SF_NLC_15_K0_5_0'
    IM_4SF_T1T3_KO_5_0 = np.loadtxt('CSS/' + InputCSSFile + '_IM.txt')
    SDR_max_4SF_T1T3_KO_5_0 = np.loadtxt('CSS/' + InputCSSFile + '_SDR_max.txt')
    Sa_max_4SF_T1T3_KO_5_0 = np.loadtxt('CSS/' + InputCSSFile + '_Sa_max.txt')
    RDR_max_4SF_T1T3_KO_5_0 = np.loadtxt('CSS/' + InputCSSFile + '_RDR_max.txt')
    RA_max_4SF_T1T3_KO_5_0 = np.loadtxt('CSS/' + InputCSSFile + '_RA_max.txt')
    VuVn_max_4SF_T1T3_KO_5_0 = np.loadtxt('CSS/' + InputCSSFile + '_VuVn_max.txt')
    x, y = SDR_max_4SF_ASCE_KO_1_0*100, SDR_max_4SF_ASCE_KO_5_0*100
    ind = np.where(x <= 6)[0]
    x, y = x[ind], y[ind]
    ind = np.where(y <= 6)[0]
    x, y = x[ind], y[ind]

    x1 = x[:,np.newaxis]
    a, _, _, _ = np.linalg.lstsq(x1, y)

    plt.plot(x1, a*x, 'r-')

    plt.plot(x, x, '-k')
    plt.scatter(x, y, alpha=0.5)
    plt.xlabel(r'$SDR_{max}$ $\zeta =$ 1.0% [%]')
    plt.ylabel(r'$SDR_{max}$ $\zeta =$ 5.0% [%]')
    plt.title(r'4- story building ASCE $K_{initial}$')
    plt.axis('square')
    plt.xlim([0, 6])
    plt.ylim([0, 6])
    plt.grid(True)
    plt.show()
    
# 8 PISOS K0
    InputCSSFile = 'ASCE_CSS_8SF_NLC_35_K0_1_0'
    IM_8SF_ASCE_KO_1_0 = np.loadtxt('CSS/' + InputCSSFile + '_IM.txt')
    SDR_max_8SF_ASCE_KO_1_0 = np.loadtxt('CSS/' + InputCSSFile + '_SDR_max.txt')
    Sa_max_8SF_ASCE_KO_1_0 = np.loadtxt('CSS/' + InputCSSFile + '_Sa_max.txt')
    RDR_max_8SF_ASCE_KO_1_0 = np.loadtxt('CSS/' + InputCSSFile + '_RDR_max.txt')
    RA_max_8SF_ASCE_KO_1_0 = np.loadtxt('CSS/' + InputCSSFile + '_RA_max.txt')
    VuVn_max_8SF_ASCE_KO_1_0 = np.loadtxt('CSS/' + InputCSSFile + '_VuVn_max.txt')
    InputCSSFile = 'T1T3_CSS_8SF_NLC_35_K0_1_0'
    IM_8SF_T1T3_KO_1_0 = np.loadtxt('CSS/' + InputCSSFile + '_IM.txt')
    SDR_max_8SF_T1T3_KO_1_0 = np.loadtxt('CSS/' + InputCSSFile + '_SDR_max.txt')
    Sa_max_8SF_T1T3_KO_1_0 = np.loadtxt('CSS/' + InputCSSFile + '_Sa_max.txt')
    RDR_max_8SF_T1T3_KO_1_0 = np.loadtxt('CSS/' + InputCSSFile + '_RDR_max.txt')
    RA_max_8SF_T1T3_KO_1_0 = np.loadtxt('CSS/' + InputCSSFile + '_RA_max.txt')
    VuVn_max_8SF_T1T3_KO_1_0 = np.loadtxt('CSS/' + InputCSSFile + '_VuVn_max.txt')
    # x, y = SDR_max_8SF_T1T3_KO_1_0*100, SDR_max_8SF_ASCE_KO_1_0*100
    # ind = np.where(x <= 6)[0]
    # x, y = x[ind], y[ind]
    # ind = np.where(y <= 6)[0]
    # x, y = x[ind], y[ind]
    # m, b = np.polyfit(x, y, 1)
    # plt.plot(x, b+m*x, '--r')
    # plt.scatter(x, y, alpha=0.5)
    # plt.xlabel(r'$SDR_{max}$ $T_{1}T_{3}$ [%]')
    # plt.ylabel(r'$SDR_{max}$ ASCE [%]')
    # plt.title(r'8 building floors $\zeta =$ 1.0% $K_{initial}$')
    # plt.axis('square')
    # plt.xlim([0, 6])
    # plt.ylim([0, 6])
    # plt.grid(True)
    # plt.show()
    
    InputCSSFile = 'ASCE_CSS_8SF_NLC_35_K0_2_5'
    IM_8SF_ASCE_KO_2_5 = np.loadtxt('CSS/' + InputCSSFile + '_IM.txt')
    SDR_max_8SF_ASCE_KO_2_5 = np.loadtxt('CSS/' + InputCSSFile + '_SDR_max.txt')
    Sa_max_8SF_ASCE_KO_2_5 = np.loadtxt('CSS/' + InputCSSFile + '_Sa_max.txt')
    RDR_max_8SF_ASCE_KO_2_5 = np.loadtxt('CSS/' + InputCSSFile + '_RDR_max.txt')
    RA_max_8SF_ASCE_KO_2_5 = np.loadtxt('CSS/' + InputCSSFile + '_RA_max.txt')
    VuVn_max_8SF_ASCE_KO_2_5 = np.loadtxt('CSS/' + InputCSSFile + '_VuVn_max.txt')
    InputCSSFile = 'T1T3_CSS_8SF_NLC_35_K0_2_5'
    IM_8SF_T1T3_KO_2_5 = np.loadtxt('CSS/' + InputCSSFile + '_IM.txt')
    SDR_max_8SF_T1T3_KO_2_5 = np.loadtxt('CSS/' + InputCSSFile + '_SDR_max.txt')
    Sa_max_8SF_T1T3_KO_2_5 = np.loadtxt('CSS/' + InputCSSFile + '_Sa_max.txt')
    RDR_max_8SF_T1T3_KO_2_5 = np.loadtxt('CSS/' + InputCSSFile + '_RDR_max.txt')
    RA_max_8SF_T1T3_KO_2_5 = np.loadtxt('CSS/' + InputCSSFile + '_RA_max.txt')
    VuVn_max_8SF_T1T3_KO_2_5 = np.loadtxt('CSS/' + InputCSSFile + '_VuVn_max.txt')
    # x, y = SDR_max_8SF_T1T3_KO_2_5 * 100, SDR_max_8SF_ASCE_KO_2_5 * 100
    # ind = np.where(x <= 6)[0]
    # x, y = x[ind], y[ind]
    # ind = np.where(y <= 6)[0]
    # x, y = x[ind], y[ind]
    # m, b = np.polyfit(x, y, 1)
    # plt.plot(x, b + m * x, '--r')
    # plt.scatter(x, y, alpha=0.5)
    # plt.xlabel(r'$SDR_{max}$ $T_{1}T_{3}$ [%]')
    # plt.ylabel(r'$SDR_{max}$ ASCE [%]')
    # plt.title(r'8 building floors $\zeta =$ 2.5% $K_{initial}$')
    # plt.axis('square')
    # plt.xlim([0, 6])
    # plt.ylim([0, 6])
    # plt.grid(True)
    # plt.show()

    InputCSSFile = 'ASCE_CSS_8SF_NLC_35_K0_5_0'
    IM_8SF_ASCE_KO_5_0 = np.loadtxt('CSS/' + InputCSSFile + '_IM.txt')
    SDR_max_8SF_ASCE_KO_5_0 = np.loadtxt('CSS/' + InputCSSFile + '_SDR_max.txt')
    Sa_max_8SF_ASCE_KO_5_0 = np.loadtxt('CSS/' + InputCSSFile + '_Sa_max.txt')
    RDR_max_8SF_ASCE_KO_5_0 = np.loadtxt('CSS/' + InputCSSFile + '_RDR_max.txt')
    RA_max_8SF_ASCE_KO_5_0 = np.loadtxt('CSS/' + InputCSSFile + '_RA_max.txt')
    VuVn_max_8SF_ASCE_KO_5_0 = np.loadtxt('CSS/' + InputCSSFile + '_VuVn_max.txt')
    InputCSSFile = 'T1T3_CSS_8SF_NLC_35_K0_5_0'
    IM_8SF_T1T3_KO_5_0 = np.loadtxt('CSS/' + InputCSSFile + '_IM.txt')
    SDR_max_8SF_T1T3_KO_5_0 = np.loadtxt('CSS/' + InputCSSFile + '_SDR_max.txt')
    Sa_max_8SF_T1T3_KO_5_0 = np.loadtxt('CSS/' + InputCSSFile + '_Sa_max.txt')
    RDR_max_8SF_T1T3_KO_5_0 = np.loadtxt('CSS/' + InputCSSFile + '_RDR_max.txt')
    RA_max_8SF_T1T3_KO_5_0 = np.loadtxt('CSS/' + InputCSSFile + '_RA_max.txt')
    VuVn_max_8SF_T1T3_KO_5_0 = np.loadtxt('CSS/' + InputCSSFile + '_VuVn_max.txt')
    x, y = SDR_max_8SF_ASCE_KO_1_0*100, SDR_max_8SF_ASCE_KO_5_0*100
    plt.plot(x, x, '--r')
    plt.scatter(x, y, alpha=0.5)
    plt.xlabel(r'$SDR_{max}$ $\zeta =$ 1.0% [%]')
    plt.ylabel(r'$SDR_{max}$ $\zeta =$ 5.0% [%]')
    plt.title(r'8-story building ASCE $K_{initial}$')
    plt.axis('square')
    plt.xlim([0, 6])
    plt.ylim([0, 6])
    plt.grid(True)
    plt.show()

# 16 PISOS K0
    InputCSSFile = 'ASCE_CSS_16SF_NLC_83_K0_1_0'
    IM_16SF_ASCE_KO_1_0 = np.loadtxt('CSS/' + InputCSSFile + '_IM.txt')
    SDR_max_16SF_ASCE_KO_1_0 = np.loadtxt('CSS/' + InputCSSFile + '_SDR_max.txt')
    Sa_max_16SF_ASCE_KO_1_0 = np.loadtxt('CSS/' + InputCSSFile + '_Sa_max.txt')
    RDR_max_16SF_ASCE_KO_1_0 = np.loadtxt('CSS/' + InputCSSFile + '_RDR_max.txt')
    RA_max_16SF_ASCE_KO_1_0 = np.loadtxt('CSS/' + InputCSSFile + '_RA_max.txt')
    VuVn_max_16SF_ASCE_KO_1_0 = np.loadtxt('CSS/' + InputCSSFile + '_VuVn_max.txt')
    InputCSSFile = 'T1T3_CSS_16SF_NLC_83_K0_1_0'
    IM_16SF_T1T3_KO_1_0 = np.loadtxt('CSS/' + InputCSSFile + '_IM.txt')
    SDR_max_16SF_T1T3_KO_1_0 = np.loadtxt('CSS/' + InputCSSFile + '_SDR_max.txt')
    Sa_max_16SF_T1T3_KO_1_0 = np.loadtxt('CSS/' + InputCSSFile + '_Sa_max.txt')
    RDR_max_16SF_T1T3_KO_1_0 = np.loadtxt('CSS/' + InputCSSFile + '_RDR_max.txt')
    RA_max_16SF_T1T3_KO_1_0 = np.loadtxt('CSS/' + InputCSSFile + '_RA_max.txt')
    VuVn_max_16SF_T1T3_KO_1_0 = np.loadtxt('CSS/' + InputCSSFile + '_VuVn_max.txt')
    # x, y = SDR_max_16SF_T1T3_KO_1_0 * 100, SDR_max_16SF_ASCE_KO_1_0 * 100
    # ind = np.where(x <= 6)[0]
    # x, y = x[ind], y[ind]
    # ind = np.where(y <= 6)[0]
    # x, y = x[ind], y[ind]
    # m, b = np.polyfit(x, y, 1)
    # plt.plot(x, b + m * x, '--r')
    # plt.scatter(x, y, alpha=0.5)
    # plt.xlabel(r'$SDR_{max}$ $T_{1}T_{3}$ [%]')
    # plt.ylabel(r'$SDR_{max}$ ASCE [%]')
    # plt.title(r'16 building floors $\zeta =$ 1.0% $K_{initial}$')
    # plt.axis('square')
    # plt.xlim([0, 6])
    # plt.ylim([0, 6])
    # plt.grid(True)
    # plt.show()

    InputCSSFile = 'ASCE_CSS_16SF_NLC_83_K0_2_5'
    IM_16SF_ASCE_KO_2_5 = np.loadtxt('CSS/' + InputCSSFile + '_IM.txt')
    SDR_max_16SF_ASCE_KO_2_5 = np.loadtxt('CSS/' + InputCSSFile + '_SDR_max.txt')
    Sa_max_16SF_ASCE_KO_2_5 = np.loadtxt('CSS/' + InputCSSFile + '_Sa_max.txt')
    RDR_max_16SF_ASCE_KO_2_5 = np.loadtxt('CSS/' + InputCSSFile + '_RDR_max.txt')
    RA_max_16SF_ASCE_KO_2_5 = np.loadtxt('CSS/' + InputCSSFile + '_RA_max.txt')
    VuVn_max_16SF_ASCE_KO_2_5 = np.loadtxt('CSS/' + InputCSSFile + '_VuVn_max.txt')
    InputCSSFile = 'T1T3_CSS_16SF_NLC_83_K0_2_5'
    IM_16SF_T1T3_KO_2_5 = np.loadtxt('CSS/' + InputCSSFile + '_IM.txt')
    SDR_max_16SF_T1T3_KO_2_5 = np.loadtxt('CSS/' + InputCSSFile + '_SDR_max.txt')
    Sa_max_16SF_T1T3_KO_2_5 = np.loadtxt('CSS/' + InputCSSFile + '_Sa_max.txt')
    RDR_max_16SF_T1T3_KO_2_5 = np.loadtxt('CSS/' + InputCSSFile + '_RDR_max.txt')
    RA_max_16SF_T1T3_KO_2_5 = np.loadtxt('CSS/' + InputCSSFile + '_RA_max.txt')
    VuVn_max_16SF_T1T3_KO_2_5 = np.loadtxt('CSS/' + InputCSSFile + '_VuVn_max.txt')
    # x, y = SDR_max_16SF_T1T3_KO_2_5 * 100, SDR_max_16SF_ASCE_KO_2_5 * 100
    # ind = np.where(x <= 6)[0]
    # x, y = x[ind], y[ind]
    # ind = np.where(y <= 6)[0]
    # x, y = x[ind], y[ind]
    # m, b = np.polyfit(x, y, 1)
    # plt.plot(x, b + m * x, '--r')
    # plt.scatter(x, y, alpha=0.5)
    # plt.xlabel(r'$SDR_{max}$ $T_{1}T_{3}$ [%]')
    # plt.ylabel(r'$SDR_{max}$ ASCE [%]')
    # plt.title(r'16 building floors $\zeta =$ 2.5% $K_{initial}$')
    # plt.axis('square')
    # plt.xlim([0, 6])
    # plt.ylim([0, 6])
    # plt.grid(True)
    # plt.show()

    InputCSSFile = 'ASCE_CSS_16SF_NLC_83_K0_5_0'
    IM_16SF_ASCE_KO_5_0 = np.loadtxt('CSS/' + InputCSSFile + '_IM.txt')
    SDR_max_16SF_ASCE_KO_5_0 = np.loadtxt('CSS/' + InputCSSFile + '_SDR_max.txt')
    Sa_max_16SF_ASCE_KO_5_0 = np.loadtxt('CSS/' + InputCSSFile + '_Sa_max.txt')
    RDR_max_16SF_ASCE_KO_5_0 = np.loadtxt('CSS/' + InputCSSFile + '_RDR_max.txt')
    RA_max_16SF_ASCE_KO_5_0 = np.loadtxt('CSS/' + InputCSSFile + '_RA_max.txt')
    VuVn_max_16SF_ASCE_KO_5_0 = np.loadtxt('CSS/' + InputCSSFile + '_VuVn_max.txt')
    InputCSSFile = 'T1T3_CSS_16SF_NLC_83_K0_5_0'
    IM_16SF_T1T3_KO_5_0 = np.loadtxt('CSS/' + InputCSSFile + '_IM.txt')
    SDR_max_16SF_T1T3_KO_5_0 = np.loadtxt('CSS/' + InputCSSFile + '_SDR_max.txt')
    Sa_max_16SF_T1T3_KO_5_0 = np.loadtxt('CSS/' + InputCSSFile + '_Sa_max.txt')
    RDR_max_16SF_T1T3_KO_5_0 = np.loadtxt('CSS/' + InputCSSFile + '_RDR_max.txt')
    RA_max_16SF_T1T3_KO_5_0 = np.loadtxt('CSS/' + InputCSSFile + '_RA_max.txt')
    VuVn_max_16SF_T1T3_KO_5_0 = np.loadtxt('CSS/' + InputCSSFile + '_VuVn_max.txt')
    x, y = SDR_max_16SF_ASCE_KO_1_0*100, SDR_max_16SF_ASCE_KO_5_0*100
    plt.plot(x, x, '--r')
    plt.scatter(x, y, alpha=0.5)
    plt.xlabel(r'$SDR_{max}$ $\zeta =$ 1.0% [%]')
    plt.ylabel(r'$SDR_{max}$ $\zeta =$ 5.0% [%]')
    plt.title(r'16-story building ASCE $K_{initial}$')
    plt.axis('square')
    plt.xlim([0, 6])
    plt.ylim([0, 6])
    plt.grid(True)
    plt.show()

# 4 PISOS KCommited
    InputCSSFile = 'CSS_4SF_NLC_Utah_11'
    IM_4SF_ASCE_KC_2_5 = np.loadtxt('CSS/' + InputCSSFile + '_IM.txt')
    SDR_max_4SF_ASCE_KC_2_5 = np.loadtxt('CSS/' + InputCSSFile + '_SDR_max.txt')
    Sa_max_4SF_ASCE_KC_2_5 = np.loadtxt('CSS/' + InputCSSFile + '_Sa_max.txt')
    RDR_max_4SF_ASCE_KC_2_5 = np.loadtxt('CSS/' + InputCSSFile + '_RDR_max.txt')
    RA_max_4SF_ASCE_KC_2_5 = np.loadtxt('CSS/' + InputCSSFile + '_RA_max.txt')
    VuVn_max_4SF_ASCE_KC_2_5 = np.loadtxt('CSS/' + InputCSSFile + '_VuVn_max.txt')
    InputCSSFile = 'RayleighCSS_4SF_NLC_Utah_11'
    IM_4SF_T1T3_KC_2_5 = np.loadtxt('CSS/' + InputCSSFile + '_IM.txt')
    SDR_max_4SF_T1T3_KC_2_5 = np.loadtxt('CSS/' + InputCSSFile + '_SDR_max.txt')
    Sa_max_4SF_T1T3_KC_2_5 = np.loadtxt('CSS/' + InputCSSFile + '_Sa_max.txt')
    RDR_max_4SF_T1T3_KC_2_5 = np.loadtxt('CSS/' + InputCSSFile + '_RDR_max.txt')
    RA_max_4SF_T1T3_KC_2_5 = np.loadtxt('CSS/' + InputCSSFile + '_RA_max.txt')
    VuVn_max_4SF_T1T3_KC_2_5 = np.loadtxt('CSS/' + InputCSSFile + '_VuVn_max.txt')
    x, y = SDR_max_4SF_ASCE_KC_2_5*100, SDR_max_4SF_ASCE_KO_2_5*100
    plt.plot(x, x, '--r')
    plt.scatter(x, y, alpha=0.5)
    plt.xlabel(r'$SDR_{max}$ $K_{committed}$ [%]')
    plt.ylabel(r'$SDR_{max}$ $K_{initial}$ [%]')
    plt.title(r'4-story building ASCE $\zeta =$ 2.5%')
    plt.axis('square')
    plt.xlim([0, 6])
    plt.ylim([0, 6])
    plt.grid(True)
    plt.show()

# 8 PISOS KCommited
    InputCSSFile = 'CSS_8SF_NLC_Utah_35'
    IM_8SF_ASCE_KC_2_5 = np.loadtxt('CSS/' + InputCSSFile + '_IM.txt')
    SDR_max_8SF_ASCE_KC_2_5 = np.loadtxt('CSS/' + InputCSSFile + '_SDR_max.txt')
    Sa_max_8SF_ASCE_KC_2_5 = np.loadtxt('CSS/' + InputCSSFile + '_Sa_max.txt')
    RDR_max_8SF_ASCE_KC_2_5 = np.loadtxt('CSS/' + InputCSSFile + '_RDR_max.txt')
    RA_max_8SF_ASCE_KC_2_5 = np.loadtxt('CSS/' + InputCSSFile + '_RA_max.txt')
    VuVn_max_8SF_ASCE_KC_2_5 = np.loadtxt('CSS/' + InputCSSFile + '_VuVn_max.txt')
    InputCSSFile = 'RayleighCSS_8SF_NLC_Utah_35'
    IM_8SF_T1T3_KC_2_5 = np.loadtxt('CSS/' + InputCSSFile + '_IM.txt')
    SDR_max_8SF_T1T3_KC_2_5 = np.loadtxt('CSS/' + InputCSSFile + '_SDR_max.txt')
    Sa_max_8SF_T1T3_KC_2_5 = np.loadtxt('CSS/' + InputCSSFile + '_Sa_max.txt')
    RDR_max_8SF_T1T3_KC_2_5 = np.loadtxt('CSS/' + InputCSSFile + '_RDR_max.txt')
    RA_max_8SF_T1T3_KC_2_5 = np.loadtxt('CSS/' + InputCSSFile + '_RA_max.txt')
    VuVn_max_8SF_T1T3_KC_2_5 = np.loadtxt('CSS/' + InputCSSFile + '_VuVn_max.txt')
    x, y = SDR_max_8SF_ASCE_KC_2_5*100, SDR_max_8SF_ASCE_KO_2_5*100
    plt.plot(x, x, '--r')
    plt.scatter(x, y, alpha=0.5)
    plt.xlabel(r'$SDR_{max}$ $K_{committed}$ [%]')
    plt.ylabel(r'$SDR_{max}$ $K_{initial}$ [%]')
    plt.title(r'8-story building ASCE $\zeta =$ 2.5%')
    plt.axis('square')
    plt.xlim([0, 6])
    plt.ylim([0, 6])
    plt.grid(True)
    plt.show()

# 16 PISOS KCommited
    InputCSSFile = 'CSS_16SF_NLC_Utah_83'
    IM_16SF_ASCE_KC_2_5 = np.loadtxt('CSS/' + InputCSSFile + '_IM.txt')
    SDR_max_16SF_ASCE_KC_2_5 = np.loadtxt('CSS/' + InputCSSFile + '_SDR_max.txt')
    Sa_max_16SF_ASCE_KC_2_5 = np.loadtxt('CSS/' + InputCSSFile + '_Sa_max.txt')
    RDR_max_16SF_ASCE_KC_2_5 = np.loadtxt('CSS/' + InputCSSFile + '_RDR_max.txt')
    RA_max_16SF_ASCE_KC_2_5 = np.loadtxt('CSS/' + InputCSSFile + '_RA_max.txt')
    VuVn_max_16SF_ASCE_KC_2_5 = np.loadtxt('CSS/' + InputCSSFile + '_VuVn_max.txt')
    InputCSSFile = 'RayleighCSS_16SF_NLC_Utah_83'
    IM_16SF_T1T3_KC_2_5 = np.loadtxt('CSS/' + InputCSSFile + '_IM.txt')
    SDR_max_16SF_T1T3_KC_2_5 = np.loadtxt('CSS/' + InputCSSFile + '_SDR_max.txt')
    Sa_max_16SF_T1T3_KC_2_5 = np.loadtxt('CSS/' + InputCSSFile + '_Sa_max.txt')
    RDR_max_16SF_T1T3_KC_2_5 = np.loadtxt('CSS/' + InputCSSFile + '_RDR_max.txt')
    RA_max_16SF_T1T3_KC_2_5 = np.loadtxt('CSS/' + InputCSSFile + '_RA_max.txt')
    VuVn_max_16SF_T1T3_KC_2_5 = np.loadtxt('CSS/' + InputCSSFile + '_VuVn_max.txt')
    x, y = SDR_max_16SF_ASCE_KC_2_5*100, SDR_max_16SF_ASCE_KO_2_5*100
    plt.plot(x, x, '--r')
    plt.scatter(x, y, alpha=0.5)
    plt.xlabel(r'$SDR_{max}$ $K_{committed}$ [%]')
    plt.ylabel(r'$SDR_{max}$ $K_{initial}$ [%]')
    plt.title(r'16-story building ASCE $\zeta =$ 2.5%')
    plt.axis('square')
    plt.xlim([0, 6])
    plt.ylim([0, 6])
    plt.grid(True)
    plt.show()

    CSS_SDRCurves = np.linspace(0.2, 6, 20)/100
    clrs_v = ['r', 'r', 'r', 'b', 'b', 'b', 'g', 'k']  # list of basic colors
    styl_v = ['-', '--', ':', '-', '--', ':', '--', '--']  # list of basic linestyles

    IM_v = np.array([IM_4SF_ASCE_KO_1_0, IM_4SF_ASCE_KO_2_5, IM_4SF_ASCE_KO_5_0, IM_4SF_T1T3_KO_1_0,
                     IM_4SF_T1T3_KO_2_5, IM_4SF_T1T3_KO_5_0, IM_4SF_ASCE_KC_2_5, IM_4SF_T1T3_KC_2_5])
    SDR_v = np.array([SDR_max_4SF_ASCE_KO_1_0, SDR_max_4SF_ASCE_KO_2_5, SDR_max_4SF_ASCE_KO_5_0,
                      SDR_max_4SF_T1T3_KO_1_0, SDR_max_4SF_T1T3_KO_2_5, SDR_max_4SF_T1T3_KO_5_0,
                      SDR_max_4SF_ASCE_KC_2_5, SDR_max_4SF_T1T3_KC_2_5])
    Ident_v = np.array(['ASCE $\zeta =$ 1.0% $K_{initial}$', 'ASCE $\zeta =$ 2.5% $K_{initial}$',
                        'ASCE $\zeta =$ 5.0% $K_{initial}$', '$T_{1}T_{3}$  $\zeta =$ 1.0% $K_{initial}$',
                        '$T_{1}T_{3}$  $\zeta =$ 2.5% $K_{initial}$', '$T_{1}T_{3}$  $\zeta =$ 5.0% $K_{initial}$',
                        'ASCE $\zeta =$ 2.5% $K_{committed}$', '$T_{1}T_{3}$  $\zeta =$ 2.5% $K_{committed}$'])

    for (IM, SDR, Ident, clrs, styl) in zip(IM_v, SDR_v, Ident_v, clrs_v, styl_v):
        hist_tot, bin_edges_tot = np.histogram(IM, bins=np.geomspace(0.0001, 5, 50))
        bin_avg_tot = (bin_edges_tot[0:-1] + bin_edges_tot[1:]) / 2
        ind_t = np.where(hist_tot > 0)[0]
        hist_tot = hist_tot[ind_t]
        bin_avg_tot = bin_avg_tot[ind_t]
        theta_v, beta_v = [], []
        for SDRCurves in CSS_SDRCurves:
            ind = np.where(SDR >= SDRCurves)[0]
            hist, bin_edges = np.histogram(IM[ind], bins=np.geomspace(0.0001, 5, 50))
            hist = hist[ind_t]
            bin_avg = (bin_edges[0:-1] + bin_edges[1:]) / 2
            bin_avg = bin_avg[ind_t]
            theta_hat_sse, beta_hat_sse = fn_sse_pc(bin_avg, hist_tot, hist)
            theta_v.append(theta_hat_sse)
            beta_v.append(beta_hat_sse)
            # print(SDRCurves, theta_hat_sse, beta_hat_sse)
        theta_v = np.array(theta_v)
        beta_v = np.array(beta_v)
        # print(CSS_SDRCurves)
        # print(theta_v)
        # print(beta_v)
        plt.plot(CSS_SDRCurves*100, theta_v, label=Ident, color=clrs, ls=styl)
    plt.grid(True)
    plt.legend(loc='upper left', fontsize=8, fancybox=True, shadow=True, ncol=1)
    plt.ylabel(r'$\theta$')
    plt.xlabel(r'j = SDR [%]')
    plt.title('4-story building')
    plt.xlim(0.2, 6)
    plt.ylim(bottom=0)
    plt.show()

    for (IM, SDR, Ident, clrs, styl) in zip(IM_v, SDR_v, Ident_v, clrs_v, styl_v):
        hist_tot, bin_edges_tot = np.histogram(IM, bins=np.geomspace(0.0001, 5, 50))
        bin_avg_tot = (bin_edges_tot[0:-1] + bin_edges_tot[1:]) / 2
        ind_t = np.where(hist_tot > 0)[0]
        hist_tot = hist_tot[ind_t]
        bin_avg_tot = bin_avg_tot[ind_t]
        theta_v, beta_v = [], []
        for SDRCurves in CSS_SDRCurves:
            ind = np.where(SDR >= SDRCurves)[0]
            hist, bin_edges = np.histogram(IM[ind], bins=np.geomspace(0.0001, 5, 50))
            hist = hist[ind_t]
            bin_avg = (bin_edges[0:-1] + bin_edges[1:]) / 2
            bin_avg = bin_avg[ind_t]
            theta_hat_sse, beta_hat_sse = fn_sse_pc(bin_avg, hist_tot, hist)
            theta_v.append(theta_hat_sse)
            beta_v.append(beta_hat_sse)
            # print(SDRCurves, theta_hat_sse, beta_hat_sse)
        theta_v = np.array(theta_v)
        beta_v = np.array(beta_v)
        # print(CSS_SDRCurves)
        # print(theta_v)
        # print(beta_v)
        plt.plot(CSS_SDRCurves * 100, beta_v, label=Ident, color=clrs, ls=styl)
    plt.grid(True)
    plt.legend(loc='upper left', fontsize=8, fancybox=True, shadow=True, ncol=1)
    plt.ylabel(r'$\beta$')
    plt.xlabel(r'j = SDR [%]')
    plt.title('4-story building')
    plt.xlim(0.2, 6)
    plt.ylim(bottom=0)
    plt.show()

#8 Pisos
CSS_SDRCurves = np.linspace(0.2, 6, 20) / 100
clrs_v = ['r', 'r', 'r', 'b', 'b', 'b', 'g', 'k']  # list of basic colors
styl_v = ['-', '--', ':', '-', '--', ':', '--', '--']  # list of basic linestyles

IM_v = np.array([IM_8SF_ASCE_KO_1_0, IM_8SF_ASCE_KO_2_5, IM_8SF_ASCE_KO_5_0, IM_8SF_T1T3_KO_1_0,
                 IM_8SF_T1T3_KO_2_5, IM_8SF_T1T3_KO_5_0, IM_8SF_ASCE_KC_2_5, IM_8SF_T1T3_KC_2_5])
SDR_v = np.array([SDR_max_8SF_ASCE_KO_1_0, SDR_max_8SF_ASCE_KO_2_5, SDR_max_8SF_ASCE_KO_5_0,
                  SDR_max_8SF_T1T3_KO_1_0, SDR_max_8SF_T1T3_KO_2_5, SDR_max_8SF_T1T3_KO_5_0,
                  SDR_max_8SF_ASCE_KC_2_5, SDR_max_8SF_T1T3_KC_2_5])
Ident_v = np.array(['ASCE $\zeta =$ 1.0% $K_{initial}$', 'ASCE $\zeta =$ 2.5% $K_{initial}$',
                    'ASCE $\zeta =$ 5.0% $K_{initial}$', '$T_{1}T_{3}$  $\zeta =$ 1.0% $K_{initial}$',
                    '$T_{1}T_{3}$  $\zeta =$ 2.5% $K_{initial}$', '$T_{1}T_{3}$  $\zeta =$ 5.0% $K_{initial}$',
                    'ASCE $\zeta =$ 2.5% $K_{committed}$', '$T_{1}T_{3}$  $\zeta =$ 2.5% $K_{committed}$'])

for (IM, SDR, Ident, clrs, styl) in zip(IM_v, SDR_v, Ident_v, clrs_v, styl_v):
    hist_tot, bin_edges_tot = np.histogram(IM, bins=np.geomspace(0.0001, 5, 50))
    bin_avg_tot = (bin_edges_tot[0:-1] + bin_edges_tot[1:]) / 2
    ind_t = np.where(hist_tot > 0)[0]
    hist_tot = hist_tot[ind_t]
    bin_avg_tot = bin_avg_tot[ind_t]
    theta_v, beta_v = [], []
    for SDRCurves in CSS_SDRCurves:
        ind = np.where(SDR >= SDRCurves)[0]
        hist, bin_edges = np.histogram(IM[ind], bins=np.geomspace(0.0001, 5, 50))
        hist = hist[ind_t]
        bin_avg = (bin_edges[0:-1] + bin_edges[1:]) / 2
        bin_avg = bin_avg[ind_t]
        theta_hat_sse, beta_hat_sse = fn_sse_pc(bin_avg, hist_tot, hist)
        theta_v.append(theta_hat_sse)
        beta_v.append(beta_hat_sse)
        # print(SDRCurves, theta_hat_sse, beta_hat_sse)
    theta_v = np.array(theta_v)
    beta_v = np.array(beta_v)
    # print(CSS_SDRCurves)
    # print(theta_v)
    # print(beta_v)
    plt.plot(CSS_SDRCurves * 100, theta_v, label=Ident, color=clrs, ls=styl)
plt.grid(True)
plt.legend(loc='upper left', fontsize=8, fancybox=True, shadow=True, ncol=1)
plt.ylabel(r'$\theta$')
plt.xlabel(r'j = SDR [%]')
plt.title('8-story building')
plt.xlim(0.2, 6)
plt.ylim(bottom=0)
plt.show()

for (IM, SDR, Ident, clrs, styl) in zip(IM_v, SDR_v, Ident_v, clrs_v, styl_v):
    hist_tot, bin_edges_tot = np.histogram(IM, bins=np.geomspace(0.0001, 5, 50))
    bin_avg_tot = (bin_edges_tot[0:-1] + bin_edges_tot[1:]) / 2
    ind_t = np.where(hist_tot > 0)[0]
    hist_tot = hist_tot[ind_t]
    bin_avg_tot = bin_avg_tot[ind_t]
    theta_v, beta_v = [], []
    for SDRCurves in CSS_SDRCurves:
        ind = np.where(SDR >= SDRCurves)[0]
        hist, bin_edges = np.histogram(IM[ind], bins=np.geomspace(0.0001, 5, 50))
        hist = hist[ind_t]
        bin_avg = (bin_edges[0:-1] + bin_edges[1:]) / 2
        bin_avg = bin_avg[ind_t]
        theta_hat_sse, beta_hat_sse = fn_sse_pc(bin_avg, hist_tot, hist)
        theta_v.append(theta_hat_sse)
        beta_v.append(beta_hat_sse)
        # print(SDRCurves, theta_hat_sse, beta_hat_sse)
    theta_v = np.array(theta_v)
    beta_v = np.array(beta_v)
    # print(CSS_SDRCurves)
    # print(theta_v)
    # print(beta_v)
    plt.plot(CSS_SDRCurves * 100, beta_v, label=Ident, color=clrs, ls=styl)
plt.grid(True)
plt.legend(loc='upper left', fontsize=8, fancybox=True, shadow=True, ncol=1)
plt.ylabel(r'$\beta$')
plt.xlabel(r'j = SDR [%]')
plt.title('8-story building')
plt.xlim(0.2, 6)
plt.ylim(bottom=0)
plt.show()

#16 Pisos
CSS_SDRCurves = np.linspace(0.2, 4, 20) / 100
clrs_v = ['r', 'r', 'r', 'b', 'b', 'b', 'g', 'k']  # list of basic colors
styl_v = ['-', '--', ':', '-', '--', ':', '--', '--']  # list of basic linestyles

IM_v = np.array([IM_16SF_ASCE_KO_1_0, IM_16SF_ASCE_KO_2_5, IM_16SF_ASCE_KO_5_0, IM_16SF_T1T3_KO_1_0,
                 IM_16SF_T1T3_KO_2_5, IM_16SF_T1T3_KO_5_0, IM_16SF_ASCE_KC_2_5, IM_16SF_T1T3_KC_2_5])
SDR_v = np.array([SDR_max_16SF_ASCE_KO_1_0, SDR_max_16SF_ASCE_KO_2_5, SDR_max_16SF_ASCE_KO_5_0,
                  SDR_max_16SF_T1T3_KO_1_0, SDR_max_16SF_T1T3_KO_2_5, SDR_max_16SF_T1T3_KO_5_0,
                  SDR_max_16SF_ASCE_KC_2_5, SDR_max_16SF_T1T3_KC_2_5])
Ident_v = np.array(['ASCE $\zeta =$ 1.0% $K_{initial}$', 'ASCE $\zeta =$ 2.5% $K_{initial}$',
                    'ASCE $\zeta =$ 5.0% $K_{initial}$', '$T_{1}T_{3}$  $\zeta =$ 1.0% $K_{initial}$',
                    '$T_{1}T_{3}$  $\zeta =$ 2.5% $K_{initial}$', '$T_{1}T_{3}$  $\zeta =$ 5.0% $K_{initial}$',
                    'ASCE $\zeta =$ 2.5% $K_{committed}$', '$T_{1}T_{3}$  $\zeta =$ 2.5% $K_{committed}$'])

for (IM, SDR, Ident, clrs, styl) in zip(IM_v, SDR_v, Ident_v, clrs_v, styl_v):
    hist_tot, bin_edges_tot = np.histogram(IM, bins=np.geomspace(0.0001, 5, 50))
    bin_avg_tot = (bin_edges_tot[0:-1] + bin_edges_tot[1:]) / 2
    ind_t = np.where(hist_tot > 0)[0]
    hist_tot = hist_tot[ind_t]
    bin_avg_tot = bin_avg_tot[ind_t]
    theta_v, beta_v = [], []
    for SDRCurves in CSS_SDRCurves:
        ind = np.where(SDR >= SDRCurves)[0]
        hist, bin_edges = np.histogram(IM[ind], bins=np.geomspace(0.0001, 5, 50))
        hist = hist[ind_t]
        bin_avg = (bin_edges[0:-1] + bin_edges[1:]) / 2
        bin_avg = bin_avg[ind_t]
        theta_hat_sse, beta_hat_sse = fn_sse_pc(bin_avg, hist_tot, hist)
        theta_v.append(theta_hat_sse)
        beta_v.append(beta_hat_sse)
        # print(SDRCurves, theta_hat_sse, beta_hat_sse)
    theta_v = np.array(theta_v)
    beta_v = np.array(beta_v)
    # print(CSS_SDRCurves)
    # print(theta_v)
    # print(beta_v)
    plt.plot(CSS_SDRCurves * 100, theta_v, label=Ident, color=clrs, ls=styl)
plt.grid(True)
plt.legend(loc='upper left', fontsize=8, fancybox=True, shadow=True, ncol=1)
plt.ylabel(r'$\theta$')
plt.xlabel(r'j = SDR [%]')
plt.title('16-story building')
plt.xlim(0.2, 4)
plt.ylim(bottom=0)
plt.show()

for (IM, SDR, Ident, clrs, styl) in zip(IM_v, SDR_v, Ident_v, clrs_v, styl_v):
    hist_tot, bin_edges_tot = np.histogram(IM, bins=np.geomspace(0.0001, 5, 50))
    bin_avg_tot = (bin_edges_tot[0:-1] + bin_edges_tot[1:]) / 2
    ind_t = np.where(hist_tot > 0)[0]
    hist_tot = hist_tot[ind_t]
    bin_avg_tot = bin_avg_tot[ind_t]
    theta_v, beta_v = [], []
    for SDRCurves in CSS_SDRCurves:
        ind = np.where(SDR >= SDRCurves)[0]
        hist, bin_edges = np.histogram(IM[ind], bins=np.geomspace(0.0001, 5, 50))
        hist = hist[ind_t]
        bin_avg = (bin_edges[0:-1] + bin_edges[1:]) / 2
        bin_avg = bin_avg[ind_t]
        theta_hat_sse, beta_hat_sse = fn_sse_pc(bin_avg, hist_tot, hist)
        theta_v.append(theta_hat_sse)
        beta_v.append(beta_hat_sse)
        # print(SDRCurves, theta_hat_sse, beta_hat_sse)
    theta_v = np.array(theta_v)
    beta_v = np.array(beta_v)
    # print(CSS_SDRCurves)
    # print(theta_v)
    # print(beta_v)
    plt.plot(CSS_SDRCurves * 100, beta_v, label=Ident, color=clrs, ls=styl)
plt.grid(True)
plt.legend(loc='upper left', fontsize=8, fancybox=True, shadow=True, ncol=1)
plt.ylabel(r'$\beta$')
plt.xlabel(r'j = SDR [%]')
plt.title('16-story building')
plt.xlim(0.2, 4)
plt.ylim(bottom=0)
plt.show()


