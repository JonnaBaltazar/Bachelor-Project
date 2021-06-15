# Show plots inline, and load main getdist plot module and samples class

import numpy as np
from pylab import *
from scipy.interpolate import interp1d
from matplotlib import *
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.patheffects as fx
from matplotlib import patches
from matplotlib import text as mtext
from scipy.signal import savgol_filter
import os

conv_eV_to_inv_s = 1.0/(6.582110e-16) # convert from eV to s^{-1}

def decay_rate_ij(mnui, mnuj, gij, gpij):

    x = mnui/mnuj
    f = x/2. + 2.0 + 2./x*log(x) - 2.0/x/x - 1.0/(2.0*x*x*x)
    h = x/2.0 -2.0 + 2./x*log(x) + 2.0/x/x - 1.0/(2.0*x*x*x)
    k = x/2.0 - 2./x*log(x) - 1.0/(2.0*x*x*x)

    Gamma = mnuj/16./np.pi * ( gij*gij*(f+k) + gpij*gpij*(h+k) ) 
    Gamma = conv_eV_to_inv_s*Gamma
   
    return Gamma

Dm21sq = 7.42e-5  # [eV^2]
Dm31sq = 2.517e-3 # [eV^2]

fig = plt.figure(figsize=[9,9])
ax = fig.add_subplot(1,1,1)

x_min = -4.0
x_max = -1.55129 # log10(0.0281) = -1.55129 
y_min = -4.0
y_max = 0.0 

################################################################################
# IceCube, 8yr
################################################################################

# Lower limit on tau2 (from nu2 --> nu1)
filename_in = '../../data/stat_analysis/ic_8yr/fixed_log10mnu1/limits_log10g21_vs_log10mnu1.dat'
log10mnu1, \
    log10g21_median, log10g21_sigma, \
    log10g21_1sigma_lo, log10g21_1sigma_hi, \
    log10g21_2sigma_lo, log10g21_2sigma_hi, \
    log10g21_3sigma_lo, log10g21_3sigma_hi, \
    log10g21_5sigma_lo, log10g21_5sigma_hi, \
    log10g21_90CL_lo, log10g21_90CL_hi \
    = np.loadtxt(filename_in, unpack=True)
Gamma21 = [decay_rate_ij(   sqrt(10.**log10mnu1[k]+Dm21sq), # = mnu2
                            10.**log10mnu1[k],
                            10.**log10g21_2sigma_hi[k],
                            10.**log10g21_2sigma_hi[k])
            for k in range(len(log10mnu1))] # [s^{-1}]
log10tau2 = [log10(1.0/x) for x in Gamma21] # [s]

# Lower limit on tau3 (from nu3 --> nu1 and nu3 --> nu2)
# nu3 --> nu1
filename_in = '../../data/stat_analysis/ic_8yr/fixed_log10mnu1/limits_log10g31_vs_log10mnu1.dat'
log10mnu1, \
    log10g31_median, log10g31_sigma, \
    log10g31_1sigma_lo, log10g31_1sigma_hi, \
    log10g31_2sigma_lo, log10g31_2sigma_hi, \
    log10g31_3sigma_lo, log10g31_3sigma_hi, \
    log10g31_5sigma_lo, log10g31_5sigma_hi, \
    log10g31_90CL_lo, log10g31_90CL_hi \
    = np.loadtxt(filename_in, unpack=True)
Gamma31 = [decay_rate_ij(   sqrt(10.**log10mnu1[k]+Dm31sq), # = mnu3
                            10.**log10mnu1[k],
                            10.**log10g31_2sigma_hi[k],
                            10.**log10g31_2sigma_hi[k])
            for k in range(len(log10mnu1))] # [s^{-1}]
# nu3 --> nu2
filename_in = '../../data/stat_analysis/ic_8yr/fixed_log10mnu1/limits_log10g32_vs_log10mnu1.dat'
log10mnu1, \
    log10g32_median, log10g32_sigma, \
    log10g32_1sigma_lo, log10g32_1sigma_hi, \
    log10g32_2sigma_lo, log10g32_2sigma_hi, \
    log10g32_3sigma_lo, log10g32_3sigma_hi, \
    log10g32_5sigma_lo, log10g32_5sigma_hi, \
    log10g32_90CL_lo, log10g32_90CL_hi \
    = np.loadtxt(filename_in, unpack=True)
Gamma32 = [decay_rate_ij(   sqrt(10.**log10mnu1[k]+Dm31sq), # = mnu3
                            10.**log10mnu1[k],
                            10.**log10g32_2sigma_hi[k],
                            10.**log10g32_2sigma_hi[k])
            for k in range(len(log10mnu1))] # [s^{-1}]
# Total
Gamma3 = [Gamma31[k]+Gamma32[k] for k in range(len(Gamma31))] # [s^{-1}]
log10tau3 = [log10(1.0/x) for x in Gamma3] # [s]

window_length = 13 #25 
polyorder = 1 #3
log10tau2 = savgol_filter(log10tau2,
    window_length, polyorder, deriv=0, mode='interp')

window_length = 13 #25 
polyorder = 1 #3
log10tau3 = savgol_filter(log10tau3,
    window_length, polyorder, deriv=0, mode='interp')

ax.plot(log10mnu1, log10tau2,
    c='firebrick', alpha=1.0, ls='-', lw=3.0, zorder=20,
    label=r'$\tau_2$, IceCube 8 yr')

ax.plot(log10mnu1, log10tau3,
    c='firebrick', alpha=1.0, ls='--', lw=3.0, zorder=20,
    label=r'$\tau_3$, IceCube 8 yr')


# ################################################################################
# # IceCube, 15yr
# ################################################################################

# # Lower limit on tau2 (from nu2 --> nu1)
# filename_in = '../data/out/stat_analysis/ic_15yr/fixed_log10mnu1/limits_log10g21_vs_log10mnu1.dat'
# log10mnu1, \
#     log10g21_median, log10g21_sigma, \
#     log10g21_1sigma_lo, log10g21_1sigma_hi, \
#     log10g21_2sigma_lo, log10g21_2sigma_hi, \
#     log10g21_3sigma_lo, log10g21_3sigma_hi, \
#     log10g21_5sigma_lo, log10g21_5sigma_hi, \
#     log10g21_90CL_lo, log10g21_90CL_hi \
#     = np.loadtxt(filename_in, unpack=True)
# Gamma21 = [decay_rate_ij(   sqrt(10.**log10mnu1[k]+Dm21sq), # = mnu2
#                             10.**log10mnu1[k],
#                             10.**log10g21_2sigma_hi[k],
#                             10.**log10g21_2sigma_hi[k])
#             for k in range(len(log10mnu1))] # [s^{-1}]
# log10tau2 = [log10(1.0/x) for x in Gamma21] # [s]

# # Lower limit on tau3 (from nu3 --> nu1 and nu3 --> nu2)
# # nu3 --> nu1
# filename_in = '../../data/out/stat_analysis/ic_15yr/fixed_log10mnu1/limits_log10g31_vs_log10mnu1.dat'
# log10mnu1, \
#     log10g31_median, log10g31_sigma, \
#     log10g31_1sigma_lo, log10g31_1sigma_hi, \
#     log10g31_2sigma_lo, log10g31_2sigma_hi, \
#     log10g31_3sigma_lo, log10g31_3sigma_hi, \
#     log10g31_5sigma_lo, log10g31_5sigma_hi, \
#     log10g31_90CL_lo, log10g31_90CL_hi \
#     = np.loadtxt(filename_in, unpack=True)
# Gamma31 = [decay_rate_ij(   sqrt(10.**log10mnu1[k]+Dm31sq), # = mnu3
#                             10.**log10mnu1[k],
#                             10.**log10g31_2sigma_hi[k],
#                             10.**log10g31_2sigma_hi[k])
#             for k in range(len(log10mnu1))] # [s^{-1}]
# # nu3 --> nu2
# filename_in = '../../data/out/stat_analysis/ic_15yr/fixed_log10mnu1/limits_log10g32_vs_log10mnu1.dat'
# log10mnu1, \
#     log10g32_median, log10g32_sigma, \
#     log10g32_1sigma_lo, log10g32_1sigma_hi, \
#     log10g32_2sigma_lo, log10g32_2sigma_hi, \
#     log10g32_3sigma_lo, log10g32_3sigma_hi, \
#     log10g32_5sigma_lo, log10g32_5sigma_hi, \
#     log10g32_90CL_lo, log10g32_90CL_hi \
#     = np.loadtxt(filename_in, unpack=True)
# Gamma32 = [decay_rate_ij(   sqrt(10.**log10mnu1[k]+Dm31sq), # = mnu3
#                             10.**log10mnu1[k],
#                             10.**log10g32_2sigma_hi[k],
#                             10.**log10g32_2sigma_hi[k])
#             for k in range(len(log10mnu1))] # [s^{-1}]
# # Total
# Gamma3 = [Gamma31[k]+Gamma32[k] for k in range(len(Gamma31))] # [s^{-1}]
# log10tau3 = [log10(1.0/x) for x in Gamma3] # [s]

# window_length = 13 #25 
# polyorder = 1 #3
# log10tau2 = savgol_filter(log10tau2,
#     window_length, polyorder, deriv=0, mode='interp')

# window_length = 13 #25 
# polyorder = 1 #3
# log10tau3 = savgol_filter(log10tau3,
#     window_length, polyorder, deriv=0, mode='interp')

# ax.plot(log10mnu1, log10tau2,
#     c='gold', alpha=1.0, ls='-', lw=3.0, zorder=20,
#     label=r'$\tau_2$, IceCube 15 yr')

# ax.plot(log10mnu1, log10tau3,
#     c='gold', alpha=1.0, ls='--', lw=3.0, zorder=20,
#     label=r'$\tau_3$, IceCube 15 yr')


# ################################################################################
# # IceCube 15yr + IceCube-Gen2 10yr
# ################################################################################

# # Lower limit on tau2 (from nu2 --> nu1)
# filename_in = '../data/out/stat_analysis/ic_icg2/fixed_log10mnu1/limits_log10g21_vs_log10mnu1.dat'
# log10mnu1, \
#     log10g21_median, log10g21_sigma, \
#     log10g21_1sigma_lo, log10g21_1sigma_hi, \
#     log10g21_2sigma_lo, log10g21_2sigma_hi, \
#     log10g21_3sigma_lo, log10g21_3sigma_hi, \
#     log10g21_5sigma_lo, log10g21_5sigma_hi, \
#     log10g21_90CL_lo, log10g21_90CL_hi \
#     = np.loadtxt(filename_in, unpack=True)
# Gamma21 = [decay_rate_ij(   sqrt(10.**log10mnu1[k]+Dm21sq), # = mnu2
#                             10.**log10mnu1[k],
#                             10.**log10g21_2sigma_hi[k],
#                             10.**log10g21_2sigma_hi[k])
#             for k in range(len(log10mnu1))] # [s^{-1}]
# log10tau2 = [log10(1.0/x) for x in Gamma21] # [s]

# # Lower limit on tau3 (from nu3 --> nu1 and nu3 --> nu2)
# # nu3 --> nu1
# filename_in = '../data/out/stat_analysis/ic_icg2/fixed_log10mnu1/limits_log10g31_vs_log10mnu1.dat'
# log10mnu1, \
#     log10g31_median, log10g31_sigma, \
#     log10g31_1sigma_lo, log10g31_1sigma_hi, \
#     log10g31_2sigma_lo, log10g31_2sigma_hi, \
#     log10g31_3sigma_lo, log10g31_3sigma_hi, \
#     log10g31_5sigma_lo, log10g31_5sigma_hi, \
#     log10g31_90CL_lo, log10g31_90CL_hi \
#     = np.loadtxt(filename_in, unpack=True)
# Gamma31 = [decay_rate_ij(   sqrt(10.**log10mnu1[k]+Dm31sq), # = mnu3
#                             10.**log10mnu1[k],
#                             10.**log10g31_2sigma_hi[k],
#                             10.**log10g31_2sigma_hi[k])
#             for k in range(len(log10mnu1))] # [s^{-1}]
# # nu3 --> nu2
# filename_in = '../data/out/stat_analysis/ic_icg2/fixed_log10mnu1/limits_log10g32_vs_log10mnu1.dat'
# log10mnu1, \
#     log10g32_median, log10g32_sigma, \
#     log10g32_1sigma_lo, log10g32_1sigma_hi, \
#     log10g32_2sigma_lo, log10g32_2sigma_hi, \
#     log10g32_3sigma_lo, log10g32_3sigma_hi, \
#     log10g32_5sigma_lo, log10g32_5sigma_hi, \
#     log10g32_90CL_lo, log10g32_90CL_hi \
#     = np.loadtxt(filename_in, unpack=True)
# Gamma32 = [decay_rate_ij(   sqrt(10.**log10mnu1[k]+Dm31sq), # = mnu3
#                             10.**log10mnu1[k],
#                             10.**log10g32_2sigma_hi[k],
#                             10.**log10g32_2sigma_hi[k])
#             for k in range(len(log10mnu1))] # [s^{-1}]
# # Total
# Gamma3 = [Gamma31[k]+Gamma32[k] for k in range(len(Gamma31))] # [s^{-1}]
# log10tau3 = [log10(1.0/x) for x in Gamma3] # [s]

# window_length = 13 #25 
# polyorder = 1 #3
# log10tau2 = savgol_filter(log10tau2,
#     window_length, polyorder, deriv=0, mode='interp')

# window_length = 13 #25 
# polyorder = 1 #3
# log10tau3 = savgol_filter(log10tau3,
#     window_length, polyorder, deriv=0, mode='interp')

# ax.plot(log10mnu1, log10tau2,
#     c='C0', alpha=1.0, ls='-', lw=3.0, zorder=20,
#     label=r'$\tau_2$, IC 15 yr + IC-Gen2 10 yr')

# ax.plot(log10mnu1, log10tau3,
#     c='C0', alpha=1.0, ls='--', lw=3.0, zorder=20,
#     label=r'$\tau_3$, IC 15 yr + IC-Gen2 10 yr')


# ################################################################################
# # All experiments combined
# ################################################################################

# # Lower limit on tau2 (from nu2 --> nu1)
# filename_in = '../data/out/stat_analysis/comb/fixed_log10mnu1/limits_log10g21_vs_log10mnu1.dat'
# log10mnu1, \
#     log10g21_median, log10g21_sigma, \
#     log10g21_1sigma_lo, log10g21_1sigma_hi, \
#     log10g21_2sigma_lo, log10g21_2sigma_hi, \
#     log10g21_3sigma_lo, log10g21_3sigma_hi, \
#     log10g21_5sigma_lo, log10g21_5sigma_hi, \
#     log10g21_90CL_lo, log10g21_90CL_hi \
#     = np.loadtxt(filename_in, unpack=True)
# Gamma21 = [decay_rate_ij(   sqrt(10.**log10mnu1[k]+Dm21sq), # = mnu2
#                             10.**log10mnu1[k],
#                             10.**log10g21_2sigma_hi[k],
#                             10.**log10g21_2sigma_hi[k])
#             for k in range(len(log10mnu1))] # [s^{-1}]
# log10tau2 = [log10(1.0/x) for x in Gamma21] # [s]

# # Lower limit on tau3 (from nu3 --> nu1 and nu3 --> nu2)
# # nu3 --> nu1
# filename_in = '../data/out/stat_analysis/comb/fixed_log10mnu1/limits_log10g31_vs_log10mnu1.dat'
# log10mnu1, \
#     log10g31_median, log10g31_sigma, \
#     log10g31_1sigma_lo, log10g31_1sigma_hi, \
#     log10g31_2sigma_lo, log10g31_2sigma_hi, \
#     log10g31_3sigma_lo, log10g31_3sigma_hi, \
#     log10g31_5sigma_lo, log10g31_5sigma_hi, \
#     log10g31_90CL_lo, log10g31_90CL_hi \
#     = np.loadtxt(filename_in, unpack=True)
# Gamma31 = [decay_rate_ij(   sqrt(10.**log10mnu1[k]+Dm31sq), # = mnu3
#                             10.**log10mnu1[k],
#                             10.**log10g31_2sigma_hi[k],
#                             10.**log10g31_2sigma_hi[k])
#             for k in range(len(log10mnu1))] # [s^{-1}]
# # nu3 --> nu2
# filename_in = '../data/out/stat_analysis/comb/fixed_log10mnu1/limits_log10g32_vs_log10mnu1.dat'
# log10mnu1, \
#     log10g32_median, log10g32_sigma, \
#     log10g32_1sigma_lo, log10g32_1sigma_hi, \
#     log10g32_2sigma_lo, log10g32_2sigma_hi, \
#     log10g32_3sigma_lo, log10g32_3sigma_hi, \
#     log10g32_5sigma_lo, log10g32_5sigma_hi, \
#     log10g32_90CL_lo, log10g32_90CL_hi \
#     = np.loadtxt(filename_in, unpack=True)
# Gamma32 = [decay_rate_ij(   sqrt(10.**log10mnu1[k]+Dm31sq), # = mnu3
#                             10.**log10mnu1[k],
#                             10.**log10g32_2sigma_hi[k],
#                             10.**log10g32_2sigma_hi[k])
#             for k in range(len(log10mnu1))] # [s^{-1}]
# # Total
# Gamma3 = [Gamma31[k]+Gamma32[k] for k in range(len(Gamma31))] # [s^{-1}]
# log10tau3 = [log10(1.0/x) for x in Gamma3] # [s]

# window_length = 13 #25 
# polyorder = 1 #3
# log10tau2 = savgol_filter(log10tau2,
#     window_length, polyorder, deriv=0, mode='interp')

# window_length = 13 #25 
# polyorder = 1 #3
# log10tau3 = savgol_filter(log10tau3,
#     window_length, polyorder, deriv=0, mode='interp')

# ax.plot(log10mnu1, log10tau2,
#     c='C1', alpha=1.0, ls='-', lw=3.0, zorder=20,
#     label=r'$\tau_2$, all future experiments combined')

# ax.plot(log10mnu1, log10tau3,
#     c='C1', alpha=1.0, ls='--', lw=3.0, zorder=20,
#     label=r'$\tau_3$, all future experiments combined')


####################################################################
# Formatting
####################################################################

ax.set_xlabel(r'Lightest neutrino mass, $\log_{10}(m_{\nu,1}/{\rm eV})$',
    fontsize=25)
ax.set_ylabel(r'Lower limit on $\nu$ lifetime (95\% C.L.), $\log_{10}(\tau_i/{\rm s})$',
    fontsize=25)

# ax.set_xscale('log')
# ax.set_yscale('log')

ax.tick_params('both', length=10, width=2, which='major')
ax.tick_params('both', length=5, width=1, which='minor')
ax.tick_params(axis='both', which='major', pad=10, direction='in')
ax.tick_params(axis='both', which='minor', pad=10, direction='in')
ax.tick_params(axis='x', which='minor', bottom=True, top=True)
ax.tick_params(axis='y', which='minor', left=True, right=True)
ax.tick_params(bottom=True, top=True, left=True, right=True)

ax.xaxis.set_major_locator(MultipleLocator(1.0))
ax.xaxis.set_minor_locator(MultipleLocator(0.2))
ax.yaxis.set_major_locator(MultipleLocator(1.0))
ax.yaxis.set_minor_locator(MultipleLocator(0.2))

ax.set_xlim([x_min, x_max])
# ax.set_ylim([y_min, y_max])

ax.xaxis.set_zorder(100)
ax.yaxis.set_zorder(100)

leg = ax.legend(loc='center right', ncol=1)
plt.show()
pylab.savefig('limits_log10tau_vs_log10mnu1.pdf',
    bbox_inches='tight', dpi=300)