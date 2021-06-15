import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
import os
from pylab import *
from matplotlib import *
import matplotlib.patheffects as path_effects
from matplotlib.ticker import LogLocator


flavor_texs = ["e", r"\mu", r"\tau"]
colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
fs = 15
coupling = r'$g_{21} = 10^{-6}, m_1=10\mathrm{eV}$'

filename_g21 = "../../data/Mass_g21.txt"
filename_g31 = "../../data/Mass_g31.txt"
filename_g32 = "../../data/Mass_g32.txt"
filename_All = "../../data/Mass_All.txt"

m1, fe_g21, fmu_g21, ftau_g21, fe_g21_SM, fmu_g21_SM, ftau_g21_SM = np.loadtxt(filename_g21, unpack=True, dtype='float', skiprows=0)

_,  fe_g31, fmu_g31, ftau_g31, fe_g31_SM, fmu_g31_SM, ftau_g31_SM = np.loadtxt(filename_g31, unpack=True, dtype='float', skiprows=0)

_,  fe_g32, fmu_g32, ftau_g32, fe_g32_SM, fmu_g32_SM, ftau_g32_SM = np.loadtxt(filename_g32, unpack=True, dtype='float', skiprows=0)

_,  fe_All, fmu_All, ftau_All, fe_All_SM, fmu_All_SM, ftau_All_SM = np.loadtxt(filename_All, unpack=True, dtype='float', skiprows=0)

v = list(plt.axis())
v[0] = m1[0]
v[1] = m1[-1]
v[2] = 0
v[3] = 35

fc = "lightblue"
alpha = 0.4

fig, ((ax1, ax2), (ax3, ax4))= plt.subplots(nrows=2, ncols=2, figsize=(10, 6))
ax1.plot(m1, fe_g21, label=r"$\nu_e$", color="royalblue")
ax1.plot(m1, fmu_g21, label=r"$\nu_\mu$", color="firebrick")
ax1.plot(m1, ftau_g21, label=r"$\nu_\tau$", color="forestgreen")

ax1.plot(m1, fe_g21_SM, linestyle="dashed",  color="royalblue", alpha=0.7)
ax1.plot(m1, fmu_g21_SM, linestyle="dashed", color="firebrick", alpha=0.7)
ax1.plot(m1, ftau_g21_SM,  linestyle="dashed", color="forestgreen", alpha=0.7)
ax1.hlines(10, 0, 5, colors = 'k', linestyle = 'dashed', label="SM")

ax1.set(ylim=(0.11,0.64), xlim=(10**(-5),4), xscale="log",title=r"$g^{(')}_{21}=10^{-6}$")
ax1.set_ylabel(r"$f_{\alpha, \oplus}$", fontsize=15)
ax1.set_xlabel(r"$\log (m_{\nu_1})$", fontsize=13)
ax1.legend(loc="upper left", fontsize=9)
ax1.fill_between([0,0.0281], [v[2], v[2]], [v[3], v[3]], facecolor = fc, edgecolor = None, alpha = alpha, zorder = -10)
ax1.fill_between(v[:2], [v[3] + 1, v[3] + 1], [v[3] + 2, v[3] + 2], facecolor = fc, edgecolor = "k", lw = 0.5, alpha = alpha)

ax2.plot(m1, fe_g31, label=r"$\nu_e$", color="royalblue")
ax2.plot(m1, fmu_g31, label=r"$\nu_\mu$", color="firebrick")
ax2.plot(m1, ftau_g31, label=r"$\nu_\tau$", color="forestgreen")

ax2.plot(m1, fe_g31_SM, linestyle="dashed", color="royalblue", alpha=0.7)
ax2.plot(m1, fmu_g31_SM, linestyle="dashed", color="firebrick", alpha=0.7)
ax2.plot(m1, ftau_g31_SM,  linestyle="dashed", color="forestgreen", alpha=0.7)
ax2.hlines(10, 0, 5, colors = 'k', linestyle = 'dashed', label="SM")

ax2.set(ylim=(0.11,0.64), xlim=(10**(-5),4), xscale="log", title=r"$g^{(')}_{31}=10^{-6}$")
ax2.set_ylabel(r"$f_{\alpha, \oplus}$", fontsize=15)
ax2.set_xlabel(r"$\log (m_{\nu_1})$", fontsize=13)
ax2.legend(loc="upper left", fontsize=9)
ax2.fill_between([0,0.0281], [v[2], v[2]], [v[3], v[3]], facecolor = fc, edgecolor = None, alpha = alpha, zorder = -10)
ax2.fill_between(v[:2], [v[3] + 1, v[3] + 1], [v[3] + 2, v[3] + 2], facecolor = fc, edgecolor = "k", lw = 0.5, alpha = alpha)

ax3.plot(m1, fe_g32, label=r"$\nu_e$", color="royalblue")
ax3.plot(m1, fmu_g32, label=r"$\nu_\mu$", color="firebrick")
ax3.plot(m1, ftau_g32, label=r"$\nu_\tau$", color="forestgreen")

ax3.plot(m1, fe_g32_SM, linestyle="dashed", color="royalblue", alpha=0.7)
ax3.plot(m1, fmu_g32_SM, linestyle="dashed", color="firebrick", alpha=0.7)
ax3.plot(m1, ftau_g32_SM,  linestyle="dashed", color="forestgreen", alpha=0.7)
ax3.hlines(10, 0, 5, colors = 'k', linestyle = 'dashed', label="SM")

ax3.set(ylim=(0.11,0.64), xlim=(10**(-5),4),xscale="log", title=r"$g^{(')}_{32}=10^{-6}$")
ax3.set_ylabel(r"$f_{\alpha, \oplus}$", fontsize=15)
ax3.set_xlabel(r"$\log (m_{\nu_1})$", fontsize=13)
ax3.legend(loc="upper left", fontsize=9)
ax3.fill_between([0,0.0281], [v[2], v[2]], [v[3], v[3]], facecolor = fc, edgecolor = None, alpha = alpha, zorder = -10)
ax3.fill_between(v[:2], [v[3] + 1, v[3] + 1], [v[3] + 2, v[3] + 2], facecolor = fc, edgecolor = "k", lw = 0.5, alpha = alpha)

ax4.plot(m1, fe_All, label=r"$\nu_e$", color="royalblue")
ax4.plot(m1, fmu_All, label=r"$\nu_\mu$", color="firebrick")
ax4.plot(m1, ftau_All, label=r"$\nu_\tau$", color="forestgreen")

ax4.plot(m1, fe_All_SM, linestyle="dashed", color="royalblue", alpha=0.7)
ax4.plot(m1, fmu_All_SM, linestyle="dashed", color="firebrick", alpha=0.7)
ax4.plot(m1, ftau_All_SM,  linestyle="dashed", color="forestgreen", alpha=0.7)
ax4.hlines(10, 0, 5, colors = 'k', linestyle = 'dashed', label="SM")

ax4.set(ylim=(0.11,0.64), xlim=(10**(-5),4), xscale="log", title=r"$g^{(')}_{21}=g^{(')}_{31}=g^{(')}_{32}=10^{-6}$")
ax4.set_ylabel(r"$f_{\alpha, \oplus}$", fontsize=15)
ax4.set_xlabel(r"$\log (m_{\nu_1})$", fontsize=13)
ax4.legend(loc="upper left", fontsize=9)
ax4.fill_between([0,0.0281], [v[2], v[2]], [v[3], v[3]], facecolor = fc, edgecolor = None, alpha = alpha, zorder = -10)
ax4.fill_between(v[:2], [v[3] + 1, v[3] + 1], [v[3] + 2, v[3] + 2], facecolor = fc, edgecolor = "k", lw = 0.5, alpha = alpha)

fig.tight_layout()
fig.savefig('../../fig/fig_Jonathan/Mass_plot.pdf')
plt.show()
