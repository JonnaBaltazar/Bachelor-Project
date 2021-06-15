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
filename = '../../data/Visible_coupling_g32.txt'

dataf = open(filename, "r")
l = [float(x) for x in dataf.readline().split()]
g21 = l[0]
g31 = l[1]
g32 = l[2]
x = l[3]
m1 = l[4]
n_gamma = int(l[5])

E_f, fe_SM, fmu_SM, fe_1, fmu_1, fe_2, fmu_2, fe_3, fmu_3, fe_4, fmu_4 = np.loadtxt(filename, unpack=True, dtype='float', skiprows=2)

ftau_SM = 1 - fe_SM - fmu_SM

ftau_1 = 1 - fe_1 - fmu_1
ftau_2 = 1 - fe_2 - fmu_2
ftau_3 = 1 - fe_3 - fmu_3
ftau_4 = 1 - fe_4 - fmu_4 

print("Parameters")
print("g21: ", g21)
print("g31: ", g31)
print("g32: ", g32)
print("z:   ", x)
print("m1:  ", m1)

v = list(plt.axis())
v[0] = E_f[0]
v[1] = E_f[-1]
v[2] = 0
v[3] = 35

fc = "lightblue"
alpha = 0.4

fig, ((ax1, ax2), (ax3, ax4))= plt.subplots(nrows=2, ncols=2, figsize=(10, 6))
ax1.plot(E_f, fe_1, label=r"$\nu_e$", color="royalblue")
ax1.plot(E_f, fmu_1, label=r"$\nu_\mu$", color="firebrick")
ax1.plot(E_f, ftau_1, label=r"$\nu_\tau$", color="forestgreen")

ax1.plot(E_f, fe_SM, linestyle="dashed", color="royalblue")
ax1.plot(E_f, fmu_SM, linestyle="dashed", color="firebrick")
ax1.plot(E_f, ftau_SM,  linestyle="dashed", color="forestgreen")
ax1.hlines(-1, 10**6, 5*10**6, colors = 'k', linestyle = 'dashed', label="SM")

ax1.set(ylim=(0.2,0.5), xscale="log", title=r"$g^{(')}_{32}=10^{-10}$")
ax1.set_ylabel(r"$f_{\alpha, \oplus}$", fontsize=15)
ax1.set_xlabel(r"$E_f{\rm\ [GeV]}$", fontsize=13)
ax1.legend(loc="upper left", fontsize=9)

ax1.fill_between([6*1e4, 3*1e6], [v[2], v[2]], [v[3], v[3]], facecolor = fc, edgecolor = None, alpha = alpha, zorder = -10)
ax1.fill_between(v[:2], [v[3] + 1, v[3] + 1], [v[3] + 2, v[3] + 2], facecolor = fc, edgecolor = "k", lw = 0.5, alpha = alpha)

ax2.plot(E_f, fe_2, label=r"$\nu_e$", color="royalblue")
ax2.plot(E_f, fmu_2, label=r"$\nu_\mu$", color="firebrick")
ax2.plot(E_f, ftau_2, label=r"$\nu_\tau$", color="forestgreen")

ax2.plot(E_f, fe_SM, linestyle="dashed", color="royalblue")
ax2.plot(E_f, fmu_SM, linestyle="dashed", color="firebrick")
ax2.plot(E_f, ftau_SM,  linestyle="dashed", color="forestgreen")
ax2.hlines(-1, 10**6, 5*10**6, colors = 'k', linestyle = 'dashed', label="SM")

ax2.set(ylim=(0.2,0.5), xscale="log", title=r"$g^{(')}_{32}=10^{-6}$")
ax2.set_ylabel(r"$f_{\alpha, \oplus}$", fontsize=15)
ax2.set_xlabel(r"$E_f{\rm\ [GeV]}$", fontsize=13)
ax2.legend(loc="upper left", fontsize=9)


ax2.fill_between([6*1e4, 3*1e6], [v[2], v[2]], [v[3], v[3]], facecolor = fc, edgecolor = None, alpha = alpha, zorder = -10)
ax2.fill_between(v[:2], [v[3] + 1, v[3] + 1], [v[3] + 2, v[3] + 2], facecolor = fc, edgecolor = "k", lw = 0.5, alpha = alpha)

ax3.plot(E_f, fe_3, label=r"$\nu_e$", color="royalblue")
ax3.plot(E_f, fmu_3, label=r"$\nu_\mu$", color="firebrick")
ax3.plot(E_f, ftau_3, label=r"$\nu_\tau$", color="forestgreen")

ax3.plot(E_f, fe_SM, linestyle="dashed", color="royalblue")
ax3.plot(E_f, fmu_SM, linestyle="dashed", color="firebrick")
ax3.plot(E_f, ftau_SM,  linestyle="dashed", color="forestgreen")
ax3.hlines(-1, 10**6, 5*10**6, colors = 'k', linestyle = 'dashed', label="SM")

ax3.set(ylim=(0.2,0.5), xscale="log", title=r"$g^{(')}_{32}=10^{-5.5}$")
ax3.set_ylabel(r"$f_{\alpha, \oplus}$", fontsize=15)
ax3.set_xlabel(r"$E_f{\rm\ [GeV]}$", fontsize=13)
ax3.legend(loc="upper left", fontsize=9)

ax3.fill_between([6*1e4, 3*1e6], [v[2], v[2]], [v[3], v[3]], facecolor = fc, edgecolor = None, alpha = alpha, zorder = -10)
ax3.fill_between(v[:2], [v[3] + 1, v[3] + 1], [v[3] + 2, v[3] + 2], facecolor = fc, edgecolor = "k", lw = 0.5, alpha = alpha)

ax4.plot(E_f, fe_4, label=r"$\nu_e$", color="royalblue")
ax4.plot(E_f, fmu_4, label=r"$\nu_\mu$", color="firebrick")
ax4.plot(E_f, ftau_4, label=r"$\nu_\tau$", color="forestgreen")

ax4.plot(E_f, fe_SM, linestyle="dashed", color="royalblue")
ax4.plot(E_f, fmu_SM, linestyle="dashed", color="firebrick")
ax4.plot(E_f, ftau_SM,  linestyle="dashed", color="forestgreen")
ax4.hlines(-1, 10**6, 5*10**6, colors = 'k', linestyle = 'dashed', label="SM")

ax4.set(ylim=(0.2,0.5), xscale="log", title=r"$g^{(')}_{32}=10^{-2.5}$")
ax4.set_ylabel(r"$f_{\alpha, \oplus}$", fontsize=15)
ax4.set_xlabel(r"$E_f{\rm\ [GeV]}$", fontsize=13)
ax4.legend(loc="upper left", fontsize=9)

ax4.fill_between([6*1e4, 3*1e6], [v[2], v[2]], [v[3], v[3]], facecolor = fc, edgecolor = None, alpha = alpha, zorder = -10)
ax4.fill_between(v[:2], [v[3] + 1, v[3] + 1], [v[3] + 2, v[3] + 2], facecolor = fc, edgecolor = "k", lw = 0.5, alpha = alpha)


fig.tight_layout()
fig.savefig('../../fig/fig_Jonathan/Visible_coupling_g32.pdf')
plt.show()