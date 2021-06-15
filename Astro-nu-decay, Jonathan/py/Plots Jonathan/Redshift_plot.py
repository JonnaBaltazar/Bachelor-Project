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

filename = '../../data/Visible_Redshift.txt'

dataf = open(filename, "r")
l = [float(x) for x in dataf.readline().split()]
a = l[0]
b = l[1]
g31 = l[2]
m1 = l[3]
gamma = l[4]

E_f, P_SM, P_d5, P_d1, P_d15, P_ykbh = np.loadtxt(filename, unpack=True, dtype='float', skiprows=1)

pivot = 42
print("Pivot is at E =", E_f[pivot])

v = list(plt.axis())
v[0] = E_f[0]
v[1] = E_f[-1]
v[2] = 0.0
v[3] = 2.00

fig, ax = plt.subplots(figsize=(12,6))
ax.plot(E_f, P_SM/P_SM[pivot], label=r"SM", linestyle="dashed", color="k")
ax.plot(E_f, P_d5/P_d5[pivot], label=r"$\delta(z-0.5)$", color="forestgreen")
ax.plot(E_f, P_d1/P_d1[pivot], label=r"$\delta(z-1)$", color="orange")
ax.plot(E_f, P_d15/P_d15[pivot], label=r"$\delta(z-1.5)$", color="darkorchid")
ax.plot(E_f, P_ykbh/P_ykbh[pivot], label="YKBH", color="crimson", linewidth=2)

ax.hlines(1.4, 400, 6*1e4, colors = 'grey', linestyle = 'solid', linewidth = 1)
ax.vlines(6*1e4, 1.4, 1.6, colors = 'grey', linestyle = 'solid', linewidth = 1)

fc = "lightblue"
alpha = 0.4
ax.fill_between([6*1e4, 3*1e6], [v[2], v[2]], [v[3], v[3]], facecolor = fc, edgecolor = None, alpha = alpha, zorder = -10)
ax.fill_between(v[:2], [v[3] + 1, v[3] + 1], [v[3] + 2, v[3] + 2], facecolor = fc, edgecolor = "k", lw = 0.5, alpha = alpha)


ax.set(ylim=(0.7,1.6), xlim=(0.4*10**3, 2*10**9), xlabel=r"$E_f$", ylabel=r"$\bar{P}^{vis}_{\mu e}$ (Rescaled)")

ax.set(xscale="log")
ax.legend(loc="center left")

fc = "lightblue"

s = "Transition probabillity"
ax.text(0.02, 0.97, s, weight='heavy', fontsize=20, transform=ax.transAxes, 
        verticalalignment='top', color='k')

ax.text(0.02, 0.9, "Benchmark parameters:", weight='heavy', fontsize=15, transform=ax.transAxes, 
        verticalalignment='top', color='k')

d = r"$m_1=0.1 \mathrm{eV}$,  $\gamma=2$, $g^{(')}_{ij}=10^{-6}$"
ax.text(0.02, 0.85, d, fontsize=15, transform=ax.transAxes, 
        verticalalignment='top', color='k')

plt.savefig("../../fig/fig_Jonathan/Redshift_plot.pdf")
plt.show()