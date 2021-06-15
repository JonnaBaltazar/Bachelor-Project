import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import LogLocator


flavor_texs = ["e", r"\mu", r"\tau"]
colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
fs = 15

filename_0 = '../../data/nu_flux_ratios_block_0_g21.dat'
filename_1 = '../../data/nu_flux_ratios_block_1_g21.dat'
filename_2 = '../../data/nu_flux_ratios_block_2_g21.dat'

fe_0, fmu_0 = np.loadtxt(filename_0, unpack=True, dtype='float', skiprows=0)
fe_1, fmu_1 = np.loadtxt(filename_1, unpack=True, dtype='float', skiprows=0)
fe_2, fmu_2 = np.loadtxt(filename_2, unpack=True, dtype='float', skiprows=0)

fe   = np.concatenate((fe_0, fe_1, fe_2))
fmu  = np.concatenate((fmu_0, fmu_1, fmu_2))

ftau = 1 - fe - fmu

filename_SM = "../../data/Flavor_Ef.txt"

_, _, _, _, fe_SM, fmu_SM, ftau_SM = np.loadtxt(filename_SM, unpack=True, dtype='float', skiprows=0)

f_tot_SM = fe_SM+fmu_SM+ftau_SM
fe_SM = fe_SM/f_tot_SM
fmu_SM = fmu_SM/f_tot_SM
ftau_SM = ftau_SM/f_tot_SM

g21 = np.linspace(-8, 0, 150)

fig, ax = plt.subplots(figsize=(12,8))

ax.plot(g21, fe, label=r"$\nu_e$", color="royalblue", zorder=2, linewidth=2.1)
ax.plot(g21, fmu, label=r"$\nu_\mu$", color="firebrick", zorder=2, linewidth=2.1)
ax.plot(g21, ftau, label=r"$\nu_\tau$", color="forestgreen", zorder=2, linewidth=2.1)

ax.plot([g21[0], g21[-1]], [fe_SM[0], fe_SM[0]], linestyle="dashed", color="royalblue", zorder=1, linewidth=2.1)
ax.plot([g21[0], g21[-1]], [fmu_SM[0], fmu_SM[0]], linestyle="dashed", color="firebrick", zorder=1, linewidth=2.1)
ax.plot([g21[0], g21[-1]], [ftau_SM[0], ftau_SM[0]],  linestyle="dashed", color="forestgreen", zorder=1, linewidth=2.1)

ax.hlines(-1, -8, 0, colors = 'k', linestyle = 'dashed', label="SM")
ax.hlines(0.36, -8.5, -3.5, colors = 'grey', linestyle = 'solid', linewidth = 1)
ax.vlines(-3.5, 0.36, 0.38, colors = 'grey', linestyle = 'solid', linewidth = 1)

plt.ylim(0.296, 0.38)
plt.xlim(-8.1,0.1)
plt.xticks(fontsize=13)
plt.yticks(fontsize=13)

plt.xlabel(r"$\log(g_{21})$", fontsize=16)
plt.ylabel(r"$\langle f_{\alpha, \oplus} \rangle$", fontsize=27)

# IceCube ROI
fc = "lightblue"
alpha = 0.4


s = "Final Ratio from (1,2,0)"
ax.text(0.02, 0.97, s, weight='heavy', fontsize=20, transform=ax.transAxes, 
        verticalalignment='top', color='k')

ax.text(0.02, 0.9, "Benchmark parameters:", weight='heavy', fontsize=15, transform=ax.transAxes, 
        verticalalignment='top', color='k')

d = r"$m_1=0.1\mathrm{eV}, z=1,  \gamma=2$ and $g_{31}=g_{32}=10^{-10}$"

ax.text(0.02, 0.85, d, fontsize=15, transform=ax.transAxes, 
        verticalalignment='top', color='k')
plt.legend(loc="upper right", fontsize=20)
plt.savefig("../../fig/fig_Jonathan/g21_plot.pdf")
plt.show()