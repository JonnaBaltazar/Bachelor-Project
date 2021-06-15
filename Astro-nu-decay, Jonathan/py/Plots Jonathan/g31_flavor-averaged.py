import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import LogLocator


flavor_texs = ["e", r"\mu", r"\tau"]
colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
fs = 15

filename_0 = '../../data/nu_flux_ratios_block_0_g31.dat'
filename_1 = '../../data/nu_flux_ratios_block_1_g31.dat'
filename_2 = '../../data/nu_flux_ratios_block_2_g31.dat'


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

g32 = np.linspace(-8, 0, 150)

fig, ax = plt.subplots(figsize=(12,8))

ax.plot(g32, fe, label=r"$\nu_e$", color="royalblue", zorder=2)
ax.plot(g32, fmu, label=r"$\nu_\mu$", color="firebrick", zorder=2)
ax.plot(g32, ftau, label=r"$\nu_\tau$", color="forestgreen", zorder=2)

ax.plot([g32[0], g32[-1]], [fe_SM[0], fe_SM[0]], linestyle="dashed", color="royalblue", zorder=1)
ax.plot([g32[0], g32[-1]], [fmu_SM[0], fmu_SM[0]], linestyle="dashed", color="firebrick", zorder=1)
ax.plot([g32[0], g32[-1]], [ftau_SM[0], ftau_SM[0]],  linestyle="dashed", color="forestgreen", zorder=1)

ax.hlines(-1, -8, 0, colors = 'k', linestyle = 'dashed', label="SM")
ax.hlines(0.52, -8.5, -2.5, colors = 'grey', linestyle = 'solid', linewidth = 1)
ax.vlines(-2.5, 0.52, 0.6, colors = 'grey', linestyle = 'solid', linewidth = 1)

plt.ylim(0.2, 0.6)
plt.xlim(-8,0)

plt.xticks(fontsize=16)
plt.yticks(fontsize=16)

plt.xlabel(r"$\log(g_{32})$", fontsize=16)
plt.ylabel(r"$\langle f_{\alpha \oplus} \rangle$", fontsize=25)

fc = "lightblue"
alpha = 0.4

s = r"Final Ratio from $\frac{1}{3}(1:2:0).$"
ax.text(0.02, 0.97, s, weight='heavy', fontsize=25, transform=ax.transAxes, 
        verticalalignment='top', color='k')

d = r"$m_1=0.1\mathrm{eV}, z=1,  \gamma=2$ and $g_{21}=g_{31}=10^{-10}.$"
ax.text(0.02, 0.88, d, fontsize=20, transform=ax.transAxes, 
        verticalalignment='top', color='k')

plt.legend(loc="upper right", fontsize=19)
plt.savefig("../../fig/fig_Jonathan/g31_plot.pdf")
plt.show()