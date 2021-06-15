import os
from numpy import *
import numpy as np
from pylab import *
from matplotlib import *
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.use('macOsX')
import matplotlib.patheffects as path_effects
import ternary

name = "icecube_8yr"
output_format = 'png'


filename = '../../data/'+name+'.txt'

fe, fmu, chi2 = np.loadtxt(filename, unpack=True, dtype='float', skiprows=0)
ftau = 1 - fe - fmu

LH = -2*np.log(chi2+0.1) 
print(LH)

mpl.rcParams['xtick.labelsize']=23
mpl.rcParams['ytick.labelsize']=23
mpl.rcParams['legend.fontsize']=23
mpl.rcParams['legend.borderpad']=0.4
mpl.rcParams['axes.labelpad']=10
mpl.rcParams['ps.fonttype']=42
mpl.rcParams['pdf.fonttype']=42

colors = ['#2a4d69', '#4b86b4', '#adcbe3', '#D3E4F4', '#63ace5']

fig = ternary.plt.figure(figsize=(11, 6))

fig, tax = ternary.figure(scale=1.0)
tax.ax.axis("off")
fig.set_facecolor('w')

fontsize = 15
label_left = r'Fraction of $\nu_\tau$'
label_right = r'Fraction of $\nu_\mu$'
label_bottom = r'Fraction of $\nu_e$'

value_min = np.min(LH)
value_max = np.max(LH)
print(value_min, value_max)
norm = matplotlib.colors.Normalize(vmin=value_min, vmax=value_max)
cmap = matplotlib.cm.get_cmap('Blues')

print(fe[np.argmax(fmu)], fmu[np.argmax(fmu)])
print(1-fe[np.argmax(fmu)]-fmu[np.argmax(fmu)])
print(max(fe), max(fmu), max(ftau), len(LH))

tax.boundary(linewidth=1.0)
tax.gridlines(color="gray", multiple=0.1, linewidth=0.5, ls='-', alpha=0.5)
tax.left_axis_label(label_left, fontsize=fontsize, offset=0.16)
tax.right_axis_label(label_right, fontsize=fontsize, offset=0.16)
tax.bottom_axis_label(label_bottom, fontsize=fontsize, offset=0.08)
for i in range(len(fe)):
        color = cmap(norm(LH[i]))
        tax.scatter([[fe[i], fmu[i], ftau[i]]], \
            marker='o', color=color, edgecolor=color,
            s=8.0, linewidths=0.5, zorder=-1, alpha=1)
tick_thing = np.linspace(min(LH), max(LH), 6)
ax1 = fig.add_axes([0.08, 0.87, 0.25, 0.03])
cb1 = mpl.colorbar.ColorbarBase(ax1, cmap=cmap, norm=norm, \
    orientation='horizontal',
    ticks=[tick_thing[0], tick_thing[1], tick_thing[2], tick_thing[3], tick_thing[4], tick_thing[5]])
cb1.set_label(r'Value', size=10, labelpad=-40)
cb1.ax.tick_params(labelsize=8)
print(tick_thing)
cb1.ax.set_xticklabels([r'-14', r'-10', r'-7', r'-3', r'1', r'5'])

tax.clear_matplotlib_ticks()
tax._redraw_labels()
tax.legend(fontsize=12)
ternary.plt.tight_layout()
ternary.plt.savefig('../../fig/fig_Jonathan/'+name+'_Likelihood_function.png', dpi=500)
print("Saved plot!")