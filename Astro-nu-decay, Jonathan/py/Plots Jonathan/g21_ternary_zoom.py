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

name = "Wrapper_Zoom"
output_format = 'png'


filename_0 = '../../data/extra_data/nu_flux_ratios_block_0_g21.dat'
filename_1 = '../../data/extra_data/nu_flux_ratios_block_1_g21.dat'
filename_2 = '../../data/extra_data/nu_flux_ratios_block_2_g21.dat'

fe_0, fmu_0 = np.loadtxt(filename_0, unpack=True, dtype='float', skiprows=0)
fe_1, fmu_1 = np.loadtxt(filename_1, unpack=True, dtype='float', skiprows=0)
fe_2, fmu_2 = np.loadtxt(filename_2, unpack=True, dtype='float', skiprows=0)

fe = np.concatenate((fe_0, fe_1, fe_2))
fmu = np.concatenate((fmu_0, fmu_1, fmu_2))

ftau = 1 - fe - fmu

mpl.rcParams['xtick.labelsize']=23
mpl.rcParams['ytick.labelsize']=23
mpl.rcParams['legend.fontsize']=23
mpl.rcParams['legend.borderpad']=0.4
mpl.rcParams['axes.labelpad']=10
mpl.rcParams['ps.fonttype']=42
mpl.rcParams['pdf.fonttype']=42

colors = ['#2a4d69', '#4b86b4', '#adcbe3', '#D3E4F4', '#63ace5']

fig = ternary.plt.figure(figsize=(11, 6))
ax1 = fig.add_subplot(2, 1, 1)
ax2 = fig.add_subplot(2, 1, 2)

tax1 = ternary.TernaryAxesSubplot(ax=ax1, scale=1)

tax1.boundary(linewidth=1.0)
tax1.gridlines(color="gray", multiple=0.1, linewidth=0.5, ls='-', alpha=0.5)
tax1.clear_matplotlib_ticks()
tax1._redraw_labels()
tax1.ax.axis("off")
tax1.ticks(axis='blr', linewidth=0.5, multiple=0.1, offset=0.022, \
    clockwise=False, tick_formats="%.1f")


tax2 = ternary.TernaryAxesSubplot(ax=ax2,scale=30)
axes_colors = {'b': 'r', 'r': 'r', 'l': 'r'}
tax2.boundary(linewidth=1.0, axes_colors=axes_colors)
tax2.gridlines(color="r", multiple=5, linewidth=1, ls='-')
tax2.ax.axis("equal")
tax2.ax.axis("off")

fontsize = 16
label_left = r'Fraction of $\nu_\tau$'
label_right = r'Fraction of $\nu_\mu$'
label_bottom = r'Fraction of $\nu_e$'

tax1.left_axis_label(label_left, fontsize=fontsize, offset=0.16)
tax1.right_axis_label(label_right, fontsize=fontsize, offset=0.16)
tax1.bottom_axis_label(label_bottom, fontsize=fontsize, offset=0.08)
tax1.set_title("Entire region",color='k', fontsize=20)

tax2.left_axis_label(label_left, fontsize=fontsize, offset=0.16, color='r')
tax2.right_axis_label(label_right, fontsize=fontsize, offset=0.16, color='r')
tax2.bottom_axis_label(label_bottom, fontsize=fontsize, offset=0.08, color='r')

#tax2.set_axis_limits({'b': [0.3, 0.4], 'l': [0.3, 0.4], 'r': [0.3, 0.4]})

tax2.set_axis_limits({'b': [0.3, 0.37], 'l': [0.315, 0.385], 'r': [0.315, 0.385]})
tax2.get_ticks_from_axis_limits(multiple=10)
tax2.set_title("Zoomed region",color='r', fontsize=20)



tick_formats = "%.3f"
tax2.set_custom_ticks(fontsize=10, offset=0.025, multiple=10, alpha=0.8,
                      axes_colors=axes_colors, tick_formats=tick_formats)

ternary.plt.tight_layout()

cmap = matplotlib.cm.get_cmap('autumn')
value_min = -8
value_max = 0
norm = matplotlib.colors.Normalize( vmin=value_min, 
                                    vmax=value_max)

g_21 = np.linspace(0, 1, len(fe));
color = cmap(g_21)

tax1.scatter([[fe, fmu, ftau]], \
    marker='o', 
    color=color,
    edgecolor=color,
    s=25, 
    label = "(1,2,0)",
    linewidths=0.5, 
    zorder=4, 
    alpha=1.)

# tax1.scatter([[fe_g31, fmu_g31, ftau_g31]], \
#     marker='o', 
#     color=color,
#     edgecolor=color,
#     s=25, 
#     label = "(1,2,0)",
#     linewidths=0.5, 
#     zorder=4, 
#     alpha=1.)

ax3 = fig.add_axes([0.36, 0.83, 0.25, 0.03])
ticks_linspace = np.linspace(min(fe), max(fe), 5)

cb1 = mpl.colorbar.ColorbarBase(ax3, cmap=cmap, norm=norm, \
    orientation='horizontal',
    ticks=[-8, -6, -4, -2, 0])

cb1.set_label(r"Value of $\log_{10}(g_{21}^{(')})$", size=12, labelpad=-50)
cb1.ax.tick_params(labelsize=8)
cb1.ax.set_xticklabels([r'-8', r'-6', r'-4', r'-2', r'0'])

points = [(fe, fmu, ftau)]
points_c = tax2.convert_coordinates(points, axisorder='brl')
tax2.scatter(points_c, marker='o', s=25, c=color, edgecolor=color, zorder=4)


# draw the zoom region on the first plot
tax1.line((0.3, 0.3, 0.4), (0.4, 0.3, 0.3), color='r', lw=2.0, zorder=4)
tax1.line((0.4, 0.3, 0.3), (0.3, 0.4, 0.3), color='r', lw=2.0, zorder=4)
tax1.line((0.3, 0.4, 0.3), (0.3, 0.3, 0.4), color='r', lw=2.0, zorder=4)

fig.set_facecolor("w")

tax1.ax.set_position([0.01, 0.05, 0.46, 0.8])
tax2.ax.set_position([0.50, 0.05, 0.46, 0.8])

tax1.resize_drawing_canvas()
tax2.resize_drawing_canvas()

tax1.clear_matplotlib_ticks()
tax2.clear_matplotlib_ticks()

tax1._redraw_labels()
tax2._redraw_labels()

ternary.plt.savefig("../../fig/fig_Jonathan/Wrapper_Zoom.png", dpi=300)
ternary.plt.show()

