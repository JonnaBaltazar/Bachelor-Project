import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import LogLocator

fs = 15

fname = "../data/Rtc.txt"
dataf = open(fname, "r")
g, m1, z, gamma = [float(x) for x in dataf.readline().split()]
dataf.close()

dts = ["Ef", "SM", "Inv", "Vis"]
dt = [(d, "f") for d in dts]
data = np.loadtxt(fname, dtype = dt, skiprows = 1)

for i in range(1, len(dts)):
	d = dts[i]
	label = r"${\rm %s}$" % d
	plt.plot(data["Ef"], data[d], label = label)

plt.xscale("log")

v = list(plt.axis())
v[0] = data["Ef"][0]
v[1] = data["Ef"][-1]
v[2] = 0
v[3] = 0.8
plt.axis(v)

plt.gca().xaxis.set_minor_locator(LogLocator(subs = "all", numticks = 1e5))
plt.gca().set_yticks(np.arange(v[2], v[3] + 1e-5, 0.05), minor = True)

plt.xlabel(r"$E_f{\rm\ [GeV]}$")
plt.ylabel(r"$R_{\rm tc}$")

s = ""
if m1 == 0:
	s += r"$m_1=0$"
else:
	s += r"$m_1=%g{\rm\ eV}$" % m1
s += "$,$ "
s += r"$z=%g$" % z
s += "\n"
s += r"$\gamma=%g$" % gamma
s += "$,$ "
s += r"$g^{(\prime)}_{ij}="
p = np.log10(g)
if int(p) == p:
	s += r"10^{%i}" % p
else:
	s += r"%.1g\times10^{%i}" % (g / (10 ** p), p)
s += "$"
plt.text(0.99, 0.98, s, ha = "right", va = "top", fontsize = fs, transform = plt.gca().transAxes)

# IceCube ROI
fc = "lightblue"
alpha = 0.4
plt.fill_between([1e4, 1e7], [v[2], v[2]], [v[3], v[3]], facecolor = fc, edgecolor = None, alpha = alpha, zorder = -10)
plt.fill_between(v[:2], [v[3] + 1, v[3] + 1], [v[3] + 2, v[3] + 2], facecolor = fc, edgecolor = "k", lw = 0.5, alpha = alpha)

plt.legend(fontsize = fs, loc = 4)

plt.savefig("../fig/Rtc.pdf")

