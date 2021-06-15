import numpy as np
import matplotlib.pyplot as plt 

gamma = np.array([1,2,3,4,5])
fmu=np.array([0.396, 0.469, 0.492, 0.499, 0.502])
ftau = np.array([0.531, 0.584, 0.601, 0.606, 0.608])
fe = np.ones_like(fmu)

ftot = fmu+ftau+fe

fe = fe/ftot
fmu = fmu/ftot
fatu = ftau/ftot

print(fe)
print(fmu)
print(ftau)

fig = plt.figure(1)
plt.scatter(gamma, fe, label="e")
plt.scatter(gamma, fmu, label="mu")
plt.scatter(gamma, ftau, label="tau")
plt.legend(loc="center right")
plt.xlabel(r"\gamma")
plt.ylabel("Ratio")
plt.show()