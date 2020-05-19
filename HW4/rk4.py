import matplotlib.pyplot as plt
from matplotlib import rc
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

data = np.loadtxt("harmonic_osc.dat")
data2 = np.loadtxt("jaynes_cummings_h1.dat")
data3 = np.loadtxt("jaynes_cummings_h0p1.dat")
data4 = np.loadtxt("jaynes_cummings_h0p01.dat")
data5 = np.loadtxt("jaynes_cummings_h0p001.dat")
data6 = np.loadtxt("jaynes_cummings.dat")


plt.rc('text', usetex=True)
plt.rc('font', family='serif')

plt.plot(data[:, 0], data[:, 1], linestyle='-')
plt.plot(data[:, 0], np.cos(data[:, 0]), linestyle='--')
plt.legend(("RK4", "Numpy Cosine"))
plt.xlim(0, 100)
plt.xlabel("$t$")
plt.ylabel("$\cos(\omega t)$")
plt.savefig('harmoinc_osc.png', dpi=300)
plt.close()

plt.plot(data[:, 1], data[:, 2])
plt.xlim(-1.5, 1.5)
plt.ylim(-1.5, 1.5)
plt.gca().set_aspect('equal', adjustable='box')
plt.xlabel("$x(t)$")
plt.ylabel("$y(t)$")
plt.savefig('harmonic_osc_phase.png', dpi=300)
plt.close()

plt.plot(data2[:, 0], data2[:, 1], linestyle='-.')
plt.plot(data3[:, 0], data3[:, 1], linestyle=':')
plt.plot(data4[:, 0], data4[:, 1], linestyle='--')
plt.plot(data5[:, 0], data5[:, 1], linestyle='-')
plt.xlabel("$t$")
plt.ylabel("$x(t)$")
plt.legend(("$h=1.0$", "$h=0.1$", "$h=0.01$", "$h=0.001$"))
plt.xlim(0, 200)
plt.savefig('chaotic_sys_xt_h.png', dpi=300)
plt.close()

plt.plot(data6[:, 0], data6[:, 1], linestyle='-')
plt.xlabel("$t$")
plt.ylabel("$x(t)$")
plt.xlim(0, 200)
plt.savefig('chaotic_sys_xt.png', dpi=300)
plt.close()

plt.plot(data6[:, 0], data6[:, 3])
plt.xlim(0, 200)
plt.xlabel("$t$")
plt.ylabel("$z(t)$")
plt.savefig('chaotic_sys_zt.png', dpi=300)
plt.close()

plt.plot(data6[:, 1], data6[:, 2])
plt.gca().set_aspect('equal', adjustable='box')
plt.xlabel("$x(t)$")
plt.ylabel("$y(t)$")
plt.savefig('chaotic_sys_xy.png', dpi=300)
plt.close()

ax = plt.gca(projection='3d')
plt.plot(data6[:, 1], data6[:, 2], data6[:, 3])
ax.set_xlabel("$x(t)$")
ax.set_ylabel("$y(t)$")
ax.set_zlabel("$z(t)$")
plt.gca().set_aspect('equal', adjustable='box')
plt.savefig('chaotic_sys_xyz.png', dpi=300)
plt.close()
