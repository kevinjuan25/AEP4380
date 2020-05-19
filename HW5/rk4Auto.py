import matplotlib.pyplot as plt
from matplotlib import rc
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

data = np.loadtxt("harmonic_osc.dat")
data2 = np.loadtxt("3_body_planets.dat")

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

# plt.plot(data[:, 0], data[:, 2], linestyle='-')
# plt.plot(data[:, 0], np.cos(data[:, 0]), linestyle='--')
# plt.legend(("RK4 Auto", "Numpy Cosine"))
# plt.xlim(0, 100)
# plt.xlabel("$t$")
# plt.ylabel("$\cos(\omega t)$")
# plt.savefig('harmonic_osc.png', dpi=300)
# plt.close()

# plt.plot(data[:, 2], data[:, 3])
# plt.xlim(-1.5, 1.5)
# plt.ylim(-1.5, 1.5)
# plt.gca().set_aspect('equal', adjustable='box')
# plt.xlabel("$x(t)$")
# plt.ylabel("$y(t)$")
# plt.savefig('harmonic_osc_phase.png', dpi=300)
# plt.close()

# plt.plot(data[:, 0], data[:, 1], linestyle='-')
# plt.xlabel("$t$")
# plt.ylabel("$h$")
# plt.xlim(0, 10)
# plt.savefig('harmonic_osc_step.png', dpi=300)
# plt.close()

# plt.plot(data2[:, 2], data2[:, 5], linestyle='-')
# plt.plot(data2[:, 3], data2[:, 6], linestyle=':')
# plt.plot(data2[:, 4], data2[:, 7], linestyle='-.')
# plt.xlabel("$x$")
# plt.ylabel("$y$")
# plt.legend(("Earth", "Moon", "Asteroid"))
# plt.savefig('trajectory.png', dpi=300)
# plt.close()

# plt.plot(data2[:, 0], data2[:, 1], linestyle='-')
# plt.xlabel("$t$")
# plt.ylabel("$h$")
# plt.savefig('ht.png', dpi=300)
# plt.close()

plt.plot(data2[:, 2], data2[:, 5], linestyle='-')
plt.xlabel("$x$")
plt.ylabel("$y$")
plt.savefig('trajectory_earth.png', dpi=300)
plt.close()
