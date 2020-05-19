import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from mpl_toolkits import mplot3d

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

rng_data = np.loadtxt('RNG_data.dat')
crystal_data = np.loadtxt('crystal_data.dat')
heightwidth_data = np.loadtxt('heightwidth_data.dat')

# plt.hist(rng_data, 100)
# plt.xlim(0, 1)
# plt.xlabel("RNG Value")
# plt.ylabel("Count")
# plt.savefig('rng_hist.png', dpi=300)
# plt.show()

# plt.clf()

# plt.scatter(rng_data[0:int(len(rng_data) / 10) - 1],
#             rng_data[1:int(len(rng_data) / 10)], s=1)
# plt.xlim(0, 1)
# plt.ylim(0, 1)
# plt.xlabel("$x_i$")
# plt.ylabel("$x_{i + 1}$")
# plt.savefig('rng_pairs.png', dpi=300)
# plt.show()

plt.clf()

xVals = np.array([])
zVals = np.array([])
count = 0
for i in range(len(crystal_data)):
    if (crystal_data[i, 1] == 13):
        xVals = np.append(xVals, crystal_data[i, 0])
        zVals = np.append(zVals, crystal_data[i, 2])

plt.subplot(111, aspect='equal')
plt.scatter(xVals, zVals)
plt.savefig('crystal_xsec.png', dpi=300)
plt.show()

plt.clf()

ax = plt.axes(projection='3d')
ax.scatter3D(crystal_data[:, 0],
             crystal_data[:, 1], crystal_data[:, 2])
ax.set_xlabel("$x$")
ax.set_ylabel("$y$")
ax.set_zlabel("$z$")
plt.savefig('3d_crystal.png', dpi=300)
plt.show()

plt.clf()

plt.plot(heightwidth_data[:, 0],
         heightwidth_data[:, 1], heightwidth_data[:, 0], heightwidth_data[:, 2])
plt.xlabel("Time Step")
plt.ylabel("Atom Lengths")
plt.legend(("Average Height", "Width"))
plt.savefig('hw.png', dpi=300)
plt.show()
