import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

orig_data = np.loadtxt('original_data.dat')
fit7_data = np.loadtxt('fit7_data.dat')
fit5_data = np.loadtxt('fit5_data.dat')
fitPoly_data = np.loadtxt('fitPoly_data.dat')


def fit7(t):
    return fit7_data[0] + fit7_data[1] * t + fit7_data[2] * t * t + fit7_data[3] * np.sin(np.pi * t / 6) + fit7_data[4] * np.sin(np.pi * t / 3) + fit7_data[5] * np.cos(np.pi * t / 6) + fit7_data[6] * np.cos(np.pi * t / 3)


def fit5(t):
    return fit5_data[0] + fit5_data[1] * t + fit5_data[2] * t * t + fit5_data[3] * np.sin(np.pi * t / 6) + fit7_data[4] * np.cos(np.pi * t / 6)


def fitPoly(t):
    return fitPoly_data[0] + fitPoly_data[1] * t + fitPoly_data[2] * t * t


plt.clf()

plt.plot(orig_data[:, 0], orig_data[:, 1])
plt.xlim(0, 500)
plt.ylim(320, 390)
plt.xlabel("$t$ in Months")
plt.ylabel("CO$_2$ Concentration in ppm")
plt.savefig('orig_data.png', dpi=300)
plt.show()

plt.clf()

plt.plot(orig_data[:, 0], orig_data[:, 1])
plt.plot(orig_data[:, 0], fit7(orig_data[:, 0]), linestyle='-.')
plt.plot(orig_data[:, 0], fitPoly(orig_data[:, 0]), linestyle=':')
plt.legend(("Original Data", "7 Parameter Fit", "Polynomial Fit"))
plt.xlim(0, 500)
plt.ylim(320, 390)
plt.xlabel("$t$ in Months")
plt.ylabel("CO$_2$ Concentration in ppm")
plt.savefig('data1.png', dpi=300)
plt.show()

plt.clf()

plt.plot(orig_data[:, 0], orig_data[:, 1])
plt.plot(orig_data[:, 0], fit5(orig_data[:, 0]), linestyle='--')
plt.plot(orig_data[:, 0], fitPoly(orig_data[:, 0]), linestyle=':')
plt.legend(("Original Data", "5 Parameter Fit", "Polynomial Fit"))
plt.xlim(0, 500)
plt.ylim(320, 390)
plt.xlabel("$t$ in Months")
plt.ylabel("CO$_2$ Concentration in ppm")
plt.savefig('data2.png', dpi=300)
plt.show()

plt.clf()

plt.plot(orig_data[:, 0], abs(fit7(orig_data[:, 0]) - orig_data[:, 1]))
plt.plot(orig_data[:, 0], abs(
    fit5(orig_data[:, 0]) - orig_data[:, 1]), linestyle=':')
plt.legend(("7 Parameter Fit Residuals", "5 Parameter Fit Residuals"))
plt.xlim(0, 500)
plt.xlabel("$t$ in Months")
plt.ylabel("Residual CO$_2$ Concentration in ppm")
plt.savefig('residual.png', dpi=300)
plt.show()
