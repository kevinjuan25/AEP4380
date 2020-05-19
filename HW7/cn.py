import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc


psi_re_t0_data = np.loadtxt('PsiRe_t0_vals.dat')
psi_im_t0_data = np.loadtxt('PsiIm_t0_vals.dat')
x_data = np.loadtxt('x_vals.dat')
PsiSq_vals_data = np.loadtxt('PsiSq_vals.dat')
potential_data = np.loadtxt('poten_vals.dat')


plt.rc('text', usetex=True)
plt.rc('font', family='serif')

# plt.plot(x_data, psi_re_t0_data, x_data, psi_im_t0_data)
# plt.xlabel("$x$ in \AA")
# plt.ylabel("$\psi$")
# plt.legend(("$\psi_{Re}$", "$\psi_{Im}$"))
# plt.xlim(100, 200)
# plt.savefig('Psi_t0.png', dpi=300)
# plt.show()

# plt.clf()

# plt.plot(x_data, potential_data)
# plt.xlabel("$x$ in \AA")
# plt.ylabel("$V(x)$ in eV")
# plt.xlim(200, 350)
# plt.ylim(0, 4.5)
# plt.savefig('potential.png', dpi=300)
# plt.show()

plt.clf()

plt.plot(x_data, PsiSq_vals_data[200, :])
plt.xlabel("$x$ in \AA")
plt.ylabel("$|\psi|^2$")
plt.xlim(0, 500)
plt.savefig('PsiSq_0p5.png', dpi=300)
plt.show()

plt.clf()

plt.plot(x_data, PsiSq_vals_data[400, :])
plt.xlabel("$x$ in \AA")
plt.ylabel("$|\psi|^2$")
plt.xlim(0, 500)
plt.savefig('PsiSq_1.png', dpi=300)
plt.show()

plt.clf()

plt.plot(x_data, PsiSq_vals_data[600, :])
plt.xlabel("$x$ in \AA")
plt.ylabel("$|\psi|^2$")
plt.xlim(0, 500)
plt.savefig('PsiSq_1p5.png', dpi=300)
plt.show()

plt.clf()

plt.plot(x_data, PsiSq_vals_data[800, :])
plt.xlabel("$x$ in \AA")
plt.ylabel("$|\psi|^2$")
plt.xlim(0, 500)
plt.tight_layout()
plt.savefig('PsiSq_2.png', dpi=300)
plt.show()

plt.clf()

plt.plot(x_data, PsiSq_vals_data[1000, :])
plt.xlabel("$x$ in \AA")
plt.ylabel("$|\psi|^2$")
plt.xlim(0, 500)
plt.savefig('PsiSq_2p5.png', dpi=300)
plt.show()
