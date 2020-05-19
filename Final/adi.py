import numpy as np
from pylab import*
from matplotlib import rc
from matplotlib.animation import FuncAnimation


n000 = np.loadtxt('000.dat')
x = np.loadtxt('x_vals.dat')
y = np.loadtxt('y_vals.dat')
t_data = np.loadtxt('t_vals.dat')
n100 = np.loadtxt('1000.dat')
n200 = np.loadtxt('2000.dat')
n300 = np.loadtxt('3000.dat')
n400 = np.loadtxt('4000.dat')
n500 = np.loadtxt('5000.dat')
n600 = np.loadtxt('6000.dat')
n700 = np.loadtxt('7000.dat')
n800 = np.loadtxt('8000.dat')
n900 = np.loadtxt('9000.dat')
n1000 = np.loadtxt('10000.dat')

data = np.zeros([11 * 101, 101, 11])

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

# clf()

# contourf(x, y, n000, 100, cmap='jet', vmin=0, vmax=200)
# colorbar()
# xlabel("$x$ (m)")
# ylabel("$y$ (m)")
# savefig('tempInit.png', dpi=300)
# show()

clf()

contourf(x, y, n500, 100, cmap='jet', vmin=0, vmax=200)
colorbar()
xlabel("$x$ (m)")
ylabel("$y$ (m)")
savefig('temp.png', dpi=300)
show()

fig = figure()
ax = axes(xlim=(0, 5), ylim=(0, 5))
xlabel("$x$ (m)")
ylabel("$y$ (m)")


for i in range(0, 11):
    if (i == 0):
        filename = '000.dat'
    else:
        filename = str(i) + '000.dat'
    data[i * 101:(i + 1) * 101, :, i] = np.loadtxt(filename)


def animate(i):
    clf()
    title("$t =$ " + str(i * 1000) + " s")
    contour = contourf(
        x, y, data[i * 101:(i + 1) * 101, :, i], 100, cmap='jet', vmin=0, vmax=200)
    maxT = np.amax(data[i * 101:(i + 1) * 101, :, i])
    colorbar(ticks=linspace(0, maxT, 5))
    return contour


anim = FuncAnimation(fig, animate, frames=11, interval=1000)
anim.save('profile.mp4', writer='ffmpeg')
show()
