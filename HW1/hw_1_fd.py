from pylab import*

# read .dat file

data = loadtxt("hw1_fd.dat", 'float')
plot(data[:, 0], data[:, 1], data[:, 0], data[:, 2],
     data[:, 0], data[:, 3], data[:, 0], data[:, 4])
xlim(-7, 7)
xlabel("x")
ylabel("y")
title("Numerical First Derivative for h = 0.5")
legend(("f(x)", "f'(x)_bd", "f'(x)_fd", "f'(x)_cd"))
savefig('hw1_first_derivative.png')
