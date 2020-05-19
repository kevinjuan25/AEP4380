from pylab import*

# read .dat file

data = loadtxt("hw1_sd.dat", 'float')
plot(data[:, 0], data[:, 1], data[:, 0], data[:, 4])
xlim(-7, 7)
xlabel("x")
ylabel("y")
title("Numerical Second Derivative for h = 0.05")
legend(("f(x)", "f''(x)"))
savefig('hw1_second_derivative.png')
