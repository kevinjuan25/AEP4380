from pylab import*

data = loadtxt("hw2.dat", 'float')

plot(data[2:103, 0], data[2:103, 1])
xlim(0, 1)
xlabel("x")
ylabel("K(x)")
title("K(x) Calculated Using Trapezoid Rule")
savefig('hw2_trapezoid.png')
