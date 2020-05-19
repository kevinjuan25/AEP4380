from pylab import*

data = loadtxt("hw2.dat")

plot(data[0:101, 0], data[0:101, 1])
xlim(0, 1)
xlabel("x")
ylabel("K(x)")
savefig('hw2_trapezoid.png')
