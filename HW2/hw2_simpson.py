from pylab import*

data = loadtxt("hw2.dat")

plot(data[102:203, 0], data[102:203, 1])
xlim(0, 1)
xlabel("x")
ylabel("K(x)")
savefig('hw2_simpson.png')
