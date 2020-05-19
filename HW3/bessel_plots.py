from pylab import*

dataj = loadtxt("besselj_func.dat")
datay = loadtxt("bessely_func.dat")

plot(dataj[:, 0], dataj[:, 1], dataj[:, 0], dataj[:, 2],
     datay[:, 0], datay[:, 1], datay[:, 0], datay[:, 2])
xlim(0, 20)
xlabel("x")
ylabel("y")
grid(True)
legend(("J_0(x)", "J_1(x)", "Y_0(x)", "Y_1(x)"))
savefig('hw3_bessel_functions.png', dpi=300)

plot(datay[:, 0], datay[:, 3])
xlim(0, 20)
xlabel("x")
ylabel("y")
grid(True)
savefig('hw3_eq_to_solve.png', dpi=300)
