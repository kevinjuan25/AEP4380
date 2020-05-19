/* AEP 4380 Homework #2

Test numerical integration
Trapezoid and Simpson's methods are tested

Run on a core i7 using clang 902.0.39.2

Kevin Juan 12 September 2018
*/

#include <assert.h>
#include <cmath>
#include <cstdlib> // plain C

#include <fstream>  // stream file IO
#include <iomanip>  // to format the output
#include <iostream> // stream IO

using namespace std;

/* Function K(x) to be integrated */
double Kx(double x, double theta) {
    return (1 / sqrt(1 - x * sin(theta) * sin(theta)));
}

/* Returns the value of the integral using Trapezoid Rule from a to b using N
 * intervals for the function *func that takes in 2 paramaters one of which is
 * x. Precondition: N must be positive
 */
double trapezoid(double (*func)(double, double), double a, double b, int N,
                 double x) {
    double h, S = 0.5 * (func(x, a) + func(x, b));
    double IOld, INew = (b - a) * S;
    int i, NInt = 1; // NInt is the minimum number of allowable intervals

    assert(N > 0); // check that the number of intervals is positive

    while (NInt <= N) {
        NInt = 2 * NInt;
        h = (b - a) / NInt;
        IOld = INew;
        for (i = 1; i <= NInt - 1; i += 2) {
            S = S + func(x, a + i * h);
        }
        INew = S * h;
    }
    return INew;
}

/* Returns the value of the integral using Simpson's Rule from a to b using N
 * intervals for the function *func that takes in 2 paramaters. Precondition: N
 * must be greater than 1.
 */
double simpson(double (*func)(double, double), double a, double b, int N,
               double x) {
    double h, S1 = func(x, a) + func(x, b), S2 = 0.0;
    double S4 = func(x, 0.5 * (a + b)), IOld;
    double INew = 0.5 * (b - a) * (S1 + 4 * S4) / 3;
    int i, NInt = 2; // NInt is the minimum number of allowable intervals

    assert(N > 1); // check that number of intervals is greater than 1

    while (NInt <= N) {
        NInt = 2 * NInt;
        h = (b - a) / NInt;
        S2 = S2 + S4;
        IOld = INew;
        S4 = 0; // set S4 to 0 to reset its value
        for (i = 1; i <= NInt - 1; i += 2) {
            // sum f(thetaA + i * dTheta) over 1, 3, 5, ... N - 1
            S4 = func(x, a + i * h) + S4;
        }
        INew = h * (S1 + 2 * S2 + 4 * S4) / 3;
    }
    return INew;
}

int main() {
    double thetaB = 4.0 * atan(1.0) / 2, thetaA = 0.0, tol = pow(10, -9);
    double xmin = 0.0, xmax = 0.9999, x;
    int i, j, NMax = 4096, pts = 100;

    ofstream fp;        // output file using streams
    fp.open("hw2.dat"); // open new file for output
    fp.precision(9);    // select 9 digits

    if (fp.fail()) {
        // or fp.bad()
        cout << "cannot open file" << endl;
        return (EXIT_SUCCESS);
    }
    fp << setw(45) << "Trapezoid Rule Values" << endl;
    fp << setw(20) << "x" << setw(20) << "Integral Value" << setw(20) << "N"
       << endl;

    /* This loop evaluates the integral for 100 equally spaced values of x for
     * K(x) using the Trapezoid method. The inner loop uses only enough
     * intervals to calculate K(x) to a tolerance of 10^-9.
     */
    for (x = xmin; x <= xmax; x += (xmax - xmin) / pts) {
        i = 2; // second fewest number of intervals to check error against tol
        while (i <= NMax &&
               abs(trapezoid(Kx, thetaA, thetaB, i, x) -
                   trapezoid(Kx, thetaA, thetaB, i / 2, x)) > tol) {
            i = 2 * i;
        }
        fp << setw(20) << x << setw(20) << trapezoid(Kx, thetaA, thetaB, i, x)
           << setw(20) << i << endl;
    }

    fp << setw(45) << "Simpson's Rule Values" << endl;

    /* This loop evaluates the integral for 100 equally spaced values of x for
     * K(x) using the Simpson's method. The inner loop uses only enough
     * intervals to calculate K(x) to a tolerance of 10^-9.
     */
    for (x = xmin; x <= xmax; x += (xmax - xmin) / pts) {
        j = 4; // second fewest number of intervals to check error against tol
        while (j <= NMax && abs(simpson(Kx, thetaA, thetaB, j, x) -
                                simpson(Kx, thetaA, thetaB, j / 2, x)) > tol) {
            j = 2 * j;
        }
        fp << setw(20) << x << setw(20) << simpson(Kx, thetaA, thetaB, j, x)
           << setw(20) << j << endl;
    }

    fp.close();
    return (EXIT_SUCCESS);
}
