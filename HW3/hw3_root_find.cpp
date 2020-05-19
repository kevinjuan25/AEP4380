/* AEP 4380 Homework #3

Test root finding
Bisection and False Position methods are tested

Run on a core i7 using clang 902.0.39.2

Kevin Juan 19 September 2018
*/

#include "bessel.h" // use bessel function file
#include "nr3.h"    // use nr3 file
#include <cmath>    // use math package
#include <cstdlib>  // plain C

#include <fstream>  // stream file IO
#include <iomanip>  // to format the output
#include <iostream> // stream IO

using namespace std;

Bessjy bessFunc; // create instance of the Bessel function

/* Employs the bisection method for a function f within bracket [xa, xb] to a
tolerance level of tol
*/
template <class T, T (*f)(T)> T bisect(T xa, T xb, T tol) {
    double f3, x1 = xa, x2 = xb, x3 = 0.5 * (x1 + x2);
    double f1 = f(x1), f2 = f(x2);
    int iter = 0;
    if ((f1 > 0 && f2 < 0) || (f1 < 0 && f2 > 0)) { // check function signs
        do {                                        // begin bisection algorithm
            x3 = 0.5 * (x1 + x2);
            f3 = f(x3);
            if ((f3 > 0 && f1 > 0) || (f3 < 0 && f1 < 0)) {
                x1 = x3;
                f1 = f3;
            } else {
                x2 = x3;
                f2 = f3;
            }
            iter++;
        } while (abs(f3) >= tol);
        cout << "The root is x = " << x3 << " found in " << iter
             << " iterations at f(x) = " << f3 << endl;
    } else { // function signs are not opposite
        cout << "f(x1) and f(x2) must have opposite signs" << endl;
    }
    return x3;
} // end bisect()

/* Employs the false position method for a function f within bracket [xa, xb]
to a tolerance level of tol
*/
template <class T, T (*f)(T)> T falsePosition(T xa, T xb, T tol) {
    double x1 = xa, x2 = xb, f1 = f(x1), f2 = f(x2);
    double x3 = x1 - f1 * (x2 - x1) / (f2 - f1), f3;
    int iter = 0;
    if ((f1 > 0 && f2 < 0) || (f1 < 0 && f2 > 0)) { // check function signs
        do { // begin false position algorithm
            x3 = x1 - f1 * (x2 - x1) / (f2 - f1);
            f3 = f(x3);
            if ((f3 > 0 && f1 > 0) || (f3 < 0 && f1 < 0)) {
                x1 = x3;
                f1 = f3;
            } else {
                x2 = x3;
                f2 = f3;
            }
            iter++;
        } while (abs(f3) >= tol);
        cout << "The root is x = " << x3 << " found in " << iter
             << " iterations at f(x) = " << f3 << endl;
    } else { // function signs are not opposite
        cout << "f(x1) and f(x2) must have opposite signs" << endl;
    }
    return x3;
} // end falsePosition()

double testBessel(double x) { // Bessel function to be tested
    return bessFunc.y0(x) * bessFunc.y1(x) * bessFunc.y1(x) * x * x -
           bessFunc.j0(x) * bessFunc.j1(x);
}

int main() {
    double x, tol = 0.0000001;
    int i;
    cout.precision(9);

    ofstream fp1;                 // output file using streams
    fp1.open("besselj_func.dat"); // open new file for output
    fp1.precision(9);             // select 9 digits

    if (fp1.fail()) {
        // or fp.bad()
        cout << "cannot open file" << endl;
        return (EXIT_SUCCESS);
    }

    for (i = 0; i <= 2000; i++) { // generate points for Bessel function plot
        x = i * 0.01;
        fp1 << setw(20) << x << setw(20) << bessFunc.j0(x) << setw(20)
            << bessFunc.j1(x) << setw(20) << endl;
    }
    fp1.close();

    ofstream fp2;                 // output file using streams
    fp2.open("bessely_func.dat"); // open new file for output
    fp2.precision(9);             // select 9 digits

    if (fp2.fail()) {
        // or fp.bad()
        cout << "cannot open file" << endl;
        return (EXIT_SUCCESS);
    }

    for (i = 75; i <= 2000; i++) { // generate points for Bessel function plot
        x = i * 0.01;
        fp2 << setw(20) << x << setw(20) << bessFunc.y0(x) << setw(20)
            << bessFunc.y1(x) << setw(20) << testBessel(x) << endl;
    }
    fp2.close();

    /* Use root finding methods for brackets obtained from plot to obtain the
     * first five roots
     */
    bisect<double, testBessel>(2, 3, tol);
    bisect<double, testBessel>(3, 4, tol);
    bisect<double, testBessel>(6, 7.4, tol);
    bisect<double, testBessel>(8.4, 8.6, tol);
    bisect<double, testBessel>(8.6, 8.7, tol);
    falsePosition<double, testBessel>(2, 3, tol);
    falsePosition<double, testBessel>(3, 4, tol);
    falsePosition<double, testBessel>(6, 7.4, tol);
    falsePosition<double, testBessel>(8.4, 8.6, tol);
    falsePosition<double, testBessel>(8.6, 8.7, tol);

    return (EXIT_SUCCESS);
}
