/* AEP 4380 Homework #7

Crank-Nicholson method is used
Tested on time-dependent Schrodinger Equation

Run on a core i7 using clang 1000.11.45.2 on macOS Mojave

Kevin Juan 31 October 2018
*/

#include <cmath>    // use math package
#include <cstdlib>  // plain C
#include <fstream>  // stream file IO
#include <iomanip>  // to format the output
#include <iostream> // stream IO
#define ARRAYT_BOUNDS_CHECK
#include "arrayt.hpp" // to use arrays
#include <complex>

using namespace std;

typedef complex<double> I;
typedef arrayt<I> arrayc;

template <class T>
void tridiag(arrayc &a, arrayc &b, arrayc &c, arrayc &psi, arrayc &d) {
    int j, jMax = b.n1();
    I bac;

    c(0) = c(0) / b(0);
    d(0) = d(0) / b(0);
    for (j = 1; j < jMax - 1; j++) {
        bac = b(j) - a(j) * c(j - 1);
        if (j < jMax - 2) {
            c(j) = c(j) / bac;
        }
        d(j) = (d(j) - a(j) * d(j - 1)) / bac;
    }

    psi(jMax - 1) = d(jMax - 1);
    for (j = jMax - 2; j >= 0; j--) {
        psi(j) = d(j) - c(j) * psi(j + 1);
    }
}

double Vx(double x, double x1, double x2, double x3, double w1, double w2,
          double w3, double V1, double V2, double V3) {
    double x1Sq = (x - x1) * (x - x1) / (w1 * w1);
    double x2Sq = (x - x2) * (x - x2) / (w2 * w2);
    double x3Sq = (x - x3) * (x - x3) / (w3 * w3);
    double Vx = V1 * exp(-x1Sq) + V2 * exp(-x2Sq) + V3 * exp(-x3Sq);

    return Vx;
}

I psiInit(double x, double L, double s, double k0) {
    I psiInit = exp(I(-(x - 0.3 * L) * (x - 0.3 * L) / (s * s), x * k0));
    return psiInit;
}

int main() {
    const static double hbar = 6.5821e-16, hBarSq2m = 3.801;
    const static double V1 = 3.8, V2 = 4.0, V3 = 4.2;
    const static double L = 500.0, x1 = 0.5 * L, x2 = 0.55 * L;
    const static double x3 = 0.6 * L, w1 = 10.0, w2 = 10.0, w3 = 10.0;
    const static double s = 15.0, k0 = 1.0, t = 2.5e-14;
    int jMax = 5001, nMax = 1001, j, n;
    double xDelta = L / (double)(jMax - 1), tDelta = t / (double)(nMax - 1);
    double poten, w = 2.0 * xDelta * xDelta / tDelta;

    arrayc a(jMax);
    arrayc b(jMax);
    arrayc c(jMax);
    arrayc d(jMax);
    arrayc psi(jMax);

    // output data for x
    ofstream fpx;
    fpx.open("x_vals.dat");
    if (fpx.fail()) {
        cout << "cannot open file" << endl;
        return (EXIT_SUCCESS);
    }
    for (j = 0; j < jMax; j++) {
        fpx << xDelta * j << endl;
    }
    fpx.close();

    // output data for t
    ofstream fpt;
    fpt.open("t_vals.dat");
    if (fpt.fail()) {
        cout << "cannot open file" << endl;
        return (EXIT_SUCCESS);
    }
    for (n = 0; n < nMax; n++) {
        fpt << tDelta * n << endl;
    }
    fpt.close();

    // output data for potential
    ofstream fpPoten;
    fpPoten.open("poten_vals.dat");
    if (fpPoten.fail()) {
        cout << "cannot open file" << endl;
        return (EXIT_SUCCESS);
    }
    for (j = 0; j < jMax; j++) {
        poten = Vx(j * xDelta, x1, x2, x3, w1, w2, w3, V1, V2, V3);
        fpPoten << poten << endl;
    }
    fpPoten.close();

    // initialize a, b, c, and d at t = 0
    for (j = 0; j < jMax - 1; j++) {
        poten = Vx(j * xDelta, x1, x2, x3, w1, w2, w3, V1, V2, V3);
        a(j) = 1.0;
        b(j) =
            I(-2.0 - xDelta * xDelta * poten / hBarSq2m, hbar * w / hBarSq2m);
        c(j) = 1.0;
        d(j) =
            -psiInit((j - 1) * xDelta, L, s, k0) +
            I(2.0 + xDelta * xDelta * poten / hBarSq2m, hbar * w / hBarSq2m) *
                psiInit(j * xDelta, L, s, k0) -
            psiInit((1 + j) * xDelta, L, s, k0);
        psi(j) = 0.0;
    }

    // output data for t = 0
    ofstream fpPsiReInit;
    fpPsiReInit.open("PsiRe_t0_vals.dat");
    if (fpPsiReInit.fail()) {
        cout << "cannot open file" << endl;
        return (EXIT_SUCCESS);
    }
    for (j = 0; j < jMax; j++) {
        fpPsiReInit << psiInit(j * xDelta, L, s, k0).real() << endl;
    }
    fpPsiReInit.close();

    ofstream fpPsiImInit;
    fpPsiImInit.open("PsiIm_t0_vals.dat");
    if (fpPsiImInit.fail()) {
        cout << "cannot open file" << endl;
        return (EXIT_SUCCESS);
    }
    for (j = 0; j < jMax; j++) {
        fpPsiImInit << psiInit(j * xDelta, L, s, k0).imag() << endl;
    }
    fpPsiImInit.close();

    ofstream fpPsiSqInit;
    fpPsiSqInit.open("PsiSq_t0_vals.dat");
    if (fpPsiSqInit.fail()) {
        cout << "cannot open file" << endl;
        return (EXIT_SUCCESS);
    }
    for (j = 0; j < jMax; j++) {
        fpPsiSqInit << (psiInit(j * xDelta, L, s, k0) *
                        conj(psiInit(j * xDelta, L, s, k0)))
                           .real()
                    << endl;
    }
    fpPsiSqInit.close();

    // output data for different times
    ofstream fpPsiSqVals;
    fpPsiSqVals.open("PsiSq_vals.dat");
    if (fpPsiSqVals.fail()) {
        cout << "cannot open file" << endl;
        return (EXIT_SUCCESS);
    }
    for (n = 0; n < nMax; n++) {
        tridiag<arrayc>(a, b, c, psi, d);
        for (j = 0; j < jMax; j++) {
            psi(j) = d(j);
            fpPsiSqVals << setw(15) << (psi(j) * conj(psi(j))).real();
        }
        fpPsiSqVals << endl;
    }
    fpPsiSqVals.close();
    return (EXIT_SUCCESS);
}