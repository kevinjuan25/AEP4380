/* AEP 4380 Final Project

ADI Method tested on the 2D Heat Equation

Run on a core i7 using clang 1000.11.45.2 on macOS Mojave

Kevin Juan 12 December 2018
*/

#include <cmath>    // use math package
#include <cstdlib>  // plain C
#include <fstream>  // stream file IO
#include <iomanip>  // to format the output
#include <iostream> // stream IO
#define ARRAYT_BOUNDS_CHECK
#include "bigarrayt.hpp" // to use arrays
#define _USE_MATH_DEFINES

using namespace std;

const static double Lx = 5.0; // in m
const static double Ly = 5.0; // in m
const static double pi = M_PI;

template <class T>
void tridiag(arrayt<T> &a, arrayt<T> &b, arrayt<T> &c, arrayt<T> &d) {
    int k, kMax = b.n1();
    T bac;

    c(0) /= b(0);
    d(0) /= b(0);
    for (k = 1; k < kMax - 1; k++) {
        bac = b(k) - a(k) * c(k - 1);
        c(k) /= bac;
        d(k) = (d(k) - a(k) * d(k - 1)) / bac;
    }

    d(kMax - 1) = (d(kMax - 1) - a(kMax - 1) * d(kMax - 2)) /
                  (b(kMax - 1) - a(kMax - 1) * c(kMax - 2));

    for (k = kMax - 2; k >= 0; k--) {
        d(k) -= c(k) * d(k + 1);
    }
    return;
}

double TempInit(double x, double y, int initProfile, int material) {
    double TempInit;
    double xo = Lx / 2.0;
    double yo = Ly / 2.0;
    double dx = x - xo;
    double dy = y - yo;
    double A; // A should not exceed material melting point
    if (material == 1) {
        A = 500.0;
    } else if (material == 2) {
        A = 200.0;
    }
    if (initProfile == 1) {
        TempInit = A * exp(-(dx * dx / 2.0 + dy * dy / 2.0));
    } else if (initProfile == 2) {
        TempInit = A * abs(sin(5.0 * x) * cos(5.0 * y));
    } else if (initProfile == 3) {
        TempInit = A * abs(sin(dx * dx + dy * dy));
    } else if (initProfile == 4) {
        TempInit = A * exp(-(x * x / 10.0 + dy * dy / 2.0));
    } else {
        TempInit = A * exp(-x / xo);
    }
    return TempInit;
}

double S(double t, int sourceType, double tMax, double qMax) {
    // source is a heat flow in J/s
    // if qMax is positive, heat is added to the system
    // if qMax is negative, heat is removed from the system
    double source;
    if (sourceType == 1) {
        source = -qMax * sin(2.0 * pi * t);
    } else if (sourceType == 2) {
        source = -qMax * exp(-t / 10.0) * sin(2 * pi * t);
    } else if (sourceType == 3) {
        source = -qMax;
    } else if (sourceType == 4) {
        source = -qMax / (1 + exp(-(t - tMax / 2.0)));
    } else {
        source = 0.0;
    }
    return source / 2.0;
}

int main() {
    int i, j, n, iMax = 101, jMax = 101, nMax;
    double dx = Lx / (iMax - 1), dy = Ly / (jMax - 1);
    double currTime, midTime;
    int initProfile, sourceType, material;
    double x, y;
    double alpha; // in m^2/s
    double Cp;    // heat capacity in J/g*C
    double rho;   // density in g/m^3
    double t;     // time in seconds
    double qMax;  // max heat flow for the source in J/s
    int nSaves, tData = 0;

    // asks user for initial profile, source type, and simulation time
    cout << "What material should be used? (1 = Aluminum, 2 = PVC)" << endl;
    cin >> material;
    cout << "What initial profile should be used? (1 = Gaussian, 2 = Bump, 3 = "
            "Ripple, 4 = Elliptic Half-Gaussian, Exponential Otherwise)"
         << endl;
    cin >> initProfile;
    cout << "What source type should be used? (1 = Oscillatory, 2 = Decaying "
            "Oscillatory, 3 = Constant, 4 = Logistic, No Source Otherwise)"
         << endl;
    cin >> sourceType;
    cout << "What value should be use for maximum heat flow in J/s?" << endl;
    cin >> qMax;
    cout << "What simulation time should be used? (Input as Seconds)" << endl;
    cin >> t;
    cout << "How many time steps should be used? (Add 1 to account for t = 0)"
         << endl;
    cin >> nMax;
    cout << "How many profiles should be saved? (Excluding t = 0)" << endl;
    cin >> nSaves;

    double dt = t / (nMax - 1);

    if (material == 1) {
        // Aluminum
        alpha = 9.7e-5; // in m^2/s
        Cp = 0.9;       // heat capacity in J/g*C
        rho = 2.71e6;   // density in g/m^3
    } else if (material == 2) {
        // PVC
        alpha = 8.0e-8; // in m^2/s
        Cp = 0.9;       // heat capacity in J/g*C
        rho = 1.38e6;   // density in g/m^3
    }

    double rx = alpha * dt / (2.0 * dx * dx);
    double ry = alpha * dt / (2.0 * dy * dy);

    // we only need one set of vectors for an equal size area
    arrayt<double> Temp(iMax, jMax);
    arrayt<double> a(iMax);
    arrayt<double> b(iMax);
    arrayt<double> c(iMax);
    arrayt<double> d(iMax);
    arrayt<int> saves(nSaves);

    // initialize array of iterations to save profiles
    for (n = 0; n < nSaves; n++) {
        saves(n) = (n + 1) * t;
    }

    // output data for x
    ofstream fpx;
    fpx.open("x_vals.dat");
    if (fpx.fail()) {
        cout << "cannot open file" << endl;
        return (EXIT_SUCCESS);
    }
    for (i = 0; i < iMax; i++) {
        fpx << dx * i << endl;
    }
    fpx.close();

    // output data for y
    ofstream fpy;
    fpy.open("y_vals.dat");
    if (fpy.fail()) {
        cout << "cannot open file" << endl;
        return (EXIT_SUCCESS);
    }
    for (j = 0; j < jMax; j++) {
        fpy << dy * j << endl;
    }
    fpy.close();

    // output data for t = 0
    ofstream fpTempInit;
    fpTempInit.open("000.dat");
    if (fpTempInit.fail()) {
        cout << "cannot open file" << endl;
        return (EXIT_SUCCESS);
    }
    for (i = 0; i < iMax; i++) {
        for (j = 0; j < jMax; j++) {
            fpTempInit << setw(15)
                       << TempInit(i * dx, j * dy, initProfile, material);
            Temp(i, j) = TempInit(i * dx, j * dy, initProfile, material);
        }
        fpTempInit << endl;
    }

    cout << "Iteration: 0, Time: 0, Temp @ Center: "
         << Temp((iMax - 1) / 2, (jMax - 1) / 2) << endl;

    for (n = 1; n < nMax; n++) {

        currTime = n * dt;
        midTime = (n - 0.5) * dt;

        // begin column-wise sweep
        for (i = 0; i < iMax; i++) {
            for (j = 0; j < jMax; j++) {
                a(j) = -rx;
                b(j) = 1.0 + 2.0 * rx;
                c(j) = -rx;
                if (j == 0) {
                    d(j) = (1.0 - 2.0 * ry) * Temp(i, j) - ry * Temp(i, j + 1) -
                           dt * S(midTime, sourceType, t, qMax) / (rho * Cp);
                } else if (j == jMax - 1) {
                    d(j) = -ry * Temp(i, j - 1) +
                           (1.0 - 2.0 * ry) * Temp(i, j) -
                           dt * S(midTime, sourceType, t, qMax) / (rho * Cp);
                } else {
                    d(j) = -ry * Temp(i, j - 1) +
                           (1.0 - 2.0 * ry) * Temp(i, j) - ry * Temp(i, j + 1) -
                           dt * S(midTime, sourceType, t, qMax) / (rho * Cp);
                }
            }

            // call to tridiag
            tridiag<double>(a, b, c, d);
            for (j = 0; j < jMax; j++) {
                Temp(i, j) = d(j);
            }
        }

        // begin row-wise sweep
        for (j = 0; j < jMax; j++) {
            for (i = 0; i < iMax; i++) {
                a(i) = -ry;
                b(i) = 1.0 + 2.0 * ry;
                c(i) = -ry;
                if (i == 0) {
                    d(i) = (1.0 - 2.0 * rx) * Temp(i, j) - rx * Temp(i + 1, j) -
                           dt * S(currTime, sourceType, t, qMax) / (rho * Cp);
                } else if (i == iMax - 1) {
                    d(i) = -rx * Temp(i - 1, j) +
                           (1.0 - 2.0 * rx) * Temp(i, j) -
                           dt * S(currTime, sourceType, t, qMax) / (rho * Cp);
                } else {
                    d(i) = -rx * Temp(i - 1, j) +
                           (1.0 - 2.0 * rx) * Temp(i, j) - rx * Temp(i + 1, j) -
                           dt * S(currTime, sourceType, t, qMax) / (rho * Cp);
                }
            }

            // call to tridiag
            tridiag<double>(a, b, c, d);
            for (i = 0; i < iMax; i++) {
                Temp(i, j) = d(i);
            }
        }
        if (n == saves(tData)) {
            cout << "Iteration: " << n << ", Time: " << currTime
                 << ", Temp @ Center: " << Temp((iMax - 1) / 2, (jMax - 1) / 2)
                 << endl;
            // output data at different times to be used for animation
            ofstream fpTemp;
            fpTemp.open(to_string(n) + ".dat");
            if (fpTemp.fail()) {
                cout << "cannot open file" << endl;
                return (EXIT_SUCCESS);
            }

            for (i = 0; i < iMax; i++) {
                for (j = 0; j < jMax; j++) {
                    fpTemp << setw(15) << Temp(i, j);
                }
                fpTemp << endl;
            }
            tData += 1;
        }
    }
    return (EXIT_SUCCESS);
}