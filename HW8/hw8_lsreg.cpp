/* AEP 4380 Homework #8

Least squares curve fitting using Gauss-Jordan is tested
Data source: http://cdiac.ornl.gov/ftp/trends/co2/ljo.dat

Run on a core i7 using clang 1000.11.45.2 on macOS Mojave

Kevin Juan 7 November 2018
*/

#include <cmath>
#include <cstdlib>
#include <fstream>  //  stream file IO
#include <iomanip>  //  to format the output
#include <iostream> //  stream IO
#include <string>   // STD strings
#include <vector>   // STD vector class
#define ARRAYT_BOUNDS_CHECK
#include "arrayt.hpp" // to use arrays
#define SWAP(a, b)                                                             \
    {                                                                          \
        double temp(a);                                                        \
        (a) = (b);                                                             \
        (b) = temp;                                                            \
    }

using namespace std;

// seven parameter model
double model7(double t, int i) {
    const static double pi = 4.0 * atan(1.0);
    if (i == 0) {
        return 1.0;
    } else if (i == 1) {
        return t;
    } else if (i == 2) {
        return t * t;
    } else if (i == 3) {
        return sin(2.0 * pi * t / 12.0);
    } else if (i == 4) {
        return sin(2.0 * pi * t / 6.0);
    } else if (i == 5) {
        return cos(2.0 * pi * t / 12.0);
    } else {
        return cos(2.0 * pi * t / 6.0);
    }
}

// five parameter model
double model5(double t, int i) {
    const static double pi = 4.0 * atan(1.0);
    if (i == 0) {
        return 1.0;
    } else if (i == 1) {
        return t;
    } else if (i == 2) {
        return t * t;
    } else if (i == 3) {
        return sin(2.0 * pi * t / 12.0);
    } else {
        return cos(2.0 * pi * t / 12.0);
    }
}

// polynomial model
double modelPoly(double t, int i) {
    if (i == 0) {
        return 1.0;
    } else if (i == 1) {
        return t;
    } else {
        return t * t;
    }
}

// nr3 gaussj
template <class T> void gaussj(arrayt<T> &a, arrayt<T> &b) {
    int i, icol, irow, j, k, l, ll, n = a.n1(), m = b.n2();
    double big, dum, pivinv;
    arrayt<T> indxc(n), indxr(n), ipiv(n);
    for (j = 0; j < n; j++)
        ipiv(j) = 0;
    for (i = 0; i < n; i++) {
        big = 0.0;
        for (j = 0; j < n; j++)
            if (ipiv(j) != 1)
                for (k = 0; k < n; k++) {
                    if (ipiv(k) == 0) {
                        if (abs(a(j, k)) >= big) {
                            big = abs(a(j, k));
                            irow = j;
                            icol = k;
                        }
                    }
                }
        ++(ipiv(icol));
        if (irow != icol) {
            for (l = 0; l < n; l++)
                SWAP(a(irow, l), a(icol, l));
            for (l = 0; l < m; l++)
                SWAP(b(irow, l), b(icol, l));
        }
        indxr(i) = irow;
        indxc(i) = icol;
        if (a(icol, icol) == 0.0)
            throw("gaussj: Singular Matrix");
        pivinv = 1.0 / a(icol, icol);
        a(icol, icol) = 1.0;
        for (l = 0; l < n; l++)
            a(icol, l) *= pivinv;
        for (l = 0; l < m; l++)
            b(icol, l) *= pivinv;
        for (ll = 0; ll < n; ll++)
            if (ll != icol) {
                dum = a(ll, icol);
                a(ll, icol) = 0.0;
                for (l = 0; l < n; l++)
                    a(ll, l) -= a(icol, l) * dum;
                for (l = 0; l < m; l++)
                    b(ll, l) -= b(icol, l) * dum;
            }
    }
    for (l = n - 1; l >= 0; l--) {
        if (indxr(l) != indxc(l))
            for (k = 0; k < n; k++)
                SWAP(a(k, indxr(l)), a(k, indxc(l)));
    }
}

// nr3 gauss j
template <class T> void gaussj(arrayt<T> &a) {
    arrayt<T> b(a.n1(), 0);
    gaussj(a, b);
}

int main() {
    int i, j, npts, year, t, nval, l, k, nEqs7 = 7, nEqs5 = 5, nEqs3 = 3;
    double co2, ymin, ymax, sumF, sumB, sumChi, fit5ChiSq, fit7ChiSq;
    double fit5Val, fit7Val;
    // use dynamically sized container classes
    string cline;
    vector<double> x, y;

    ifstream fp; // input file stream
    fp.open("ljo.dat");
    if (fp.fail()) {
        cout << "Can't open file." << endl;
        exit(0);
    }
    ofstream fpOrig;
    fpOrig.open("original_data.dat");
    if (fpOrig.fail()) {
        cout << "cannot open file" << endl;
        return (EXIT_SUCCESS);
    }

    //--------- read data from file in complicted format -------
    // skip first 16 lines
    for (i = 0; i < 16; i++)
        getline(fp, cline); // read a whole line
    t = 0;                  // time in months
    npts = 0;               // number of data points
    ymin = 1000.0;
    ymax = -ymin;
    for (i = 0; i < 70; i++) {
        fp >> year;
        if (0 == i)
            nval = 11;
        else
            nval = 12; // line is short(?)
        for (j = 0; j < nval; j++) {
            fp >> co2;
            if (co2 > 0.0) {
                x.push_back(t);   // use auto sizing because we don't
                y.push_back(co2); // know how many elements there will be
                if (y[npts] > ymax)
                    ymax = y[npts]; // x and y index like an array
                if (y[npts] < ymin)
                    ymin = y[npts]; // could use co2 here also
                npts++;
            }
            fpOrig << setw(5) << t << setw(10) << co2 << endl;
            t += 1;
        }
        if (year >= 2007)
            break;          // end of file
        getline(fp, cline); // read rest of line
    }

    // initialize arrays
    arrayt<double> sigmaSq(npts);
    arrayt<double> F7(nEqs7, nEqs7);
    arrayt<double> b7(nEqs7, nEqs7);
    arrayt<double> F5(nEqs5, nEqs5);
    arrayt<double> b5(nEqs5, nEqs5);
    arrayt<double> F3(nEqs3, nEqs3);
    arrayt<double> b3(nEqs3, nEqs3);

    // calculate error values
    for (i = 0; i < npts; i++) {
        sigmaSq(i) = 0.002 * 0.002 * y[i] * y[i];
    }

    // initialize F_lk for 7 parameter
    for (l = 0; l < nEqs7; l++) {
        for (k = 0; k < nEqs7; k++) {
            sumF = 0.0;
            sumB = 0.0;
            for (i = 0; i < npts; i++) {
                sumF += model7(x[i], l) * model7(x[i], k) / sigmaSq(i);
                sumB += y[i] * model7(x[i], l) / sigmaSq(i);
            }
            F7(l, k) = sumF;
            if (k == 0) {
                b7(l, k) = sumB;
            } else {
                b7(l, k) = 0.0;
            }
        }
    }

    // initialize F_lk for 5 parameter
    for (l = 0; l < nEqs5; l++) {
        for (k = 0; k < nEqs5; k++) {
            sumF = 0.0;
            sumB = 0.0;
            for (i = 0; i < npts; i++) {
                sumF += model5(x[i], l) * model5(x[i], k) / sigmaSq(i);
                sumB += y[i] * model5(x[i], l) / sigmaSq(i);
            }
            F5(l, k) = sumF;
            if (k == 0) {
                b5(l, k) = sumB;
            } else {
                b5(l, k) = 0.0;
            }
        }
    }

    // initialize F_lk for polynomial fit
    for (l = 0; l < nEqs3; l++) {
        for (k = 0; k < nEqs3; k++) {
            sumF = 0.0;
            sumB = 0.0;
            for (i = 0; i < npts; i++) {
                sumF += modelPoly(x[i], l) * modelPoly(x[i], k) / sigmaSq(i);
                sumB += y[i] * modelPoly(x[i], l) / sigmaSq(i);
            }
            F3(l, k) = sumF;
            if (k == 0) {
                b3(l, k) = sumB;
            } else {
                b3(l, k) = 0.0;
            }
        }
    }

    // Gauss-Jordan elimination
    gaussj(F7, b7);
    gaussj(F5, b5);
    gaussj(F3, b3);

    // output 7 parameter error
    for (i = 0; i < nEqs7; i++) {
        for (j = 0; j < nEqs7; j++) {
            if (i == j) {
                cout << "Error in parameter a_" << i + 1
                     << " for the 7 parameter fit is " << F7(i, j) << endl;
            }
        }
    }

    // output 5 parameter error
    for (i = 0; i < nEqs5; i++) {
        for (j = 0; j < nEqs5; j++) {
            if (i == j) {
                cout << "Error in parameter a_" << i + 1
                     << " for the 5 parameter fit is " << F5(i, j) << endl;
            }
        }
    }

    // output 5 parameter error
    for (i = 0; i < nEqs3; i++) {
        for (j = 0; j < nEqs3; j++) {
            if (i == j) {
                cout << "Error in parameter a_" << i + 1
                     << " for the polynomial fit is " << F3(i, j) << endl;
            }
        }
    }

    // output 7 parameter fit
    ofstream fpFit7;
    fpFit7.open("fit7_data.dat");
    if (fpFit7.fail()) {
        cout << "cannot open file" << endl;
        return (EXIT_SUCCESS);
    }
    for (i = 0; i < nEqs7; i++) {
        fpFit7 << b7(i, 0) << endl;
    }

    // output 5 parameter fit
    ofstream fpFit5;
    fpFit5.open("fit5_data.dat");
    if (fpFit5.fail()) {
        cout << "cannot open file" << endl;
        return (EXIT_SUCCESS);
    }
    for (i = 0; i < nEqs5; i++) {
        fpFit5 << b5(i, 0) << endl;
    }

    // output polynomial fit
    ofstream fpFitPoly;
    fpFitPoly.open("fitPoly_data.dat");
    if (fpFitPoly.fail()) {
        cout << "cannot open file" << endl;
        return (EXIT_SUCCESS);
    }
    for (i = 0; i < nEqs3; i++) {
        fpFitPoly << b3(i, 0) << endl;
    }

    // chi-square calculation for 7 parameter
    sumChi = 0.0;
    for (i = 0; i < npts; i++) {
        fit7Val = b7(0, 0) * model7(i, 0) + b7(1, 0) * model7(i, 1) +
                  b7(2, 0) * model7(i, 2) + b7(3, 0) * model7(i, 3) +
                  b7(4, 0) * model7(i, 4) + b7(5, 0) * model7(i, 5) +
                  b7(6, 0) * model7(i, 6);
        sumChi += (y[i] - fit7Val) * (y[i] - fit7Val) / sigmaSq(i);
    }
    sumChi /= (npts - nEqs7);
    cout << "Chi-Squared for 7 parameter fit: " << sumChi << endl;

    // chi-square calculation for 7 parameter
    sumChi = 0.0;
    for (i = 0; i < npts; i++) {
        fit5Val = b5(0, 0) * model5(i, 0) + b5(1, 0) * model5(i, 1) +
                  b5(2, 0) * model5(i, 2) + b5(3, 0) * model5(i, 3) +
                  b5(4, 0) * model5(i, 4);
        sumChi += (y[i] - fit5Val) * (y[i] - fit5Val) / sigmaSq(i);
    }
    sumChi /= (npts - nEqs5);
    cout << "Chi-Squared for 5 parameter fit: " << sumChi << endl;

    return (EXIT_SUCCESS);
}