/* AEP 4380 Homework #6

Finite difference method using successive over relaxation
Gauss-Seidel is used
Tested on a potential field with electrodes

Run on a core i7 using clang 1000.11.45.2 on macOS Mojave

Kevin Juan 24 October 2018
*/

#include <cmath>    // use math package
#include <cstdlib>  // plain C
#include <fstream>  // stream file IO
#include <iomanip>  // to format the output
#include <iostream> // stream IO
#define ARRAYT_BOUNDS_CHECK
#include "arrayt.hpp" // to use arrays

using namespace std;

int main() {
    float h, w, phiFD, delVMax = 1000.0, tol, dPhi, V2 = 1000.0, V1 = 0.0;
    float totalE = 0.0, ErSq, EzSq, C, VC = V2 - V1;
    int zBound = 50, rBound = 20, i, j, iter = 0, maxIter = 100000;
    const static float e0 = 8.8542 * 1.0e-3;

    // aks user for inputs for w, h, and tolerance
    cout << "What relxation parameter should be used?" << endl;
    cin >> w;
    cout << "What h value should be used?" << endl;
    cin >> h;
    cout << "What level of tolerance should be used?" << endl;
    cin >> tol;

    int zMax = (int)(zBound / h + 1.5), rMax = (int)(rBound / h + 1.5);
    int i25mm = (int)(25 / h + 0.5), j3mm = (int)(3 / h + 0.5);
    int j7mm = (int)(7 / h + 0.5), j5mm = (int)(5 / h + 0.5);
    int i17mm = (int)(17 / h + 0.5), i33mm = (int)(33 / h + 0.5);
    int i10mm = (int)(10 / h + 0.5), i40mm = (int)(40 / h + 0.5);
    int i21p5mm = (int)(215 / (10 * h) + 0.5);
    int i28p5mm = (int)(285 / (10 * h) + 0.5);
    arrayt<float> phi(zMax, rMax);
    arrayt<char> flag(zMax, rMax);
    arrayt<float> Ez(zMax, rMax);
    arrayt<float> Er(zMax, rMax);

    ofstream fpz;
    fpz.open("z_vals.dat");
    if (fpz.fail()) {
        cout << "cannot open file" << endl;
        return (EXIT_SUCCESS);
    }
    ofstream fpr;
    fpr.open("r_vals.dat");

    if (fpr.fail()) {
        cout << "cannot open file" << endl;
        return (EXIT_SUCCESS);
    }

    // write z and r values to file
    for (i = 0; i < zMax; i++) {
        fpz << setw(15) << i * h - zBound / 2 << endl;
    }
    for (j = 0; j < rMax; j++) {
        fpr << setw(15) << j * h << endl;
    }
    fpz.close();
    fpr.close();

    // initialize flag array and grid with boundaries and electrodes
    for (i = 0; i < zMax; i++) {
        for (j = 0; j < rMax; j++) {
            phi(i, j) = 0.0;
            flag(i, j) = '0';
            if (i == i25mm && (j >= j3mm && j <= j7mm)) {
                phi(i, j) = V2;
                flag(i, j) = '1';
            } else if (j == j7mm && (i >= i17mm && i <= i33mm)) {
                phi(i, j) = V2;
                flag(i, j) = '1';
            } else if (j == j3mm && (i >= i10mm && i <= i21p5mm)) {
                phi(i, j) = V1;
                flag(i, j) = '1';
            } else if (j == j5mm && (i >= i28p5mm && i <= i40mm)) {
                phi(i, j) = V1;
                flag(i, j) = '1';
            }
        }
    } // end phi and flag initialization

    // begin finite difference method with successive over-relaxation
    do {
        delVMax = -1.0;
        for (i = 1; i < zMax - 1; i++) {
            for (j = 0; j < rMax - 1; j++) {
                if (flag(i, j) == '0') {
                    if (j == 0) {
                        phiFD =
                            (4.0 * phi(i, 1) + phi(i + 1, 0) + phi(i - 1, 0)) /
                            6.0;
                        dPhi = phiFD - phi(i, j);
                        if (abs(dPhi) > delVMax) {
                            delVMax = abs(dPhi);
                        }
                        phi(i, j) = phi(i, j) + w * dPhi;
                    } else {
                        phiFD = (phi(i, j + 1) + phi(i, j - 1) + phi(i + 1, j) +
                                 phi(i - 1, j)) /
                                    4.0 +
                                (phi(i, j + 1) - phi(i, j - 1)) / (8.0 * j);
                        dPhi = phiFD - phi(i, j);
                        if (abs(dPhi) > delVMax) {
                            delVMax = abs(dPhi);
                        }
                        phi(i, j) = phi(i, j) + w * dPhi;
                    }
                }
            }
        }
        iter++;
        cout << "Iterations: " << iter << " w: " << w << " delVMax: " << delVMax
             << " V(0,0): " << phi(i25mm, 0) << endl;
    } while (iter < maxIter && delVMax > tol); // end finite difference method

    // write phi values to file
    ofstream fpDat;
    fpDat.open("phi_vals.dat");
    if (fpDat.fail()) {
        cout << "cannot open file" << endl;
        return (EXIT_SUCCESS);
    }
    for (i = 0; i < zMax; i++) {
        for (j = 0; j < rMax; j++) {
            fpDat << setw(15) << phi(i, j);
        }
        fpDat << endl;
    }
    fpDat.close();

    ofstream fpEz;
    fpEz.open("Ez_vals.dat");
    if (fpEz.fail()) {
        cout << "cannot open file" << endl;
        return (EXIT_SUCCESS);
    }

    ofstream fpEr;
    fpEr.open("Er_vals.dat");
    if (fpDat.fail()) {
        cout << "cannot open file" << endl;
        return (EXIT_SUCCESS);
    }

    // begin electric field calculations noting that edges are lost
    for (j = 0; j < rMax; j++) {
        Ez(0, j) = 0;
        Ez(zMax - 1, j) = 0;
        for (i = 1; i < zMax - 1; i++) {
            Ez(i, j) = -(phi(i + 1, j) - phi(i - 1, j)) / (2.0 * h);
        }
    }
    for (i = 0; i < zMax; i++) {
        Er(i, 0) = 0;
        Er(i, rMax - 1) = 0;
        for (j = 1; j < rMax - 1; j++) {
            Er(i, j) = -(phi(i, j + 1) - phi(i, j - 1)) / (2.0 * h);
        }
    }

    // write electric field values to files
    for (i = 0; i < zMax; i++) {
        for (j = 0; j < rMax; j++) {
            fpEz << setw(15) << Ez(i, j);
        }
        fpEz << endl;
    }
    for (i = 0; i < zMax; i++) {
        for (j = 0; j < rMax; j++) {
            fpEr << setw(15) << Er(i, j);
        }
        fpEr << endl;
    }
    fpEz.close();
    fpEr.close();

    // calculate capacitance
    for (i = 0; i < zMax; i++) {
        for (j = 0; j < rMax; j++) {
            EzSq = Ez(i, j) * Ez(i, j);
            ErSq = Er(i, j) * Er(i, j);
            totalE += j * (ErSq + EzSq);
        }
    }
    cout << totalE << endl;
    C = 2 * h * h * h * M_PI * totalE * e0 / (VC * VC);
    cout << "Capacitance in picoFarads: " << C << " pF" << endl;
}