/* AEP 4380 Homework #9

Monte Carlo method tested on crystal growth

Run on a core i7 using clang 1000.11.45.2 on macOS Mojave

Kevin Juan 14 November 2018
*/

#include "nr3.h"
#include "ran.h"
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <fstream>  //  stream file IO
#include <iomanip>  //  to format the output
#include <iostream> //  stream IO
#define ARRAYT_BOUNDS_CHECK
#include "bigarrayt.hpp"

using namespace std;

double height(arrayt<int> &L) {
    double h = 0;
    int i, j, k, Nx = L.n1(), Ny = L.n2(), Nz = L.n3();
    int lastZ = 1;

    for (i = 0; i < Nx; i++) {
        for (j = 0; j < Ny; j++) {
            for (k = 0; k < Nz; k++) {
                if (L(i, j, k) == 1 && k + 1 > lastZ) {
                    lastZ = k + 1;
                }
            }
            h += (double)lastZ;
            lastZ = 1;
        }
    }
    return h / (Nx * Ny);
}

double width(arrayt<int> &L, double h) {
    double wSq, hSq = 0;
    int i, j, k, Nx = L.n1(), Ny = L.n2(), Nz = L.n3();
    int lastZ = 1;

    for (i = 0; i < Nx; i++) {
        for (j = 0; j < Ny; j++) {
            for (k = 0; k < Nz; k++) {
                if (L(i, j, k) == 1 && k + 1 > lastZ) {
                    lastZ = k + 1;
                }
            }
            hSq += (double)(lastZ * lastZ);
            lastZ = 1;
        }
    }
    hSq /= (double)(Nx * Ny);
    wSq = hSq - h * h;
    return sqrt(wSq);
}

inline int periodic(int i, int n) { return (i + n) % n; }

double PTotal(arrayt<int> &L, bool dep, int x, int y, int z) {
    const static double PMax = 19.1041;
    int i, j, k, ix, jy;
    double PTotal = 0.0, px, p0 = 0.5, p1 = 0.5;

    if (dep == true) {
        px = p0;
    } else {
        px = p1;
    }
    for (i = -1; i <= 1; i++) {
        ix = periodic(x + i, L.n1());
        for (j = -1; j <= 1; j++) {
            jy = periodic(y + j, L.n2());
            for (k = -1; k <= 1; k++) {
                if (i != 0 || j != 0 || k != 0) {
                    PTotal += L(ix, jy, z + k) * px /
                              (PMax * sqrt(i * i + j * j + k * k));
                }
            }
        }
    }
    return PTotal;
}

// return true if point is a perimeter, false otherwise
bool isPerimeter(arrayt<int> &L, int x, int y, int z) {
    int i, j, k, ix, jy;
    bool perimeter = false;

    if (L(x, y, z) == 0) {
        for (i = -1; i <= 1; i += 2) {
            ix = periodic(x + i, L.n1());
            for (j = -1; j <= 1; j += 2) {
                jy = periodic(y + j, L.n2());
                for (k = -1; k <= 1; k += 2) {
                    if (L(ix, jy, z + k) == 1) {
                        perimeter = true;
                    }
                }
            }
        }
    }
    return perimeter;
}

// return 1 or -1 to get random neighbor for x, y, or z
int randNN(double randVal) {
    if (randVal > 0.5) {
        return 1;
    } else {
        return -1;
    }
}

int main() {
    unsigned int seed1, seed2;
    const static int Nx = 50, Ny = 25, Nz = 25, Na = 10;
    const static double PDep = 0.5;
    bool dep;
    double h, w;
    int i, j, k, tIter = 0, numDep = Nx * Ny, count;
    arrayt<int> x(Nx * Ny * Nz); // initialize array for x perimeters
    arrayt<int> y(Nx * Ny * Nz); // initialize array for y perimeters
    arrayt<int> z(Nx * Ny * Nz); // initialize array for z perimeters

    seed1 = 1;
    seed2 = time(NULL);
    struct Ranq1 randNum1(seed1);
    struct Ranq1 randNum2(seed2);

    // output RNG data
    ofstream fpRNG;
    fpRNG.open("RNG_data.dat");
    if (fpRNG.fail()) {
        cout << "cannot open file" << endl;
        return (EXIT_SUCCESS);
    }
    // initialize seed
    i = 0;
    // generate new random number from different seed
    while (i < 100000) {
        fpRNG << randNum1.doub() << endl;
        i++;
    }

    // initialize lattice points
    arrayt<int> L(Nx, Ny, Nz);
    for (i = 0; i < Nx; i++) {
        for (j = 0; j < Ny; j++) {
            for (k = 0; k < Nz; k++) {
                if (k == 0) {
                    L(i, j, k) = 1;
                } else {
                    L(i, j, k) = 0;
                }
            }
        }
    }

    // output height and width data
    ofstream fpHW;
    fpHW.open("heightwidth_data.dat");
    if (fpHW.fail()) {
        cout << "cannot open file" << endl;
        return (EXIT_SUCCESS);
    }

    // initial height and width
    h = height(L);
    w = width(L, h);
    fpHW << setw(10) << tIter << setw(10) << h << setw(10) << w << endl;

    while (numDep <= Nx * Ny * Na) {
        double R1 = randNum2.doub();
        double R2 = randNum2.doub();
        double R3 = randNum2.doub();
        double randVal = randNum2.doub();
        double xRand = randNum2.doub(); // random number for x neighbor
        double yRand = randNum2.doub(); // random number for y neighbor
        double zRand = randNum2.doub(); // random number for z neighbor
        int xVal, yVal, zVal, xInc, yInc, zInc;
        int xNbr, yNbr, zNbr;

        count = 0;
        // determine perimeter sites
        for (i = 0; i < Nx; i++) {
            for (j = 0; j < Ny; j++) {
                for (k = 1; k < Nz - 1; k++) {
                    if (isPerimeter(L, i, j, k) == true) {
                        x(count) = i;
                        y(count) = j;
                        z(count) = k;
                        count++;
                    }
                }
            }
        }

        // randomize deposition or diffusion
        if (R1 < PDep) {
            dep = true;
        } else {
            dep = false;
        }
        int position = randVal * count;
        xVal = x(position);
        yVal = y(position);
        zVal = z(position);
        xInc = randNN(xRand);
        yInc = randNN(yRand);
        zInc = randNN(zRand);
        double perimPTotal = PTotal(L, dep, xVal, yVal, zVal);

        // deposition
        if (dep == true) {
            if (perimPTotal > R2) {
                L(xVal, yVal, zVal) = 1;
                numDep++;
            }
        } else {
            xNbr = periodic(xVal + xInc, Nx);
            yNbr = periodic(yVal + yInc, Ny);
            zNbr = zVal + zInc;

            // diffusion
            if (L(xNbr, yNbr, zNbr) == 1 && zNbr != 0) {
                double neighborPTotal = PTotal(L, dep, xNbr, yNbr, zNbr);
                if (neighborPTotal < perimPTotal || perimPTotal > R3) {
                    L(xNbr, yNbr, zNbr) = 0;
                    L(xVal, yVal, zVal) = 1;

                    // move back if PTotal is 0 at new site
                    if (PTotal(L, dep, xVal, yVal, zVal) == 0.0) {
                        L(xNbr, yNbr, zNbr) = 1;
                        L(xVal, yVal, zVal) = 0;
                    }
                }
            }
        }
        h = height(L);
        w = width(L, h);
        tIter++;
        fpHW << setw(10) << tIter << setw(10) << h << setw(10) << w << endl;
    }
    fpHW.close();

    // output crystal data
    ofstream fp;
    fp.open("crystal_data.dat");
    if (fp.fail()) {
        cout << "cannot open file" << endl;
        return (EXIT_SUCCESS);
    }
    for (i = 0; i < Nx; i++) {
        for (j = 0; j < Ny; j++) {
            for (k = 0; k < Nz; k++) {
                if (L(i, j, k) == 1) {
                    fp << setw(5) << i << setw(5) << j << setw(5) << k << endl;
                }
            }
        }
    }
    fp.close();
}