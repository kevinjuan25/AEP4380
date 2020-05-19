/* AEP 4380 Homework #9

Fast Fourier Transform using fftw3
Tested on a 2D Wave Equation

Run on a core i7 using clang 1000.11.45.2 on macOS Mojave

Kevin Juan 28 November 2018
*/

#include "fftw3.h"
#include <cmath>
#include <cstdlib>
#include <fstream>  //  stream file IO
#include <iomanip>  //  to format the output
#include <iostream> //  stream IO
#define _USE_MATH_DEFINES

using namespace std;

int main() {
    int Nx = 512, Ny = 512, i, j;
    double k, t = 0.1, Lx = 500.0, Ly = 500.0;
    double xi[] = {0.4 * Lx, 0.5 * Lx, 0.6 * Lx};
    double yi[] = {0.4 * Ly, 0.6 * Ly, 0.4 * Ly};
    double si[] = {10.0, 20.0, 10.0};
    double Ai[] = {1.0, 2.0, 1.0};
    double kx[Ny];
    double ky[Ny];
    double v = 343.0, x, y, init1, init2, init3;

    fftw_complex *psi;
    fftw_plan forward, inverse;
    psi = (fftw_complex *)fftw_malloc(Nx * Ny * sizeof(fftw_complex));

    if (psi == NULL) {
        cout << "Cannot allocate array y" << endl;
        exit(EXIT_FAILURE);
    }
    ofstream fpx;
    fpx.open("x_data.dat");
    if (fpx.fail()) {
        cout << "cannot open file" << endl;
        return (EXIT_SUCCESS);
    }
    ofstream fpy;
    fpy.open("y_data.dat");
    if (fpy.fail()) {
        cout << "cannot open file" << endl;
        return (EXIT_SUCCESS);
    }
    ofstream fpReInit;
    fpReInit.open("wave_Re_init_data.dat");
    if (fpReInit.fail()) {
        cout << "cannot open file" << endl;
        return (EXIT_SUCCESS);
    }

    // define plans
    forward = fftw_plan_dft_2d(Nx, Ny, psi, psi, FFTW_FORWARD, FFTW_ESTIMATE);
    inverse = fftw_plan_dft_2d(Nx, Ny, psi, psi, FFTW_BACKWARD, FFTW_ESTIMATE);

    // initialize wave vectors
    for (i = 0; i < Nx; i++) {
        if (i >= Nx / 2) {
            kx[i] = (i - Nx) / Lx;
            ky[i] = (i - Ny) / Ly;
        } else {
            kx[i] = i / Lx;
            ky[i] = i / Ly;
        }
    }

    // initialize waves
    for (j = Ny - 1; j >= 0; j--) {
        y = j * Ly / Ny;
        fpy << y << endl;
        for (i = 0; i < Nx; i++) {
            x = i * Lx / Nx;
            if (j == 0) {
                fpx << x << endl;
            }
            init1 =
                Ai[0] *
                exp(-((x - xi[0]) * (x - xi[0]) + (y - yi[0]) * (y - yi[0])) /
                    (si[0] * si[0]));
            init2 =
                Ai[1] *
                exp(-((x - xi[1]) * (x - xi[1]) + (y - yi[1]) * (y - yi[1])) /
                    (si[1] * si[1]));
            init3 =
                Ai[2] *
                exp(-((x - xi[2]) * (x - xi[2]) + (y - yi[2]) * (y - yi[2])) /
                    (si[2] * si[2]));
            psi[j + i * Ny][0] = init1 + init2 + init3;
            psi[j + i * Ny][1] = 0.0;
            fpReInit << setw(15) << psi[j + i * Ny][0] << setw(15);
        }
        fpReInit << endl;
    }

    fftw_execute_dft(forward, psi, psi); // forward fft

    ofstream fpRe;
    fpRe.open("wave_Re_data.dat");
    if (fpRe.fail()) {
        cout << "cannot open file" << endl;
        return (EXIT_SUCCESS);
    }

    // multiply dij by cosine
    for (i = 0; i < Nx; i++) {
        for (j = 0; j < Ny; j++) {
            k = sqrt(kx[i] * kx[i] + ky[j] * ky[j]);
            psi[j + i * Ny][0] =
                psi[j + i * Ny][0] * cos(2.0 * M_PI * v * k * t);
            psi[j + i * Ny][1] =
                psi[j + i * Ny][1] * cos(2.0 * M_PI * v * k * t);
        }
    }

    fftw_execute_dft(inverse, psi, psi); // inverse fft

    // write new wave to file
    for (j = Ny - 1; j >= 0; j--) {
        for (i = 0; i < Nx; i++) {
            fpRe << setw(15) << psi[j + i * Ny][0] / (Nx * Ny) << setw(15);
        }
        fpRe << endl;
    }
}