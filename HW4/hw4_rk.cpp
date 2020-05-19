/* AEP 4380 Homework #4

Test ode solver
4th order Runge-Kutta is tested on a harmonic oscillator and chaotic system

Run on a core i7 using clang 1000.11.45.2 on macOS Mojave

Kevin Juan 3 October 2018
*/

#include <cmath>   // use math package
#include <cstdlib> // plain C

#include <fstream>  // stream file IO
#include <iomanip>  // to format the output
#include <iostream> // stream IO

using namespace std;

void rk4(double yold[], double ynew[], double h, double t, int N,
         void frhs(double[], double, double[])) {
    int i;
    double *k1, *k2, *k3, *k4, *temp;

    k1 = new double[5 * N];
    k2 = k1 + N;
    k3 = k2 + N;
    k4 = k3 + N;
    temp = k4 + N;

    frhs(yold, t, k1);
    for (i = 0; i < N; i++) {
        temp[i] = yold[i] + 0.5 * h * k1[i];
    }

    frhs(temp, t + 0.5 * h, k2);
    for (i = 0; i < N; i++) {
        temp[i] = yold[i] + 0.5 * h * k2[i];
    }

    frhs(temp, t + 0.5 * h, k3);
    for (i = 0; i < N; i++) {
        temp[i] = yold[i] + h * k3[i];
    }

    frhs(temp, t + h, k4);
    for (i = 0; i < N; i++) {
        ynew[i] = yold[i] + h * (k1[i] + 2.0 * (k2[i] + k3[i]) + k4[i]) / 6.0;
    }

    delete[] k1;
    return;
}

void harmonicosc(double y[], double t, double f[]) // harmonic oscillator rhs
{
    f[0] = y[1];  // dy_o/dt
    f[1] = -y[0]; // dy_1/dt
    return;
}

void rhsJC(double y[], double t, double f[]) // Jaynes-Cummings rhs
{
    f[0] = -y[1];                         // dx/dt
    f[1] = y[0] + y[2] * y[3];            // dy/dt
    f[2] = -y[1] * y[3];                  // dz/dt
    f[3] = y[4];                          // de_o/dt
    f[4] = 1.0 * f[1] - 1.0 * 1.0 * y[3]; // de_1/dt
    return;
}

int main() {
    int istep, nstep = 200000, N = 5, i;
    double yold[N], ynew[N], t, tinit = 0.0, h = 0.001;

    // // Runge-Kutta for harmonic oscillator
    // yold[0] = 1.0;
    // yold[1] = 0.0;

    // ofstream fp1;                 // output file using streams
    // fp1.open("harmonic_osc.dat"); // open new file for output
    // fp1.precision(9);             // select 9 digits

    // if (fp1.fail())
    // {
    //     // or fp.bad()
    //     cout << "cannot open file" << endl;
    //     return (EXIT_SUCCESS);
    // }

    // for (istep = 0; istep <= nstep; istep++)
    // {
    //     t = tinit + istep * h;
    //     rk4(yold, ynew, h, t, N, harmonicosc);
    //     fp1 << setw(20) << t << setw(20) << ynew[0] << setw(20) << ynew[1] <<
    //     endl; for (i = 0; i < N; i++)
    //     {
    //         yold[i] = ynew[i];
    //     }
    // }
    // fp1.close();

    // Runge-Kutta for Jaynes-Cummings
    yold[0] = 0.0;      // x
    yold[1] = 0.0;      // y
    yold[2] = 1.0;      // z
    yold[3] = 0.000001; // e
    yold[4] = 0.0;      // e'

    ofstream fp2;                    // output file using streams
    fp2.open("jaynes_cummings.dat"); // open new file for output
    fp2.precision(9);                // select 9 digits

    if (fp2.fail()) {
        // or fp.bad()
        cout << "cannot open file" << endl;
        return (EXIT_SUCCESS);
    }

    for (istep = 0; istep <= nstep; istep++) {
        t = tinit + istep * h;
        rk4(yold, ynew, h, t, N, rhsJC);
        fp2 << setw(10) << t << setw(20) << ynew[0] << setw(20) << ynew[1]
            << setw(20) << ynew[2] << setw(20) << ynew[3] << setw(20) << ynew[4]
            << endl;
        for (i = 0; i < N; i++) {
            yold[i] = ynew[i];
        }
    }
    fp2.close();
    return (EXIT_SUCCESS);
}