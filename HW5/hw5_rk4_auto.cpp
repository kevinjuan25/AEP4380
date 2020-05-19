/* AEP 4380 Homework #5

Test automatic step size ode solver
Runge-Kutta using automatic step sizing is tested on a harmonic oscillator and
three body problem

Run on a core i7 using clang 1000.11.45.2 on macOS Mojave

Kevin Juan 12 October 2018
*/

#include <cmath>   // use math package
#include <cstdlib> // plain C

#include <fstream>  // stream file IO
#include <iomanip>  // to format the output
#include <iostream> // stream IO

using namespace std;

void rkAdapt(double yold[], double ynew[], double h, int nEqns,
             void frhs(double[], double, double[]), ofstream &fp, double tMin,
             double tMax, double tol) {
    int i;
    double *k1, *k2, *k3, *k4, *k5, *k6, *temp, *scale, *yerr, *ynewSt, *delta;
    const static double c2 = 0.2, c3 = 0.3, c4 = 0.8, c5 = 8.0 / 9.0, c6 = 1.0;
    const static double a21 = 0.2, a31 = 3.0 / 40.0, a32 = 9.0 / 40.0;
    const static double a41 = 44.0 / 45.0, a42 = -56.0 / 15.0, a43 = 32.0 / 9.0;
    const static double a51 = 19372.0 / 6561.0, a52 = -25360.0 / 2187.0,
                        a53 = 64448.0 / 6561.0, a54 = -212.0 / 729.0;
    const static double a61 = 9017.0 / 3168.0, a62 = -355.0 / 33.0,
                        a63 = 46732.0 / 5247.0, a64 = 49.0 / 176.0,
                        a65 = -5103.0 / 18656.0;
    const static double b1 = 35.0 / 384.0, b2 = 0.0, b3 = 500.0 / 1113.0,
                        b4 = 125.0 / 192.0, b5 = -2187.0 / 6784.0,
                        b6 = 11.0 / 84.0;
    const static double bSt1 = 5179.0 / 57600.0, bSt2 = 0.0,
                        bSt3 = 7571.0 / 16695.0, bSt4 = 393.0 / 640.0,
                        bSt5 = -92097.0 / 339200.0, bSt6 = 187.0 / 2100.0,
                        bSt7 = 1.0 / 40.0;
    double deltaMax, t = tMin;
    const static double smallVal = 1.0e-30;

    k1 = new double[11 * nEqns];
    k2 = k1 + nEqns;
    k3 = k2 + nEqns;
    k4 = k3 + nEqns;
    k5 = k4 + nEqns;
    k6 = k5 + nEqns;
    temp = k6 + nEqns;
    ynewSt = temp + nEqns;
    scale = ynewSt + nEqns;
    yerr = scale + nEqns;
    delta = yerr + nEqns;

    // output initial values to file
    fp << setw(17) << t << setw(17) << h;
    for (i = 0; i < nEqns; i++) {
        fp << setw(17) << yold[i];
    }
    fp << endl;

    while (t <= tMax) {
        frhs(yold, t, k1);
        for (i = 0; i < nEqns; i++) {
            temp[i] = yold[i] + h * a21 * k1[i];
        }

        frhs(temp, t + c2 * h, k2);
        for (i = 0; i < nEqns; i++) {
            temp[i] = yold[i] + h * (a31 * k1[i] + a32 * k2[i]);
        }

        frhs(temp, t + c3 * h, k3);
        for (i = 0; i < nEqns; i++) {
            temp[i] = yold[i] + h * (a41 * k1[i] + a42 * k2[i] + a43 * k3[i]);
        }

        frhs(temp, t + c4 * h, k4);
        for (i = 0; i < nEqns; i++) {
            temp[i] = yold[i] + h * (a51 * k1[i] + a52 * k2[i] + a53 * k3[i] +
                                     a54 * k4[i]);
        }

        frhs(temp, t + c5 * h, k5);
        for (i = 0; i < nEqns; i++) {
            temp[i] = yold[i] + h * (a61 * k1[i] + a62 * k2[i] + a63 * k3[i] +
                                     a64 * k4[i] + a65 * k5[i]);
        }

        frhs(temp, t + c6 * h, k6);
        for (i = 0; i < nEqns; i++) {
            ynew[i] = yold[i] + h * (b1 * k1[i] + b2 * k2[i] + b3 * k3[i] +
                                     b4 * k4[i] + b5 * k5[i] + b6 * k6[i]);
        }

        frhs(ynew, t + h, temp);
        for (i = 0; i < nEqns; i++) {
            ynewSt[i] =
                yold[i] +
                h * (bSt1 * k1[i] + bSt2 * k2[i] + bSt3 * k3[i] + bSt4 * k4[i] +
                     bSt5 * k5[i] + bSt6 * k6[i] + bSt7 * temp[i]);
        }

        deltaMax = -1000.0; // initialize an impossible deltaMax
        for (i = 0; i < nEqns; i++) {
            yerr[i] = abs(ynew[i] - ynewSt[i]);
            scale[i] = abs(ynew[i]) + abs(h * k1[i]) + smallVal;
            delta[i] = yerr[i] / scale[i];
            if (delta[i] > deltaMax) {
                deltaMax = delta[i];
            }
        }

        if (deltaMax < tol) {
            t += h;
            fp << setw(17) << t << setw(17) << h;
            for (i = 0; i < nEqns; i++) {
                fp << setw(17) << ynew[i];
                yold[i] = ynew[i];
            }
            fp << endl;
            if (deltaMax < tol / 2.0) {
                h = h * pow(abs(tol / deltaMax), 1.0 / 5.0);
            }
        } else {
            h = h / 5.0;
        }
        cout << "t: " << t << endl;
    }
    delete[] k1;
    return;
}

void harmonicosc(double y[], double t, double f[]) // harmonic oscillator rhs
{
    f[0] = y[1];  // dy_0/dt
    f[1] = -y[0]; // dy_1/dt
    return;
}

void planets(double y[], double t, double f[]) {
    int N = 3, xo = 0, yo = N, vxo = 2 * N, vyo = 3 * N, i, j;
    double distCube, dx, dy, dist;
    const static double G = 6.674e-11;
    const static double mass[] = {5.976e24, 0.0123 * 5.976e24,
                                  0.2 * 0.0123 * 5.976e24}; // mass in kg
    const static double radii[] = {6378e3, 3476e3,
                                   0.5 * 3476e3}; // radii in meters

    for (i = 0; i < N; i++) {
        f[i + vxo] = 0.0; // initialize rhs for dv_x/dt
        f[i + vyo] = 0.0; // initialize rhs for dv_y/dt
        for (j = 0; j < N; j++) {
            if (i != j) {
                dx = y[j + xo] - y[i + xo];
                dy = y[j + yo] - y[i + yo];
                dist = sqrt(dx * dx + dy * dy);
                distCube = dist * dist * dist;
                if (dist <= radii[i] + radii[j]) {
                    cout << "A collision has occurred" << endl;
                    exit(1);
                }
                f[i + vxo] += G * mass[j] * dx / distCube; // dv_x/dt
                f[i + vyo] += G * mass[j] * dy / distCube; // dv_y/dt
            }
        }
        f[i + xo] = y[i + vxo]; // dx_i/dt
        f[i + yo] = y[i + vyo]; // dy_i/dt
    }
    return;
}

int main() {
    int nEqns = 12;
    double yold[nEqns], ynew[nEqns], h = 0.001;
    const static double convertDay =
        60.0 * 60.0 * 24.0; // day to second conversion

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

    // rkAdapt(yold, ynew, h, nEqns, harmonicosc, fp1, 0.0, 100.0, 10e-7);

    // fp1.close();

    // Runge-Kutta for 3-body motion
    yold[0] = 0.0;     // x_earth
    yold[1] = 0.0;     // x_moon
    yold[2] = -4.97e8; // x_moon-2
    yold[3] = 0.0;     // y_earth
    yold[4] = 3.84e8;  // y_moon
    yold[5] = 0.0;     // y_moon2
    yold[6] = -12.593; // vx_eartg
    yold[7] = 1019.0;  // vx_moon
    yold[8] = 965.0;   // vx_moon-2
    yold[9] = 0.0;     // vy_earth
    yold[10] = 0.0;    // vy_moon
    yold[11] = 820.0;  // vy_moon-2

    ofstream fp2;
    fp2.open("3_body_planets.dat");
    fp2.precision(9);

    if (fp2.fail()) {
        cout << "cannot open file" << endl;
        return (EXIT_SUCCESS);
    }

    rkAdapt(yold, ynew, h, nEqns, planets, fp2, 0.0, 200.0 * convertDay,
            10e-10);

    fp2.close();
}