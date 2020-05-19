/* AEP 4380 Homework #1

Test numerical derivatives
Forward, backward, and central difference methods are tests

Run on a core i7 using clang 902.0.39.2

Kevin Juan 2 September 2018
*/

#include <cmath>
#include <cstdlib> // plain C

#include <fstream>  // stream file IO
#include <iomanip>  // to format the output
#include <iostream> // stream IO

using namespace std;

int main() {
  int i, n = 200;
  double h = 0.5, xmin = -7.0, xmax = +7.0, x, dx, f1, fpfd, fpcd, fpbd;
  double feval(double);

  ofstream fp; // output file using streams

  fp.open("hw1_fd.dat"); // open new file for output
  if (fp.fail()) {       // or fp.bad()
    cout << "cannot open file" << endl;
    return (EXIT_SUCCESS);
  }

  dx = (xmax - xmin) / (n - 1);
  for (i = 0; i < n; i++) {
    x = xmin + i * dx;
    f1 = feval(x);
    fpbd = (f1 - feval(x - h)) / h;
    fpfd = (feval(x + h) - f1) / h;
    fpcd = (feval(x + h) - feval(x - h)) / (2 * h);

    // data file for python
    fp << setw(15) << x << setw(15) << f1 << setw(15) << fpbd << setw(15)
       << fpfd << setw(15) << fpcd << endl;

    // print to screen
    cout << setw(15) << x << setw(15) << f1 << setw(15) << fpbd << setw(15)
         << fpfd << setw(15) << fpcd << endl;
  }

  fp.close();

  return (EXIT_SUCCESS);
}

double feval(double x) { return (sin(x) * exp(-0.04 * x * x)); }
