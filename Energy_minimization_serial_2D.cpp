// Include library
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <vector>
using namespace std;

// define required variables
//#define epsilon 0.0103 // epsilon for LJ potential
//#define sigma 3.3      // sigma for LJ potential

#define epsilon 1.0 // epsilon for LJ potential
#define sigma 1.0   // sigma for LJ potential
#define Tol 1e-2    // tolerance for steepest descent method
#define h 1e-4      // step size

/* Define all required function */
// Function to calculate potential energy
double V_lj(double r) {
  double V =
      (4 * epsilon * (pow(sigma, 12) / pow(r, 12) - pow(sigma, 6) / pow(r, 6)));
  return V;
}
// function to calculate first derivative of Potential energy
double d_V_lj_dr(double r) {
  double dV = (-48 * epsilon * pow(sigma, 12) / pow(r, 13) +
               24 * epsilon * pow(sigma, 6) / pow(r, 7));
  return dV;
}

// find mode of given list
double mod(vector<double> vect) {

  double sum = 0;
  for (double x : vect)
    sum = sum + x * x;
  // cout << sum;
  return sqrt(sum);
}

// calculate total energy of system
double Enrgy(const vector<double> &q_s) {
  double sum = 0;
  double r_ij;
  int N = q_s.size();
  int j;
#pragma omp parallel for private(j)
  for (int i = 0; i < N; i++) {
    for (j = 0; j < N; j++) {
      if (i < j) {
        r_ij = abs(q_s[j] - q_s[i]);
        sum = V_lj(r_ij) + sum;
      }
    }
  }
  return sum;
}

// Find total force experiences by i-th atom due to neighbour atoms

vector<double> Force(const vector<double> &q_sx,
                     const vector<double> &q_sy) // ds
{
  double F_x, F_y, F_r;
  double r_ij;
  double a = 1.0; // lattice distance
  double r_cut = 3 * a;
  int N = q_sx.size();
  int M = q_sy.size();
  vector<double> X(N);
  vector<double> Y(N);
  vector<double> Z(N);

  int j;
#pragma omp parallel for private(j)
  for (int i = 0; i < N; i++) {
    F_x = 0;
    F_y = 0;

    for (j = i; j < N; j++) {
      if (i != j) {
        r_ij = abs(sqrt((q_sx[i] - q_sx[j]) * (q_sx[i] - q_sx[j]) +
                        ((q_sy[i] - q_sy[j]) * (q_sy[i] - q_sy[j]))));
        // Include
        // cout << r_ij << endl;
        if (r_ij < r_cut && r_ij != 0) {
          // cout << r_ij << endl;

          F_x = F_x + (-d_V_lj_dr(r_ij) * ((q_sx[i] - q_sx[j]) / r_ij));
          F_y = F_y + (-d_V_lj_dr(r_ij) * ((q_sy[i] - q_sy[j]) / r_ij));
        }
      }
    }
    // Total force experiences by i-th atom due to neighbour atoms
    X[i] = F_x;
    Y[i] = F_y;
    Z[i] = sqrt(F_x * F_x + F_y * F_y);
    // cout << Z[i] << endl;
  }
  return X, Y, Z;
}

/* Main program*/
int main() {
  // Define variables
  int N = 4;
  int M = 4;              // number of atoms
  vector<double> q_sx(N); // atom postion
  vector<double> q_sy(M); // atom postion

  vector<double> f_sx(N); // force on each atom
  vector<double> f_sy(M); // force on each atom
  vector<double> f_R(N);

  vector<double> fx(N);
  vector<double> fy(M);

  // Give starting value for function
  double a = 1.0; // lattice distance
  double b = 1.0; // lattice distance

  double err;

  // Define initial postion of each atom (q_s_0)

  for (int i = 0; i < N; i++) {
    q_sx[i] = i * a;

    q_sy[i] = i * b;
    // Find force vector for N atoms
    f_sx, f_sy, f_R = Force(q_sx, q_sy);
    // calculate error

    err = abs(mod(f_R));
    cout << err << endl;

    // to write output data
    ofstream force;
    force.open("force.txt");

    ofstream energy;
    energy.open("energy.txt");

    ofstream pos;
    pos.open("pos.txt");

    // To keep counting number of iterations
    int count = 0;
    // Gradient descent loop
    while (err > Tol) {
      count += 1;
      // Dirichlet boundary condition
      q_sx[0] = 0;
      q_sy[0] = 0;
      for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {

          q_sx[i] = (q_sx[j] - q_sx[i]) + h * (f_sx[j] - f_sx[i]);
          q_sy[j] = (q_sy[j] - q_sy[i]) + h * (f_sy[j] - f_sy[i]);

          pos << "(" << q_sx[i] << "," << q_sy[j] << ")" << endl;
          force << "(" << f_sx[j] << "," << f_sy[j] << ")" << endl;
        }
        pos << endl;
        force << endl;
        f_sx, f_sy, f_R = Force(q_sx, q_sy);
        err = abs(mod(f_R));
        cout << err << endl;
        vector<double> r;
        for (int i = 0; i < N; i++) {

          for (int j = i; j < M; j++) {
            if (i != j) {
              double r_ij = (sqrt((q_sx[j] - q_sx[i]) * (q_sx[j] - q_sx[i]) +
                                  ((q_sy[j] - q_sy[i]) * (q_sy[j] - q_sy[i]))));
              if (r_ij != 0) {
                r.push_back(r_ij);
                // cout<<r_ij<<endl;
              }
            }
          }
        }
        double P = Enrgy(r);
        energy << P << endl;
      }
      force.close();
      pos.close();
      energy.close();
    }
  }

  return 0;
}
