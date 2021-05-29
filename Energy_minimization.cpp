#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <utility>
#include <vector>
// Include library

using namespace std;

// define required variables
//#define epsilon 0.0103 // epsilon for LJ potential
//#define sigma 3.3      // sigma for LJ potential

#define epsilon 1.0 // epsilon for LJ potential
#define sigma 1.0   // sigma for LJ potential
#define Tol 1e-2    // tolerance for steepest descent method
#define h 1e-5      // step size
#define a 1.0
#define b 1.0
#define n_rows 10
#define n_cols 10
#define n_atoms 10 * 10

/* Define all required function */
// Function to calculate potential energy
double V_lj(double r) {
  double V =
      (4 * epsilon * (pow(sigma, 12) / pow(r, 12) - pow(sigma, 6) / pow(r, 6)));
  return V;
}

// function to calculate first derivative of Potential energy
void d_V_lj_dr(double r_x, double r_y) {
  double dVx = (-48 * epsilon * pow(sigma, 12) / pow(r_x, 13) +
                24 * epsilon * pow(sigma, 6) / pow(r_x, 7));
  double dVx = (-48 * epsilon * pow(sigma, 12) / pow(r_y, 13) +
                24 * epsilon * pow(sigma, 6) / pow(r_y, 7));
}

/* For energy*/
double Enrgy(const vector<double> &q_s) {
  double sum = 0;
  //   double r_ij;
  //   int N = q_s.size();

  for (int i = 0; i < n_atoms; i++) {
    if (q_s[i] > 0) {
      // for (int j = 0; j < N; j++) {
      //   if (i < j) {
      // r_ij = abs(q_s[j] - q_s[i]);
      sum = V_lj(q_s[i]) + sum;
    }
  }

  return sum;
}

/* Main program*/
int main(int argc, char **argv) {
  /* To build 2 D lattice system with each position as one atomic postion*/
  pair<double, double> apair;          // Pair inside the vector
  vector<pair<double, double>> v_temp; // Temporary vector to store each pair
  vector<vector<pair<double, double>>> pair2dvector; // Required 2D vectors

  vector<double> dist(int(n_rows * n_cols));
  for (int i = 0; i < n_rows; i++) {
    for (int j = 0; j < n_cols; j++) {
      apair.first = a * i;
      apair.second = a * j;
      v_temp.push_back(apair);
    }
    pair2dvector.push_back(v_temp);
    v_temp.clear();
  }
  /* Find r_ij distance for all possible pair of atoms*/

  for (vector<vector<pair<double, double>>>::iterator it = pair2dvector.begin();
       it != pair2dvector.end(); ++it) {
    v_temp = *it;
    for (vector<pair<double, double>>::iterator it2 = v_temp.begin();
         it2 != v_temp.end(); ++it2) {
      apair = *it2;

      double r_ijx = apair.first;
      double r_ijy = apair.second;
      double d_ij = sqrt(r_ijx * r_ijx + r_ijy * r_ijy);
      auto it = dist.insert(dist.begin(), d_ij);

      double err;
      d_V_lj_dr(r_ijx,r_ijy);
      

      //   cout << "(" << r_ijx << "," << r_ijy << ") ; ";
    }
    // cout << endl;
  }

  // int n_atoms = n_rows * n_cols;

  // for (int k = 0; k < n_atoms; k++) {
  //   cout << dist[k] << endl;
  // }

  //   }

  // vector<double> f_s(N); // force on each atom
  // vector<double> f(N);

  // Give starting value for function
  // double a = 1.0; // lattice distance

  // Define initial postion of each atom (q_s_0)
  // for (int i = 0; i < N; i++) {
  //   q_s[i] = i * a;
   f_s = Force(q_s);
  // calc
  // Find force vector for N atoms
  //ulate error
  // err = mod(f_s);

  double P = Enrgy(dist);

  cout << P << endl;
  // // number of atoms
  // vector<double> q_s(N); // atom postion
  // vector<double> f_s(N); // force on each atom
  // vector<double> f(N);

  // // Give starting value for function
  // //double a = 1.0; // lattice distance
  // double err;

  // // Define initial postion of each atom (q_s_0)
  // for (int i = 0; i < N; i++) {
  //   q_s[i] = i * a;

  //   // Find force vector for N atoms
  //   f_s = Force(q_s);
  //   // calculate error
  //   err = mod(f_s);
  //   // cout << err << endl;

  //   // to write output data
  //   ofstream force;
  //   force.open("force.txt");

  //   ofstream energy;
  //   energy.open("energy.txt");

  //   ofstream pos;
  //   pos.open("pos.txt");
  //   // To keep counting number of iterations
  //   int count = 0;
  //   // Gradient descent loop
  //   while (err > Tol) {
  //     count += 1;
  //     for (int i = 0; i < N; i++) {
  //       q_s[i] = q_s[i] + h * f_s[i];
  //       pos << q_s[i] << " ";
  //       force << f_s[i] << " ";
  //     }
  //     pos << endl;
  //     force << endl;
  //     f_s = Force(q_s);
  //     err = mod(f_s);
  //     double P = Enrgy(q_s);
  //     energy << P << endl;
  //   }
  //   force.close();
  //   pos.close();
  //   energy.close();
  // }
  return 0;
}