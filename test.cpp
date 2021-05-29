
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
  double dVy = (-48 * epsilon * pow(sigma, 12) / pow(r_y, 13) +
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
int main(int argc, char **argv) {
  int M, N;
  M = 10;
  N = 10;

  pair<double, double> apair;
  vector<pair<double, double>> v_temp;
  vector<vector<pair<double, double>>> pair2dvector;
  //   vector<double> dist(int(N * M));
  for (int i = 0; i < M; i++) {
    for (int j = 0; j < N; j++) {
      apair.first = a * i;
      apair.second = a * j;
      v_temp.push_back(apair);
    }
    pair2dvector.push_back(v_temp);
    v_temp.clear();
  }
  int count_i = 0;
  int count_j = 0;
  vector<double> dist(int(n_rows * n_cols));
  double old_ref_x = 0.0;
  double old_ref_y = 0.0;
  for (vector<vector<pair<double, double>>>::iterator it = pair2dvector.begin();
       it != pair2dvector.end(); ++it) {
    v_temp = *it;
    // R_IJ[count_i][count_j]

    double new_ref_y, new_ref_x, d_ij;
    for (vector<pair<double, double>>::iterator it2 = v_temp.begin();
         it2 != v_temp.end(); ++it2) {
      apair = *it2;
      // cout << "(" << count_j << "," << count_i << ")" << endl;
      new_ref_x = apair.first;
      new_ref_y = apair.second;

      double r_ijx = new_ref_x - old_ref_x;
      double r_ijy = new_ref_y - old_ref_y;


      d_ij = sqrt(r_ijx * r_ijx + r_ijy * r_ijy);
      auto it = dist.insert(dist.begin(), d_ij);
      cout << d_ij << endl;

      // auto it = dist.insert(dist.begin(), d_ij);
      //   cout << "(" << r_ijx << "," << r_ijy << ") ; ";
      count_j = count_j + 1;
    }
    old_ref_x = new_ref_x;
    old_ref_y = new_ref_y;
    count_i = count_i + 1;
    cout << endl;
  }
  for (int k = 0; k < 100; k++) {
    cout << dist[k] << "          ";
    cout << endl;
  }
  // cout << count_i << endl;
  return 0;
}

/* Define all required function */
// // Function to calculate potential energy
// double V_lj(double r) {
//   double V =
//       (4 * epsilon * (pow(sigma, 12) / pow(r, 12) - pow(sigma, 6) / pow(r,
//       6)));
//   return V;
// }
// // function to calculate first derivative of Potential energy
// double d_V_lj_dr(double r) {
//   double dV = (-48 * epsilon * pow(sigma, 12) / pow(r, 13) +
//                24 * epsilon * pow(sigma, 6) / pow(r, 7));
//   return dV;
// }

// // find mode of given list
// double mod(vector<double> vect) {
//   double sum = 0;
//   for (double x : vect)
//     sum = sum + x * x;
//   return sqrt(sum);
// }

// calculate total energy of system
// double Enrgy(const vector<double> &q_s) {
//   double sum = 0;
//   double r_ij;
//   int N = q_s.size();

//   for (int i = 0; i < N; i++) {
//     for (int j = 0; j < N; j++) {
//       if (i < j) {
//         r_ij = abs(q_s[j] - q_s[i]);
//         sum = V_lj(r_ij) + sum;
//       }
//     }
//   }
//   return sum;
// }

// // Find total force experiences by i-th atom due to neighbour atoms

// vector<double> Force(const vector<double> &q_s) // ds
// {
//   double F;
//   double r_ij;
//   //  double a = 1.0; // lattice distance
//   double r_cut = 3 * a;
//   int N = q_s.size();
//   vector<double> X(N);
//   for (int i = 0; i < N; i++) {
//     F = 0;
//     for (int j = 0; j < N; j++) {
//       if (i != j) {
//         r_ij = abs(q_s[j] - q_s[i]);
//         // Include cutoff
//         if (r_ij < r_cut) {
//           F = F + (-d_V_lj_dr(r_ij) * ((q_s[i] - q_s[j]) / r_ij));
//         }
//       }
//     }
//     // Total force experiences by i-th atom due to neighbour atoms
//     X[i] = F;
//     // cout << X[i] << endl;
//   }
//   return X;
// }
