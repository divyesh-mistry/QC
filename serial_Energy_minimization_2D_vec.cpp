
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
double d_V_lj_dr(double r) {
  double dV;
  if (r != 0) {

    dV = (-48 * epsilon * pow(sigma, 12) / pow(r, 13) +
          24 * epsilon * pow(sigma, 6) / pow(r, 7));
  } else {
    dV = 0;
  }
  return dV;
}

/* For energy*/
double Enrgy(const vector<double> &q_s) {
  double sum = 0;
  //   double r_ij;
  //   int N = q_s.size();

  for (int i = 0; i < n_atoms; i++) {
    if (q_s[i] != 0) {
      // for (int j = 0; j < N; j++) {
      //   if (i < j) {
      // r_ij = abs(q_s[j] - q_s[i]);
      sum = V_lj(q_s[i]) + sum;
    }
  }

  return sum;
}
double mod(vector<double> vect) {
  double sum = 0;
  for (double x : vect)
    sum = sum + x * x;
  return sqrt(sum);
}
int main(int argc, char **argv) {

  pair<double, double> apair;
  vector<pair<double, double>> v_temp;
  vector<vector<pair<double, double>>> pair2dvector;
  //   vector<double> dist(int(N * M));
  for (int i = 0; i < n_rows; i++) {
    for (int j = 0; j < n_cols; j++) {
      apair.first = a * i;
      apair.second = a * j;
      v_temp.push_back(apair);
    }
    pair2dvector.push_back(v_temp);
    v_temp.clear();
  }
  int count_i = 0;
  int count_j = 0;

  double old_ref_x = 0.0;
  double old_ref_y = 0.0;

  /* To calculate force on each atoms*/
  // pair<double, double> fpair;          // Pair inside the vector
  pair<double, double> fpair;
  vector<pair<double, double>> F_I; // Temporary vector to store each pair

  // FOr ENERGY calculations
  vector<double> enrgy;
  vector<double> dist;
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

      d_ij = (sqrt((r_ijx * r_ijx) + (r_ijy * r_ijy)));

      dist.insert(dist.begin(), d_ij);

      //  cout << d_ij << endl;

      double fx = (-d_V_lj_dr(r_ijx));
      double fy = (-d_V_lj_dr(r_ijy));

      fpair.first = fx;
      fpair.second = fy;

      count_j = count_j + 1;
      F_I.push_back(fpair);
    }
    old_ref_x = new_ref_x;
    old_ref_y = new_ref_y;
    enrgy.push_back(Enrgy(dist));
    count_i = count_i + 1;
    for (int k = 0; k < n_atoms; k++) {

      cout << enrgy[k] << endl;
    }
    dist.clear();
  }
  double Fr[n_atoms];
  for (int i = 0; i < n_atoms; i++) {
    // "first" and "second" are used to access
    double fx_i = F_I[i].first;
    double fy_i = F_I[i].second;

    Fr[i] = Fr[i] + sqrt((fx_i * fx_i) + (fy_i * fy_i));

    // 1st and 2nd element of pair respectively
    // cout << "(" << F_I[i].first << ", " << F_I[i].second << ")" << endl;
    cout << Fr[i] << endl;
    cout << enrgy[i] << endl; //<< ", " << F_I[i].second << ")" <<

    // // Gradient descent loop
    // while (err > Tol) {
    //   count += 1;
    //   for (int i = 0; i < N; i++) {
    //     q_s[i] = q_s[i] + h * F[i];
    //     pos << q_s[i] << " ";
    //     force << f_s[i] << " ";
    //   }
    //   pos << endl;
    //   force << endl;
    //   f_s = Force(q_s);
    //   err = mod(f_s);
    //   double P = Enrgy(q_s);
    //   energy << P << endl;
  }

  return 0;
}
