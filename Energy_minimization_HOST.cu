// Include library
#include "cuda_runtime.h"
#include <fstream>
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <time.h>
#include <vector>
using namespace std;

// define required variables
//#define epsilon 0.0103 // epsilon for LJ potential
//#define sigma 3.3      // sigma for LJ potential

#define epsilon 1.0 // epsilon for LJ potential
#define sigma 1.0   // sigma for LJ potential
#define Tol 1e-1    // tolerance for steepest descent method
#define h 1e-5      // step size
#define N_atoms 200

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
  return sqrt(sum);
}

// calculate total energy of system
__global__ void Enrgy(const vector<double> &q_s) {
  // Write CUDA kernal here
  int thread = threadIdx.x;
  int block = blockIdx.x;
  int dimen = blockDim.x;

  int i = thread + block * dimen;

  double sum = 0;
  double r_ij;
  //   int N_atoms = q_s.size();
  int j;

  if ((block != 0) && (block != 201)) {
    if ((thread != 0) && (thread != 201)) {
      {
        r_ij = abs(q_s[j] - q_s[i]);
        sum = V_lj(r_ij) + sum;
      }
    }
  }
  return sum;
}

// Find total force experiences by i-th atom due to neighbour atoms

__global__ void Force(const vector<double> &q_s) // ds
{
  double F;
  double r_ij;
  double a = 1.0; // lattice distance
  double r_cut = 3 * a;
  //   int N_atoms = q_s.size();
  vector<double> X(N_atoms);
  int thread = threadIdx.x;
  int block = blockIdx.x;
  int dimen = blockDim.x;

  int i = thread + block * dimen;
  int j = thread + (block - 1) * dimen;
  if ((block != 0) && (block != 201)) {
    F = 0;
    if ((thread != 0) && (thread != 201)) {
      {

        r_ij = abs(q_s[j] - q_s[i]);
        // Include cutoff
        if (r_ij < r_cut) {
          F = F + (-d_V_lj_dr(r_ij) * ((q_s[i] - q_s[j]) / r_ij));
        }
      }
      // Total force experiences by i-th atom due to neighbour atoms
      X[i] = F;
      // cout << X[i] << endl;
    }
    return X;
  }

  /* Main program*/
  int main(int argc, char **argv) {
    // Define variables
    // int N_atoms = 31;            // number of atoms
    vector<double> q_s(N_atoms); // atom postion
    vector<double> f_s(N_atoms); // force on each atom
    vector<double> f(N_atoms);

    // Give starting value for function
    double a = 1.0; // lattice distance
    double err;
    float dt;
    clock_t start, end;
    int i, j, iter;
    double elapsed;
    double *dev_q_s, *dev_f_s;
    int size = (N_atoms + 1) * sizeof(float);

    start = clock();
    cudaMalloc((void **)&dev_q_s, size);
    cudaMalloc((void **)&dev_f_s, size);
    // Define initial postion of each atom (q_s_0)

    int j;
    for (int i = 0; i < N_atoms; i++) {
      q_s[i] = i * a;

      // Find force vector for N_atoms atoms
      Energy<<<200, 1>>>(dev_q_s);
      Force<<<200, 1>>>(dev_q_s);
      f_s = Force(q_s);
      // calculate error
      err = mod(f_s);
      // cout << err << endl;
      cudaMemcpy(dev_told, &q_s, size, cudaMemcpyHostToDevice);
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
        for (j = 0; j < N_atoms; j++) {

          q_s[j] = q_s[j] + h * f_s[j];

          // Dirichlet boundary condition
          q_s[0] = 0;
          pos << q_s[j] << " ";
          force << f_s[j] << " ";
        }
        pos << endl;
        force << endl;
        Energy<<<200, 1>>>(dev_q_s);
        Force<<<200, 1>>>(dev_q_s);
        err = mod(f_s);
        //   double P = Enrgy(q_s);
        //   energy << P << endl;
      }
      force.close();
      pos.close();
      energy.close();
    }
    end = clock();
    elapsed = ((double)(end - start)) / CLOCKS_PER_SEC;
    printf(" \n Time taken is %f \n", elapsed);

    return 0;
  }
