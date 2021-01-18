// Include library
#include <cstdio>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
using namespace std;

// define required variables
//#define epsilon 0.0103 // epsilon for LJ potential
//#define sigma 3.3      // sigma for LJ potential

#define epsilon 1.0 // epsilon for LJ potential
#define sigma 1.0   // sigma for LJ potential
#define Tol 1e-2    // tolerance for steepest descent method
#define h 1e-5      // step size

/* Define all required function */
// Function to calculate potential energy
double V_lj(double r)
{
    double V = (4 * epsilon * (pow(sigma, 12) / pow(r, 12) - pow(sigma, 6) / pow(r, 6)));
    return V;
}
// function to calculate first derivative of Potential energy
double d_V_lj_dr(double r)
{
    double dV = (-48 * epsilon * pow(sigma, 12) / pow(r, 13) + 24 * epsilon * pow(sigma, 6) / pow(r, 7));
    return dV;
}

//find mode of given list
double mod(vector<double> vect)
{
    double sum = 0;
    for (double x : vect)
        sum = sum + x * x;
    return sqrt(sum);
}

// calculate total energy of system
double Enrgy(vector<double> q_s)
{
    double sum = 0;
    double r_ij;
    int N = q_s.size();

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            if (i < j)
            {
                r_ij = abs(q_s[j] - q_s[i]);
                sum = V_lj(r_ij) + sum;
            }
        }
    }
    return sum;
}

// Find total force experiences by i-th atom due to neighbour atoms

vector<double> Force(const vector<double> &q_s) // ds
{
    double F;
    double r_ij;
    double a = 1.0; // lattice distance
    double r_cut = 3 * a;
    int N = q_s.size();
    vector<double> X(N);
    for (int i = 0; i < N; i++)
    {
        F = 0;
        for (int j = 0; j < N; j++)
        {
            if (i != j)
            {
                r_ij = abs(q_s[j] - q_s[i]);
                // Include cutoff
                if (r_ij < r_cut)
                {
                    F = F + (-d_V_lj_dr(r_ij) * ((q_s[i] - q_s[j]) / r_ij));
                }
            }
        }
        //Total force experiences by i-th atom due to neighbour atoms
        X[i] = F;
        //cout << X[i] << endl;
    }
    return X;
}

/* Main program*/
int main()
{
    // Define variables
    int N = 11;            // number of atoms
    vector<double> q_s(N); // atom postion
    vector<double> f_s(N); // force on each atom
    vector<double> f(N);

    // Give starting value for function
    double a = 1.0; // lattice distance
    double err;

    // Define initial postion of each atom (q_s_0)
    for (int i = 0; i < N; i++)
    {
        q_s[i] = i * a;

        // Find force vector for N atoms
        f_s = Force(q_s);
        // calculate error
        err = mod(f_s);
        //cout << err << endl;

        //to write output data
        ofstream force;
        force.open("force.txt");

        ofstream energy;
        energy.open("energy.txt");

        ofstream pos;
        pos.open("pos.txt");
        // To keep counting number of iterations
        int count = 0;
        // Gradient descent loop
        while (err > Tol)
        {
            count += 1;
            for (int i = 0; i < N; i++)
            {
                q_s[i] = q_s[i] + h * f_s[i];
                pos << q_s[i] << " ";
                force << f_s[i] << " ";
            }
            pos << endl;
            force << endl;
            f_s = Force(q_s);
            err = mod(f_s);
            double P = Enrgy(q_s);
            energy << P << endl;
        }
        force.close();
        pos.close();
        energy.close();
    }
    return 0;
}