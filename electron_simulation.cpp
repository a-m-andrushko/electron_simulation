#include <iostream>
#include <string>
#include <cmath>
#include <fstream>
#include <sstream>
#include <iomanip>

using namespace std;

// PARAMETERS FROM 'parameters.txt'
int N, n, steps;
double kappa, omega, tau, Delta_tau;

struct Particle
{
    double* Psi_R;
    double* Psi_I;
    double* H_R;
    double* H_I;
};

void readParametersFromFile(const string& filename)
{
    ifstream file(filename);
    if (!file)
    {
        cerr << "Error: Unable to open " << filename << "!" << endl;
        cerr << "Default parameters are inserted." << endl;
        N = 100;
        n = 1;
        kappa = 0.0;
        omega = 0.0;
        tau = 0.0;
        Delta_tau = 0.0001;
        steps = 1000;
    }

    file >> N >> n >> kappa >> omega >> tau >> Delta_tau >> steps;

    file.close();

    cout << "Parameters:\n"
            << "N = " << N << "\n"
            << "n = " << n << "\n"
            << "kappa = " << kappa << "\n"
            << "omega = " << omega << "\n"
            << "tau = " << tau << "\n"
            << "Delta_tau = " << Delta_tau << "\n"
            << "steps = " << steps << "\n";
}

void initialiseWaveFunction(Particle* electron)
{
    double Delta_x = 1.0 / N;
    for (int k = 0; k <= N; k++)
    {
        double x_k = k * Delta_x;
        electron -> Psi_R[k] = sqrt(2.0) * sin(n * M_PI * x_k);
        electron -> Psi_I[k] = 0.0;
    }
}

void calculateHamiltonian(Particle* electron)
{
    double Delta_x = 1.0 / N;

    for (int k = 1; k < N; k++)
    {
        double x_k = k * Delta_x;

        double kinetic_R = -0.5 * (electron -> Psi_R[k + 1] + electron -> Psi_R[k - 1] - 2 * electron -> Psi_R[k]) / (Delta_x * Delta_x);
        double kinetic_I = -0.5 * (electron -> Psi_I[k + 1] + electron -> Psi_I[k - 1] - 2 * electron -> Psi_I[k]) / (Delta_x * Delta_x);

        double potential = kappa * (x_k - 0.5) * sin(omega * tau);

        electron -> H_R[k] = kinetic_R + potential * electron -> Psi_R[k];
        electron -> H_I[k] = kinetic_I + potential * electron -> Psi_I[k];
    }

    electron -> H_R[0] = electron -> H_R[N] = electron -> H_I[0] = electron -> H_I[N] = 0.0;
}

void integrate(Particle* electron)
{
    for (int k = 1; k < N; k++)
        electron -> Psi_R[k] += electron -> H_I[k] * (Delta_tau / 2.0);

    calculateHamiltonian(electron);

    for (int k = 1; k < N; k++)
        electron -> Psi_I[k] -= electron -> H_R[k] * Delta_tau;

    calculateHamiltonian(electron);

    for (int k = 1; k < N; k++)
        electron -> Psi_R[k] += electron -> H_I[k] * (Delta_tau / 2.0);
}

void CalculateQuantities(Particle* electron, double* norm, double* x_expect, double* energy)
{
    double Delta_x = 1.0 / N;
    *norm = 0.0;
    *x_expect = 0.0;
    *energy = 0.0;

    for (int k = 0; k <= N; k += 2)
    {
        double x_k = k * Delta_x;
        double Psi_sq = (electron -> Psi_R[k] * electron -> Psi_R[k]) + (electron -> Psi_I[k] * electron -> Psi_I[k]);

        *norm += Psi_sq;
        *x_expect += x_k * Psi_sq;
        *energy += electron -> Psi_R[k] * electron -> H_R[k] + electron -> Psi_I[k] * electron -> H_I[k];
    }

    *norm *= Delta_x;
    *x_expect *= Delta_x;
    *energy *= Delta_x;
}

int main()
{
    readParametersFromFile("parameters.txt");

    // DYNAMIC MEMORY ALLOCATION
    Particle electron;
    electron.Psi_R = new double[N + 1];
    electron.Psi_I = new double[N + 1];
    electron.H_R = new double[N + 1];
    electron.H_I = new double[N + 1];

    initialiseWaveFunction(&electron);

    // FILES TO STORE DATA
    ofstream outfile1("wave_function.txt");
    if (!outfile1)
        cerr << "Unable to open wave_function.txt!";

    ofstream outfile2("quantities.txt");
    if (!outfile2)
        cerr << "Unable to open quantities.txt!";

    ofstream outfile3("density.txt");
    if (!outfile3)
        cerr << "Unable to open quantities.txt!";

    outfile1 << fixed << setprecision(5);

    outfile2 << fixed << setprecision(5);

    outfile3 << fixed << setprecision(5);

    double norm = 0.0;
    double x_expect = 0.0;
    double energy = 0.0;
    double rho = 0.0;

    // EVOLUTION LOOP
    for (int step = 0; step < steps; ++step)
    {
        calculateHamiltonian(&electron);

        integrate(&electron);

        CalculateQuantities(&electron, &norm, &x_expect, &energy);        

        for (int k = 0; k <= N; k++)
        {
            outfile1 << setw(12) << k << setw(12) << electron.Psi_R[k] << setw(12) << electron.Psi_I[k] << "\n";

            if (k % 2 == 0)
            {
                rho = (electron.Psi_R[k] * electron.Psi_R[k]) + (electron.Psi_I[k] * electron.Psi_I[k]);
                outfile3 << k << setw(12) << rho << "\n";
            }
        }

        outfile1 << endl << endl;

        outfile3 << endl << endl;

        outfile2 << setw(12) << tau << setw(12) << norm << setw(12) << x_expect << setw(12) << energy << "\n";

        tau += Delta_tau;
    }

    outfile1.close();
    outfile2.close();
    outfile3.close();

    delete[] electron.Psi_R;
    delete[] electron.Psi_I;
    delete[] electron.H_R;
    delete[] electron.H_I;

    return 0;
}
