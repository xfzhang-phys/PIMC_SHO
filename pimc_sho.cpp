#include <iostream>
#include <string>
#include "path.h"
#include "estimator.h"


int main(int argc, char** argv) {
    int nbeads = 16;
    int nsteps = 1e6;
    int nthermal = 2e5;
    int nbins = 1e4;
    double temp = 0.5;
    double beta = 1.0 / temp;
    double max_disp = 0.5;
    string method = "PCV";

    // read parameters from command line
    if (argc > 1) {
        cout << "You should type arguments like: pimc.x [temp] [nbeads] [nsteps] [nthermal] [method]" << endl;
        temp = stof(string(argv[1])); beta = 1.0 / temp;
        nbeads = stoi(string(argv[2]));
        nsteps = stoi(string(argv[3]));  nbins = nsteps / 100; 
        nthermal = stoi(string(argv[4]));
        method = string(argv[5]);
    }

    // random number generator
    random_device r;
    mt19937_64 gen(r());

    Path path(nbeads, beta);
    Estimator estimator(nsteps, nbins, method);

    // thermalization
    for (int istep = 0; istep < nthermal; istep++) {
        if (nbeads > 1) {
            path.move_bisection(&gen);
        }
        path.move_com(max_disp, &gen);

        // auto update max_disp
        if (istep % 100 == 0 && istep > 0) {
            path.update_com_max_disp(&max_disp);
        }
    }

    // statistics
    for (int istep = 0; istep < nsteps; istep++) {
        if (nbeads > 1) {
            path.move_bisection(&gen);
        }
        path.move_com(max_disp, &gen);

        // accumulator
        estimator.accumulate(path);

        // auto update max_disp
        if (istep % 100 == 0 && istep > 0) {
            path.update_com_max_disp(&max_disp);
        }
    }

    estimator.output();
}