#include <iostream>
#include <string>
//#include <chrono>
#include "path.h"
#include "estimator.h"
#include "input.h"

//using namespace std::chrono;

int main(int argc, char** argv) {

    // timer
    //auto start = steady_clock::now();

    Input input_params(32, 1e6, 2e5, 1e4, 0.5, 0.5, "T");
    // read parameters from command line
    if (argc > 1) {
        input_params.get_input_params(argc, argv);
    }

    // random number generator
    random_device r;
    mt19937_64 gen(r());

    Path path(input_params.nbeads, input_params.beta);
    Estimator estimator(input_params.nsteps, input_params.nbins, input_params.method);

    // thermalization
    for (int istep = 0; istep < input_params.nthermal; istep++) {
        if (input_params.nbeads > 1) {
            path.move_bisection(&gen);
        }
        path.move_com(input_params.max_disp, &gen);

        // auto update max_disp
        if (istep % 100 == 0 && istep > 0) {
            path.update_com_max_disp(&(input_params.max_disp));
        }
    }

    // statistics
    for (int istep = 0; istep < input_params.nsteps; istep++) {
        if (input_params.nbeads > 1) {
            path.move_bisection(&gen);
        }
        path.move_com(input_params.max_disp, &gen);

        // accumulator
        estimator.accumulate(path);

        // auto update max_disp
        if (istep % 100 == 0 && istep > 0) {
            path.update_com_max_disp(&(input_params.max_disp));
        }
    }

    estimator.output();

    // timer
    //auto end = steady_clock::now();
    //auto elapsed_time = duration_cast<seconds> (end - start);
    //cout << "Elapsed time: " << elapsed_time.count() << "s." << endl;
    int tmp;
    cin >> tmp;
}