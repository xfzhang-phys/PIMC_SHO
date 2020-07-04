#include <iostream>
#include <string>
#include <chrono>
#include "path.h"
#include "estimator.h"
#include "input.h"
#include "histogram.h"

using namespace std::chrono;
using std::random_device;


int main(int argc, char** argv) {
    
    // timer
    auto start = steady_clock::now();

    Input input_params(64, 5e6, 1e6, 1e4, 1.0, 0.05, -2.0, 2.0, 0.01, "PCV");
    // read parameters from command line
    if (argc > 1) {
        input_params.get_input_params(argc, argv);
    }
    // random number generator
    random_device r;
    mt19937_64 gen(r());

    Path path(input_params.nbeads, input_params.beta);
    Estimator estimator(input_params.nsteps, input_params.nbins, input_params.nbeads, input_params.method);
    Histogram histogram(input_params.hmin, input_params.hmax, input_params.hstep); 

    // thermalization
    for (int istep = 0; istep < input_params.nthermal; istep++) {
        if (input_params.nbeads > 1) {
            path.move_bisection(&gen);
        }
        // path.move_com(input_params.max_disp, &gen);

        // auto update max_disp
        if (istep % 100 == 0 && istep > 0) {
            // path.update_com_max_disp(&(input_params.max_disp));
        }
    }

    // statistics
    for (int istep = 0; istep < input_params.nsteps; istep++) {
        if (input_params.nbeads > 1) {
            path.move_bisection(&gen);
        }
        // path.move_com(input_params.max_disp, &gen);

        // accumulator
        estimator.accumulate(path);
        estimator.accumulate_ctau(path);

        // histogram
        histogram.get_hist(path);

        // auto update max_disp
        if (istep % 100 == 0 && istep > 0) {
            // path.update_com_max_disp(&(input_params.max_disp));
        }
    }

    estimator.output();
    histogram.output();
    
    // timer
    auto end = steady_clock::now();
    auto elapsed_time = duration<double> (end - start);
    std::cout << "Elapsed time: " << elapsed_time.count() << " s." << std::endl;
}