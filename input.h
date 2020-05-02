#ifndef _INPUT_H_
#define _INPUT_H_

#include <iostream>
#include <string>

using std::string;
using std::stoi;
using std::stod;


class Input {
public:
    int nbeads;
    int nsteps;
    int nthermal;
    int nbins;
    double temp;
    double beta;
    double max_disp;
    string method;

    Input(int _nbeads = 16, int _nsteps = 1e6, int _nthermal = 2e5, int _nbins = 1e4,
        double _temp = 0.5, double _max_disp = 0.5, string _method = "T");
    ~Input();

    void get_input_params(int argc, char** argv);
};

#endif
