#include "input.h"

Input::Input(int _nbeads, int _nsteps, int _nthermal, int _nbins, double _temp, double _max_disp, string _method) :
    nbeads(_nbeads), nsteps(_nsteps), nthermal(_nthermal), nbins(_nbins), temp(_temp), max_disp(_max_disp), method(_method) {

    beta = 1.0 / temp;
}

Input::~Input() {
}

void Input::get_input_params(int argc, char** argv) {
    int iarg = 1;
    while (iarg < argc)
    {
        if (string(argv[iarg]) == "-n") {
            nsteps = stoi(argv[iarg + 1]);
            iarg += 2;
        }
        else if (string(argv[iarg]) == "-w") {
            nthermal = stoi(argv[iarg + 1]);
            iarg += 2;
        }
        else if (string(argv[iarg]) == "-p") {
            nbeads = stoi(argv[iarg + 1]);
            iarg += 2;
        }
        else if (string(argv[iarg]) == "-b") {
            nbins = stoi(argv[iarg + 1]);
            iarg += 2;
        }
        else if (string(argv[iarg]) == "-t") {
            temp = stod(argv[iarg + 1]);
            beta = 1.0 / temp;
            iarg += 2;
        }
        else if (string(argv[iarg]) == "-d") {
            max_disp = stod(argv[iarg + 1]);
            iarg += 2;
        }
        else if (string(argv[iarg]) == "-m") {
            method = string(argv[iarg + 1]);
            iarg += 2;
        }
        else if (string(argv[iarg]) == "-h") {
            cout << "-n :        total Monte Carlo sweeps" << endl;
            cout << "-w :        thermalizatoin sweeps" << endl;
            cout << "-d :        max step length" << endl;
            cout << "-t :        temperature in Kelvin" << endl;
            cout << "-b :        number of bins for error analysis" << endl;
            cout << "-p :        number of time-slices" << endl;
            cout << "-m :        choose method for estimators (T, PT, CV and PCV)" << endl;
            cout << "-h :        help" << endl;
            exit(0);
        }
        else
        {
            cout << "Type -h for help." << endl;
            exit(1);
        }
    }
    // for classical monte carlo, estimators are forced to the primitive methods
    if (nbeads == 1) method = "T";
}