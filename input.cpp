#include "input.h"

Input::Input(int _nbeads, int _nsteps, int _nthermal, int _nbins, double _temp, double _max_disp, string _method) :
    nbeads(_nbeads), nsteps(_nsteps), nthermal(_nthermal), nbins(_nbeads), temp(_temp), max_disp(_max_disp), method(_method) {

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
            printf("-n :        total Monte Carlo sweeps\n");
            printf("-w :        thermalizatoin sweeps\n");
            printf("-d :        max step length\n");
            printf("-t :        temperature in Kelvin\n");
            printf("-b :        number of bins for error analysis\n");
            printf("-p :        number of time-slices\n");
            printf("-m :        choose method for estimators (T, PT, CV and PCV)\n");
            printf("-h :        help\n");
            exit(0);
        }
        else
        {
            printf("Type -h for help.\n");
            exit(1);
        }
    }
}