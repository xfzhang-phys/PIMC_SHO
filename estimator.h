#ifndef _ESTIMATOR_H_
#define _ESTIMATOR_H_

#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <utility>
#include "path.h"

using std::string;
using std::vector;
using std::pair;
using std::map;


class Estimator {
public:
    Estimator(int _nsteps, int _nbins, int _nbeads, string _method = "T");
    ~Estimator();

    void accumulate(Path _path);
    void accumulate_ctau(Path _path);
    void output();

private:
    int nbins;  // bin size
    int len;    // number of binned data
    int icount; // counter for data accumulation
    int ccount; // counter for correlation function accumulation
    int idx;
    int cidx;
    // choose estimators: T (thermodynamic), CV (centroid viral) and their projected version PT, PCV
    string method;

    map<string, double> accmltr;    // accumulator
    map<string, vector<double>> estimators;

    vector<double> accmltr_ctau;    // accumulator for imaginary-time correlation function
    vector<vector<double>> ctau;    // imaginary-time correlation function

    double calc_pot(double _pos);
    double calc_kcv(Path _path);
    pair<double, double> calc_avg_and_err_jackknife(vector<double> _data);
};

#endif