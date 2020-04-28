#ifndef _ESTIMATOR_H_
#define _ESTIMATOR_H_

#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <utility>
#include "path.h"

using namespace std;


class Estimator {
public:
    Estimator(int _nsteps, int _nbins);
    ~Estimator();

    void accumulate(Path _path);
    void output();

private:
    int nbins;  // bin size
    int len;    // number of binned data
    int icount; // counter for data accumulation
    int idx;

    map<string, double> tmp_en;
    map<string, vector<double>> estimators;

    pair<double, double> calc_avg_and_err_jackknife(vector<double> _data);
};

#endif