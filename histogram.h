#ifndef _HISTOGRAM_H_
#define _HISTOGRAM_H_

#include <string>
#include <map>
#include <vector>
#include "path.h"

using std::string;
using std::map;
using std::vector;

class Histogram {
public:
    double hmin;
    double hmax;
    double hstep;

    map<string, vector<double>> hist;
    map<string, vector<double>> bins;

    Histogram(double _hmin, double _hmax, double _hstep);
    ~Histogram();

    void get_hist(Path _path);
    void output();

private:
    int hbins;
    double hcount;
};

#endif
