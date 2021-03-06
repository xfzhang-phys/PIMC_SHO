#ifndef _PATH_H_
#define _PATH_H_

#include <iostream>
#include <vector>
#include <random>
#include <string>
#include <map>

using std::mt19937_64;
using std::uniform_int_distribution;
using std::uniform_real_distribution;
using std::normal_distribution;
using std::vector;


class Path {
public:
    int nbeads;     // number of time-slices
    double beta;    // inversion temperature
    double pos_com; // position of center of mass
    std::vector<double> pos;
    std::map<std::string, double> en;

    Path(int _nbeads, double _beta);
    ~Path();

    void move_bisection(mt19937_64* _gen);
    void move_com(double _max_disp, mt19937_64* _gen);
    void update_com_max_disp(double* _max_disp);


private:
    int lmax;
    int accepted_com_trials;
    int all_com_trials;

    void calc_toten();
    double calc_kin();
    double calc_pot(double _pos);
};

#endif

