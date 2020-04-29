#include "path.h"


Path::Path(int _nbeads, double _beta) : nbeads(_nbeads), beta(_beta) {
   
    // lmax setting for bisection
    if (nbeads >= 4 && nbeads < 8) lmax = 2;
    else if (nbeads >= 2 && nbeads < 4) lmax = 1;
    else lmax = 3;
    // initialize configuration
    pos.resize(nbeads);
    for (int ibead = 0; ibead < nbeads; ibead++) {
        pos.at(ibead) = 0.0;
        pos_com += pos.at(ibead);
    }
    pos_com /= nbeads;
    // energy estimators
    en["esprng"] = 0.0; en["ekin"] = 0.0; en["epot"] = 0.0; en["toten"] = 0.0;
    calc_toten();
    // initialize max_disp updator
    accepted_com_trials = 0;
    all_com_trials = 0;
}

Path::~Path() {
}

void Path::move_bisection(mt19937_64* gen) {
    int acc = 1;
    double en_dpot = 0.0;
    double slt = sqrt(beta / nbeads);
    vector<double> tmp_pos;

    // random distribution
    uniform_int_distribution<int> uniform_dist_int(0, nbeads-1);
    uniform_real_distribution<double> uniform_dist(0, 1);
    normal_distribution<double> normal_dist(0, 1);

    tmp_pos.resize(nbeads);
    for (int ibead = 0; ibead < nbeads; ibead++) {
        tmp_pos.at(ibead) = pos.at(ibead);
    }

    int bs_ibead = uniform_dist_int(*gen);
    // Loop for bisection level
    double en_dpot_prev = 0.0;
    for (int il = lmax; il > 0; il--) {
        int lbead = 1 << (il - 1);
        double max_disp = slt * sqrt(lbead);

        // Loop for atoms should be moved at the current level
        en_dpot = 0.0;
        for (int ilb = lbead; ilb < (1<<lmax); ilb += (1<<il)) {
            int ibead = (ilb + bs_ibead) % nbeads;
            int prev_bead = ((ilb - lbead) + bs_ibead) % nbeads;
            int next_bead = ((ilb + lbead) + bs_ibead) % nbeads;
            tmp_pos.at(ibead) = 0.5 * (tmp_pos.at(prev_bead) + tmp_pos.at(next_bead))
                + max_disp * normal_dist(*gen);
            en_dpot += calc_pot(tmp_pos.at(ibead)) - calc_pot(pos.at(ibead));
        }
        en_dpot *= lbead;
        en_dpot += 0.5 * en_dpot_prev;

        if (uniform_dist(*gen) < exp(-1.0*beta*(en_dpot-en_dpot_prev))) {
            en_dpot_prev = en_dpot;
        }
        else {
            acc = 0;
            break;
        }
    }

    if (acc) {
        // update configuration
        pos_com = 0.0;
        for (int ibead = 0; ibead < nbeads; ibead++) {
            pos.at(ibead) = tmp_pos.at(ibead);
            pos_com += pos.at(ibead);
        }
        pos_com /= nbeads;
        // update energy
        en.at("esprng") = calc_kin();
        en.at("ekin") = 1.0 * nbeads / (2.0 * beta) - en.at("esprng");
        en.at("epot") += en_dpot;
        en.at("toten") = en.at("ekin") + en.at("epot");
    }
}

void Path::move_com(double max_disp, mt19937_64* gen) {
    double en_dpot = 0.0;

    // random distribution
    uniform_real_distribution<double> uniform_dist(0, 1);

    double disp = max_disp * (1.0 - 2.0 * uniform_dist(*gen));
    for (int ibead = 0; ibead < nbeads; ibead++) {
        en_dpot += calc_pot(pos.at(ibead)+disp) - calc_pot(pos.at(ibead));
    }

    if (uniform_dist(*gen) < exp(-1.0*beta*en_dpot)) {
        // update configuration
        for (int ibead = 0; ibead < nbeads; ibead++) {
            pos.at(ibead) += disp;
        }
        pos_com += disp;
        // update energy
        en.at("epot") += en_dpot;
        en.at("toten") += en_dpot;
        accepted_com_trials++;
    }
    all_com_trials++;
}

double Path::calc_pot(double _pos) {
    // let hbar = omega = k_B = (hbar^2/2m) = 1
    // potential of harmonic oscillator: V = 1/2 m w^2 x^2
    return 0.25 * _pos * _pos / nbeads;
}

double Path::calc_kin() {
    // let hbar = omega = k_B = (hbar^2/2m) = 1
    double omega_p2 = nbeads / (beta * beta);
    double en_spring = 0.0;

    for (int ibead = 0; ibead < nbeads; ibead++) {
        int next_bead = (ibead + 1) % nbeads;
        en_spring += (pos.at(next_bead) - pos.at(ibead)) * (pos.at(next_bead) - pos.at(ibead));
    }

    return 0.25 * omega_p2 * en_spring;
}

void Path::calc_toten() {

    for (int ibead = 0; ibead < nbeads; ibead++) {
        en.at("epot") += calc_pot(pos.at(ibead));
    }
    en.at("esprng") = calc_kin();
    en.at("ekin") = 1.0 * nbeads / (2.0 * beta) - en.at("esprng");
    en.at("toten") = en.at("ekin") + en.at("epot");
}

void Path::update_com_max_disp(double* max_disp) {
    (*max_disp) *= (1.0 * accepted_com_trials / all_com_trials) / 0.5;
    accepted_com_trials = 0;
    all_com_trials = 0;
}