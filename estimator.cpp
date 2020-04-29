#include "estimator.h"


Estimator::Estimator(int _nsteps, int _nbins, bool _improved) : nbins(_nbins), improved(_improved) {
    // initialize local variables
    len = _nsteps / nbins;
    icount = 0; idx = 0;

    // initialize accumulators for binning
    accmltr["en"] = 0.0;

    // initialize accumulators for improved method to calculate energy and specific heat.
    // Ref: J.Chem.Phys.116,5951(2002); J.Chem.Phys.117,3020(2002).
    if (improved == true) {
        accmltr["f_e"] = 0.0; accmltr["f_e^2"] = 0.0; accmltr["ef_e"] = 0.0;
        accmltr["f_cv"] = 0.0; accmltr["f_cv^2"] = 0.0; accmltr["cvpf_cv"] = 0.0; accmltr["cvp"] = 0.0;
        accmltr["e^-bdv"] = 0.0;
    }
    else {
        accmltr["en2"] = 0.0; accmltr["cv_kin"] = 0.0;
    }

    // initialize binned estimators for statistics
    estimators["<E>"].resize(len, 0.0); estimators["<Cv>"].resize(len, 0.0);
}

Estimator::~Estimator() {
}

void Estimator::accumulate(Path path) {
    int nbeads = path.nbeads;
    double beta = path.beta;
    double toten = path.en.at("toten");
    double esprng = path.en.at("esprng");
    // for classical MC, the estimators are forced to primitive method
    if (nbeads == 1) improved = false;

    accmltr.at("en") += toten;
    // improved accumulator
    if (improved == true) {
        double Ge = 0.5 / beta;
        double Gcv = 0.75 / (beta * beta);
        double ekin = path.en.at("ekin");
        double epot = path.en.at("epot");
        double dV = calc_pot(path.pos_com) - epot;
        double ebdv = exp(-1.0 * beta * dV);
        double cvp = toten * toten + (nbeads / (2.0 * beta) - 2.0 * esprng) / beta;
        double f_e = ekin * ebdv;
        double f_cv = (ekin * ekin + (nbeads / (2.0 * beta) - 2.0 * esprng) / beta) * ebdv;

        accmltr.at("f_e") += f_e;
        accmltr.at("f_e^2") += f_e * f_e;
        accmltr.at("ef_e") += toten * f_e;
        accmltr.at("f_cv") += f_cv;
        accmltr.at("f_cv^2") += f_cv * f_cv;
        accmltr.at("cvp") += cvp;
        accmltr.at("cvpf_cv") += cvp * f_cv;
        accmltr.at("e^-bdv") += ebdv;
        accmltr.at("e^-bvc") += ebvc;
        icount++;

        if (icount % nbins == 0) {
            accmltr.at("en") /= nbins;
            accmltr.at("f_e") /= nbins; accmltr.at("f_e^2") /= nbins; accmltr.at("ef_e") /= nbins;
            accmltr.at("f_cv") /= nbins; accmltr.at("f_cv^2") /= nbins; accmltr.at("cvp") /= nbins;
            accmltr.at("cvpf_cv") /= nbins; accmltr.at("e^-bdv") /= nbins;

            double alpha_e = (accmltr.at("ef_e") - accmltr.at("en") * accmltr.at("f_e"))
                / (accmltr.at("f_e^2") - accmltr.at("f_e") * accmltr.at("f_e"));
            double alpha_cv = (accmltr.at("cvpf_cv") - accmltr.at("cvp") * accmltr.at("f_cv"))
                / (accmltr.at("f_cv^2") - accmltr.at("f_cv") * accmltr.at("f_cv"));

            estimators.at("<E>").at(idx) = accmltr.at("en") - alpha_e * accmltr.at("f_e") + alpha_e * Ge * accmltr.at("e^-bdv");
            estimators.at("<Cv>").at(idx) = (accmltr.at("cvp") - alpha_cv * accmltr.at("f_cv") + alpha_cv * Gcv * accmltr.at("e^-bdv")
                - estimators.at("<E>").at(idx) * estimators.at("<E>").at(idx)) * beta * beta;

            accmltr.at("en") = 0.0;
            accmltr.at("f_e") = 0.0; accmltr.at("f_e^2") = 0.0; accmltr.at("ef_e") = 0.0;
            accmltr.at("f_cv") = 0.0; accmltr.at("f_cv^2") = 0.0; accmltr.at("cvp") = 0.0;
            accmltr.at("cvpf_cv") = 0.0; accmltr.at("e^-bdv") = 0.0;
            idx++;
        }
    }
    // primitive accumulators
    else {
        accmltr.at("en2") += toten * toten;
        accmltr.at("cv_kin") += (nbeads / (2.0 * beta) - 2.0 * esprng) / beta;
        icount++;

        if (icount % nbins == 0) {
            accmltr.at("en") /= nbins; accmltr.at("en2") /= nbins; accmltr.at("cv_kin") /= nbins;

            estimators.at("<E>").at(idx) = accmltr.at("en");
            estimators.at("<Cv>").at(idx) = beta * beta
                * (accmltr.at("en2") - accmltr.at("en") * accmltr.at("en") + accmltr.at("cv_kin"));

            accmltr.at("en") = 0.0; accmltr.at("en2") = 0.0; accmltr.at("cv_kin") = 0.0;
            idx++;
        }
    }
}

void Estimator::output() {
    pair<double, double> en, cv;
    FILE* fp;

    fp = fopen("estimator.dat", "w");
    fprintf(fp, "# <E>                <Cv>\n");
    for (int i = 0; i < len; i++) {
        fprintf(fp, "%.12lf    %.12lf\n", estimators.at("<E>").at(i), estimators.at("<Cv>").at(i));
    }
    fclose(fp);

    // Simple output (error analysis by jackknife method)
    en = calc_avg_and_err_jackknife(estimators.at("<E>"));
    cv = calc_avg_and_err_jackknife(estimators.at("<Cv>"));
    
    cout << "E = " << en.first << " +/- " << en.second << endl;
    cout << "Cv = " << cv.first << " +/- " << cv.second << endl;
}

double Estimator::calc_pot(double pos) {
    return 0.25 * pos * pos;
}

pair<double, double> Estimator::calc_avg_and_err_jackknife(vector<double> data) {
    double avg = 0.0;

    for (int i = 0; i < len; i++) {
        avg += data.at(i);
    }
    avg /= len;

    double tmp_bar = 0.0;
    vector<double> tmp;
    tmp.resize(len);
    for (int i = 0; i < len; i++) {
        tmp.at(i) = (len * avg - data.at(i)) / (len - 1);
        tmp_bar += tmp.at(i);
    }
    tmp_bar /= len;

    double var = 0.0;
    for (int i = 0; i < len; i++) {
        var += (tmp.at(i) - tmp_bar) * (tmp.at(i) - tmp_bar);
    }
    double err = sqrt(((len - 1.0) / len) * var);

    return { avg, err };
}
