#include "estimator.h"


Estimator::Estimator(int _nsteps, int _nbins, int _nbeads, string _method) : nbins(_nbins), method(_method) {
    // initialize local variables
    len = _nsteps / nbins;
    icount = 0; idx = 0;
    ccount = 0; cidx = 0;

    // initialize accumulators for binning
    accmltr["en"] = 0.0;

    // initialize accumulators for improved method to calculate energy and specific heat.
    // Ref: J.Chem.Phys.116,5951(2002); J.Chem.Phys.117,3020(2002).
    if (method == "PT") {
        accmltr["f_e"] = 0.0; accmltr["f_e^2"] = 0.0; accmltr["ef_e"] = 0.0;
        accmltr["f_cv"] = 0.0; accmltr["f_cv^2"] = 0.0; accmltr["cvpf_cv"] = 0.0; accmltr["cvp"] = 0.0;
        accmltr["e^-bdv"] = 0.0;
    }
    else if (method == "PCV") {
        accmltr["et"] = 0.0;
        accmltr["f_cv"] = 0.0; accmltr["f_cv^2"] = 0.0; accmltr["cvpf_cv"] = 0.0; accmltr["cvp"] = 0.0;
        accmltr["e^-bdv"] = 0.0;
    }
    else if (method == "CV") {
        accmltr["et"] = 0.0;
        accmltr["etecv"] = 0.0; 
    }
    else {
        accmltr["en2"] = 0.0; accmltr["cv_kin"] = 0.0;
    }

    // initialize binned estimators for statistics
    estimators["<E>"].resize(len, 0.0); estimators["<Cv>"].resize(len, 0.0);

    // initialize accumulator for correlation function
    accmltr_ctau.resize(_nbeads, 0.0);
    ctau.resize(len, vector<double>(_nbeads, 0.0));
}

Estimator::~Estimator() {
}

void Estimator::accumulate(Path path) {
    int nbeads = path.nbeads;
    double beta = path.beta;
    double toten = path.en["toten"];
    double esprng = path.en["esprng"];
    // for classical MC, the estimators are forced to primitive method
    if (nbeads == 1) method = "T";

    // projected thermodynamic accumulator
    if (method == "PT") {
        double Ge = 0.5 / beta;
        double Gcv = 0.75 / (beta * beta);
        double ekin = path.en["ekin"];
        double epot = path.en["epot"];
        double dV = calc_pot(path.pos_com) - epot;
        double ebdv = exp(-1.0 * beta * dV);
        double cvp = toten * toten + (nbeads / (2.0 * beta) - 2.0 * esprng) / beta;
        double f_e = ekin * ebdv;
        double f_cv = (ekin * ekin + (nbeads / (2.0 * beta) - 2.0 * esprng) / beta) * ebdv;
   
        accmltr["en"] += toten;
        accmltr["f_e"] += f_e;
        accmltr["f_e^2"] += f_e * f_e;
        accmltr["ef_e"] += toten * f_e;
        accmltr["f_cv"] += f_cv;
        accmltr["f_cv^2"] += f_cv * f_cv;
        accmltr["cvp"] += cvp;
        accmltr["cvpf_cv"] += cvp * f_cv;
        accmltr["e^-bdv"] += ebdv;
        icount++;

        if (icount % nbins == 0) {
            accmltr["en"] /= nbins;
            accmltr["f_e"] /= nbins; accmltr["f_e^2"] /= nbins; accmltr["ef_e"] /= nbins;
            accmltr["f_cv"] /= nbins; accmltr["f_cv^2"] /= nbins; accmltr["cvp"] /= nbins;
            accmltr["cvpf_cv"] /= nbins; accmltr["e^-bdv"] /= nbins;

            double alpha_e = (accmltr["ef_e"] - accmltr["en"] * accmltr["f_e"])
                / (accmltr["f_e^2"] - accmltr["f_e"] * accmltr["f_e"]);
            double alpha_cv = (accmltr["cvpf_cv"] - accmltr["cvp"] * accmltr["f_cv"])
                / (accmltr["f_cv^2"] - accmltr["f_cv"] * accmltr["f_cv"]);

            estimators["<E>"][idx] = accmltr["en"] - alpha_e * accmltr["f_e"] + alpha_e * Ge * accmltr["e^-bdv"];
            estimators["<Cv>"][idx] = (accmltr["cvp"] - alpha_cv * accmltr["f_cv"] + alpha_cv * Gcv * accmltr["e^-bdv"]
                - estimators["<E>"][idx] * estimators["<E>"][idx]) * beta * beta;

            accmltr["en"] = 0.0;
            accmltr["f_e"] = 0.0; accmltr["f_e^2"] = 0.0; accmltr["ef_e"] = 0.0;
            accmltr["f_cv"] = 0.0; accmltr["f_cv^2"] = 0.0; accmltr["cvp"] = 0.0;
            accmltr["cvpf_cv"] = 0.0; accmltr["e^-bdv"] = 0.0;
            idx++;
        }
    }
    // projected centroid viral estimators
    else if (method == "PCV") {
        double Gcv = 0.75 / (beta * beta);
        double ekin = path.en["ekin"];
        double epot = path.en["epot"];
        double ekcv = calc_kcv(path);
        double ecv = ekcv + epot;
        double dV = calc_pot(path.pos_com) - epot;
        double ebdv = exp(-1.0 * beta * dV);
        double cvp = toten * ecv + 1.0 / (2.0 * beta * beta);
        double f_cv = (ekin * (0.5 / beta) + 0.5 / (beta * beta)) * ebdv;

        accmltr["en"] += ecv;
        accmltr["et"] += toten;
        accmltr["f_cv"] += f_cv;
        accmltr["f_cv^2"] += f_cv * f_cv;
        accmltr["cvp"] += cvp;
        accmltr["cvpf_cv"] += cvp * f_cv;
        accmltr["e^-bdv"] += ebdv;
        icount++;

        if (icount % nbins == 0) {
            accmltr["en"] /= nbins; accmltr["et"] /= nbins;
            accmltr["f_cv"] /= nbins; accmltr["f_cv^2"] /= nbins; accmltr["cvp"] /= nbins;
            accmltr["cvpf_cv"] /= nbins; accmltr["e^-bdv"] /= nbins;

            double alpha_cv = (accmltr["cvpf_cv"] - accmltr["cvp"] * accmltr["f_cv"])
                / (accmltr["f_cv^2"] - accmltr["f_cv"] * accmltr["f_cv"]);

            estimators["<E>"][idx] = accmltr["en"];
            estimators["<Cv>"][idx] = (accmltr["cvp"] - alpha_cv * accmltr["f_cv"] + alpha_cv * Gcv * accmltr["e^-bdv"]
                - accmltr["en"] * accmltr["et"]) * beta * beta;

            accmltr["en"] = 0.0; accmltr["et"] = 0.0;
            accmltr["f_cv"] = 0.0; accmltr["f_cv^2"] = 0.0; accmltr["cvp"] = 0.0;
            accmltr["cvpf_cv"] = 0.0; accmltr["e^-bdv"] = 0.0;
            idx++;
        }
    }
    // centroid viral estimators
    else if (method == "CV") {
        double epot = path.en["epot"];
        double ekcv = calc_kcv(path);
        double ecv = ekcv + epot;

        accmltr["en"] += ecv;
        accmltr["et"] += toten;
        accmltr["etecv"] += toten * ecv;
        icount++;

        if (icount % nbins == 0) {
            accmltr["en"] /= nbins;
            accmltr["et"] /= nbins;
            accmltr["etecv"] /= nbins;

            estimators["<E>"][idx] = accmltr["en"];
            estimators["<Cv>"][idx] = (accmltr["etecv"] - accmltr["en"] * accmltr["et"] + 0.5 / (beta*beta)) * beta * beta;

            accmltr["en"] = 0.0;
            accmltr["et"] = 0.0;
            accmltr["etecv"] = 0.0;
            idx++;
        }
    }
    // primitive (thermodynamic) accumulators
    else {
        accmltr["en"] += toten;
        accmltr["en2"] += toten * toten;
        accmltr["cv_kin"] += (nbeads / (2.0 * beta) - 2.0 * esprng) / beta;
        icount++;

        if (icount % nbins == 0) {
            accmltr["en"] /= nbins; accmltr["en2"] /= nbins; accmltr["cv_kin"] /= nbins;

            estimators["<E>"][idx] = accmltr["en"];
            estimators["<Cv>"][idx] = beta * beta
                * (accmltr["en2"] - accmltr["en"] * accmltr["en"] + accmltr["cv_kin"]);

            accmltr["en"] = 0.0; accmltr["en2"] = 0.0; accmltr["cv_kin"] = 0.0;
            idx++;
        }
    }
}

void Estimator::accumulate_ctau(Path path) {
    int taup;

    for (int dtau = 0; dtau < path.nbeads; dtau++) {
        for (int tau = 0; tau < path.nbeads; tau++) {
            taup = (tau + dtau) % path.nbeads;
            accmltr_ctau[dtau] += path.pos[tau] * path.pos[taup];
        }
        accmltr_ctau[dtau] /= path.nbeads;
    }
    ccount++;

    if (ccount % nbins == 0) {
        for (int dtau = 0; dtau < path.nbeads; dtau++) {
            ctau[cidx][dtau] = accmltr_ctau[dtau] / nbins;
            accmltr_ctau[dtau] = 0.0;
        }
        cidx++;
    }
}

void Estimator::output() {
    pair<double, double> en, cv;
    FILE* fp;
    
    // output estimators
    fp = fopen("estimator.dat", "w");
    fprintf(fp, "# <E>                <Cv>\n");
    for (int i = 0; i < len; i++) {
         fprintf(fp, "%.12lf    %.12lf\n", estimators["<E>"][i], estimators["<Cv>"][i]);
    }
    fclose(fp);

    // Simple output (error analysis by jackknife method)
    en = calc_avg_and_err_jackknife(estimators["<E>"]);
    cv = calc_avg_and_err_jackknife(estimators["<Cv>"]);
    
    std::cout << "E = " << en.first << " +/- " << en.second << std::endl;
    std::cout << "Cv = " << cv.first << " +/- " << cv.second << std::endl;

    // output imaginary-time correlation function
    fp = fopen("D:\\ctau.dat", "w");
    for (auto c : ctau) {
        for (auto ct : c) {
            fprintf(fp, "%20.12lf", ct);
        }
        fprintf(fp, "%20.12lf\n", c.front());
    }
    fclose(fp);
}

double Estimator::calc_pot(double pos) {
    return 0.25 * pos * pos;
}

double Estimator::calc_kcv(Path path) {
    double ycp = 0.0;
    for (int i = 0; i < path.nbeads; i++) {
        ycp += (path.pos[i] - path.pos_com) * 0.5 * path.pos[i];
    }
    ycp /= 2 * path.nbeads;

    return ycp + 0.5 / path.beta;
}

pair<double, double> Estimator::calc_avg_and_err_jackknife(vector<double> data) {
    double avg = 0.0;

    for (int i = 0; i < len; i++) {
        avg += data[i];
    }
    avg /= len;

    double tmp_bar = 0.0;
    vector<double> tmp;
    tmp.resize(len);
    for (int i = 0; i < len; i++) {
        tmp[i] = (len * avg - data[i]) / (len - 1);
        tmp_bar += tmp[i];
    }
    tmp_bar /= len;

    double var = 0.0;
    for (int i = 0; i < len; i++) {
        var += (tmp[i] - tmp_bar) * (tmp[i] - tmp_bar);
    }
    double err = sqrt(((len - 1.0) / len) * var);

    return { avg, err };
}
