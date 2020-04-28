#include "estimator.h"


Estimator::Estimator(int _nsteps, int _nbins) : nbins(_nbins) {
    // initialize local variables
    len = _nsteps / nbins;
    icount = 0; idx = 0;

    // initialize accumulators for binning
    tmp_en["en"] = 0.0; tmp_en["en2"] = 0.0; tmp_en["cv_kin"] = 0.0;

    // initialize binned estimators for statistics
    estimators["<E>"].resize(len, 0.0); estimators["<E^2>"].resize(len, 0.0);
    estimators["<Cv_kin>"].resize(len, 0.0); estimators["<Cv>"].resize(len, 0.0);
}

Estimator::~Estimator() {
}

void Estimator::accumulate(Path path) {
    double toten = path.en.at("toten");

    tmp_en.at("en") += toten;
    tmp_en.at("en2") += toten * toten;
    tmp_en.at("cv_kin") += (path.nbeads / (2.0 * path.beta) - 2.0 * path.en.at("esprng")) / path.beta;
    icount++;

    if (icount % nbins == 0) {
        tmp_en.at("en") /= nbins; tmp_en.at("en2") /= nbins; tmp_en.at("cv_kin") /= nbins;

        estimators.at("<E>").at(idx) = tmp_en.at("en");
        estimators.at("<E^2>").at(idx) = tmp_en.at("en2");
        estimators.at("<Cv_kin>").at(idx) = tmp_en.at("cv_kin");
        estimators.at("<Cv>").at(idx) = path.beta * path.beta
            * (tmp_en.at("en2") - tmp_en.at("en") * tmp_en.at("en") + tmp_en.at("cv_kin"));
 
        tmp_en["en"] = 0.0; tmp_en["en2"] = 0.0; tmp_en["cv_kin"] = 0.0;
        idx++;
    }
}

void Estimator::output() {
    pair<double, double> en, cv;
    FILE* fp;

    fp = fopen("estimator.dat", "w");
    fprintf(fp, "# <E>                <E^2>                <Cv_kin>                <Cv>\n");
    for (int i = 0; i < len; i++) {
        fprintf(fp, "%.12lf    %.12lf    %.12lf    %.12lf\n", estimators.at("<E>").at(i), estimators.at("<E^2>").at(i),
            estimators.at("<Cv_kin>").at(i), estimators.at("<Cv>").at(i));
    }
    fclose(fp);

    // Simple output (error analysis by jackknife method)
    en = calc_avg_and_err_jackknife(estimators.at("<E>"));
    cv = calc_avg_and_err_jackknife(estimators.at("<Cv>"));
    
    cout << "E = " << en.first << " +/- " << en.second << endl;
    cout << "Cv = " << cv.first << " +/- " << cv.second << endl;
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