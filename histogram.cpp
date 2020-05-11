#include "histogram.h"

Histogram::Histogram(double _hmin, double _hmax, double _hstep) : hmin(_hmin), hmax(_hmax), hstep(_hstep) {
    hbins = static_cast<int>(ceil((hmax - hmin) / hstep));
    hcount = 0.0;

    hist["x"].resize(hbins, 0.0);
    bins["x"].resize(hbins);

    for (int i = 0; i < hbins; i++) {
        bins["x"][i] = (i + 0.5) * hstep + hmin;
    }
}

Histogram::~Histogram() {
}

void Histogram::get_hist(Path _path) {
    for (int ibead = 0; ibead < _path.nbeads; ibead++) {
        int idx = static_cast<int>(floor((_path.pos[ibead] - hmin) / hstep));
        if (idx >= 0 && idx < hbins) {
            hist["x"][idx] += 1.0;
            hcount += 1.0;
        }
    }
}

void Histogram::output() {
    FILE* fp;

    fp = fopen("histogram.dat", "w");
    for (int ibin = 0; ibin < hbins; ibin++) {
        fprintf(fp, "%8.3lf    %16.12lf\n", bins["x"][ibin], hist["x"][ibin]/hcount);
    }
    fclose(fp);
}