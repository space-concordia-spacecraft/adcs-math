#pragma once

#include "Utils.h"

namespace adcs::teme {
    typedef struct iau80data
    {
        int    iar80[107][6];
        double rar80[107][5];
    } iau80data;

    vec3 predictSunPosition(double jdtdb, double jdtdbF);

    mat3 precess(double julianCenturies);

    void nutation(double ttt, double ddpsi, double ddeps,
                  const iau80data &iau80arr, eOpt opt,
                  double& deltapsi, double& deltaeps, double& trueeps, double& meaneps, double& omega,
                  mat3 &nut);

    void fundarg(double ttt, eOpt opt,
                 double& l, double& l1, double& f, double& d, double& omega,
                 double& lonmer, double& lonven, double& lonear, double& lonmar,
                 double& lonjup, double& lonsat, double& lonurn, double& lonnep,
                 double& precrate);

}