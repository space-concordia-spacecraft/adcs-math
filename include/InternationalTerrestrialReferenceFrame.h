#pragma once

#include "Utils.h"

namespace adcs::itrf {
    void temeToItrf(vec3 steme, vec3 rteme, vec3 vteme, vec3 ateme,
                    double julianCenturies, double jdut1,
                    vec3 &secef, vec3 &recef, vec3 &vecef, vec3 &aecef, mat3 &DCM_ET,
                    double lod, double xp, double yp, int eqeterms);

    mat3 polarm(double xp, double yp, double ttt, eOpt opt);

    void xyz_ell3(vec3 position, double &latitude, double &longitude, double &altitude, double &h);

    void igrfs(double latitude, double longitude, double altitude, vec3 &position);

    void lg2ct(vec3 magneticField, double latitude, double longitude, vec3 &coordinateDifferences);

}