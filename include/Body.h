#pragma once

#include "Utils.h"

namespace adcs::body {
    void itrfToBody(vec3 be, vec3 ve, vec3 pe, vec3 se, mat3 dcm_be, vec3 &b, vec3 &v, vec3 &p, vec3 &s);
    void disturbance_torque(vec3 dipoleMoment, vec3 b, vec3 p, vec3 inertiaMatrix, vec3 s, vec3 v, double h, vec3 &torqueDisturbance);
};