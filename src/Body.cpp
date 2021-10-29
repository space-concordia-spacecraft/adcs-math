#include "adcs/Body.h"

namespace adcs::body {
    void itrfToBody(vec3 be, vec3 ve, vec3 pe, vec3 se, mat3 dcm_be, vec3 &b, vec3 &v, vec3 &p, vec3 &s) {
        double angle = acos(glm::dot(se, pe)/(glm::length(se)*glm::length(pe))) * 180.0 / CONST_PI;

        b = be * dcm_be;
        v = ve * dcm_be;
        p = pe * dcm_be;

        if(angle > 109.79)
            s = vec3(0,0,0);
        else
            s = se * dcm_be;
    }

    void disturbance_torque(vec3 dipoleMoment, vec3 b, vec3 p, vec3 inertiaMatrix, vec3 s, vec3 v, double h, vec3 &torqueDisturbance){
        // magnetic field
        vec3 magneticField = glm::cross(dipoleMoment, (b * 1e-9));

        // gravity gradient
        p *= 1000;

        double tempdot = glm::dot(p,p);
        double thirdTerm = pow(tempdot, 2.5);

        vec3 temp = p * inertiaMatrix;
        vec3 fourthTerm = glm::cross(temp, p);

        vec3 gravityGradient = vec3(1,1,1);
        gravityGradient *= CONST_GRAVITATIONAL_CONSTANT;
        gravityGradient *= 3;
        gravityGradient *= thirdTerm;
        gravityGradient *= fourthTerm;

        // solar radiation pressure

        // atmospheric drag

    }
}