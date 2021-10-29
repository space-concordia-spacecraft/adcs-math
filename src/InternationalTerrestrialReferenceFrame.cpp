#include "adcs/InternationalTerrestrialReferenceFrame.h"

namespace adcs::itrf {
    void temeToItrf(vec3 steme, vec3 rteme, vec3 vteme, vec3 ateme,
                    double julianCenturies, double jdut1,
                    vec3 &secef, vec3 &recef, vec3 &vecef, vec3 &aecef, mat3 &DCM_ET,
                    double lod, double xp, double yp, int eqeterms) {
        double omega, gmstg, thetasa, gmst;

        vec3 omegaearth, rpef, vpef, apef, spef, omgxr, omgxomgxr, omgxv, tempvec1, tempvec, temp;

        mat3 st, stdot, pm, pmp, stp;

        gmst = gstime(jdut1);

        // find omeage from nutation theory
        omega = 125.04452222 + (-6962890.5390 * julianCenturies + 7.455 * julianCenturies * julianCenturies + 0.008 * julianCenturies * julianCenturies * julianCenturies) / 3600.0;
        omega = fmod(omega, 360.0) * CONST_DEG_TO_RAD;

        // teme does not include the geometric terms here
        // after 1997, kinematic terms apply
        if ((jdut1 > 2450449.5) && (eqeterms > 0))
        {
            gmstg = gmst
                    + 0.00264*CONST_PI / (3600.0 * 180.0)*sin(omega)
                    + 0.000063*CONST_PI / (3600.0 * 180.0)*sin(2.0 *omega);
        }
        else
            gmstg = gmst;
        gmstg = fmod(gmstg, 2.0*CONST_PI);

        thetasa = 7.29211514670698e-05 * (1.0 - lod / 86400.0);
        omegaearth[0] = 0.0;
        omegaearth[1] = 0.0;
        omegaearth[2] = thetasa;

        st[0][0] = cos(gmstg);
        st[0][1] = -sin(gmstg);
        st[0][2] = 0.0;
        st[1][0] = sin(gmstg);
        st[1][1] = cos(gmstg);
        st[1][2] = 0.0;
        st[2][0] = 0.0;
        st[2][1] = 0.0;
        st[2][2] = 1.0;

        pm = polarm(xp, yp, julianCenturies, e80);

        rpef  = transpose(st) * rteme;
        recef = transpose(pm) * rpef;

        spef  = transpose(st) * steme;
        secef = transpose(pm) * spef;

        thetasa = 7.29211514670698e-05 * (1.0  - lod/86400.0 );
        omegaearth = vec3(0, 0, thetasa);

        vpef  = (transpose(st) * vteme) - cross(omegaearth, rpef);
        vecef = transpose(pm) * vpef;

        temp  = cross(omegaearth,rpef);

        aecef = transpose(pm)*(transpose(st)*ateme - cross(omegaearth,temp) - (2.0f * cross(omegaearth,vpef)));
        DCM_ET = transpose(pm) * st;
    }

    mat3 polarm(double xp, double yp, double julianCenturies, eOpt opt) {
        double cosxp, cosyp, sinxp, sinyp, sp, cossp, sinsp;

        mat3 pm = mat3();

        cosxp = cos(xp);
        sinxp = sin(xp);
        cosyp = cos(yp);
        sinyp = sin(yp);

        if ((opt == e80) | (opt == e96))
        {
            pm[0][0] = cosxp;
            pm[0][1] = 0.0;
            pm[0][2] = -sinxp;
            pm[1][0] = sinxp * sinyp;
            pm[1][1] = cosyp;
            pm[1][2] = cosxp * sinyp;
            pm[2][0] = sinxp * cosyp;
            pm[2][1] = -sinyp;
            pm[2][2] = cosxp * cosyp;
        }
        else
        {
            // approximate sp value in rad
            sp = -47.0e-6 * julianCenturies * CONST_PI / (180.0 * 3600.0);
            cossp = cos(sp);
            sinsp = sin(sp);

            // form the matrix
            pm[0][0] = cosxp * cossp;
            pm[0][1] = -cosyp * sinsp + sinyp * sinxp * cossp;
            pm[0][2] = -sinyp * sinsp - cosyp * sinxp * cossp;
            pm[1][0] = cosxp * sinsp;
            pm[1][1] = cosyp * cossp + sinyp * sinxp * sinsp;
            pm[1][2] = sinyp * cossp - cosyp * sinxp * sinsp;
            pm[2][0] = sinxp;
            pm[2][1] = -sinyp * cosxp;
            pm[2][2] = cosyp * cosxp;
        }

        return pm;
    }  //  polarm

    void xyz_ell3(vec3 position, double &latitude, double &longitude, double &altitude, double &h) {
        double a = 6378137.0;
        double finv = 298.257222101;

        double f = 1 / finv;
        double b = a * (1 - f);
        double e2 = 1 - pow((1 - f), 2);

        longitude = atan2(position.y,position.x); //longitude
        double e = e2*pow((a/b),2);
        double p = sqrt((pow(position.x, 2))+pow(position.y, 2));
        double r = sqrt(pow(p, 2)+pow(position.z, 2));
        double u = atan(b*position.z*(1+e*b/r)/(a*p));
        latitude = atan((position.z+(e*b*pow(sin(u),3) ))/(p-(e2*a*pow(cos(u),3)))); // latitude
        double v = a/sqrt(1-(e2*pow(sin(latitude),2) ));
        h = p*cos(latitude)+position.z*sin(latitude)-a*a/v;
        altitude = r; // altitude
    }

    void igrfs(double latitude, double longitude, double altitude, vec3 &magneticField) {
        double costheta = cos(((CONST_PI / 2.0)- latitude));
        double sintheta = sin(((CONST_PI / 2.0)- latitude));

        double r = altitude;
        double phi = longitude;

        double gh [195] = {-29442, -1501, 4797.1, -2445.1, 3012.9, -2845.6, 1676.7, -641.9, 1350.7, -2352.3, -115.3, 1225.6, 244.9, 582, -538.4, 907.6, 813.7, 283.3, 120.4, -188.7, -334.9, 180.9, 70.4, -329.5, -232.6, 360.1, 47.3, 192.4, 197, -140.9, -119.3, -157.5, 16, 4.1, 100.2, 70, 67.7, -20.8, 72.7, 33.2, -129.9, 58.9, -28.9, -66.7, 13.2, 7.3, -70.9, 62.6, 81.6, -76.1, -54.1, -6.8, -19.5, 51.8, 5.7, 15, 24.4, 9.4, 3.4, -2.8, -27.4, 6.8, -2.2, 24.2, 8.8, 10.1, -16.9, -18.3, -3.2, 13.3, -20.6, -14.6, 13.4, 16.2, 11.7, 5.7, -15.9, -9.1, -2, 2.1, 5.4, 8.8, -21.6, 3.1, 10.8, -3.3, 11.8, 0.7, -6.8, -13.3, -6.9, -0.1, 7.8, 8.7, 1, -9.1, -4, -10.5, 8.4, -1.9, -6.3, 3.2, 0.1, -0.4, 0.5, 4.6, -0.5, 4.4, 1.8, -7.9, -0.7, -0.6, 2.1, -4.2, 2.4, -2.8, -1.8, -1.2, -3.6, -8.7, 3.1, -1.5, -0.1, -2.3, 2, 2, -0.7, -0.8, -1.1, 0.6, 0.8, -0.7, -0.2, 0.2, -2.2, 1.7, -1.4, -0.2, -2.5, 0.4, -2, 3.5, -2.4, -1.9, -0.2, -1.1, 0.4, 0.4, 1.2, 1.9, -0.8, -2.2, 0.9, 0.3, 0.1, 0.7, 0.5, -0.1, -0.3, 0.3, -0.4, 0.2, 0.2, -0.9, -0.9, -0.1, 0, 0.7, 0, -0.9, -0.9, 0.4, 0.4, 0.5, 1.6, -0.5, -0.5, 1, -1.2, -0.2, -0.1, 0.8, 0.4, -0.1, -0.1, 0.3, 0.4, 0.1, 0.5, 0.5, -0.3, -0.4, -0.4, -0.3, -0.8};

        int nmax = 13;

        float cosphi[nmax + 1];
        float sinphi[nmax + 1];

        for (int i = 1; i <= nmax; i++) {
            cosphi[i] = cos(i * phi);
            sinphi[i] = sin(i * phi);
        }


        int Pmax = (nmax + 1) * (nmax + 2) / 2;

        double Br = 0;
        double Bt = 0;
        double Bp = 0;

        double P[Pmax + 1];

        P[1] = 1;
        P[3] = sintheta;

        double dP[Pmax + 1];

        dP[1] = 0;
        dP[3] = costheta;


        int m = 1;
        int n = 0;
        int coefindex = 1;

        double a_r = pow((CONST_EARTH_RADIUS / r), 2);

        for (int i = 2; i <= Pmax; i++) {
            if (n < m) {
                m = 0;
                n++;
                a_r = a_r * (CONST_EARTH_RADIUS / r);
            }
            if (m < n && i != 3) {
                int last1n = i - n;
                int last2n = i - 2 * n + 1;
                // do magn field
                P[i] = (2 * n - 1) / sqrt(pow(n, 2) - pow(m, 2)) * costheta * P[last1n] -  sqrt((pow((n - 1), 2) - pow(m, 2)) / (pow(n, 2) - pow(m, 2))) * P[last2n];
                dP[i] = (2 * n - 1) / sqrt(pow(n, 2) - pow(m, 2)) * (costheta * dP[last1n] -
                                                                     sintheta * P[last1n]) - sqrt((pow((n - 1), 2) - pow(m, 2)) / (pow(n, 2) - pow(m, 2))) *
                                                                                             dP[last2n];
            }
            else if(i != 3) {
                int lastn = i - n - 1;
                P[i] = sqrt(1 - 1 / (2 * m)) * sintheta * P[lastn];
                dP[i] = sqrt(1 - 1 / (2 * m)) * (sintheta * dP[lastn] + costheta * P[lastn]);
            }
            if (m == 0) {
                double coef = a_r * gh[coefindex];
                Br = Br + (n + 1) * coef * P[i];
                Bt = Bt - coef * dP[i];
                coefindex++;
            }
            else {
                double coef = a_r * (gh[coefindex] * cosphi[m] + gh[coefindex + 1] * sinphi[m]);
                Br = Br + (n + 1) * coef * P[i];
                Bt = Bt - coef * dP[i];

                if (sintheta == 0) {
                    Bp = Bp - costheta * a_r * (-gh[coefindex] * sinphi[m] + gh[coefindex + 1] * cosphi[m]) * dP[i];
                }
                else {
                    Bp = Bp - 1 / sintheta * a_r * m * (-gh[coefindex] * sinphi[m] +
                                                        gh[coefindex + 1] * cosphi[m]) * P[i];
                }
                coefindex = coefindex + 2;
            }

            m++;

        }

        //Bn
        magneticField[0] = -Bt;
        //Be
        magneticField[1] = Bp;
        //Bd
        magneticField[2] = Br;

    }

    void lg2ct(vec3 position, double latitude, double longitude, vec3 &coordinateDifferences) {

    }
}