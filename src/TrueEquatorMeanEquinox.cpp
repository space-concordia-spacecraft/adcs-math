#include "adcs/TrueEquatorMeanEquinox.h"

namespace adcs::teme {
    // simple function to store help levels that can be set on/off during execution
    void sethelp(char& iauhelp, char iauopt)
    {
        static char iaustore;

        if (iauopt != ' ')
        {
            iauhelp = iauopt;
            iaustore = iauopt;
        }
        else
            iauhelp = iaustore;
    }

    vec3 predictSunPosition(double jdtdb, double jdtdbF) {
        // Julian Centuries of Terrestrial time
        double julianCenturies = (jdtdb + jdtdbF - 2451545.0) / 36525.0;

        // ----- START SUNALMANAC ----- //

        vec3 rsun(0, 0, 0), reci(0,0,0), rteme(0, 0, 0);
        double meanlong, meananomaly, eclplong, obliquity, magr, rtasc, decl;

        // -------------------------  implementation   -----------------
        // -------------------  initialize values   --------------------
        meanlong = 280.460 + 36000.77 * julianCenturies;
        meanlong = fmod(meanlong, 360.0);  //deg

        meananomaly = 357.5277233 + 35999.05034 * julianCenturies;
        meananomaly = fmod(meananomaly * CONST_DEG_TO_RAD, CONST_TWO_PI);  //rad
        if (meananomaly < 0.0)
        {
            meananomaly = CONST_TWO_PI + meananomaly;
        }
        eclplong = meanlong + 1.914666471 * sin(meananomaly)
                   + 0.019994643 * sin(2.0 * meananomaly); //deg
        obliquity = 23.439291 - 0.0130042 * julianCenturies;  //deg
        meanlong = meanlong * CONST_DEG_TO_RAD;
        if (meanlong < 0.0)
        {
            meanlong = CONST_TWO_PI + meanlong;
        }
        eclplong = eclplong * CONST_DEG_TO_RAD;
        obliquity = obliquity * CONST_DEG_TO_RAD;

        // --------- find magnitude of sun vector, ) components ------
        magr = 1.000140612 - 0.016708617 * cos(meananomaly)
               - 0.000139589 * cos(2.0 *meananomaly);    // in au's

        rsun[0] = magr * cos(eclplong);
        rsun[1] = magr * cos(obliquity) * sin(eclplong);
        rsun[2] = magr * sin(obliquity) * sin(eclplong);

        rtasc = atan(cos(obliquity) * tan(eclplong));

        // --- check that rtasc is in the same quadrant as eclplong ----
        if (eclplong < 0.0)
        {
            eclplong = eclplong + CONST_TWO_PI;    // make sure it's in 0 to 2pi range
        }
        if (fabs(eclplong - rtasc) > CONST_PI * 0.5)
        {
            rtasc = rtasc + 0.5 * CONST_PI * round((eclplong - rtasc) / (0.5 * CONST_PI));
        }
        decl = asin(sin(obliquity) * sin(eclplong));

        // ----- DONE SUNALMANAC ----- //

        // ----- START MOD2ECI ----- //

        mat3 prec = precess(julianCenturies);
        reci = rsun * prec;

        // ----- DONE MOD2ECI ----- //

        // ----- START ECI2TEME ----- //

        iau80data iau;

        double deltapsi = 0, deltaeps = 0, trueeps = 0, meaneps = 0, omega = 0;
        mat3 nut;

        nutation(julianCenturies, 0, 0,
                 iau, e80,
                 deltapsi, deltaeps, trueeps, meaneps, omega,
                 nut);

        double eqeg = deltapsi * cos(meaneps);
        eqeg = fmod(eqeg, 2.0*CONST_PI);

        mat3 eqe;

        eqe[0][0] = cos(eqeg);
        eqe[0][1] = sin(eqeg);
        eqe[0][2] = 0.0;
        eqe[1][0] = -sin(eqeg);
        eqe[1][1] = cos(eqeg);
        eqe[1][2] = 0.0;
        eqe[2][0] = 0.0;
        eqe[2][1] = 0.0;
        eqe[2][2] = 1.0;

        mat3 tm = eqe * transpose(nut) * transpose(prec);

        rteme = tm * reci;

        // ----- DONE ECI2TEME ----- //

        return rteme;

    }

    mat3 precess(double julianCenturies) {
        double convertToRad = CONST_PI / (180.0 * 3600.0);
        double julianCenturiesSquared = pow(julianCenturies, 2);
        double julianCenturiesCubed = pow(julianCenturies, 3);

        double psia  =             5038.7784*julianCenturies - 1.07259*julianCenturiesSquared - 0.001147*julianCenturiesCubed;
        double wa    = 84381.448                             + 0.05127*julianCenturiesSquared - 0.007726*julianCenturiesCubed;
        double ea    = 84381.448 -   46.8150*julianCenturies - 0.00059*julianCenturiesSquared + 0.001813*julianCenturiesCubed;
        double xa    =               10.5526*julianCenturies - 2.38064*julianCenturiesSquared - 0.001125*julianCenturiesCubed;
        double zeta  =             2306.2181*julianCenturies + 0.30188*julianCenturiesSquared + 0.017998*julianCenturiesCubed;
        double theta =             2004.3109*julianCenturies - 0.42665*julianCenturiesSquared - 0.041833*julianCenturiesCubed;
        double z     =             2306.2181*julianCenturies + 1.09468*julianCenturiesSquared + 0.018203*julianCenturiesCubed;

        psia  = psia  * convertToRad;
        wa    = wa    * convertToRad;
        ea    = ea    * convertToRad;
        xa    = xa    * convertToRad;
        zeta  = zeta  * convertToRad;
        theta = theta * convertToRad;
        z     = z     * convertToRad;

        double coszeta  = cos(zeta);
        double sinzeta  = sin(zeta);
        double costheta = cos(theta);
        double sintheta = sin(theta);
        double cosz     = cos(z);
        double sinz     = sin(z);

        double row1col1 = (coszeta * costheta * cosz) - (sinzeta * sinz);
        double row1col2 = (coszeta * costheta * sinz) + (sinzeta * cosz);
        double row1col3 =  coszeta * sintheta;
        double row2col1 = (-sinzeta * costheta * cosz) - (coszeta * sinz);
        double row2col2 = (-sinzeta * costheta * sinz) + (coszeta * cosz);
        double row2col3 = -sinzeta * sintheta;
        double row3col1 = -sintheta * cosz;
        double row3col2 = -sintheta * sinz;
        double row3col3 =  costheta;

        return mat3(row1col1, row2col1, row3col1, row1col2, row2col2, row3col2,row1col3, row2col3, row3col3);
    }

    /* -----------------------------------------------------------------------------
    *
    *                           function nutation
    *
    *  this function calulates the transformation matrix that accounts for the
    *    effects of nutation.
    *
    *  author        : david vallado                  719-573-2600   27 jun 2002
    *
    *  revisions
    *    vallado     - consolidate with iau 2000                     14 feb 2005
    *    vallado     - conversion to c++                             21 feb 2005
    *    vallado     - conversion to c#                              16 Nov 2011
    *
    *  inputs          description                                 range / units
    *    ttt         - julian centuries of tt
    *    ddpsi       - delta psi correction to gcrf                      rad
    *    ddeps       - delta eps correction to gcrf                      rad
    *    iau80arr    - record containing the iau80 constants rad
    *    opt         - method option                               e00cio, e00a, e96, e80
    *
    *  outputs       :
    *    deltapsi    - nutation in longitude angle                       rad
    *    trueeps     - true obliquity of the ecliptic                    rad
    *    meaneps     - mean obliquity of the ecliptic                    rad
    *    raan        -                                                   rad
    *    nut         - transform matrix for tod - mod
    *
    *  locals        :
    *    iar80       - integers for fk5 1980
    *    rar80       - reals for fk5 1980                                rad
    *    l           -                                                   rad
    *    ll          -                                                   rad
    *    f           -                                                   rad
    *    d           -                                                   rad
    *    deltaeps    - change in obliquity                               rad
    *
    *  coupling      :
    *    fundarg     - find fundamental arguments
    *    fmod      - modulus division
    *
    *  references    :
    *    vallado       2013, 213, 224
    * --------------------------------------------------------------------------- */
    void nutation
            (
                    double ttt, double ddpsi, double ddeps,
                    const iau80data &iau80arr, eOpt opt,
                    double& deltapsi, double& deltaeps, double& trueeps, double& meaneps, double& omega,
                    mat3 &nut
            )
    {
        // locals
        double  l, l1, f, d,
                lonmer, lonven, lonear, lonmar, lonjup, lonsat, lonurn, lonnep, precrate,
                cospsi, sinpsi, coseps, sineps, costrueeps, sintrueeps;
        int  i;
        double tempval;

        char iauhelp;
        sethelp(iauhelp, ' ');

        // ---- determine coefficients for iau 1980 nutation theory ----
        meaneps = ((0.001813  * ttt - 0.00059) * ttt - 46.8150) * ttt + 84381.448;
        meaneps = fmod(meaneps / 3600.0, 360.0);
        meaneps = meaneps * CONST_DEG_TO_RAD;

        fundarg(ttt, opt, l, l1, f, d, omega,
                lonmer, lonven, lonear, lonmar, lonjup, lonsat, lonurn, lonnep, precrate);

        deltapsi = 0.0;
        deltaeps = 0.0;
        for (i = 105; i >= 0; i--)
        {
            tempval = iau80arr.iar80[i][0] * l + iau80arr.iar80[i][1] * l1 + iau80arr.iar80[i][2] * f +
                      iau80arr.iar80[i][3] * d + iau80arr.iar80[i][4] * omega;
            deltapsi = deltapsi + (iau80arr.rar80[i][0] + iau80arr.rar80[i][1] * ttt)  *  sin(tempval);
            deltaeps = deltaeps + (iau80arr.rar80[i][2] + iau80arr.rar80[i][3] * ttt) * cos(tempval);
        }

        // --------------- find nutation parameters --------------------
        deltapsi = fmod(deltapsi + ddpsi, 2.0 * CONST_PI);
        deltaeps = fmod(deltaeps + ddeps, 2.0 * CONST_PI);

        trueeps = meaneps + deltaeps;

        cospsi = cos(deltapsi);
        sinpsi = sin(deltapsi);
        coseps = cos(meaneps);
        sineps = sin(meaneps);
        costrueeps = cos(trueeps);
        sintrueeps = sin(trueeps);

        nut[0][0] = cospsi;
        nut[0][1] = costrueeps * sinpsi;
        nut[0][2] = sintrueeps * sinpsi;
        nut[1][0] = -coseps * sinpsi;
        nut[1][1] = costrueeps * coseps * cospsi + sintrueeps * sineps;
        nut[1][2] = sintrueeps * coseps * cospsi - sineps * costrueeps;
        nut[2][0] = -sineps * sinpsi;
        nut[2][1] = costrueeps * sineps * cospsi - sintrueeps * coseps;
        nut[2][2] = sintrueeps * sineps * cospsi + costrueeps * coseps;

        // alternate approach
        //MathTimeLib::rot1mat(trueeps, n1);
        //MathTimeLib::rot3mat(deltapsi, n2);
        //MathTimeLib::rot1mat(-meaneps, n3);
        //MathTimeLib::matmult(n2, n1, tr1, 3, 3, 3);
        //MathTimeLib::matmult(n3, tr1, nut, 3, 3, 3);

        if (iauhelp == 'y')
            printf("meaneps %11.7f dp  %11.7f de  %11.7f te  %11.7f  \n", meaneps * 180 / CONST_PI, deltapsi * 180 / CONST_PI,
                   deltaeps * 180 / CONST_PI, trueeps * 180 / CONST_PI);
    }  //  nutation

    /* -----------------------------------------------------------------------------
         *
         *                           function fundarg
         *
         *  this function calulates the delauany variables and planetary values for
         *  several theories.
         *
         *  author        : david vallado                  719-573-2600   16 jul 2004
         *
         *  revisions
         *    vallado     - conversion to c++                             23 nov 2005
         *    vallado     - conversion to c#                              16 Nov 2011
         *
         *  inputs          description                                  range / units
         *    ttt         - julian centuries of tt
         *    opt         - method option                                e00cio, e00a, e96, e80
         *
         *  outputs       :
         *    l           - mean anomaly of the moon                          rad
         *    l1          - mean anomaly of the Sun                           rad
         *    f           - mean longitude of the Moon minus that of asc node rad
         *    d           - mean elongation of the Moon from the Sun          rad
         *    omega       - mean longitude of the ascending node of the Moon  rad
         *    planetary longitudes                                          rad
         *
         *  locals        :
         *
         *  coupling      :
         *    none        -
         *
         *  references    :
         *    vallado       2013, 210-211, 225
         * --------------------------------------------------------------------------- */

    void fundarg
            (
                    double ttt, eOpt opt,
                    double& l, double& l1, double& f, double& d, double& omega,
                    double& lonmer, double& lonven, double& lonear, double& lonmar,
                    double& lonjup, double& lonsat, double& lonurn, double& lonnep,
                    double& precrate
            )
    {
        double oo3600;
        char iauhelp;

        sethelp(iauhelp, ' ');
        oo3600 = 1.0 / 3600.0;
        l = l1 = f = d = omega = lonmer = lonven = lonear = lonmar = lonjup = lonsat = lonurn = lonnep = precrate = 0.0;

        // ---- determine coefficients for various iers nutation theories ----
        // ----  iau-2010 cio theory and iau-2000a theory
        if (opt == e00cio || opt == e00a)
        {
            // ------ form the delaunay fundamental arguments in ", converted to rad
            l = ((((-0.00024470 * ttt + 0.051635) * ttt + 31.8792) * ttt + 1717915923.2178) * ttt + 485868.249036) * oo3600;
            l1 = ((((-0.00001149 * ttt + 0.000136) * ttt - 0.5532) * ttt + 129596581.0481) * ttt + 1287104.793048) * oo3600;
            f = ((((+0.00000417 * ttt - 0.001037) * ttt - 12.7512) * ttt + 1739527262.8478) * ttt + 335779.526232) * oo3600;
            d = ((((-0.00003169 * ttt + 0.006593) * ttt - 6.3706) * ttt + 1602961601.2090) * ttt + 1072260.703692) * oo3600;
            omega = ((((-0.00005939 * ttt + 0.007702) * ttt + 7.4722) * ttt - 6962890.5431) * ttt + 450160.398036) * oo3600;

            // ------ form the planetary arguments in ", converted to rad
            lonmer = (908103.259872 + 538101628.688982  * ttt) * oo3600;
            lonven = (655127.283060 + 210664136.433548  * ttt) * oo3600;
            lonear = (361679.244588 + 129597742.283429  * ttt) * oo3600;
            lonmar = (1279558.798488 + 68905077.493988  * ttt) * oo3600;
            lonjup = (123665.467464 + 10925660.377991  * ttt) * oo3600;
            lonsat = (180278.799480 + 4399609.855732  * ttt) * oo3600;
            lonurn = (1130598.018396 + 1542481.193933  * ttt) * oo3600;
            lonnep = (1095655.195728 + 786550.320744  * ttt) * oo3600;
            precrate = ((1.112022 * ttt + 5028.8200) * ttt) * oo3600;
            // these are close (all in rad) - usually 1e-10, but some are as high as 1e-06
            //lonmer = (4.402608842 + 2608.7903141574 * ttt) % twopi;
            //lonven = (3.176146697 + 1021.3285546211 * ttt) % twopi;
            //lonear = (1.753470314 + 628.3075849991 * ttt) % twopi;
            //lonmar = (6.203480913 + 334.0612426700 * ttt) % twopi;
            //lonjup = (0.599546497 + 52.9690962641 * ttt) % twopi;
            //lonsat = (0.874016757 + 21.3299104960 * ttt) % twopi;
            //lonurn = (5.481293872 + 7.4781598567 * ttt) % twopi;
            //lonnep = (5.311886287 + 3.8133035638 * ttt) % twopi;
            //precrate = (0.024381750 + 0.00000538691 * ttt ) *ttt;
        }

        // ---- iau-2000b theory
        if (opt == e00b)
        {
            // ------ form the delaunay fundamental arguments in deg
            l = (1717915923.2178  * ttt + 485868.249036) * oo3600;
            l1 = (129596581.0481  * ttt + 1287104.79305) * oo3600;
            f = (1739527262.8478  * ttt + 335779.526232) * oo3600;
            d = (1602961601.2090  * ttt + 1072260.70369) * oo3600;
            omega = (-6962890.5431  * ttt + 450160.398036) * oo3600;

            // ------ form the planetary arguments in deg
            lonmer = 0.0;
            lonven = 0.0;
            lonear = 0.0;
            lonmar = 0.0;
            lonjup = 0.0;
            lonsat = 0.0;
            lonurn = 0.0;
            lonnep = 0.0;
            precrate = 0.0;
            // instead uses a constant rate
            // dplan = -0.135 * oo3600 * deg2rad;
            // deplan = 0.388 * oo3600 * deg2rad;
        }

        // ---- iau-1996 theory
        if (opt == e96)
        {
            // ------ form the delaunay fundamental arguments in deg
            l = ((((-0.00024470 * ttt + 0.051635) * ttt + 31.8792) * ttt + 1717915923.2178) * ttt) * oo3600 + 134.96340251;
            l1 = ((((-0.00001149 * ttt - 0.000136) * ttt - 0.5532) * ttt + 129596581.0481) * ttt) * oo3600 + 357.52910918;
            f = ((((+0.00000417 * ttt + 0.001037) * ttt - 12.7512) * ttt + 1739527262.8478) * ttt) * oo3600 + 93.27209062;
            d = ((((-0.00003169 * ttt + 0.006593) * ttt - 6.3706) * ttt + 1602961601.2090) * ttt) * oo3600 + 297.85019547;
            omega = ((((-0.00005939 * ttt + 0.007702) * ttt + 7.4722) * ttt - 6962890.2665) * ttt) * oo3600 + 125.04455501;
            // ------ form the planetary arguments in deg
            lonmer = 0.0;
            lonven = 181.979800853 + 58517.8156748   * ttt;
            lonear = 100.466448494 + 35999.3728521   * ttt;
            lonmar = 355.433274605 + 19140.299314    * ttt;
            lonjup = 34.351483900 + 3034.90567464  * ttt;
            lonsat = 50.0774713998 + 1222.11379404  * ttt;
            lonurn = 0.0;
            lonnep = 0.0;
            precrate = (0.0003086 * ttt + 1.39697137214) * ttt;
        }

        // ---- iau-1980 theory
        if (opt == e80)
        {
            // ------ form the delaunay fundamental arguments in deg
            l = ((((0.064) * ttt + 31.310) * ttt + 1717915922.6330) * ttt) * oo3600 + 134.96298139;
            l1 = ((((-0.012) * ttt - 0.577) * ttt + 129596581.2240) * ttt) * oo3600 + 357.52772333;
            f = ((((0.011) * ttt - 13.257) * ttt + 1739527263.1370) * ttt) * oo3600 + 93.27191028;
            d = ((((0.019) * ttt - 6.891) * ttt + 1602961601.3280) * ttt) * oo3600 + 297.85036306;
            omega = ((((0.008) * ttt + 7.455) * ttt - 6962890.5390) * ttt) * oo3600 + 125.04452222;
            // ------ form the planetary arguments in deg
            // iers tn13 shows no planetary
            // seidelmann shows these equations
            // circ 163 shows no planetary
            // ???????
            lonmer = 252.3 + 149472.0  * ttt;
            lonven = 179.9 + 58517.8  * ttt;
            lonear = 98.4 + 35999.4  * ttt;
            lonmar = 353.3 + 19140.3  * ttt;
            lonjup = 32.3 + 3034.9  * ttt;
            lonsat = 48.0 + 1222.1  * ttt;
            lonurn = 0.0;
            lonnep = 0.0;
            precrate = 0.0;
        }

        // ---- convert units from deg to rad
        l = fmod(l, 360.0)      *  CONST_DEG_TO_RAD;
        l1 = fmod(l1, 360.0)     *  CONST_DEG_TO_RAD;
        f = fmod(f, 360.0)      *  CONST_DEG_TO_RAD;
        d = fmod(d, 360.0)      *  CONST_DEG_TO_RAD;
        omega = fmod(omega, 360.0)  *  CONST_DEG_TO_RAD;

        lonmer = fmod(lonmer, 360.0) * CONST_DEG_TO_RAD;
        lonven = fmod(lonven, 360.0) * CONST_DEG_TO_RAD;
        lonear = fmod(lonear, 360.0) * CONST_DEG_TO_RAD;
        lonmar = fmod(lonmar, 360.0) * CONST_DEG_TO_RAD;
        lonjup = fmod(lonjup, 360.0) * CONST_DEG_TO_RAD;
        lonsat = fmod(lonsat, 360.0) * CONST_DEG_TO_RAD;
        lonurn = fmod(lonurn, 360.0) * CONST_DEG_TO_RAD;
        lonnep = fmod(lonnep, 360.0) * CONST_DEG_TO_RAD;
        precrate = fmod(precrate, 360.0) * CONST_DEG_TO_RAD;

        if (iauhelp == 'y')
        {
            printf("fa %11.7f  %11.7f  %11.7f  %11.7f  %11.7f deg \n", l * 180 / CONST_PI, l1 * 180 / CONST_PI, f * 180 / CONST_PI, d * 180 / CONST_PI, omega * 180 / CONST_PI);
            printf("fa %11.7f  %11.7f  %11.7f  %11.7f deg \n", lonmer * 180 / CONST_PI, lonven * 180 / CONST_PI, lonear * 180 / CONST_PI, lonmar * 180 / CONST_PI);
            printf("fa %11.7f  %11.7f  %11.7f  %11.7f deg \n", lonjup * 180 / CONST_PI, lonsat * 180 / CONST_PI, lonurn * 180 / CONST_PI, lonnep * 180 / CONST_PI);
            printf("fa %11.7f  \n", precrate * 180 / CONST_PI);
        }
    }  // procedure fundarg
}