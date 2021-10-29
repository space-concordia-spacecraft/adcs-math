#include "adcs/Utils.h"

namespace adcs {
    // older sgp4ext methods
    /* -----------------------------------------------------------------------------
    *
    *                           function gstime_SGP4
    *
    *  this function finds the greenwich sidereal time.
    *
    *  author        : david vallado                  719-573-2600    1 mar 2001
    *
    *  inputs          description                    range / units
    *    jdut1       - julian date in ut1             days from 4713 bc
    *
    *  outputs       :
    *    gstime      - greenwich sidereal time        0 to 2pi rad
    *
    *  locals        :
    *    temp        - temporary variable for doubles   rad
    *    tut1        - julian centuries from the
    *                  jan 1, 2000 12 h epoch (ut1)
    *
    *  coupling      :
    *    none
    *
    *  references    :
    *    vallado       2013, 187, eq 3-45
    * --------------------------------------------------------------------------- */

    double  gstime(double jdut1)
    {
        double       temp, tut1;

        tut1 = (jdut1 - 2451545.0) / 36525.0;
        temp = -6.2e-6* tut1 * tut1 * tut1 + 0.093104 * tut1 * tut1 +
               (876600.0 * 3600 + 8640184.812866) * tut1 + 67310.54841;  // sec
        temp = fmod(temp * CONST_DEG_TO_RAD / 240.0, CONST_TWO_PI); //360/86400 = 1/240, to deg, to rad

        // ------------------------ check quadrants ---------------------
        if (temp < 0.0)
            temp += CONST_TWO_PI;

        return temp;
    }  // gstime
}