#ifndef GRAVTEST_H
#define GRAVTEST_H

/**
 * @file gravtest.h
 * @brief Gravitational force calculations between particles with softened Changa Spline
 */

#include <vector>
#include <cmath>
#include <string>
#include "Vector3D.h"
#include "randvec.h"

class GravTest {
public:

    /**
     * Calculate gravitational spline coefficients a and b
     * Identical to Changa's SPLINE function in gravity.h
     * a is used for potential calculations, b is used for force/acceleration
     * @param r2 Squared distance between particles (r²)
     * @param a Output parameter for potential coefficient
     * @param b Output parameter for force/acceleration coefficient
     */

    void calculateSplineForce(double r2, double twoh, double& a, double& b) const {
        double r, u, dir;
        
        // Handle zero distance case (self-interaction)
        if (r2 == 0.0) {
            a = 0.0;
            b = 0.0;
            return;
        }
        
        r = std::sqrt(r2);

        if (r < twoh) {
            double dih_calc = 2.0 / twoh;  // Calculate dih same as Changa
            u = r * dih_calc;
            if (u < 1.0) {
                a = dih_calc * (7.0/5.0 
                       - 2.0/3.0*u*u 
                       + 3.0/10.0*u*u*u*u
                       - 1.0/10.0*u*u*u*u*u);
                b = dih_calc*dih_calc*dih_calc * (4.0/3.0 
                           - 6.0/5.0*u*u 
                           + 1.0/2.0*u*u*u);
            }
            else {
                dir = 1.0/r;
                a = -1.0/15.0*dir 
                    + dih_calc * (8.0/5.0 
                           - 4.0/3.0*u*u + u*u*u
                           - 3.0/10.0*u*u*u*u 
                           + 1.0/30.0*u*u*u*u*u);
                b = -1.0/15.0*dir*dir*dir 
                    + dih_calc*dih_calc*dih_calc * (8.0/3.0 - 3.0*u 
                           + 6.0/5.0*u*u 
                           - 1.0/6.0*u*u*u);
            }
        }
        else {
            a = 1.0/r;
            b = a*a*a;
        }
    }
};
#endif // GRAVTEST_H