#ifndef GRAVTEST_H
#define GRAVTEST_H

/**
 * @file gravtest.h
 * @brief Gravitational force calculations between particles with softened Changa Spline
 */

#include <vector>
#include <cmath>
#include <fstream>
#include <string>
#include "Vector3D.h"
#include "randvec.h"

// Physics Simulation Configuration

// Number of iterations - Configurable value
const int NUM_ITERATIONS = 1;

// Softening parameter - Half-distance at which forces are softened between two particles
const double dih = 0.05;

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

    void calculateSplineForce(double r2, double& a, double& b) {
        double r, u, dir;
        double twoh = 2*dih;  // dih is the softening length (half-distance)
        
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
    
    /**
     * Calculate acceleration components due to gravitational force from one particle to another
     * @param particle1 The particle experiencing the acceleration
     * @param particle2 The source particle generating the shared gravitational force
     * @param accX Output acceleration in x direction
     * @param accY Output acceleration in y direction
     * @param accZ Output acceleration in z direction
     */

    void calculateAcceleration(const Particle& particle1, const Particle& particle2,
                               double& accX, double& accY, double& accZ) {
        Vector3D<double> displacement = particle2.position - particle1.position;
        double r2 = displacement.lengthSquared();
        if (r2 > 0) {
            double r = std::sqrt(r2);
            double a, b;
            calculateSplineForce(r2, a, b);
            Vector3D<double> acc = displacement * (-particle2.mass * b);
            accX = acc.x;
            accY = acc.y;
            accZ = acc.z;
        }
        else {
            accX = 0.0;
            accY = 0.0;
            accZ = 0.0;
        }
    }
    
    /**
     * Calculate gravitational potential for a particle due to all others
     * Optimized to use r² for potential vectorization
     * @param particles Vector of all particles
     * @param particleIndex Index of the particle to calculate potential for
     * @return Gravitational potential (scalar) for the particle
     */
    double calculatePotential(const std::vector<Particle>& particles, size_t particleIndex) {
        double potential = 0.0;
        const Particle& particle1 = particles[particleIndex];
        
        for (size_t j = 0; j < particles.size(); ++j) {
            if (particleIndex != j) {
                const Particle& particle2 = particles[j];
                Vector3D<double> displacement = particle2.position - particle1.position;
                double r2 = displacement.lengthSquared();
                if (r2 > 0) {
                    double a, b;
                    calculateSplineForce(r2, a, b);
                    double r = std::sqrt(r2);
                    potential -= particle2.mass * a;
                }
            }
        }
        
        return potential;
    }
};
#endif // GRAVTEST_H