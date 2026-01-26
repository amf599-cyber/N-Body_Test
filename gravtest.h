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
#include "randvec.h"

// Physics Simulation Configuration

// Number of iterations - Configurable value
const int NUM_ITERATIONS = 100000;

// Softening parameter - Half-distance at which forces are softened between two particles
const double dih = 0.05;

// Time step per iteration - small value for accurate integration
// Total simulation time = NUM_ITERATIONS * dt
const double dt = 0.0001;
const int t = NUM_ITERATIONS;

class GravTest {
public:
    GravTest() {
    }
    
    std::vector<double> particleAccelerationsX;
    std::vector<double> particleAccelerationsY;
    std::vector<double> particleAccelerationsZ;

    /**
     * Write particle acceleration results to a CSV file for Python analysis
     * @param filename Output filename
     * @param iteration Iteration number (for time evolution tracking)
     */
    
    void writeResultsToFile(const std::string& filename, int iteration) {
        std::ofstream file;
        
        // Open in append mode so we can add iterations sequentially
        file.open(filename, std::ios::app);
        
        if (!file.is_open()) {
            std::cerr << "Error: Could not open file " << filename << " for writing." << std::endl;
            return;
        }
        
        // Write header on first iteration
        if (iteration == 0) {
            file << "iteration,particle_id,accX,accY,accZ,magnitude\n";
        }
        
        // Write acceleration data for each particle
        for (size_t i = 0; i < particleAccelerationsX.size(); ++i) {
            double accX = particleAccelerationsX[i];
            double accY = particleAccelerationsY[i];
            double accZ = particleAccelerationsZ[i];
            double magnitude = std::sqrt(accX*accX + accY*accY + accZ*accZ);
            
            file << iteration << "," << i << "," 
                 << accX << "," << accY << "," << accZ << "," 
                 << magnitude << "\n";
        }
        
        file.close();
    }
    
    /**
     * Calculate gravitational spline force based on shared distance r
     * @param r Distance between particles
     * @return Force magnitude based on softened spline kernel
     */

    double calculateSplineForce(double r) {     // Based on gravity.h Changa Code
        double u, a, b, dir;  // Local working variables

        // Cases for softened force within 2*softening length
        if (r < 2.0 * dih) {  
            u = r/dih;
        // Case for softened force within 1*softening length
            if (u < 1.0) {
                a = (dih*((7.0)/(5.0)) 
                     -((2.0)/(3.0))*u*u 
                     +((3.0)/(10.0))*u*u*u*u
                     -((1.0)/(10.0))*u*u*u*u*u);
                b = (dih*dih*dih*((4.0)/(3.0))
                     -((6.0)/(5.0))*u*u 
                     +((1.0)/(2.0))*u*u*u);
            }
        // Case for softened force between 1-2*softening length
            else {
                dir = (1.0)/r;
                a = (((-1.0)/(15.0))*dir 
                     +dih*((8.0)/(5.0)) 
                     -((4.0)/(3.0))*u*u + u*u*u
                     -((3.0)/(10.0))*u*u*u*u 
                     +((1.0)/(30.0))*u*u*u*u*u);
                b = (((-1.0)/(15.0))*dir*dir*dir 
                     +dih*dih*dih*((8.0)/(3.0)) 
                     -(3.0)*u + ((6.0)/(5.0))*u*u 
                     -((1.0)/(6.0))*u*u*u);
            }
        }
        // Distance beyond 2*softening length: standard Newtonian gravity
        else {
            a = (1.0)/r;
            b = a*a*a;
        }
        return b;
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
        // Calculate displacement vector from particle1 to particle2
        double dx = particle2.x - particle1.x;
        double dy = particle2.y - particle1.y;
        double dz = particle2.z - particle1.z;
        
        // Apply periodic boundary conditions (wrapping in 1x1x1 domain) minimum image convention
        if (dx > 0.5) dx -= 1.0;
        else if (dx < -0.5) dx += 1.0;
        
        if (dy > 0.5) dy -= 1.0;
        else if (dy < -0.5) dy += 1.0;
        
        if (dz > 0.5) dz -= 1.0;
        else if (dz < -0.5) dz += 1.0;
        
        double r = std::sqrt(dx*dx + dy*dy + dz*dz);
        
        // Calculate force from spline kernel
        double force = calculateSplineForce(r);
        
        // Gravitational acceleration: a = -m2*b*(r vector)
        // Force is attractive: points from particle2 toward particle1 (negative of displacement)
        // Guard against division by zero (same particle position)
        if (r > 0.0) {
            double factor = -particle2.mass * force / r;
            accX = factor * dx;
            accY = factor * dy;
            accZ = factor * dz;
        }
        else {
            accX = 0.0;
            accY = 0.0;
            accZ = 0.0;
        }
    }
};

#endif // GRAVTEST_H
