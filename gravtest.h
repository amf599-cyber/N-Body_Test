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
const int NUM_ITERATIONS = 1;

// Softening parameter - Half-distance at which forces are softened between two particles
const double dih = 0.05;

// Time step per iteration - small value for accurate integration
// Total simulation time = NUM_ITERATIONS * dt
// Not using time evolution in this version, but defined for future implementation
// const double dt = 0.0001;

// Setting total simulation run steps for nested loops
const int t = NUM_ITERATIONS;

class GravTest {
public:
    GravTest() {
    }

    /**
     * Write particle data to a CSV file with specific columns
     * Output format: delta_x, delta_y, delta_z, acc_x, acc_y, acc_z, pot_x, pot_y, pot_z
     * @param filename Output filename
     * @param iteration Iteration number (for time evolution tracking)
     * @param particles Reference to the particles vector for position and acceleration data
     */
    
    void writeResultsToFile(const std::string& filename, int iteration, const std::vector<Particle>& particles) {
        std::ofstream file;
        
        // Open in append mode so we can add iterations sequentially
        file.open(filename, std::ios::app);
        
        if (!file.is_open()) {
            std::cerr << "Error: Could not open file " << filename << " for writing." << std::endl;
            return;
        }
        
        // Write header on first iteration
        if (iteration == 0) {
            file << "delta_x,delta_y,delta_z,acc_x,acc_y,acc_z,pot_x,pot_y,pot_z\n";
        }
        
        // Write particle data for each particle
        for (size_t i = 0; i < particles.size(); ++i) {
            double dx = particles[i].position.x;
            double dy = particles[i].position.y;
            double dz = particles[i].position.z;
            double accX = particles[i].acceleration.x;
            double accY = particles[i].acceleration.y;
            double accZ = particles[i].acceleration.z;
            double potential = calculatePotential(particles, i);
            double accMagnitude = particles[i].acceleration.length();
            
            // Calculate potential components by decomposing along acceleration direction
            double potX = 0.0, potY = 0.0, potZ = 0.0;
            if (accMagnitude > 0.0) {
                double factor = potential / accMagnitude;
                potX = factor * accX;
                potY = factor * accY;
                potZ = factor * accZ;
            }
            
            file << dx << "," << dy << "," << dz << ","
                 << accX << "," << accY << "," << accZ << ","
                 << potX << "," << potY << "," << potZ << "\n";
        }
        
        file.close();
    }
    
    /**
     * Calculate gravitational spline force based on squared distance r²
     * Optimized to avoid computing square root outside the function
     * @param r2 Squared distance between particles (r²)
     * @return Force magnitude based on softened spline kernel
     */

    double calculateSplineForce(double r2) {     // Based on gravity.h Changa Code
        double r = std::sqrt(r2);  // Calculate r from r² only once
        double u, a, b, dir;  // Local working variables
        double dih2 = dih * dih;  // Pre-calculate dih²

        // Cases for softened force within 2*softening length
        if (r2 < 4.0 * dih2) {  // Equivalent to r < 2*dih squaring each side
            u = r / dih;
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
        double dx = particle2.position.x - particle1.position.x;
        double dy = particle2.position.y - particle1.position.y;
        double dz = particle2.position.z - particle1.position.z;
                
        // Boundary conditions applied, now compute r² and pass to spline
        double r2 = dx*dx + dy*dy + dz*dz;
        
        // Calculate force from spline kernel (function handles sqrt internally)
        double force = calculateSplineForce(r2);
        
        // Gravitational acceleration: a = -m2*b*(r vector)
        // Force is attractive: points from particle2 toward particle1 (negative of displacement)
        // Guard against division by zero (same particle position)
        if (r2 > 1e-20) {  // Avoid singularity in acceleration calculation with epsilon threshold
            double r = std::sqrt(r2);
            double factor = -particle2.mass * force / r;
            accX = factor * dx;
            accY = factor * dy;
            accZ = factor * dz;
        }
        // Case for zero vector/self-interaction of particles: no acceleration
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
            if (particleIndex != j) {  // Skip self-interaction
                const Particle& particle2 = particles[j];
                
                // Calculate displacement vector
                double dx = particle2.position.x - particle1.position.x;
                double dy = particle2.position.y - particle1.position.y;
                double dz = particle2.position.z - particle1.position.z;
                
                // Compute r² (avoid computing r until necessary)
                double r2 = dx*dx + dy*dy + dz*dz;
                
                // Calculate force from spline kernel using r²
                if (r2 > 1e-20) {  // Avoid singularity in potential calculation with epsilon threshold
                    double force = calculateSplineForce(r2);
                    double r = std::sqrt(r2);
                    // Gravitational potential: φ = -m/r (multiplied by spline kernel b)
                    // The spline function b is already dimensionless, so we divide by r to get proper potential scaling
                    potential -= particle2.mass * force / r;
                }
            }
        }
        
        return potential;
    }
};
#endif // GRAVTEST_H