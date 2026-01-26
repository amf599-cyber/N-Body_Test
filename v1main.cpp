#include <iostream>
#include <vector>
#include <cmath>
#include <memory>
#include <ctime>
#include <cstdio>

#include "gravtest.h"
#include "randvec.h"

/**
 * @file v1main.cpp
 * @brief Main entry point for the V1 particle forces simulation
 * Contains implementations of RandomVector and GravTest classes
 */

/**
 * Calculate Euclidean distance between two particles
 * @param p1 First particle (the one we're calculating forces on)
 * @param p2 Second particle (source of gravitational force)
 * @return Distance between particles
 */

// Calculate Euclidean distance between two particles in 3D space
double calculateEuclideanDistance(const Particle& p1, const Particle& p2) {
    double dx = p1.x - p2.x;
    double dy = p1.y - p2.y;
    double dz = p1.z - p2.z;
    return std::sqrt(dx*dx + dy*dy + dz*dz);
}

/**
 * Perform one iteration of gravitational force calculation and update particle positions/velocities
 * @param particles Vector of particles to update
 * @param gravityCalculator GravTest instance for calculations
 * @param iteration Current iteration number
 * @param verbose Whether to print detailed output
 */

void performForceCalculationIteration(std::vector<Particle>& particles, GravTest& gravityCalculator, 
                                      int iteration, bool verbose = false) {
    // Initialize temporary vectors to accumulate new accelerations for this iteration
    // Previous accelerations are retained and will be added to these
    std::vector<double> newAccelerationX(particles.size(), 0.0);
    std::vector<double> newAccelerationY(particles.size(), 0.0);
    std::vector<double> newAccelerationZ(particles.size(), 0.0);
    
    // Reset particle accelerations to zero at the start of each iteration
    // This ensures it only accumulate forces from this iteration, not previous ones
    for (size_t i = 0; i < particles.size(); ++i) {
        particles[i].accelerationX = 0.0;
        particles[i].accelerationY = 0.0;
        particles[i].accelerationZ = 0.0;
    }
    
    if (verbose) {
        std::cout << "\n=== Iteration " << iteration+1 << " ===" << std::endl;
    }
    
    // Nested loop: for each particle, calculate forces from all other particles
    for (size_t i = 0; i < particles.size(); ++i) {
        double totalAccX = 0.0;
        double totalAccY = 0.0;
        double totalAccZ = 0.0;
        
        if (verbose && i < 3) {  // Only show first 3 particles in detail for iteration > 0
            std::cout << "Particle " << i+1 << " - calculating accelerations from all other particles:" << std::endl;
        }
        
        for (size_t j = 0; j < particles.size(); ++j) {
            if (i != j) {  // Skip self-interaction
                // Calculate acceleration from particle j on particle i
                double accX = 0.0, accY = 0.0, accZ = 0.0;
                gravityCalculator.calculateAcceleration(particles[i], particles[j], 
                                                        accX, accY, accZ);
                
                // Sum the accelerations
                totalAccX += accX;
                totalAccY += accY;
                totalAccZ += accZ;
                
                // Accumulate new accelerations
                newAccelerationX[i] += accX;
                newAccelerationY[i] += accY;
                newAccelerationZ[i] += accZ;
                
                // Also accumulate in particle (for velocity updates below)
                particles[i].accelerationX += accX;
                particles[i].accelerationY += accY;
                particles[i].accelerationZ += accZ;
            }
        }
        
        // Store cumulative accelerations in GravTest for file output
        // This includes all accumulated forces from this and all previous iterations
        gravityCalculator.particleAccelerationsX[i] = particles[i].accelerationX;
        gravityCalculator.particleAccelerationsY[i] = particles[i].accelerationY;
        gravityCalculator.particleAccelerationsZ[i] = particles[i].accelerationZ;
        
        if (verbose && i < 3) {
            std::cout << "  Total acceleration on Particle " << i+1 << ": (" 
                      << totalAccX << ", " << totalAccY << ", " << totalAccZ << ")" << std::endl;
        }
    }
    
    // Update velocities and positions based on accelerations
    for (size_t i = 0; i < particles.size(); ++i) {
        // v = v + a*dt
        particles[i].vx += particles[i].accelerationX * dt;
        particles[i].vy += particles[i].accelerationY * dt;
        particles[i].vz += particles[i].accelerationZ * dt;
        
        // x = x + v*dt
        particles[i].x += particles[i].vx * dt;
        particles[i].y += particles[i].vy * dt;
        particles[i].z += particles[i].vz * dt;
        
        // Apply periodic boundary conditions
        particles[i].applyPeriodicBoundary();
    }
}

int main() {
    std::cout << "=== V1 Particle Forces Simulation ===" << std::endl;
    
    try {
        // Initialize random vector generator
        RandomVector randGen;
        randGen.setSeed(27);
        
        std::cout << "\nGenerating particle population..." << std::endl;
        std::vector<Particle> particles = randGen.generateParticlePopulation(NUM_PARTICLES);
        
        std::cout << "Successfully generated " << particles.size() << " particles with uniform mass " 
                  << UNIFORM_MASS << std::endl;
        
        // Display sample particles
        std::cout << "\nSample particle positions and velocities:" << std::endl;
        for (size_t i = 0; i < particles.size() && i < 10; ++i) {
            std::cout << "  Particle " << i+1 << ": pos=(" << particles[i].x << ", " 
                      << particles[i].y << ", " << particles[i].z << ")"
                      << " vel=(" << particles[i].vx << ", " << particles[i].vy << ", " 
                      << particles[i].vz << ")"
                      << " mass=" << particles[i].mass << std::endl;
        }
        
        // Initialize GravTest instance for force calculations
        GravTest gravityCalculator;
        gravityCalculator.particleAccelerationsX.resize(particles.size(), 0.0);
        gravityCalculator.particleAccelerationsY.resize(particles.size(), 0.0);
        gravityCalculator.particleAccelerationsZ.resize(particles.size(), 0.0);
        
        // Open output file and prepare for multiple iterations
        std::string outputFile = "particle_accelerations.csv";
        std::remove(outputFile.c_str());  // Clear any existing file
        
        std::cout << "\n=== Starting Time Evolution Simulation ===" << std::endl;
        std::cout << "Running " << t << " iterations of force calculation..." << std::endl;
        std::cout << "Physical parameters:" << std::endl;
        std::cout << "  dih (softening parameter) = " << dih << std::endl;
        std::cout << "  dt (time step) = " << dt << std::endl;
        std::cout << "  Total simulation time = " << (t * dt) << std::endl;
        
        // Time loop: perform force calculations for t iterations
        for (int iteration = 0; iteration < (int)t; ++iteration) {
            performForceCalculationIteration(particles, gravityCalculator, iteration, 
                                            iteration == 0);  // Only verbose for first iteration
            
            // Write results to file after each iteration
            gravityCalculator.writeResultsToFile(outputFile, iteration);
            
            if ((iteration + 1) % 10 == 0) {
                std::cout << "  Completed iteration " << (iteration + 1) << " of " << (int)t << std::endl;
            }
        }
        
        std::cout << "\n=== All iterations completed ===" << std::endl;
        std::cout << "Results written to " << outputFile << std::endl;
        std::cout << "File contains " << (int)t << " iterations x " << particles.size() << " particles\n";
        return 0;
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}

