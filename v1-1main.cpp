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
 * Perform one iteration of gravitational force calculation and update particle positions/velocities
 * @param particles Vector of particles to update
 * @param gravityCalculator GravTest instance for calculations
 * @param iteration Current iteration number
 * @param verbose Whether to print detailed output
 */

void performForceCalculationIteration(std::vector<Particle>& particles, GravTest& gravityCalculator, 
                                      int iteration, bool verbose = false) {
    // Reset particle accelerations to zero at the start of each iteration
    // This ensures we only accumulate forces from this iteration
    for (size_t i = 0; i < particles.size(); ++i) {
        particles[i].acceleration.x = 0.0;
        particles[i].acceleration.y = 0.0;
        particles[i].acceleration.z = 0.0;
    }
    
    if (verbose) {
        // std::cout << "\n=== Iteration " << iteration+1 << " ===" << std::endl;
    }
    
    // Nested loop: for each particle, calculate forces from all other particles
    for (size_t i = 0; i < particles.size(); ++i) {
        if (verbose && i < 3) {  // Only show first 3 particles in detail
            // std::cout << "Particle " << i+1 << " - calculating accelerations from all other particles:" << std::endl;
        }
        
        for (size_t j = 0; j < particles.size(); ++j) {
            if (i != j) {  // Skip self-interaction
                // Calculate acceleration from particle j on particle i
                double accX = 0.0, accY = 0.0, accZ = 0.0;
                gravityCalculator.calculateAcceleration(particles[i], particles[j], 
                                                        accX, accY, accZ);
                // Accumulate in particle
                particles[i].acceleration.x += accX;
                particles[i].acceleration.y += accY;
                particles[i].acceleration.z += accZ;
            }
        }
        
        if (verbose && i < 3) {
            // std::cout << "  Total acceleration on Particle " << i+1 << ": (" 
            //           << totalAccX << ", " << totalAccY << ", " << totalAccZ << ")" << std::endl;
        }
    }
    
    // Print particle data: positions, accelerations, and potentials
    std::cout << "\n=== Particle Data ===" << std::endl;
    std::cout << "delta_x, delta_y, delta_z, acc_x, acc_y, acc_z, pot_x, pot_y, pot_z" << std::endl;
    for (size_t i = 0; i < particles.size(); ++i) {
        double dx = particles[i].position.x;
        double dy = particles[i].position.y;
        double dz = particles[i].position.z;
        double accX = particles[i].acceleration.x;
        double accY = particles[i].acceleration.y;
        double accZ = particles[i].acceleration.z;
        double potential = gravityCalculator.calculatePotential(particles, i);
        double accMagnitude = particles[i].acceleration.length();
        
        // Calculate potential components by decomposing along acceleration direction
        double potX = 0.0, potY = 0.0, potZ = 0.0;
        if (accMagnitude > 0.0) {
            double factor = potential / accMagnitude;
            potX = factor * accX;
            potY = factor * accY;
            potZ = factor * accZ;
        }
        
        std::cout << dx << ", " << dy << ", " << dz << ", "
                  << accX << ", " << accY << ", " << accZ << ", "
                  << potX << ", " << potY << ", " << potZ << std::endl;
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
        
        // Initialize GravTest instance for force calculations
        GravTest gravityCalculator;
        
        // Open output file and prepare for multiple iterations
        std::string outputFile = "particle_accelerations.csv";
        std::remove(outputFile.c_str());  // Clear any existing file
        
        // Time loop: perform force calculations for t iterations
        for (int iteration = 0; iteration < (int)t; ++iteration) {
            performForceCalculationIteration(particles, gravityCalculator, iteration, 
                                            iteration == 0);  // Only verbose for first iteration
            
            // Write results to file after each iteration
            gravityCalculator.writeResultsToFile(outputFile, iteration, particles);
            
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

