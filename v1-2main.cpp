#include <iostream>
#include <vector>
#include <cmath>

#include "gravtest.h"
#include "randvec.h"

/**
 * @file v1main.cpp
 * @brief Main entry point for the V1 particle forces simulation
 * Contains implementations of RandomVector and GravTest classes
 * Using pre-existing Vector3D/TypeSelection classes from Changa 
 * for vector operations
 */

/**
 * Calculate gravitational potential between two particles
 * @param particle1 The first particle
 * @param particle2 The second particle
 * @param gravityCalculator GravTest instance for calculations
 * @return Gravitational potential between the two particles
 */
double calculatePotentialBetweenPair(const Particle& particle1, const Particle& particle2,
                                      GravTest& gravityCalculator) {
    // Calculate displacement vector
    Vector3D<double> displacement = particle2.position - particle1.position;
    double r2 = displacement.lengthSquared();
    
    // Avoid singularity at zero distance
    if (r2 > 0) {
        double r = std::sqrt(r2);
        double a, b;
        gravityCalculator.calculateSplineForce(r2, a, b);
        // Gravitational potential: φ = -m*a
        // a is the potential coefficient from the spline kernel
        return -particle2.mass * a;
    }
    return 0.0;
}

/**
 * Write particle pair interaction data to CSV file
 * @param file Output file stream (must be open)
 * @param distance Separation distance between particles
 * @param accMagnitude Magnitude of acceleration vector
 * @param potMagnitude Magnitude of potential vector
 */
void writeTestParticlePairToFile(std::ofstream& file, double distance,
                                 double accMagnitude, double potMagnitude) {
    file << distance << "," << accMagnitude << "," << potMagnitude << "\n";
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
        
        // Test particle pair interactions and write to file
        std::cout << "\n=== Testing Particle Pair Interactions ===" << std::endl;
        std::cout << "Writing particle pair interactions for first particle with all " << NUM_PARTICLES << " particles to CSV: particle_accelerations_test.csv" << std::endl;
        
        std::string testOutputFile = "particle_accelerations_test.csv";
        std::remove(testOutputFile.c_str());
        
        std::ofstream csvFile(testOutputFile);
        if (!csvFile.is_open()) {
            std::cerr << "Error: Could not open " << testOutputFile << " for writing." << std::endl;
            return 1;
        }
        
        // Write header
        csvFile << "distance,acc_magnitude,pot_magnitude\n";
        
        // void partForce(Particle *part, std::vector<Particles> &particles){}
        // Particles operating on, the vector of all particles
        // *part replaces particles[0] line 92 except inside function
        // updates acceleration in the Particle struct, add softening to same class
        // Initialize softening, and use in calculation of 2h
        // vectorize lines 93-103 from v1main (SIMD)

        // Test particle pair interactions including self-interaction case
        for (size_t j = 0; j < particles.size(); ++j) {
            // Calculate displacement and acceleration
            Vector3D<double> displacement = particles[j].position - particles[0].position;
            double accX = 0.0, accY = 0.0, accZ = 0.0;
            gravityCalculator.calculateAcceleration(particles[0], particles[j], accX, accY, accZ);
            
            // Write to CSV file with calculated magnitude and potential
            writeTestParticlePairToFile(csvFile, displacement.length(), 
                std::sqrt(accX * accX + accY * accY + accZ * accZ),
                calculatePotentialBetweenPair(particles[0], particles[j], gravityCalculator));
        }
        
        csvFile.close();
        
        std::cout << "\nParticle pair interactions written to " << testOutputFile << std::endl;
        
        return 0;
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}