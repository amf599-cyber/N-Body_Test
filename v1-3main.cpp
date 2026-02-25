#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <limits>

#include "gravtest.h"
#include "randvec.h"

/**
 * @file v1-3main.cpp
 * @brief Main entry point for the V1 particle forces simulation
 * Computes per-interaction gravitational forces now using partBucketForce from Changa
 * Outputs displacement, acceleration, and potential for each particle pair interaction
 */

/**
 * Format a Vector3D and scalar value as a CSV line with arbitrary precision
 * @param displacement 3D displacement vector
 * @param acceleration 3D acceleration vector
 * @param potential Scalar potential value
 * @return Formatted CSV string with full double precision
 */
std::string formatToCSV(const Vector3D<double>& displacement,
                        const Vector3D<double>& acceleration,
                        double potential) {
    std::ostringstream oss;
    oss << std::setprecision(std::numeric_limits<double>::max_digits10)
        << std::scientific
        << displacement.x << "," 
        << displacement.y << "," << displacement.z << ","
        << acceleration.x << "," << acceleration.y << "," 
        << acceleration.z << "," << potential;
    return oss.str();
}

/**
 * Calculate forces on a particle due to all particles in a collection
 * Modified partBucketForce method from gravity.h in Changa
 * Writes per-interaction data to CSV file
 * 
 * @param targetParticle The particle on which forces are calculated (modified)
 * @param particles Vector of source particles (const reference)
 * @param gravityCalculator GravTest instance for spline calculations
 * @param csvFile Output file stream for writing per-interaction data
 */
void partBucketForce(Particle& targetParticle, 
                     const std::vector<Particle>& particles,
                     const GravTest& gravityCalculator,
                     std::ofstream& csvFile) {
    double potentialCoeff, forceCoeff;
    
    // Loop through all particles in the collection
    for (size_t j = 0; j < particles.size(); ++j) {
        const Particle& sourceParticle = particles[j];
        
        // Calculate displacement vector from target to source particle using Vector3D methods
        const Vector3D<double> displacement = sourceParticle.position - targetParticle.position;
        double distanceSq = displacement.lengthSquared();
        
        // Combined softening length (uniform for now, dih = 0.5)
        double softeningLength = targetParticle.softening + sourceParticle.softening;
        
        // Calculate spline coefficients for force and potential
        gravityCalculator.calculateSplineForce(distanceSq, softeningLength, potentialCoeff, forceCoeff);
        
        // Calculate acceleration for this interaction
        Vector3D<double> interactionAccel = displacement * (forceCoeff * sourceParticle.mass);
        
        // Calculate potential for this interaction
        double interactionPotential = -sourceParticle.mass * potentialCoeff;
        
        // Write to CSV efficiently
        csvFile << formatToCSV(displacement, interactionAccel, interactionPotential) << "\n";
        
        // Accumulate forces using Vector3D operations
        targetParticle.acceleration += interactionAccel;
        targetParticle.potential += interactionPotential;
    }
}

int main() {
    std::cout << "=== V1 Particle Forces Simulation with Per-Interaction Data ===" << std::endl;
    
    try {
        // Initialize random vector generator
        RandomVector randomGenerator;
        randomGenerator.setSeed(27);
        
        std::cout << "\nGenerating particle population..." << std::endl;
        std::vector<Particle> particles = randomGenerator.generateParticlePopulation(NUM_PARTICLES);
        
        std::cout << "Successfully generated " << particles.size() << " particles with uniform mass " 
                  << UNIFORM_MASS << std::endl;
        
        // Initialize GravTest instance for force calculations
        const GravTest gravityCalculator;
        
        std::cout << "\n=== Computing Per-Interaction Forces ===" << std::endl;
        
        Particle& targetParticle = particles[0];
        
        // Open output file for per-interaction results
        const std::string outputFile = "particle_accelerations_test.csv";
        std::remove(outputFile.c_str());
        
        std::ofstream csvFile(outputFile);
        if (!csvFile.is_open()) {
            std::cerr << "Error: Could not open " << outputFile << " for writing." << std::endl;
            return 1;
        }
        
        // Write header for per-interaction data
        csvFile << "d_x,d_y,d_z,a_x,a_y,a_z,pot\n";
        
        // Calculate individual interactions and write to file using partBucketForce
        partBucketForce(targetParticle, particles, gravityCalculator, csvFile);
        
        csvFile.close();
        
        std::cout << "Computed " << particles.size() << " per-interaction calculations" << std::endl;
        std::cout << "\nInteraction data written to " << outputFile << std::endl;
        
        return 0;
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}