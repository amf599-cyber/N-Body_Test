#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <limits>

#include "gravtest.h"
#include "randvec.h"

#include <experimental/simd>
namespace stdx = std::experimental;
using doublev = stdx::rebind_simd_t<double, stdx::native_simd<float>>;

/**
 * @file v1-3main.cpp
 * @brief Main entry point for the V1 particle forces simulation
 * Computes per-interaction gravitational forces using native SIMD vectorization
 * Outputs displacement, acceleration, and potential for each particle pair interaction
 */

/**
 * Format vectorized results for CSV output
 * Extracts individual elements from SIMD vectors and formats as CSV
 * @param disp_x, disp_y, disp_z Displacement components in SIMD format
 * @param accel_x, accel_y, accel_z Acceleration components in SIMD format
 * @param pot_vec Potential values in SIMD format
 * @param idx Index within the SIMD vector to extract
 * @return Formatted CSV string with full double precision
 */
template <typename SIMDType>
std::string formatToCSVSIMD(const SIMDType& disp_x, const SIMDType& disp_y, 
                            const SIMDType& disp_z, const SIMDType& accel_x, 
                            const SIMDType& accel_y, const SIMDType& accel_z,
                            const SIMDType& pot_vec, int idx) {
    std::ostringstream oss;
    oss << std::setprecision(std::numeric_limits<double>::max_digits10)
        << std::scientific
        << disp_x[idx] << ","
        << disp_y[idx] << "," << disp_z[idx] << ","
        << accel_x[idx] << "," << accel_y[idx] << ","
        << accel_z[idx] << "," << pot_vec[idx];
    return oss.str();
}

/**
 * Calculate forces on a target particle due to all particles in SIMD SoA format
 * Processes particles in their native Structure of Arrays (SoA) SIMD representation
 * No manual AoS->SoA conversion needed - particles initialized in optimal vectorized format
 * Automatically adapts to available SIMD capabilities (x86 AVX/AVX-512, ARM NEON/SVE, scalar)
 * 
 * @param targetParticle The particle on which forces are calculated 
 * @param particleVectors Vector of SoA ParticleVector<doublev> objects (const reference)
 * @param gravityCalculator GravTest instance for spline calculations
 * @param csvFile Output file stream for writing per-interaction test data
 */
void partBucketForce(Particle& targetParticle,
                     const std::vector<ParticleVector<doublev>>& particleVectors,
                     const GravTest& gravityCalculator,
                     std::ofstream& csvFile) {
    // Store target position for displacement calculations
    const Vector3D<double> target_pos = targetParticle.position;
    
    // Accumulator for totals using Vector3D
    Vector3D<double> total_accel(0.0);
    double total_potential = 0.0;
    
    const int vector_width = doublev::size();
    
    // Process each SoA vector batch directly (no manual AoS->SoA conversion needed)
    for (const auto& pv : particleVectors) {
        // Calculate displacement vectors directly from SoA data as Vector3D<doublev>
        Vector3D<doublev> displacement(pv.pos_x - target_pos.x,
                                       pv.pos_y - target_pos.y,
                                       pv.pos_z - target_pos.z);
        
        // Calculate squared distances using Vector3D::lengthSquared()
        doublev dist_sq = displacement.lengthSquared();
        
        // Extract displacement components for later use
        doublev disp_x = displacement.x;
        doublev disp_y = displacement.y;
        doublev disp_z = displacement.z;
        
        // Calculate spline coefficients using vectorized operation
        doublev potential_coeff, force_coeff;
        gravityCalculator.calculateSplineForceSIMD<doublev>(dist_sq, dih * 2.0, potential_coeff, force_coeff);
        
        // Calculate acceleration components and potential in condensed form
        doublev accel_x = force_coeff * pv.mass_vec * disp_x;
        doublev accel_y = force_coeff * pv.mass_vec * disp_y;
        doublev accel_z = force_coeff * pv.mass_vec * disp_z;
        doublev potential_interaction = -pv.mass_vec * potential_coeff;
        
        // Extract and write per-interaction data to CSV, accumulating totals
        for (int i = 0; i < vector_width; ++i) {
            csvFile << formatToCSVSIMD(disp_x, disp_y, disp_z,
                                      accel_x, accel_y, accel_z,
                                      potential_interaction, i) << "\n";
            
            // Accumulate totals using Vector3D semantics
            total_accel.x += accel_x[i];
            total_accel.y += accel_y[i];
            total_accel.z += accel_z[i];
            total_potential += potential_interaction[i];
        }
    }
    
    // Update target particle with accumulated values
    targetParticle.acceleration = total_accel;
    targetParticle.potential = total_potential;
}

int main() {
    std::cout << "=== V1 Particle Forces Simulation with Per-Interaction Data ===" << std::endl;
    std::cout << "SIMD Vector Width: " << doublev::size() << " elements" << std::endl;
    
    try {
        // Initialize random vector generator
        RandomVector randomGenerator;
        randomGenerator.setSeed(27);
        
        std::cout << "\nGenerating particle population in SIMD SoA format..." << std::endl;
        std::vector<ParticleVector<doublev>> particleVectors = randomGenerator.generateParticlePopulationSIMD<doublev>(NUM_PARTICLES);
        
        std::cout << "Successfully generated " << NUM_PARTICLES << " particles in " 
                  << particleVectors.size() << " SIMD vectors (vector width: " << doublev::size() << ")" << std::endl;
        
        // Get first particle from SoA vectors for target calculation
        // Convert first element of first SoA vector back to AoS for target comparison
        // Eventually can be reworked to outer for loop to iterate over all particles/time steps
        Particle targetParticle = particleVectors[0].getParticle(0);
        
        // Initialize GravTest instance for force calculations
        const GravTest gravityCalculator;
        
        std::cout << "\n=== Computing Per-Interaction Forces ===" << std::endl;
        
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
        
        // Compute forces using SIMD SoA vectors directly
        // Process particles in their native vectorized format
        partBucketForce(targetParticle, particleVectors, gravityCalculator, csvFile);
        
        csvFile.close();
        
        std::cout << "Computed " << NUM_PARTICLES << " per-interaction calculations" << std::endl;
        std::cout << "\nInteraction data written to " << outputFile << std::endl;

        return 0;
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}