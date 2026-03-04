#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <limits>
#include <chrono>

#include "simd_traits.h"
#include "particle_buffer_soa.h"
#include "gravtest_simd.h"
#include "randvec_simd.h"
#include "Vector3D_simd.h"

/**
 * @file v1-4main.cpp
 * @brief SIMD-optimized particle forces simulation with Structure of Arrays layout
 * 
 * Complete rewrite of partBucketForce using:
 * - Structure of Arrays (SoA) memory layout for SIMD-friendly access
 * - Vectorized computation with boolean masking utilized in gravtest_simd.h
 * - Chunked processing based on system SIMD width (SSE2, AVX2, AVX-512)
 * - Batched CSV output to reduce I/O overhead
 * 
 * Supported architectures:
 *   - AVX-512: 8 doubles per SIMD operation (512-bit vectors)
 *   - AVX2/AVX: 4 doubles per SIMD operation (256-bit vectors)
 *   - SSE2: 2 doubles per SIMD operation (128-bit vectors)
 *   - Scalar: 1 double fallback
 * 
 * @author Alex Feucht University of Washington Physics/Astronomy Undergraduate
 * @version 1.4
 */

/**
 * @struct CSVResultBatch
 * @brief Batch buffer for CSV output to reduce I/O overhead
 * 
 * Instead of writing each interaction individually, accumulate results
 * and write in batches to reduce syscall overhead
 */
struct CSVResultBatch {
    std::vector<Vector3D<double>> displacements;
    std::vector<Vector3D<double>> accelerations;
    std::vector<double> potentials;
    
    static constexpr size_t BATCH_SIZE = 8192;  // Tune based on system
    
    /**
     * Check if batch is full and should be flushed
     */
    bool isFull() const {
        return displacements.size() >= BATCH_SIZE;
    }
    
    /**
     * Clear all batch data
     */
    void clear() {
        displacements.clear();
        accelerations.clear();
        potentials.clear();
    }
    
    /**
     * Get current batch size
     */
    size_t size() const {
        return displacements.size();
    }
};

/**
 * Format Vector3D and scalar values as CSV row with full double precision
 * Does not use vectorized SIMD formatting since this is only for testing output
 * Not performance-critical, only happens once per set of interactions
 * Scales similarly either way, could implement SIMD explicit method if needed
 * @param displacement Packed 3D displacement vector
 * @param acceleration 3D acceleration vector
 * @param potential Scalar potential value
 * @return Formatted CSV string
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
 * Write batched results to CSV file
 * @param batch Result batch to write
 * @param csvFile Output file stream
 */
void writeBatchToCSV(const CSVResultBatch& batch, std::ofstream& csvFile) {
    for (size_t i = 0; i < batch.size(); ++i) {
        csvFile << formatToCSV(batch.displacements[i], 
                              batch.accelerations[i], 
                              batch.potentials[i]) << "\n";
    }
}

/**
 * SIMD-Optimized partBucketForce implementation
 * 
 * Processes particle-particle interactions with vectorized computation.
 * Chunks source particles by SIMD vector width and applies where-based
 * masking for different softening regimes.
 * 
 * @param target_particle_idx Index of target particle in source buffer
 * @param source_buffer SoA particle buffer (source particles)
 * @param gravityCalculator GravTestSIMD instance
 * @param csvFile Output CSV file stream
 */
void partBucketForceVectorized(
    size_t target_particle_idx,
    ParticleBufferSOA& source_buffer,
    const GravTestSIMD& gravityCalculator,
    std::ofstream& csvFile)
{
    constexpr size_t VECTOR_WIDTH = DefaultSIMDTraits::VECTOR_WIDTH;
    size_t num_particles = source_buffer.size();
    
    if (num_particles == 0) return;

    // Extract target particle data (scalar broadcast to all SIMD lanes)
    double target_x = source_buffer.pos_x[target_particle_idx];
    double target_y = source_buffer.pos_y[target_particle_idx];
    double target_z = source_buffer.pos_z[target_particle_idx];
    double target_softening = source_buffer.softening[target_particle_idx];
    
    double& target_acc_x = source_buffer.acc_x[target_particle_idx];
    double& target_acc_y = source_buffer.acc_y[target_particle_idx];
    double& target_acc_z = source_buffer.acc_z[target_particle_idx];
    double& target_potential = source_buffer.potential[target_particle_idx];

    CSVResultBatch batch;
    
    // Process particles in chunks of VECTOR_WIDTH
    for (size_t chunk_start = 0; chunk_start < num_particles; chunk_start += VECTOR_WIDTH) {
        size_t chunk_size = std::min(VECTOR_WIDTH, num_particles - chunk_start);
        
        // Prepare vectors for this chunk
        std::vector<double> disp_x(chunk_size), disp_y(chunk_size), disp_z(chunk_size);
        std::vector<double> softenings(chunk_size), masses(chunk_size);
        SIMDMask mask(chunk_size);
        
        // Vectorized displacement calculation (SIMD-optimized)
        Vec3SIMD::displacement(
            target_x, target_y, target_z,
            &source_buffer.pos_x[chunk_start],
            &source_buffer.pos_y[chunk_start],
            &source_buffer.pos_z[chunk_start],
            disp_x.data(), disp_y.data(), disp_z.data(),
            chunk_size
        );
        
        // Load scalar metadata (softenings, masses) in one pass
        for (size_t lane = 0; lane < chunk_size; ++lane) {
            mask.valid[lane] = true;
            softenings[lane] = target_softening + source_buffer.softening[chunk_start + lane];
            masses[lane] = source_buffer.mass[chunk_start + lane];
        }
        
        // Vectorized force calculation with where-based masking
        std::vector<double> accel_x, accel_y, accel_z, potential;
        gravityCalculator.calculateForcesVectorized(
            disp_x, disp_y, disp_z, softenings, masses, mask,
            accel_x, accel_y, accel_z, potential, chunk_size
        );
        
        // Vectorized accumulation (no branches in hot path)
        for (size_t lane = 0; lane < chunk_size; ++lane) {
            target_acc_x += accel_x[lane];
            target_acc_y += accel_y[lane];
            target_acc_z += accel_z[lane];
            target_potential += potential[lane];
            batch.displacements.push_back(Vector3D<double>(disp_x[lane], disp_y[lane], disp_z[lane]));
            batch.accelerations.push_back(Vector3D<double>(accel_x[lane], accel_y[lane], accel_z[lane]));
            batch.potentials.push_back(potential[lane]);
        }
        
        // Flush batch only after chunk processing (reduces branch frequency)
        if (batch.isFull()) {
            writeBatchToCSV(batch, csvFile);
            batch.clear();
        }
    }
    
    // Flush remaining batch data
    if (!batch.displacements.empty()) {
        writeBatchToCSV(batch, csvFile);
    }
}

/**
 * Print SIMD capability information
 */
void printSIMDInfo() {
    std::cout << "\n=== SIMD Configuration ===" << std::endl;
    std::cout << "SIMD Level: " << DefaultSIMDTraits::SIMD_LEVEL << std::endl;
    std::cout << "Vector Width: " << DefaultSIMDTraits::VECTOR_WIDTH << " doubles" << std::endl;
    std::cout << "Byte Width: " << DefaultSIMDTraits::BYTE_WIDTH << " bytes" << std::endl;
    std::cout << "Memory Alignment: " << DefaultSIMDTraits::ALIGNMENT << " bytes" << std::endl;
    std::cout << "Vectorized: " << (DefaultSIMDTraits::VECTORIZED ? "Yes" : "No") << std::endl;
}

int main() {
    try {
        printSIMDInfo();
        
        std::cout << "Generating particles..." << std::endl;
        RandomParticleGeneratorSIMD particleGenerator;
        particleGenerator.setSeed(27);
        ParticleBufferSOA particleBuffer = particleGenerator.generatePopulationVectorized(NUM_PARTICLES);

        std::cout << "Computing forces..." << std::endl;
        const GravTestSIMD gravityCalculator;
        
        std::cout << "\n=== Computing Per-Interaction Forces (SIMD-Optimized) ===" << std::endl;
        
        // Process target particle 0
        size_t target_idx = 0;
        
        // Open output file for per-interaction results
        const std::string outputFile = "particle_accelerations_test.csv";
        std::remove(outputFile.c_str());
        
        std::ofstream csvFile(outputFile);
        if (!csvFile.is_open()) {
            std::cerr << "Error: Could not open " << outputFile << " for writing." << std::endl;
            return 1;
        }
        
        // Write CSV header
        csvFile << "d_x,d_y,d_z,a_x,a_y,a_z,pot\n";

        // Calculate forces using SIMD-optimized chunking and masking
        partBucketForceVectorized(target_idx, particleBuffer, gravityCalculator, csvFile);
        
        csvFile.close();
        
        std::cout << "Complete. Output written to " << outputFile << std::endl;
        
        return 0;
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}
