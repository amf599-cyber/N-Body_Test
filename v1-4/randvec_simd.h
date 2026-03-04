#ifndef RANDVEC_SIMD_H
#define RANDVEC_SIMD_H

/**
 * @file randvec_simd.h
 * @brief Vectorized random particle generation directly in SoA format
 * Generates entire particle population with native SIMD-friendly data layout
 */

#include <random>
#include "particle_buffer_soa.h"
#include "Vector3D_simd.h"

// Default particle population size - 2^20 particles (factor of 2 for SIMD alignment)
constexpr size_t NUM_PARTICLES = 1048576;

/**
 * @class RandomParticleGeneratorSIMD
 * @brief Generates particles directly in SoA format using vectorized initialization
 */
class RandomParticleGeneratorSIMD {
private:
    std::mt19937 generator;
    std::uniform_real_distribution<double> dist;

public:
    RandomParticleGeneratorSIMD() 
        : generator(std::random_device{}()), // Uniform distribution constructor in [0, 1)
          dist(0.0, 1.0) {}

    /**
     * Set initialization generator random seed
     * @param seed Seed value
     */
    void setSeed(unsigned int seed) {
        generator.seed(seed);
    }

    /**
     * Generate particle population with batched random initialization
     * Pre-generates batches of random numbers for better cache locality and SIMD alignment
     * 
     * @param n Number of particles to generate
     * @return ParticleBufferSOA with all particles initialized in SoA layout
     */
    ParticleBufferSOA generatePopulationVectorized(int n = NUM_PARTICLES) {
        ParticleBufferSOA buffer;
        buffer.reserve(n);

        constexpr size_t BATCH_SIZE = 1024;  // Pre-generate in batches for cache efficiency
        std::vector<double> batch(BATCH_SIZE * 3);
        
        for (int i = 0; i < n; i += BATCH_SIZE) {
            int batch_count = std::min((int)BATCH_SIZE, n - i);
            
            // Pre-generate random numbers for this batch
            for (int j = 0; j < batch_count * 3; ++j) {
                batch[j] = dist(generator);
            }
            
            // Add particles from pre-generated batch
            for (int j = 0; j < batch_count; ++j) {
                buffer.addParticle(
                    Vector3D<double>(batch[j*3], batch[j*3 + 1], batch[j*3 + 2]),
                    1.0,     // mass
                    0.05     // softening
                );
            }
        }

        return buffer;
    }
};

#endif // RANDVEC_SIMD_H
