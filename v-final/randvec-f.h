#ifndef RANDVEC_H
#define RANDVEC_H

/**
 * @file randvec.h
 * @brief Random vector generation utilities for particle initialization
 * Supports both scalar and SIMD-vectorized particle generation
 * Uses Vector3D from Vector3D.h for 3D vector calculations
 */

#include <vector>
#include <random>
#include <iostream>
#include <experimental/simd>
#include "Vector3D.h"

namespace stdx = std::experimental;

// Constants for particle simulation
const int NUM_PARTICLES = 268435456;  // Total particles to generate (2^28 for extended scale test, crashes at 2^29 on 16GB RAM)
const double UNIFORM_MASS = 1.0;
const double dih = 0.05;  // Softening parameter - Half-distance at which forces are softened

/**
 * @struct Particle
 * @brief Represents a particle with position, acceleration, and mass
 * Uses Vector3D for efficient vector operations
 */
struct Particle {
    Vector3D<double> position;   // Position in 3D Cartesian coordinates
    Vector3D<double> acceleration;  // Acceleration in 3D
    double mass;                 // Particle mass (uniform = 1.0)
    double potential;            // Gravitational potential energy
    double softening;            // Softening length for force calculations
    
    // Constructor
    Particle() 
        : position(0.0), acceleration(0.0), mass(UNIFORM_MASS), potential(0.0), softening(dih) {}
    
    // Constructor with all parameters
    Particle(const Vector3D<double>& pos, const Vector3D<double>& accel, 
             double m, double pot, double soft)
        : position(pos), acceleration(accel), mass(m), potential(pot), softening(soft) {}
};

/**
 * @struct ParticleVector
 * @brief SIMD-optimized container for vectorized particle data (Structure of Arrays - SoA)
 * Holds position components in SIMD vector format for efficient vectorized processing
 * Internally uses SoA layout but can convert to/from AoS (Array of Structures) implicitly
 */
template <typename SIMDType>
struct ParticleVector {
    using value_type = typename SIMDType::value_type;
    using simd_type = SIMDType;
    
    SIMDType pos_x;     // X positions packed in SIMD vector
    SIMDType pos_y;     // Y positions packed in SIMD vector
    SIMDType pos_z;     // Z positions packed in SIMD vector
    SIMDType accel_x;   // X accelerations packed in SIMD vector
    SIMDType accel_y;   // Y accelerations packed in SIMD vector
    SIMDType accel_z;   // Z accelerations packed in SIMD vector
    SIMDType mass_vec;  // Mass values packed in SIMD vector
    SIMDType potential_vec;  // Potential values packed in SIMD vector
    SIMDType softening_vec;  // Softening lengths packed in SIMD vector
    
    // Constructor
    ParticleVector() 
        : pos_x(0.0), pos_y(0.0), pos_z(0.0),
          accel_x(0.0), accel_y(0.0), accel_z(0.0),
          mass_vec(UNIFORM_MASS), potential_vec(0.0), softening_vec(dih) {}
    
    /**
     * Convert scalar Particle (AoS) to element in ParticleVector (SoA)
     * @param p Particle to store at index idx
     * @param idx Index within the SIMD vector
     */
    void setParticle(const Particle& p, int idx) {
        pos_x[idx] = p.position.x;
        pos_y[idx] = p.position.y;
        pos_z[idx] = p.position.z;
        accel_x[idx] = p.acceleration.x;
        accel_y[idx] = p.acceleration.y;
        accel_z[idx] = p.acceleration.z;
        mass_vec[idx] = p.mass;
        potential_vec[idx] = p.potential;
        softening_vec[idx] = p.softening;
    }
    
    /**
     * Convert element in ParticleVector (SoA) to scalar Particle (AoS)
     * @param idx Index within the SIMD vector
     * @return Particle object extracted from SoA layout
     */
    Particle getParticle(int idx) const {
        Particle p;
        p.position = Vector3D<double>(pos_x[idx], pos_y[idx], pos_z[idx]);
        p.acceleration = Vector3D<double>(accel_x[idx], accel_y[idx], accel_z[idx]);
        p.mass = mass_vec[idx];
        p.potential = potential_vec[idx];
        p.softening = softening_vec[idx];
        return p;
    }
};

class RandomVector {
public:
    RandomVector() 
        : rng(std::random_device{}()) {
    }
    
    /**
     * Generate a population of SIMD-vectorized particles (SoA format)
     * Particles are created as scalar objects then implicitly packed into Structure of Arrays format
     * @tparam SIMDType Native SIMD type (e.g., floatv, doublev)
     * @param n Number of logical particles to generate
     * @return Vector of SIMD-packed ParticleVector objects in SoA format
     * 
     * Key features:
     * - Internal generation uses scalar uniform distribution between COORD_MIN and COORD_MAX for positions
     * - Implicit AoS→SoA conversion via setParticle() method
     * - Produces ParticleVector objects ready for vectorized force calculations
     * - Handles padding for non-aligned particle counts
     */
    template <typename SIMDType>
    std::vector<ParticleVector<SIMDType>> generateParticlePopulationSIMD(int n = NUM_PARTICLES) {
        std::vector<ParticleVector<SIMDType>> particleVectors;
        const int simdWidth = SIMDType::size();
        const int numVectors = (n + simdWidth - 1) / simdWidth;
        
        particleVectors.reserve(numVectors);
        
        std::uniform_real_distribution<double> uniform_real(0.0, 1.0);
        
        for (int v = 0; v < numVectors; ++v) {
            ParticleVector<SIMDType> pv;
            const int vectorParticleCount = std::min(simdWidth, n - v * simdWidth);
            
            // Generate particles and implicitly pack into SoA vectors
            for (int i = 0; i < vectorParticleCount; ++i) {
                // Generate random particle as AoS - condensed to single operations
                Particle p(
                    Vector3D<double>(uniform_real(rng), uniform_real(rng), uniform_real(rng)),
                    Vector3D<double>(0.0, 0.0, 0.0),
                    1.0, 0.0, dih
                );
                pv.setParticle(p, i);
            }
            
            // Pad remaining SIMD lanes if n is not divisible by simdWidth
            for (int i = vectorParticleCount; i < simdWidth; ++i) {
                pv.setParticle(Particle(), i);
            }
            
            particleVectors.push_back(pv);
        }
        
        return particleVectors;
    }
    
    /**
     * Set seed for random number generator
     */
    void setSeed(unsigned int seed) {
        rng.seed(seed);
    }

private:
    std::mt19937 rng;
};

#endif // RANDVEC_H