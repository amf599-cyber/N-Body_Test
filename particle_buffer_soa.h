#ifndef PARTICLE_BUFFER_SOA_H
#define PARTICLE_BUFFER_SOA_H

/**
 * @file particle_buffer_soa.h
 * @brief Structure of Arrays (SoA) particle buffer for SIMD-efficient memory access
 * Reorganizes particle data from AoS (Array of Structures) to SoA layout to maximize
 * cache efficiency and enable vectorized operations without stride penalties
 * 
 * Memory layout (SoA):
 * [x0 x1 x2 x3 ... xn] [y0 y1 y2 y3 ... yn] [z0 z1 z2 z3 ... zn]
 * 
 * vs. AoS:
 * [x0 y0 z0 ax0 ay0 az0 ...] [x1 y1 z1 ax1 ay1 az1 ...] ...
 */

#include <vector>
#include <stdexcept>
#include "simd_traits.h"

/**
 * @class ParticleBufferSOA
 * @brief Structure of Arrays particle storage with SIMD-aligned memory
 * 
 * Stores particles in columnar format (all x-coordinates together, all y-coordinates, etc)
 * enabling efficient vectorized access patterns and cache locality optimization
 */
class ParticleBufferSOA {
public:
    // Allocator type for SIMD-aligned memory
    template <typename T>
    using AlignedVector = std::vector<T, AlignedAllocator<T, DefaultSIMDTraits::ALIGNMENT>>;

    /**
     * Particle properties stored in SoA layout
     */
    AlignedVector<double> pos_x, pos_y, pos_z;    // Position coordinates (3 arrays)
    AlignedVector<double> acc_x, acc_y, acc_z;    // Acceleration (accumulator)
    AlignedVector<double> mass;                    // Mass per particle
    AlignedVector<double> potential;               // Potential energy per particle
    AlignedVector<double> softening;               // Softening parameter per particle
    
    size_t num_particles;

    /**
     * Default constructor
     */
    ParticleBufferSOA() : num_particles(0) {}

    /**
     * Reserve capacity for particles (preallocate SIMD-aligned memory)
     * @param capacity Number of particles to reserve space for
     */
    void reserve(size_t capacity) {
        pos_x.reserve(capacity);
        pos_y.reserve(capacity);
        pos_z.reserve(capacity);
        acc_x.reserve(capacity);
        acc_y.reserve(capacity);
        acc_z.reserve(capacity);
        mass.reserve(capacity);
        potential.reserve(capacity);
        softening.reserve(capacity);
    }

    /**
     * Add a single particle to the buffer
     * Converts from AoS (Particle struct) to SoA layout
     * @param position Particle position vector
     * @param particle_mass Particle mass
     * @param softening_param Softening length parameter
     */
    void addParticle(const Vector3D<double>& position,
                     double particle_mass,
                     double softening_param) {
        pos_x.push_back(position.x);
        pos_y.push_back(position.y);
        pos_z.push_back(position.z);
        
        acc_x.push_back(0.0);
        acc_y.push_back(0.0);
        acc_z.push_back(0.0);
        
        mass.push_back(particle_mass);
        potential.push_back(0.0);
        softening.push_back(softening_param);
        
        num_particles++;
    }

    /**
     * Get number of particles in buffer
     * @return Number of particles
     */
    size_t size() const { return num_particles; }

    /**
     * Check if buffer is empty
     * @return True if num_particles == 0
     */
    bool empty() const { return num_particles == 0; }

    /**
     * Clear all particle data and reset size
     */
    void clear() {
        pos_x.clear();
        pos_y.clear();
        pos_z.clear();
        acc_x.clear();
        acc_y.clear();
        acc_z.clear();
        mass.clear();
        potential.clear();
        softening.clear();
        num_particles = 0;
    }

    /**
     * Get position of particle at index i
     * @param i Particle index
     * @return Vector3D position
     */
    Vector3D<double> getPosition(size_t i) const {
        if (i >= num_particles) throw std::out_of_range("Particle index out of range");
        return Vector3D<double>(pos_x[i], pos_y[i], pos_z[i]);
    }

    /**
     * Get acceleration of particle at index i
     * @param i Particle index
     * @return Vector3D acceleration
     */
    Vector3D<double> getAcceleration(size_t i) const {
        if (i >= num_particles) throw std::out_of_range("Particle index out of range");
        return Vector3D<double>(acc_x[i], acc_y[i], acc_z[i]);
    }

    /**
     * Get potential of particle at index i
     * @param i Particle index
     * @return Potential energy
     */
    double getPotential(size_t i) const {
        if (i >= num_particles) throw std::out_of_range("Particle index out of range");
        return potential[i];
    }

    /**
     * Add acceleration to particle at index i (for accumulation during force calculations)
     * @param i Particle index
     * @param accel Acceleration vector to add
     */
    void addAcceleration(size_t i, const Vector3D<double>& accel) {
        if (i >= num_particles) throw std::out_of_range("Particle index out of range");
        acc_x[i] += accel.x;
        acc_y[i] += accel.y;
        acc_z[i] += accel.z;
    }

    /**
     * Add potential to particle at index i
     * @param i Particle index
     * @param pot Potential energy to add
     */
    void addPotential(size_t i, double pot) {
        if (i >= num_particles) throw std::out_of_range("Particle index out of range");
        potential[i] += pot;
    }

    /**
     * Get raw aligned pointer to position X data
     * Used for direct vectorized operations
     * @return Pointer to contiguous x-coordinate data
     */
    double* getPosXPtr() { return pos_x.data(); }
    const double* getPosXPtr() const { return pos_x.data(); }

    /**
     * Get raw aligned pointer to position Y data
     * @return Pointer to contiguous y-coordinate data
     */
    double* getPosYPtr() { return pos_y.data(); }
    const double* getPosYPtr() const { return pos_y.data(); }

    /**
     * Get raw aligned pointer to position Z data
     * @return Pointer to contiguous z-coordinate data
     */
    double* getPosZPtr() { return pos_z.data(); }
    const double* getPosZPtr() const { return pos_z.data(); }

    /**
     * Get raw aligned pointer to mass data
     * @return Pointer to contiguous mass data
     */
    double* getMassPtr() { return mass.data(); }
    const double* getMassPtr() const { return mass.data(); }

    /**
     * Get raw aligned pointer to softening data
     * @return Pointer to contiguous softening data
     */
    double* getSofteningPtr() { return softening.data(); }
    const double* getSofteningPtr() const { return softening.data(); }

    /**
     * Get raw aligned pointer to acceleration X data
     * @return Pointer to contiguous acceleration x-coordinate data
     */
    double* getAccXPtr() { return acc_x.data(); }

    /**
     * Get raw aligned pointer to acceleration Y data
     * @return Pointer to contiguous acceleration y-coordinate data
     */
    double* getAccYPtr() { return acc_y.data(); }

    /**
     * Get raw aligned pointer to acceleration Z data
     * @return Pointer to contiguous acceleration z-coordinate data
     */
    double* getAccZPtr() { return acc_z.data(); }

    /**
     * Get raw aligned pointer to potential data
     * @return Pointer to contiguous potential data
     */
    double* getPotentialPtr() { return potential.data(); }

    /**
     * Reset all accumulated values (acceleration, potential)
     * Useful before computing new force interactions
     */
    void resetAccumulators() {
        std::fill(acc_x.begin(), acc_x.end(), 0.0);
        std::fill(acc_y.begin(), acc_y.end(), 0.0);
        std::fill(acc_z.begin(), acc_z.end(), 0.0);
        std::fill(potential.begin(), potential.end(), 0.0);
    }

    /**
     * Get total number of bytes allocated for this buffer
     * @return Total memory footprint in bytes
     */
    size_t getMemoryUsage() const {
        return (pos_x.capacity() + pos_y.capacity() + pos_z.capacity() +
                acc_x.capacity() + acc_y.capacity() + acc_z.capacity() +
                mass.capacity() + potential.capacity() + softening.capacity()) * sizeof(double);
    }
};

#endif // PARTICLE_BUFFER_SOA_H
