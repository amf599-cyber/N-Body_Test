#ifndef RANDVEC_H
#define RANDVEC_H

/**
 * @file randvec.h
 * @brief Random vector generation utilities for particle initialization
 * Uses Vector3D from Vector3D.h for 3D vector calculations
 */

#include <vector>
#include <random>
#include <ctime>
#include <iostream>
#include "Vector3D.h"

// Constants for particle simulation
const int NUM_PARTICLES = 1000000;  
const double UNIFORM_MASS = 1.0;
const double COORD_MIN = 0.0;
const double COORD_MAX = 1.0;
const double INIT_VELOCITY = 0.0;
const double INIT_ACCELERATION = 0.0;

/**
 * @struct Particle
 * @brief Represents a particle with position, velocity, and mass
 * Uses Vector3D for efficient vector operations
 */
struct Particle {
    Vector3D<double> position;   // Position in 3D Cartesian coordinates
    Vector3D<double> velocity;   // Velocity in 3D (initially 0)
    Vector3D<double> acceleration;  // Acceleration in 3D
    double mass;                 // Particle mass (uniform = 1.0)
    
    // Constructor
    Particle() 
        : position(0.0), velocity(0.0), acceleration(0.0), mass(UNIFORM_MASS) {}
    
    /**
     * Calculate displacement from this particle to another
     * @param other The other particle
     * @return Displacement vector
     */
    Vector3D<double> displacementTo(const Particle& other) const {
        return other.position - position;
    }
    
    /**
     * Output stream operator for file output
     * Outputs: particle_id, delta_x, delta_y, delta_z, mass
     */
    friend std::ostream& operator<<(std::ostream& os, const Particle& p) {
        os << p.position.x << "," << p.position.y << "," << p.position.z << "," << p.mass;
        return os;
    }
};

class RandomVector {
public:
    RandomVector() 
        : generator(std::time(nullptr)), 
          distribution(COORD_MIN, COORD_MAX) {
    }
    
    /**
     * Generate a population of particles with uniform mass
     * @param n Number of particles to generate (default: NUM_PARTICLES)
     * @return Vector of Particle objects randomly distributed in 1x1x1 space
     */
    std::vector<Particle> generateParticlePopulation(int n = NUM_PARTICLES) {
        std::vector<Particle> particles;
        particles.reserve(n);
        for (int i = 0; i < n; ++i) {
            Particle p;
            p.position = Vector3D<double>(distribution(generator), 
                                          distribution(generator), 
                                          distribution(generator));
            particles.push_back(p);
        }
        return particles;
    }
    
    /**
     * Set seed for random number generator
     */
    void setSeed(unsigned int seed) {
        generator.seed(seed);
    }

private:
    std::mt19937 generator;
    std::uniform_real_distribution<double> distribution;
};

#endif // RANDVEC_H
