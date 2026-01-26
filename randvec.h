#ifndef RANDVEC_H
#define RANDVEC_H

/**
 * @file randvec.h
 * @brief Random vector generation utilities for particle initialization
 */

#include <vector>
#include <random>
#include <ctime>

// Constants for particle simulation
const int NUM_PARTICLES = 100;
const double UNIFORM_MASS = 1.0;
const double COORD_MIN = 0.0;
const double COORD_MAX = 1.0;
const double INIT_VELOCITY = 0.0;
const double INIT_ACCELERATION = 0.0;

/**
 * @struct Particle
 * @brief Represents a particle with position, velocity, and mass
 */
struct Particle {
    double x, y, z;      // Position in 3D Cartesian coordinates
    double vx, vy, vz;   // Velocity in 3D (initially 0)
    double mass;         // Particle mass
    double accelerationX;
    double accelerationY;
    double accelerationZ;
    
    /**
     * Apply periodic boundary conditions (wrapping)
     * Particles that exceed bounds [0, 1] wrap to the other side
     */
    void applyPeriodicBoundary() {
        // Wrap x coordinate
        if (x < COORD_MIN) {
            x += (COORD_MAX - COORD_MIN);
        } else if (x > COORD_MAX) {
            x -= (COORD_MAX - COORD_MIN);
        }
        
        // Wrap y coordinate
        if (y < COORD_MIN) {
            y += (COORD_MAX - COORD_MIN);
        } else if (y > COORD_MAX) {
            y -= (COORD_MAX - COORD_MIN);
        }
        
        // Wrap z coordinate
        if (z < COORD_MIN) {
            z += (COORD_MAX - COORD_MIN);
        } else if (z > COORD_MAX) {
            z -= (COORD_MAX - COORD_MIN);
        }
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
        for (int i = 0; i < n; ++i) {
            Particle p;
            p.x = distribution(generator);
            p.y = distribution(generator);
            p.z = distribution(generator);
            p.vx = INIT_VELOCITY;  // velocity x-component
            p.vy = INIT_VELOCITY;  // velocity y-component
            p.vz = INIT_VELOCITY;  // velocity z-component
            p.mass = UNIFORM_MASS;
            p.accelerationX = INIT_ACCELERATION;
            p.accelerationY = INIT_ACCELERATION;
            p.accelerationZ = INIT_ACCELERATION;
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
    
    /**
     * Apply periodic boundary conditions to all particles in a population
     * Particles that exceed bounds [0, 1] wrap to the other side
     * @param particles Vector of particles to apply wrapping to
     */
    
    static void applyPeriodicBoundaryConditions(std::vector<Particle>& particles) {
        for (auto& particle : particles) {
            particle.applyPeriodicBoundary();
        }
    }

private:
    std::mt19937 generator;
    std::uniform_real_distribution<double> distribution;
};

#endif // RANDVEC_H
