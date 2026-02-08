#ifndef RANDVEC_H
#define RANDVEC_H

/**
 * @file randvec.h
 * @brief Random vector generation utilities for particle initialization
 * Incorporates Vector3D.h operations for 3D vector calculations from Changa
 */

#include <vector>
#include <random>
#include <ctime>
#include <cmath>
#include <iostream>
#include <fstream>

// Constants for particle simulation
const int NUM_PARTICLES = 1000;
const double UNIFORM_MASS = 1.0;
const double COORD_MIN = 0.0;
const double COORD_MAX = 1.0;
const double INIT_VELOCITY = 0.0;
const double INIT_ACCELERATION = 0.0;

/**
 * @class Vector3D
 * @brief A templated 3D vector class for Cartesian coordinates
 * Based on Vector3D from Changa project with essential operations
 */
template <typename T = double>
class Vector3D {
public:
	typedef T componentType;
	
	T x, y, z;

	// Constructors
	Vector3D(T a = 0) : x(a), y(a), z(a) { }
	Vector3D(T a, T b, T c) : x(a), y(b), z(c) { }

	template <typename T2>
	Vector3D(T2* arr) : x(*arr), y(*(arr + 1)), z(*(arr + 2)) { }

	template <typename T2>
	Vector3D(const Vector3D<T2>& v) : x(static_cast<T>(v.x)), y(static_cast<T>(v.y)), z(static_cast<T>(v.z)) { }

	// Vector length operations
	inline T length() const {
		return static_cast<T>(std::sqrt(static_cast<double>(x * x + y * y + z * z)));
	}

	inline T lengthSquared() const {
		return x * x + y * y + z * z;
	}

	inline Vector3D<T>& normalize() {
		return *this /= length();
	}

	// Array subscripting
	inline T& operator[](int index) {
		switch(index) {
			case 0: return x;
			case 1: return y;
			case 2: return z;
		}
		return x;
	}

	// Assignment
	template <typename T2>
	inline Vector3D<T>& operator=(const Vector3D<T2>& v) {
		x = v.x;
		y = v.y;
		z = v.z;
		return *this;
	}

	// Equality
	inline bool operator==(const Vector3D<T>& v) const { 
		return (x == v.x) && (y == v.y) && (z == v.z);
	}

	inline bool operator!=(const Vector3D<T>& v) const { 
		return (x != v.x) || (y != v.y) || (z != v.z);
	}

	// Vector addition/subtraction
	inline Vector3D<T> operator+(const Vector3D<T>& v) const {
		return Vector3D<T>(x + v.x, y + v.y, z + v.z);
	}
	
	inline Vector3D<T> operator-(const Vector3D<T>& v) const {
		return Vector3D<T>(x - v.x, y - v.y, z - v.z);
	}
	
	inline Vector3D<T>& operator+=(const Vector3D<T>& v) {
		x += v.x;
		y += v.y;
		z += v.z;
		return *this;
	}
	
	inline Vector3D<T>& operator-=(const Vector3D<T>& v) {
		x -= v.x;
		y -= v.y;
		z -= v.z;
		return *this;
	}

	// Scalar multiplication
	inline Vector3D<T> operator*(const T& s) const {
		return Vector3D<T>(x * s, y * s, z * s);
	}
	
	inline Vector3D<T>& operator*=(const T& s) {
		x *= s;
		y *= s;
		z *= s;
		return *this;
	}

	// Scalar division
	inline Vector3D<T> operator/(const T& s) const {
		return Vector3D<T>(x / s, y / s, z / s);
	}
	
	inline Vector3D<T>& operator/=(const T& s) {
		x /= s;
		y /= s;
		z /= s;
		return *this;
	}

	// Unary negation
	inline Vector3D<T> operator-() const {
		return Vector3D<T>(-x, -y, -z);
	}

	// Component-wise multiplication
	inline Vector3D<T> operator*(const Vector3D<T>& v) const {
		return Vector3D<T>(x * v.x, y * v.y, z * v.z);
	}

	inline Vector3D<T>& operator*=(const Vector3D<T>& v) {
		x *= v.x;
		y *= v.y;
		z *= v.z;
		return *this;
	}

	// Component-wise division
	inline Vector3D<T> operator/(const Vector3D<T>& v) const {
		return Vector3D<T>(x / v.x, y / v.y, z / v.z);
	}

	inline Vector3D<T>& operator/=(const Vector3D<T>& v) {
		x /= v.x;
		y /= v.y;
		z /= v.z;
		return *this;
	}

	// Output operator
	friend std::ostream& operator<<(std::ostream& os, const Vector3D<T>& v) {
		os << '(' << v.x << ' ' << v.y << ' ' << v.z << ')';
		return os;
	}
};

// Scalar multiplication on the left
template <typename T, typename T2>
inline Vector3D<T> operator*(const T2& s, const Vector3D<T>& v) {
	return Vector3D<T>(v.x * s, v.y * s, v.z * s);
}

// Dot product
template <typename T, typename T2>
inline T dot(const Vector3D<T>& a, const Vector3D<T2>& b) {
	return (a.x * b.x + a.y * b.y + a.z * b.z);
}

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
        for (int i = 0; i < n; ++i) {
            Particle p;
            p.position.x = distribution(generator);
            p.position.y = distribution(generator);
            p.position.z = distribution(generator);
            p.velocity.x = INIT_VELOCITY;
            p.velocity.y = INIT_VELOCITY;
            p.velocity.z = INIT_VELOCITY;
            p.acceleration.x = INIT_ACCELERATION;
            p.acceleration.y = INIT_ACCELERATION;
            p.acceleration.z = INIT_ACCELERATION;
            p.mass = UNIFORM_MASS;
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
