# V1 Project
Alex Feucht University of Washington Comprehensive Physics and Astronomy Undergraduate

## Overview
V1 is the first version of the N-Body Vectorization project I am working on for Professor Quinn. This project utilizes some existing code from the Changa Project, namely gravity.h and vector3D.h. This is the first attempt at implementing a better vectorization solution to these coding modules using updated versions of the C++ Standard Template Library. All unaltered and repurposed code remains attributed to original authors.

## Using partBucketForce code for inner loop (v1.3)
- Removed existing inner for loop from main in exchange for an additional method, partBucketForce from gravity.h
- Added additional CSV helper method for more efficient and less redundant output to file inside partBucketForce
- Significantly cleaned up redundant code blocks and made sure all of the main methods were being used primarily
- Validated that the same analytic outputs were achieved for both potential and acceleration vs. separation distance
- in Python plots in the greater than 2 * dih range, with softening occuring at separations below this
- Cleaned up variable definitions and access between blocks of code and includes
- Added more descriptive variable names in some cases where it could add clarity
- Overall formatting for readability and ease of following code structure

## Particle Pair Interactions only (v1.2)
- Testing individual particle pair interactions between first particle and all others, including self-interaction case
- Removed epsilon threshold when passing to gravtest.h, including radius of 0 or greater (up to max separation distance)
- Condensed multiple vector operations into smaller/more efficient lines of code using Vector3D class
- Included entire unaltered Vector3D.h and TypeSelection.h headers from Changa code base
- Using only inner loop of nested for loop to highlight softening and Newtonian analytical results for r>=2h
- Fixed csv ofstream to only include separation distances, abs. magnitude acceleration, magnitude potential
- Fixed some of the logic in the gravtest.h to more closely resemble exact code from Changa Gravity.h
- Ran Python analysis on 1,000,000 particles and their interactions with first particle to validate methods

## Recent Refactoring (v1.1)
- Integrated Vector3D template class from Changa project into randvec.h for efficient 3D vector operations
- Refactored Particle struct to use Vector3D objects (position, velocity, acceleration) instead of individual component doubles
- Implemented output stream operator for Particle class
- Streamlined file output to include only positional data (delta_x, delta_y, delta_z) and uniform mass via std:ofstream
- Commented out periodic and global boundary conditions for future implementation of existing Changa code modules
- Set NUM_ITERATIONS to 1 for validity and efficiency testing of single iteration of force/potential calculations

## Features
- **gravtest.h**: Gravitational force calculations with softened Changa Spline kernel
  - Configurable softening parameter (dih) and time step (dt)
  - Customizable NUM_ITERATIONS for simulation length
  - CSV output of particle position and mass data
- **randvec.h**: Random vector generation and Vector3D operations
  - Particle population initialization with uniform mass distribution
  - Complete Vector3D template class with standard operations
  - Vector algebra: addition, subtraction, scalar multiplication/division, dot product, length calculations
  - Periodic boundary condition support

## Building

### Prerequisites
- C++17 or later
- GCC or Clang compiler
- Make

### Compilation
```bash
make
```

To rebuild from scratch:
```bash
make rebuild
```

To clean build artifacts:
```bash
make clean
```

## Project Structure
```
v1/
├── gravtest.h      # Gravitational force testing interface
├── randvec.h       # Random vector generation interface
├── Vector3D.h      # Changa vector operations code
├── TypeSelection.h # Changa code required due to Vector3D.h include
├── Makefile        # Build configuration
├── v1-2main.cpp      # Main source file
└── README.md       # This file
```

## Usage
After building, run the application:
```bash
./v1_app
```