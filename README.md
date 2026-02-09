# V1 Project
Alex Feucht University of Washington Comprehensive Physics and Astronomy Undergraduate

## Overview
V1 is the first version of the N-Body Vectorization project I am working on for Professor Quinn. This project utilizes some existing code from the Changa Project, namely gravity.h and vector3D.h. This is the first attempt at implementing a better vectorization solution to these coding modules using updated versions of the C++ Standard Template Library.

## Recent Refactoring (v1.1)
- Integrated Vector3D template class from Changa project into randvec.h for efficient 3D vector operations
- Refactored Particle struct to use Vector3D objects (position, velocity, acceleration) instead of individual component doubles
- Implemented output stream operator for Particle class
- Streamlined file output to include only components of positional data, acceleration data, potential data, and uniform mass via std:ofstream
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
├── Makefile        # Build configuration
├── v1main.cpp      # Main source file
└── README.md       # This file
```

## Usage
After building, run the application:
```bash
./v1_app
```
