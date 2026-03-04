# V1 Project
Alex Feucht University of Washington Comprehensive Physics and Astronomy Undergraduate

## Overview
V1 is the first version of the N-Body Vectorization project I am working on for Professor Quinn. This project utilizes some existing code from the Changa Project, namely gravity.h and vector3D.h. This is the first attempt at implementing a better vectorization solution to these coding modules using updated versions of the C++ Standard Template Library. All unaltered and repurposed code remains attributed to original authors.

## SIMD Vectorization Overhaul (v1.4)
Complete rewrite of N-body force calculation engine with SIMD optimization and SoA memory layout:

### Architecture Changes
- **Structure of Arrays (SoA) Layout**: Reorganized particle data from AoS to columnar SoA format
  - Position data: separate arrays for x, y, z coordinates
  - Acceleration: separate acc_x, acc_y, acc_z arrays
  - Enables stride-1 access patterns critical for SIMD vectorization
  - Cache-efficient columnar data layout maximizes memory bandwidth utilization

### SIMD Implementation
- **Automatic Architecture Detection**: Compile-time detection of available SIMD capabilities
  - AVX-512: 8-wide double precision (512 bits)
  - AVX2: 4-wide double precision (256 bits) - default on modern x86_64
  - SSE2: 2-wide double precision (128 bits) - fallback
  - Scalar: Single-element fallback for unsupported architectures
- **Explicit SIMD Intrinsics**: Direct use of architecture-specific intrinsics (_mm256, _mm512 families)
- **Branchless Operations**: Avoids branch misprediction in hot loops using ternary operators and where-semantics
- **Batch Vectorization**: Processes multiple particles in single SIMD register for improved throughput

### Core Components
- **simd_traits.h**: Trait-based SIMD configuration providing compile-time constants for vector width and alignment
- **particle_buffer_soa.h**: SIMD-aligned memory allocator and columnar particle storage container
- **Vector3D_simd.h**: Packed vector operations (distance, dot product, cross product) on SoA data
- **gravtest_simd.h**: Vectorized gravitational force and potential calculation with softening regime classification
- **randvec_simd.h**: Batched random particle generation directly in SoA format

### Performance Features
- SIMD-aligned memory allocation (32-byte for AVX2, 64-byte for AVX-512)
- Reduced memory bandwidth requirements through cache-friendly SoA layout
- Minimized data dependencies enabling instruction-level parallelism
- Pre-computation of spline coefficients in vectorized batches
- Compiled with `-O3 -march=native` for maximum optimization

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

### SIMD Optimizations
- **Vectorized Arithmetic**: All core operations (distance, acceleration, potential) vectorized for multi-element parallelism
- **Architecture-Aware Compilation**: Automatic selection of best available SIMD instruction set
- **Intrinsic-Based Operations**: Direct SIMD intrinsics for distance calculations, dot products, and scalar operations
- **Memory Alignment**: SIMD-aware allocators ensure 32/64-byte aligned buffers for unaligned access penalties

### Gravitational Physics
- **Softening Spline Kernel**: Continuous Changa-style spline softening avoiding singularities at r=0
- **Multi-Regime Softening**: Automatic classification and handling of 4 softening regimes:
  - Self-Interaction (r = 0): No net acceleration or potentials
  - Inner kernel (0 < r < h): Full cubic spline interpolation
  - Outer kernel (h ≤ r < 2h): Modified spline for smoothness
  - Newtonian (r ≥ 2h): 1/r gravitational approximation
- **Particle Interactions**: Branchless pair-force calculations for all N² particle interactions
- **Potential & Acceleration**: Simultaneous computation of gravitational potential energy and acceleration vectors

### Data Management
- **Structure of Arrays**: Column-oriented memory layout optimized for SIMD access patterns
- **Scalable Population**: Support for millions of particles (default 2²⁰ = 1,048,576 particles)
- **CSV Output**: Efficient batched CSV writing of particle accelerations and diagnostics for testing
- **Aligned Memory**: Custom allocators ensuring cache-optimal memory alignment

### Testing & Validation
- **Per-Interaction Forces**: Full computation of pairwise gravitational forces with analytical validation done in Python
- **Softening Validation**: Verification against Newtonian predictions in different softening regimes
- **CSV Diagnostics**: Output includes separation distances, acceleration magnitudes, and potential energies
- **Performance Reporting**: Runtime display of SIMD configuration and kernel statistics

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
├── Core SIMD Implementation (v1.4)
│   ├── gravtest_simd.h         # SIMD-optimized gravitational force calculations
│   ├── randvec_simd.h          # Vectorized random particle generation in SoA format
│   ├── Vector3D_simd.h         # SIMD-accelerated 3D vector operations
│   ├── particle_buffer_soa.h   # Structure of Arrays (SoA) particle buffer layout
│   └── simd_traits.h           # SIMD configuration and trait detection
│
├── Legacy Implementation (v1.1-1.3)
│   ├── gravtest.h              # Original gravitational force testing interface
│   ├── randvec.h               # Original random vector generation interface
│   ├── Vector3D.h              # Changa vector operations code
│   └── TypeSelection.h         # Changa code required by Vector3D.h
│
├── Main Applications
│   ├── v1-4main.cpp            # Current main: SIMD with SoA layout (v1.4)
│   ├── v1-3main.cpp            # Previous: partBucketForce implementation (v1.3)
│   ├── v1-2main.cpp            # Previous: Particle pair interactions (v1.2)
│   └── v1main.cpp              # Original implementation
│
├── Build and Documentation
│   ├── Makefile                # Build configuration
│   └── README.md               # This file
│
└── Output Files
    └── particle_accelerations_test.csv  # Computed particle acceleration data
```

## Usage
After building, run the application:
```bash
./v1_app
```