# V1 Project
Alex Feucht University of Washington Comprehensive Physics and Astronomy Undergraduate

## Overview
V1 is the first version of the N-Body Vectorization project I am working on for Professor Quinn. This project utilizes some existing code from the Changa Project, namely gravity.h and vector3D. This is the first attempt at implementing a better vectorization solution to these coding modules using updated versions of the C++ Standard Template Library.

## Features
- **gravtest.h**: Gravitational force testing utilities
- **randvec.h**: Random vector generation for particle population initialization

## Building

### Prerequisites
- C++11 or later
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