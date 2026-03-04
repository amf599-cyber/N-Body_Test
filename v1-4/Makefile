# V1 Project Makefile - Version History
# V1.0-V1.3: Scalar implementation with gravtest.h, randvec.h
# V1.4: SIMD-optimized implementation with vectorized operations
#       - Architecture-aware SIMD detection (AVX2, SSE2, scalar fallback)
#       - Structure of Arrays (SoA) memory layout
#       - Branchless where-based conditional masking
#       - Explicit SIMD intrinsics for guaranteed vectorization
#       - Direct SoA particle generation

CXX = g++
CXXFLAGS = -std=c++17 -Wall -Wextra -O3 -march=native
TARGET = v1-4

# Source files for V1.4 SIMD implementation
SOURCES = v1-4main.cpp
OBJECTS = $(SOURCES:.cpp=.o)

# Header dependencies for V1.4
# simd_traits.h: Architecture detection, Vector3D, SIMDMask, AlignedAllocator
# particle_buffer_soa.h: Structure of Arrays particle storage
# gravtest_simd.h: Vectorized spline calculations with where-based masking
# Vector3D_SIMD.h: Packed vector operations with explicit SIMD intrinsics
# randvec_simd.h: Direct SoA particle generation
HEADERS = simd_traits.h particle_buffer_soa.h gravtest_simd.h Vector3D_simd.h randvec_simd.h

# Default target
all: $(TARGET)

# Link object files to create executable
$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^ -lm

# Compile source files to object files
%.o: %.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -c $<

# Clean build artifacts
clean:
	rm -f $(OBJECTS) $(TARGET) particle_accelerations_test.csv

# Rebuild from scratch
rebuild: clean all

# Run the SIMD-optimized simulation
run: $(TARGET)
	./$(TARGET)

# Phony targets
.PHONY: all clean rebuild run
