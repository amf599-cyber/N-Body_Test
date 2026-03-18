# V1 Project Makefile - SIMD N-Body Force Simulation
# Supports native SIMD vectorization (AVX, NEON, SVE, scalar fallback)
# Pure SIMD pipeline: Generation → SoA Processing → CSV Output
CXX = g++
CXXFLAGS = -std=c++17 -Wall -Wextra -O2 -march=native
TARGET = vf_app

# Source files
SOURCES = v-fmain.cpp
OBJECTS = $(SOURCES:.cpp=.o)

# Header dependencies
HEADERS = gravtest-f.h randvec-f.h TypeSelection.h Vector3D.h

# Default target
all: $(TARGET)

# Link object files to create executable
$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^

# Compile source files to object files
%.o: %.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -c $<

# Clean build artifacts
clean:
	rm -f $(OBJECTS) $(TARGET)

# Rebuild from scratch
rebuild: clean all

# Phony targets
.PHONY: all clean rebuild
