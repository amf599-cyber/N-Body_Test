# V1 Project Makefile
CXX = g++
CXXFLAGS = -std=c++17 -Wall -Wextra -O2
TARGET = v1_app

# Source files
SOURCES = v1-3main.cpp
OBJECTS = $(SOURCES:.cpp=.o)
PRECOMPILED = gravtest.o randvec.o

# Header dependencies
HEADERS = gravtest.h randvec.h

# Default target
all: $(TARGET)

# Link object files to create executable
$(TARGET): $(OBJECTS) $(PRECOMPILED)
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
