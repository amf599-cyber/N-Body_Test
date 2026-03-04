#ifndef SIMD_TRAITS_H
#define SIMD_TRAITS_H

/**
 * @file simd_traits.h
 * @brief SIMD architecture detection and vector width configuration
 * Automatically detects available SIMD capabilities and provides compile-time constants
 * Supports SSE2 (4-wide), AVX/AVX2 (8-wide), and AVX-512 (16-wide) for double precision
 */

#include <cstddef>
#include <cstdint>
#include <vector>
#include <limits>

// Detect available SIMD capabilities at compile time
#if defined(__AVX512F__) && defined(__AVX512D__)
    #define SIMD_AVX512 1
    #define SIMD_AVX2 0
    #define SIMD_SSE2 0
#elif defined(__AVX2__) || defined(__AVX__)
    #define SIMD_AVX512 0
    #define SIMD_AVX2 1
    #define SIMD_SSE2 0
#elif defined(__SSE2__) || defined(_M_X64) || defined(_M_AMD64)
    #define SIMD_AVX512 0
    #define SIMD_AVX2 0
    #define SIMD_SSE2 1
#else
    #define SIMD_AVX512 0
    #define SIMD_AVX2 0
    #define SIMD_SSE2 0
#endif

/**
 * @struct SIMDTraits
 * @brief Compile-time SIMD configuration for double precision floating point
 * Provides vector width, alignment, and feature flags based on detected architecture
 */
template <typename T = double>
struct SIMDTraits {
    // Vector width: number of scalar elements processed in parallel
    #if SIMD_AVX512
        static constexpr size_t VECTOR_WIDTH = 8;      // 8 doubles × 64 bits = 512 bits
        static constexpr size_t BYTE_WIDTH = 64;       // 512 bits = 64 bytes
        static constexpr const char* SIMD_LEVEL = "AVX-512";
    #elif SIMD_AVX2
        static constexpr size_t VECTOR_WIDTH = 4;      // 4 doubles × 64 bits = 256 bits
        static constexpr size_t BYTE_WIDTH = 32;       // 256 bits = 32 bytes
        static constexpr const char* SIMD_LEVEL = "AVX2";
    #elif SIMD_SSE2
        static constexpr size_t VECTOR_WIDTH = 2;      // 2 doubles × 64 bits = 128 bits
        static constexpr size_t BYTE_WIDTH = 16;       // 128 bits = 16 bytes
        static constexpr const char* SIMD_LEVEL = "SSE2";
    #else
        static constexpr size_t VECTOR_WIDTH = 1;      // Scalar fallback
        static constexpr size_t BYTE_WIDTH = 8;
        static constexpr const char* SIMD_LEVEL = "Scalar";
    #endif

    // Alignment requirements for efficient loading/storing
    static constexpr size_t ALIGNMENT = BYTE_WIDTH;
    
    // Feature flags
    static constexpr bool HAS_AVX512 = (SIMD_AVX512 == 1);
    static constexpr bool HAS_AVX2 = (SIMD_AVX2 == 1);
    static constexpr bool HAS_SSE2 = (SIMD_SSE2 == 1);
    static constexpr bool VECTORIZED = (VECTOR_WIDTH > 1);
};

// Convenience typedefs
using DefaultSIMDTraits = SIMDTraits<double>;

/**
 * @struct Vector3D
 * @brief Minimal 3D vector for AoS interop and batch output
 * Used in v1-4main.cpp for CSV output formatting only
 */
template <typename T = double>
struct Vector3D {
    T x, y, z;
    
    Vector3D() : x(0), y(0), z(0) {}
    Vector3D(T x_, T y_, T z_) : x(x_), y(y_), z(z_) {}
};

/**
 * @struct SIMDMask
 * @brief Boolean mask for SIMD conditional operations
 * Each lane holds a scalar boolean representing validity/condition for vectorized operations
 */
struct SIMDMask {
    enum class SofteningRegime : uint8_t {
        SELF_INTERACTION = 0,      // r = 0: Skip entirely
        INNER_KERNEL = 1,          // 0 < r < h: Use spline regime 1
        OUTER_KERNEL = 2,          // h ≤ r < 2h: Use spline regime 2
        NEWTONIAN = 3              // r ≥ 2h: Use Newtonian formula (1/r)
    };

    std::vector<SofteningRegime> regimes;  // Per-lane softening regime classification
    std::vector<bool> valid;               // Per-lane validity mask

    SIMDMask(size_t vector_width = DefaultSIMDTraits::VECTOR_WIDTH) 
        : regimes(vector_width, SofteningRegime::NEWTONIAN),
          valid(vector_width, true) {}

    void reset(size_t size) {
        regimes.assign(size, SofteningRegime::NEWTONIAN);
        valid.assign(size, true);
    }
};

/**
 * @struct AlignedAllocator
 * @brief Custom allocator for SIMD-aligned memory allocation
 * Ensures allocated memory respects SIMD alignment requirements
 * Required for efficient vectorized memory access patterns
 */
template <typename T, size_t ALIGNMENT = DefaultSIMDTraits::ALIGNMENT>
struct AlignedAllocator {
    using value_type = T;
    using pointer = T*;
    using const_pointer = const T*;
    using reference = T&;
    using const_reference = const T&;
    using size_type = std::size_t;
    using difference_type = std::ptrdiff_t;

    template <typename U>
    struct rebind { typedef AlignedAllocator<U, ALIGNMENT> other; };

    pointer allocate(size_type n) {
        if (n == 0) return nullptr;
        return static_cast<pointer>(::operator new(n * sizeof(T), 
                                                   std::align_val_t(ALIGNMENT)));
    }

    void deallocate(pointer p, size_type) {
        if (p == nullptr) return;
        ::operator delete(p, std::align_val_t(ALIGNMENT));
    }

    size_type max_size() const {
        return std::numeric_limits<size_type>::max() / sizeof(T);
    }
};

#endif // SIMD_TRAITS_H
