#ifndef VECTOR3D_SIMD_H
#define VECTOR3D_SIMD_H

/**
 * @file Vector3D_SIMD.h
 * @brief SIMD-optimized Vector3D operations for SoA particle data
 * 
 * Provides efficient packed vector operations on Structure of Arrays layout
 * with explicit SIMD intrinsics (AVX2/SSE2) and scalar fallback.
 * Architecture-aware compilation: automatically selects best available SIMD level.
 */

#include <cmath>
#include "simd_traits.h"

// Include SIMD intrinsics based on detected architecture
#if SIMD_AVX2
    #include <immintrin.h>
#elif SIMD_SSE2
    #include <emmintrin.h>
#endif

/**
 * @namespace Vec3SIMD
 * @brief Vectorized 3D vector operations for SoA-layout particle data
 * All operations work on packed arrays of coordinates
 * 
 * Uses architecture-specific intrinsics when available:
 * - AVX2: 4-wide double precision operations
 * - SSE2: 2-wide double precision operations
 * - Scalar: Single-element fallback
 */
namespace Vec3SIMD {

/**
 * Compute squared distance between packed position data
 * Optimized with SIMD intrinsics for bulk distance calculations
 */
inline void distanceSquared(
    const double* dx, const double* dy, const double* dz,
    double* r2, size_t width)
{
#if SIMD_AVX2
    // AVX2: Process 4 doubles per iteration
    constexpr size_t VWIDTH = 4;
    size_t i = 0;
    
    for (; i + VWIDTH <= width; i += VWIDTH) {
        __m256d vdx = _mm256_loadu_pd(&dx[i]);
        __m256d vdy = _mm256_loadu_pd(&dy[i]);
        __m256d vdz = _mm256_loadu_pd(&dz[i]);
        
        __m256d vr2 = _mm256_add_pd(
            _mm256_mul_pd(vdx, vdx),
            _mm256_add_pd(
                _mm256_mul_pd(vdy, vdy),
                _mm256_mul_pd(vdz, vdz)
            )
        );
        
        _mm256_storeu_pd(&r2[i], vr2);
    }
    
    // Handle remainder with scalar loop
    for (; i < width; ++i) {
        r2[i] = dx[i]*dx[i] + dy[i]*dy[i] + dz[i]*dz[i];
    }
#elif SIMD_SSE2
    // SSE2: Process 2 doubles per iteration
    constexpr size_t VWIDTH = 2;
    size_t i = 0;
    
    for (; i + VWIDTH <= width; i += VWIDTH) {
        __m128d vdx = _mm_loadu_pd(&dx[i]);
        __m128d vdy = _mm_loadu_pd(&dy[i]);
        __m128d vdz = _mm_loadu_pd(&dz[i]);
        
        __m128d vr2 = _mm_add_pd(
            _mm_mul_pd(vdx, vdx),
            _mm_add_pd(
                _mm_mul_pd(vdy, vdy),
                _mm_mul_pd(vdz, vdz)
            )
        );
        
        _mm_storeu_pd(&r2[i], vr2);
    }
    
    // Handle remainder
    for (; i < width; ++i) {
        r2[i] = dx[i]*dx[i] + dy[i]*dy[i] + dz[i]*dz[i];
    }
#else
    // Scalar fallback
    for (size_t i = 0; i < width; ++i) {
        r2[i] = dx[i]*dx[i] + dy[i]*dy[i] + dz[i]*dz[i];
    }
#endif
}

/**
 * Compute scalar product between two packed vectors
 * Optimized dot product using SIMD intrinsics
 */
inline void dotProduct(
    const double* a_x, const double* a_y, const double* a_z,
    const double* b_x, const double* b_y, const double* b_z,
    double* result, size_t width)
{
#if SIMD_AVX2
    // AVX2: Process 4 doubles per iteration
    constexpr size_t VWIDTH = 4;
    size_t i = 0;
    
    for (; i + VWIDTH <= width; i += VWIDTH) {
        __m256d vax = _mm256_loadu_pd(&a_x[i]);
        __m256d vay = _mm256_loadu_pd(&a_y[i]);
        __m256d vaz = _mm256_loadu_pd(&a_z[i]);
        
        __m256d vbx = _mm256_loadu_pd(&b_x[i]);
        __m256d vby = _mm256_loadu_pd(&b_y[i]);
        __m256d vbz = _mm256_loadu_pd(&b_z[i]);
        
        __m256d vres = _mm256_add_pd(
            _mm256_mul_pd(vax, vbx),
            _mm256_add_pd(
                _mm256_mul_pd(vay, vby),
                _mm256_mul_pd(vaz, vbz)
            )
        );
        
        _mm256_storeu_pd(&result[i], vres);
    }
    
    // Handle remainder
    for (; i < width; ++i) {
        result[i] = a_x[i]*b_x[i] + a_y[i]*b_y[i] + a_z[i]*b_z[i];
    }
#elif SIMD_SSE2
    // SSE2: Process 2 doubles per iteration
    constexpr size_t VWIDTH = 2;
    size_t i = 0;
    
    for (; i + VWIDTH <= width; i += VWIDTH) {
        __m128d vax = _mm_loadu_pd(&a_x[i]);
        __m128d vay = _mm_loadu_pd(&a_y[i]);
        __m128d vaz = _mm_loadu_pd(&a_z[i]);
        
        __m128d vbx = _mm_loadu_pd(&b_x[i]);
        __m128d vby = _mm_loadu_pd(&b_y[i]);
        __m128d vbz = _mm_loadu_pd(&b_z[i]);
        
        __m128d vres = _mm_add_pd(
            _mm_mul_pd(vax, vbx),
            _mm_add_pd(
                _mm_mul_pd(vay, vby),
                _mm_mul_pd(vaz, vbz)
            )
        );
        
        _mm_storeu_pd(&result[i], vres);
    }
    
    // Handle remainder
    for (; i < width; ++i) {
        result[i] = a_x[i]*b_x[i] + a_y[i]*b_y[i] + a_z[i]*b_z[i];
    }
#else
    // Scalar fallback
    for (size_t i = 0; i < width; ++i) {
        result[i] = a_x[i]*b_x[i] + a_y[i]*b_y[i] + a_z[i]*b_z[i];
    }
#endif
}

/**
 * Scale packed vectors by scalar coefficients
 * Optimized element-wise multiplication using SIMD
 */
inline void scale(
    const double* x, const double* y, const double* z,
    const double* scale,
    double* out_x, double* out_y, double* out_z,
    size_t width)
{
#if SIMD_AVX2
    // AVX2: Process 4 doubles per iteration
    constexpr size_t VWIDTH = 4;
    size_t i = 0;
    
    for (; i + VWIDTH <= width; i += VWIDTH) {
        __m256d vx = _mm256_loadu_pd(&x[i]);
        __m256d vy = _mm256_loadu_pd(&y[i]);
        __m256d vz = _mm256_loadu_pd(&z[i]);
        __m256d vscale = _mm256_loadu_pd(&scale[i]);
        
        _mm256_storeu_pd(&out_x[i], _mm256_mul_pd(vx, vscale));
        _mm256_storeu_pd(&out_y[i], _mm256_mul_pd(vy, vscale));
        _mm256_storeu_pd(&out_z[i], _mm256_mul_pd(vz, vscale));
    }
    
    // Handle remainder
    for (; i < width; ++i) {
        out_x[i] = x[i] * scale[i];
        out_y[i] = y[i] * scale[i];
        out_z[i] = z[i] * scale[i];
    }
#elif SIMD_SSE2
    // SSE2: Process 2 doubles per iteration
    constexpr size_t VWIDTH = 2;
    size_t i = 0;
    
    for (; i + VWIDTH <= width; i += VWIDTH) {
        __m128d vx = _mm_loadu_pd(&x[i]);
        __m128d vy = _mm_loadu_pd(&y[i]);
        __m128d vz = _mm_loadu_pd(&z[i]);
        __m128d vscale = _mm_loadu_pd(&scale[i]);
        
        _mm_storeu_pd(&out_x[i], _mm_mul_pd(vx, vscale));
        _mm_storeu_pd(&out_y[i], _mm_mul_pd(vy, vscale));
        _mm_storeu_pd(&out_z[i], _mm_mul_pd(vz, vscale));
    }
    
    // Handle remainder
    for (; i < width; ++i) {
        out_x[i] = x[i] * scale[i];
        out_y[i] = y[i] * scale[i];
        out_z[i] = z[i] * scale[i];
    }
#else
    // Scalar fallback
    for (size_t i = 0; i < width; ++i) {
        out_x[i] = x[i] * scale[i];
        out_y[i] = y[i] * scale[i];
        out_z[i] = z[i] * scale[i];
    }
#endif
}

/**
 * Add packed vectors element-wise
 * Optimized vector addition using SIMD
 */
inline void add(
    const double* a_x, const double* a_y, const double* a_z,
    const double* b_x, const double* b_y, const double* b_z,
    double* out_x, double* out_y, double* out_z,
    size_t width)
{
#if SIMD_AVX2
    // AVX2: Process 4 doubles per iteration
    constexpr size_t VWIDTH = 4;
    size_t i = 0;
    
    for (; i + VWIDTH <= width; i += VWIDTH) {
        __m256d vax = _mm256_loadu_pd(&a_x[i]);
        __m256d vay = _mm256_loadu_pd(&a_y[i]);
        __m256d vaz = _mm256_loadu_pd(&a_z[i]);
        
        __m256d vbx = _mm256_loadu_pd(&b_x[i]);
        __m256d vby = _mm256_loadu_pd(&b_y[i]);
        __m256d vbz = _mm256_loadu_pd(&b_z[i]);
        
        _mm256_storeu_pd(&out_x[i], _mm256_add_pd(vax, vbx));
        _mm256_storeu_pd(&out_y[i], _mm256_add_pd(vay, vby));
        _mm256_storeu_pd(&out_z[i], _mm256_add_pd(vaz, vbz));
    }
    
    // Handle remainder
    for (; i < width; ++i) {
        out_x[i] = a_x[i] + b_x[i];
        out_y[i] = a_y[i] + b_y[i];
        out_z[i] = a_z[i] + b_z[i];
    }
#elif SIMD_SSE2
    // SSE2: Process 2 doubles per iteration
    constexpr size_t VWIDTH = 2;
    size_t i = 0;
    
    for (; i + VWIDTH <= width; i += VWIDTH) {
        __m128d vax = _mm_loadu_pd(&a_x[i]);
        __m128d vay = _mm_loadu_pd(&a_y[i]);
        __m128d vaz = _mm_loadu_pd(&a_z[i]);
        
        __m128d vbx = _mm_loadu_pd(&b_x[i]);
        __m128d vby = _mm_loadu_pd(&b_y[i]);
        __m128d vbz = _mm_loadu_pd(&b_z[i]);
        
        _mm_storeu_pd(&out_x[i], _mm_add_pd(vax, vbx));
        _mm_storeu_pd(&out_y[i], _mm_add_pd(vay, vby));
        _mm_storeu_pd(&out_z[i], _mm_add_pd(vaz, vbz));
    }
    
    // Handle remainder
    for (; i < width; ++i) {
        out_x[i] = a_x[i] + b_x[i];
        out_y[i] = a_y[i] + b_y[i];
        out_z[i] = a_z[i] + b_z[i];
    }
#else
    // Scalar fallback
    for (size_t i = 0; i < width; ++i) {
        out_x[i] = a_x[i] + b_x[i];
        out_y[i] = a_y[i] + b_y[i];
        out_z[i] = a_z[i] + b_z[i];
    }
#endif
}

/**
 * Accumulate packed vectors (in-place addition)
 * Optimized for tight accumulation loops
 */
inline void accumulate(
    double* out_x, double* out_y, double* out_z,
    const double* a_x, const double* a_y, const double* a_z,
    size_t width)
{
#if SIMD_AVX2
    // AVX2: Process 4 doubles per iteration
    constexpr size_t VWIDTH = 4;
    size_t i = 0;
    
    for (; i + VWIDTH <= width; i += VWIDTH) {
        __m256d voutx = _mm256_loadu_pd(&out_x[i]);
        __m256d vouty = _mm256_loadu_pd(&out_y[i]);
        __m256d voutz = _mm256_loadu_pd(&out_z[i]);
        
        __m256d vax = _mm256_loadu_pd(&a_x[i]);
        __m256d vay = _mm256_loadu_pd(&a_y[i]);
        __m256d vaz = _mm256_loadu_pd(&a_z[i]);
        
        _mm256_storeu_pd(&out_x[i], _mm256_add_pd(voutx, vax));
        _mm256_storeu_pd(&out_y[i], _mm256_add_pd(vouty, vay));
        _mm256_storeu_pd(&out_z[i], _mm256_add_pd(voutz, vaz));
    }
    
    // Handle remainder
    for (; i < width; ++i) {
        out_x[i] += a_x[i];
        out_y[i] += a_y[i];
        out_z[i] += a_z[i];
    }
#elif SIMD_SSE2
    // SSE2: Process 2 doubles per iteration
    constexpr size_t VWIDTH = 2;
    size_t i = 0;
    
    for (; i + VWIDTH <= width; i += VWIDTH) {
        __m128d voutx = _mm_loadu_pd(&out_x[i]);
        __m128d vouty = _mm_loadu_pd(&out_y[i]);
        __m128d voutz = _mm_loadu_pd(&out_z[i]);
        
        __m128d vax = _mm_loadu_pd(&a_x[i]);
        __m128d vay = _mm_loadu_pd(&a_y[i]);
        __m128d vaz = _mm_loadu_pd(&a_z[i]);
        
        _mm_storeu_pd(&out_x[i], _mm_add_pd(voutx, vax));
        _mm_storeu_pd(&out_y[i], _mm_add_pd(vouty, vay));
        _mm_storeu_pd(&out_z[i], _mm_add_pd(voutz, vaz));
    }
    
    // Handle remainder
    for (; i < width; ++i) {
        out_x[i] += a_x[i];
        out_y[i] += a_y[i];
        out_z[i] += a_z[i];
    }
#else
    // Scalar fallback
    for (size_t i = 0; i < width; ++i) {
        out_x[i] += a_x[i];
        out_y[i] += a_y[i];
        out_z[i] += a_z[i];
    }
#endif
}

/**
 * Compute displacements from target to sources (SoA-optimized)
 * Uses SIMD for bulk displacement computation
 */
inline void displacement(
    double target_x, double target_y, double target_z,
    const double* src_x, const double* src_y, const double* src_z,
    double* disp_x, double* disp_y, double* disp_z,
    size_t width)
{
#if SIMD_AVX2
    // AVX2: Process 4 doubles per iteration
    constexpr size_t VWIDTH = 4;
    size_t i = 0;
    
    __m256d vtarget_x = _mm256_set1_pd(target_x);
    __m256d vtarget_y = _mm256_set1_pd(target_y);
    __m256d vtarget_z = _mm256_set1_pd(target_z);
    
    for (; i + VWIDTH <= width; i += VWIDTH) {
        __m256d vsrc_x = _mm256_loadu_pd(&src_x[i]);
        __m256d vsrc_y = _mm256_loadu_pd(&src_y[i]);
        __m256d vsrc_z = _mm256_loadu_pd(&src_z[i]);
        
        _mm256_storeu_pd(&disp_x[i], _mm256_sub_pd(vsrc_x, vtarget_x));
        _mm256_storeu_pd(&disp_y[i], _mm256_sub_pd(vsrc_y, vtarget_y));
        _mm256_storeu_pd(&disp_z[i], _mm256_sub_pd(vsrc_z, vtarget_z));
    }
    
    // Handle remainder
    for (; i < width; ++i) {
        disp_x[i] = src_x[i] - target_x;
        disp_y[i] = src_y[i] - target_y;
        disp_z[i] = src_z[i] - target_z;
    }
#elif SIMD_SSE2
    // SSE2: Process 2 doubles per iteration
    constexpr size_t VWIDTH = 2;
    size_t i = 0;
    
    __m128d vtarget_x = _mm_set1_pd(target_x);
    __m128d vtarget_y = _mm_set1_pd(target_y);
    __m128d vtarget_z = _mm_set1_pd(target_z);
    
    for (; i + VWIDTH <= width; i += VWIDTH) {
        __m128d vsrc_x = _mm_loadu_pd(&src_x[i]);
        __m128d vsrc_y = _mm_loadu_pd(&src_y[i]);
        __m128d vsrc_z = _mm_loadu_pd(&src_z[i]);
        
        _mm_storeu_pd(&disp_x[i], _mm_sub_pd(vsrc_x, vtarget_x));
        _mm_storeu_pd(&disp_y[i], _mm_sub_pd(vsrc_y, vtarget_y));
        _mm_storeu_pd(&disp_z[i], _mm_sub_pd(vsrc_z, vtarget_z));
    }
    
    // Handle remainder
    for (; i < width; ++i) {
        disp_x[i] = src_x[i] - target_x;
        disp_y[i] = src_y[i] - target_y;
        disp_z[i] = src_z[i] - target_z;
    }
#else
    // Scalar fallback
    for (size_t i = 0; i < width; ++i) {
        disp_x[i] = src_x[i] - target_x;
        disp_y[i] = src_y[i] - target_y;
        disp_z[i] = src_z[i] - target_z;
    }
#endif
}

/**
 * Multiply two packed arrays element-wise
 * Optimized multiplication using SIMD
 */
inline void multiply(
    const double* a, const double* b,
    double* out, size_t width)
{
#if SIMD_AVX2
    // AVX2: Process 4 doubles per iteration
    constexpr size_t VWIDTH = 4;
    size_t i = 0;
    
    for (; i + VWIDTH <= width; i += VWIDTH) {
        __m256d va = _mm256_loadu_pd(&a[i]);
        __m256d vb = _mm256_loadu_pd(&b[i]);
        _mm256_storeu_pd(&out[i], _mm256_mul_pd(va, vb));
    }
    
    // Handle remainder
    for (; i < width; ++i) {
        out[i] = a[i] * b[i];
    }
#elif SIMD_SSE2
    // SSE2: Process 2 doubles per iteration
    constexpr size_t VWIDTH = 2;
    size_t i = 0;
    
    for (; i + VWIDTH <= width; i += VWIDTH) {
        __m128d va = _mm_loadu_pd(&a[i]);
        __m128d vb = _mm_loadu_pd(&b[i]);
        _mm_storeu_pd(&out[i], _mm_mul_pd(va, vb));
    }
    
    // Handle remainder
    for (; i < width; ++i) {
        out[i] = a[i] * b[i];
    }
#else
    // Scalar fallback
    for (size_t i = 0; i < width; ++i) {
        out[i] = a[i] * b[i];
    }
#endif
}

/**
 * Fused multiply-add: out = a * b + c
 * Optimized using native FMA instruction when available
 */
inline void fma(
    const double* a, const double* b, const double* c,
    double* out, size_t width)
{
#if SIMD_AVX2
    // AVX2: Process 4 doubles per iteration with FMA
    constexpr size_t VWIDTH = 4;
    size_t i = 0;
    
    for (; i + VWIDTH <= width; i += VWIDTH) {
        __m256d va = _mm256_loadu_pd(&a[i]);
        __m256d vb = _mm256_loadu_pd(&b[i]);
        __m256d vc = _mm256_loadu_pd(&c[i]);
        
        // Check if FMA is available (requires AVX2 + FMA instruction set)
        #ifdef __FMA__
            _mm256_storeu_pd(&out[i], _mm256_fmadd_pd(va, vb, vc));
        #else
            // Fallback: separate multiply and add
            _mm256_storeu_pd(&out[i], 
                _mm256_add_pd(_mm256_mul_pd(va, vb), vc));
        #endif
    }
    
    // Handle remainder
    for (; i < width; ++i) {
        out[i] = a[i] * b[i] + c[i];
    }
#elif SIMD_SSE2
    // SSE2: Process 2 doubles per iteration
    constexpr size_t VWIDTH = 2;
    size_t i = 0;
    
    for (; i + VWIDTH <= width; i += VWIDTH) {
        __m128d va = _mm_loadu_pd(&a[i]);
        __m128d vb = _mm_loadu_pd(&b[i]);
        __m128d vc = _mm_loadu_pd(&c[i]);
        
        _mm_storeu_pd(&out[i], 
            _mm_add_pd(_mm_mul_pd(va, vb), vc));
    }
    
    // Handle remainder
    for (; i < width; ++i) {
        out[i] = a[i] * b[i] + c[i];
    }
#else
    // Scalar fallback
    for (size_t i = 0; i < width; ++i) {
        out[i] = a[i] * b[i] + c[i];
    }
#endif
}

/**
 * Negate packed array
 * Scalar operation (not performance-critical)
 */
inline void negate(
    const double* in, double* out, size_t width)
{
    for (size_t i = 0; i < width; ++i) {
        out[i] = -in[i];
    }
}

/**
 * Broadcast scalar to packed vector
 * Scalar operation (not performance-critical)
 */
inline void broadcast(
    double scalar,
    double* out_x, double* out_y, double* out_z,
    size_t width)
{
    for (size_t i = 0; i < width; ++i) {
        out_x[i] = scalar;
        out_y[i] = scalar;
        out_z[i] = scalar;
    }
}

/**
 * Zero-initialize packed vectors
 * Scalar operation (not performance-critical)
 */
inline void zero(
    double* out_x, double* out_y, double* out_z,
    size_t width)
{
    broadcast(0.0, out_x, out_y, out_z, width);
}

} // namespace Vec3SIMD

#endif // VECTOR3D_SIMD_H
