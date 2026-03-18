#ifndef GRAVTEST_H
#define GRAVTEST_H

/**
 * @file gravtest.h
 * @brief Gravitational force calculations between particles with softened Changa Spline
 * Supports vectorized SIMD operations
 */

#include <cmath>
#include <experimental/simd>
#include "Vector3D.h"

namespace stdx = std::experimental;

class GravTest {
public:
    /**
     * Vectorized gravitational spline force calculation using native SIMD
     * Uses implicit SIMD masks to handle different distance regimes
     * @param r2_vec Vector of squared distances (r²)
     * @param twoh Softening parameter (same for all vectors)
     * @param a_vec Output vector for potential coefficients
     * @param b_vec Output vector for force/acceleration coefficients
     */
    template <typename SIMDType>
    void calculateSplineForceSIMD(const SIMDType& r2_vec, SIMDType& twoh,
                                   SIMDType& a_vec, SIMDType& b_vec) const {
        // Initialize output vectors
        a_vec = SIMDType(0.0);
        b_vec = SIMDType(0.0);

        // Compute square root and vectors for different regimes
        SIMDType r_vec = stdx::sqrt(r2_vec);
        SIMDType dih_calc = SIMDType(2.0 / twoh);
        
        // Create masks for different regimes
        // is_zero mask for self-interaction (r2 == 0), skip all calculations
        auto is_zero = (r2_vec == SIMDType(0.0));
        auto is_short_range = (r_vec < twoh) && !is_zero;
        auto is_long_range = !is_short_range && !is_zero;
        
        // For regime 3 (r >= twoh && !is_zero), compute and apply immediately
        SIMDType dir(0.0);
        stdx::where(is_long_range, dir) = SIMDType(1.0) / r_vec;
        SIMDType dir3 = dir * dir * dir;
        
        stdx::where(is_long_range, a_vec) = dir;
        stdx::where(is_long_range, b_vec) = dir3;
        
        // For short-range non-zero interactions, compute u = r * (2/twoh)
        SIMDType u_vec(0.0);
        stdx::where(is_short_range, u_vec) = r_vec * dih_calc;
        
        // Create regime masks (already exclude zero)
        auto is_u_small = (u_vec < 1.0);
        auto is_regime_2a = is_short_range && is_u_small;
        auto is_regime_2b = is_short_range && !is_u_small;
        
        // === Compute powers of u_vec for both regimes 2a and 2b ===
        SIMDType u2 = u_vec * u_vec;
        SIMDType u3 = u2 * u_vec;
        SIMDType u4 = u2 * u2;
        SIMDType u5 = u4 * u_vec;
        
        // === Compute Regime 2a coefficients (r < twoh && u < 1 && !is_zero) ===
        // Compiler natively handles coefficients in SIMD vectorized way
        SIMDType a_2a = dih_calc * (SIMDType(7.0/5.0) - SIMDType(2.0/3.0) * u2 + SIMDType(3.0/10.0) * u4 - SIMDType(1.0/10.0) * u5);
        SIMDType b_2a = dih_calc * dih_calc * dih_calc * (SIMDType(4.0/3.0) - SIMDType(6.0/5.0) * u2 + SIMDType(1.0/2.0) * u3);
        
        // === Compute Regime 2b coefficients (r < twoh && u >= 1 && !is_zero) ===
        // Reusing u2, u3, u4, u5 for coefficient calculations
        SIMDType a_2b = SIMDType(-1.0/15.0) * dir + dih_calc * (SIMDType(8.0/5.0) - SIMDType(4.0/3.0) * u2 + u3 - SIMDType(3.0/10.0) * u4 + SIMDType(1.0/30.0) * u5);
        SIMDType b_2b = SIMDType(-1.0/15.0) * dir3 + dih_calc * dih_calc * dih_calc * (SIMDType(8.0/3.0) - SIMDType(3.0) * u_vec + SIMDType(6.0/5.0) * u2 - SIMDType(1.0/6.0) * u3);
        
        // === Apply masks and set results using conditional select ===
        // Use where(...) = value pattern which returns a where_expression
        stdx::where(is_regime_2a, a_vec) = a_2a;
        stdx::where(is_regime_2a, b_vec) = b_2a;
        
        stdx::where(is_regime_2b, a_vec) = a_2b;
        stdx::where(is_regime_2b, b_vec) = b_2b;
        
        // is_zero automatically stays at 0.0 from initialization
    }
};
#endif // GRAVTEST_H