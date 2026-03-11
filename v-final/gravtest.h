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
    void calculateSplineForceSIMD(const SIMDType& r2_vec, double twoh,
                                   SIMDType& a_vec, SIMDType& b_vec) const {
        // Initialize output vectors
        a_vec = SIMDType(0.0);
        b_vec = SIMDType(0.0);

        // Compute square root and vectors for different regimes
        SIMDType r_vec = stdx::sqrt(r2_vec);
        SIMDType dih_calc = SIMDType(2.0 / twoh);
        
        // Create masks for different regimes
        // is_zero mask for self-interaction (r2 == 0), handled by init acc/pot to 0.0
        auto is_short_range = (r_vec < twoh);
        
        // For short-range interactions, compute u = r * (2/twoh)
        SIMDType u_vec = r_vec * dih_calc;
        
        // Create regime masks
        auto is_u_small = (u_vec < 1.0);
        auto is_regime_2a = is_short_range && is_u_small;
        auto is_regime_2b = is_short_range && !is_u_small;
        auto is_long_range = !is_short_range;
        
        // === Compute Regime 2a coefficients (r < twoh && u < 1) ===
        SIMDType u2 = u_vec * u_vec;
        SIMDType u3 = u2 * u_vec;
        SIMDType u4 = u2 * u2;
        SIMDType u5 = u4 * u_vec;
        
        SIMDType a_2a = dih_calc * (7.0/5.0 - 2.0/3.0 * u2 + 3.0/10.0 * u4 - 1.0/10.0 * u5);
        SIMDType b_2a = dih_calc * dih_calc * dih_calc * (4.0/3.0 - 6.0/5.0 * u2 + 1.0/2.0 * u3);
        
        // === Compute Regime 2b coefficients (r < twoh && u >= 1) ===
        SIMDType dir = SIMDType(1.0) / r_vec;
        SIMDType dir3 = dir * dir * dir;
        
        SIMDType u2_b = u_vec * u_vec;
        SIMDType u3_b = u2_b * u_vec;
        SIMDType u4_b = u2_b * u2_b;
        SIMDType u5_b = u4_b * u_vec;
        
        SIMDType a_2b = -1.0/15.0 * dir + dih_calc * (8.0/5.0 - 4.0/3.0 * u2_b + u3_b - 3.0/10.0 * u4_b + 1.0/30.0 * u5_b);
        SIMDType b_2b = -1.0/15.0 * dir3 + dih_calc * dih_calc * dih_calc * (8.0/3.0 - 3.0 * u_vec + 6.0/5.0 * u2_b - 1.0/6.0 * u3_b);
        
        // === Compute Regime 3 coefficients (r >= twoh) ===
        SIMDType a_3 = dir;
        SIMDType b_3 = dir3;
        
        // === Apply masks and set results using conditional select ===
        // Use where(...) = value pattern which returns a where_expression
        stdx::where(is_regime_2a, a_vec) = a_2a;
        stdx::where(is_regime_2a, b_vec) = b_2a;
        
        stdx::where(is_regime_2b, a_vec) = a_2b;
        stdx::where(is_regime_2b, b_vec) = b_2b;
        
        stdx::where(is_long_range, a_vec) = a_3;
        stdx::where(is_long_range, b_vec) = b_3;
        
        // is_zero is already handled by initialization to 0.0
    }
};
#endif // GRAVTEST_H