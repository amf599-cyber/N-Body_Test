#ifndef GRAVTEST_SIMD_H
#define GRAVTEST_SIMD_H

/**
 * @file gravtest_simd.h
 * @brief Vectorized gravitational spline calculations with softening regime classification
 * 
 * Implements SIMD-friendly spline coefficient calculation with boolean masking for
 * different softening regimes:
 *   - r = 0: Self-interaction (skip)
 *   - 0 < r < h: Inner kernel regime (use full spline)
 *   - h ≤ r < 2h: Outer kernel regime (use modified spline)
 *   - r ≥ 2h: Newtonian regime (1/r approximation)
 * 
 * Uses where-like conditional operations to avoid branching within hot loops
 */

#include <vector>
#include <cmath>
#include "simd_traits.h"

/**
 * Branchless where-based conditional operations
 * Uses ternary operators instead of if/else to avoid branch misprediction
 * This mimics SIMD where semantics: compute all paths, mask-select results
 */
template<typename T>
inline T where(bool condition, T true_val, T false_val) {
    return condition ? true_val : false_val;
}

/**
 * @struct Vector3D_Batch
 * @brief Batched vector storage keeping x, y, z components together
 * 
 * Groups vector components into contiguous arrays for cleaner batch operations
 * and reduced code repetition in force calculations.
 */
struct Vector3D_Batch {
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> z;
    
    Vector3D_Batch() = default;
    
    explicit Vector3D_Batch(size_t width) : x(width), y(width), z(width) {}
    
    size_t size() const { return x.size(); }
    
    void resize(size_t width) {
        x.resize(width);
        y.resize(width);
        z.resize(width);
    }
    
    /**
     * Scale vector by coefficient array: v_out = v * coeffs
     */
    static Vector3D_Batch scale(const Vector3D_Batch& v, const std::vector<double>& coeffs) {
        Vector3D_Batch result(v.size());
        for (size_t i = 0; i < v.size(); ++i) {
            result.x[i] = v.x[i] * coeffs[i];
            result.y[i] = v.y[i] * coeffs[i];
            result.z[i] = v.z[i] * coeffs[i];
        }
        return result;
    }
    
    /**
     * Element-wise addition: result += other
     */
    static void add(Vector3D_Batch& result, const Vector3D_Batch& other, size_t width) {
        for (size_t i = 0; i < width; ++i) {
            result.x[i] += other.x[i];
            result.y[i] += other.y[i];
            result.z[i] += other.z[i];
        }
    }
};

/**
 * @class GravTestSIMD
 * @brief Vectorized gravitational spline calculations for SIMD processing
 * 
 * Processes multiple particle pairs simultaneously, classifying each pair into
 * appropriate softening regime and applying where-based conditional coefficients
 */
class GravTestSIMD {
public:
    /**
     * Default constructor
     */
    GravTestSIMD() = default;

    /**
     * Calculate spline coefficients for a batch of particle pairs
     * 
     * Uses where-based conditional logic to avoid branch misprediction penalties.
     * Each lane is classified into softening regime, then appropriate formula is applied.
     * 
     * @param r2_vec Vector of squared distances (one per SIMD lane)
     * @param twoh_vec Vector of combined softening lengths (one per SIMD lane)
     * @param mask SIMDMask classifying each lane into softening regime
     * @param a_out Output vector of potential coefficients
     * @param b_out Output vector of force/acceleration coefficients
     * @param vector_width Number of valid lanes to process
     */
    void calculateSplineForceVectorized(
        const std::vector<double>& r2_vec,
        const std::vector<double>& twoh_vec,
        SIMDMask& mask,
        std::vector<double>& a_out,
        std::vector<double>& b_out,
        size_t vector_width) const
    {
        // Ensure output vectors are properly sized
        a_out.resize(vector_width);
        b_out.resize(vector_width);

        // Process each lane independently with branchless where-based masking
        for (size_t lane = 0; lane < vector_width; ++lane) {
            double r2 = r2_vec[lane];
            double twoh = twoh_vec[lane];
            bool valid = mask.valid[lane];

            // Self-interaction regime (r^2 = 0)
            bool is_self = (r2 == 0.0);
            
            // Compute distance and normalized coordinate
            double r = std::sqrt(r2);
            double dih_calc = 2.0 / twoh;
            double u = r * dih_calc;
            
            // Regime classification with where-based masking
            bool is_inner_kernel = (r < twoh) && (u < 1.0);
            bool is_outer_kernel = (r < twoh) && (u >= 1.0);
            
            // Apply regime masks using branchless where semantics
            SIMDMask::SofteningRegime regime = 
                where(is_self, SIMDMask::SofteningRegime::SELF_INTERACTION,
                where(is_inner_kernel, SIMDMask::SofteningRegime::INNER_KERNEL,
                where(is_outer_kernel, SIMDMask::SofteningRegime::OUTER_KERNEL,
                      SIMDMask::SofteningRegime::NEWTONIAN)));
            
            // Compute all regime coefficients (masked selection avoids branching)
            // Inner kernel: 0 < r < h (u < 1.0)
            double u2_inner = u * u;
            double u3_inner = u2_inner * u;
            double u4_inner = u2_inner * u2_inner;
            double u5_inner = u4_inner * u;
            
            double a_inner = dih_calc * (7.0/5.0 
                - 2.0/3.0 * u2_inner
                + 3.0/10.0 * u4_inner
                - 1.0/10.0 * u5_inner);
            
            double b_inner = dih_calc * dih_calc * dih_calc * (4.0/3.0
                - 6.0/5.0 * u2_inner
                + 1.0/2.0 * u3_inner);
            
            // Outer kernel: h ≤ r < 2h (u >= 1.0)
            double dir = 1.0 / r;
            double u2_outer = u * u;
            double u3_outer = u2_outer * u;
            double u4_outer = u2_outer * u2_outer;
            double u5_outer = u4_outer * u;
            
            double a_outer = -1.0/15.0 * dir
                + dih_calc * (8.0/5.0
                    - 4.0/3.0 * u2_outer
                    + u3_outer
                    - 3.0/10.0 * u4_outer
                    + 1.0/30.0 * u5_outer);
            
            double b_outer = -1.0/15.0 * dir * dir * dir
                + dih_calc * dih_calc * dih_calc * (8.0/3.0
                    - 3.0 * u
                    + 6.0/5.0 * u2_outer
                    - 1.0/6.0 * u3_outer);
            
            // Newtonian: r >= 2h
            double a_newtonian = dir;
            double b_newtonian = dir * dir * dir;
            
            // Self-interaction yields zero
            double a_self = 0.0;
            double b_self = 0.0;
            
            // Where-based coefficient selection (no branches, pure masking)
            double a = where(is_self, a_self,
                    where(is_inner_kernel, a_inner,
                    where(is_outer_kernel, a_outer,
                          a_newtonian)));
            
            double b = where(is_self, b_self,
                    where(is_inner_kernel, b_inner,
                    where(is_outer_kernel, b_outer,
                          b_newtonian)));
            
            // Mask invalid lanes to zero output
            a_out[lane] = where(valid, a, 0.0);
            b_out[lane] = where(valid, b, 0.0);
            mask.regimes[lane] = regime;
        }
    }

    /**
     * Vectorized batch spline calculation for efficient processing
     * 
     * Processes multiple displacement vectors and applies where-based masks
     * for efficient data compression and regime-specific computation.
     * 
     * @param disp_x Vector of x-displacements per lane
     * @param disp_y Vector of y-displacements per lane
     * @param disp_z Vector of z-displacements per lane
     * @param softenings Vector of combined softening lengths per lane
     * @param masses Vector of source particle masses per lane
     * @param mask Output mask classifying regimes and validity
     * @param accel_x_out Output x-acceleration per lane
     * @param accel_y_out Output y-acceleration per lane
     * @param accel_z_out Output z-acceleration per lane
     * @param potential_out Output potential values per lane
     * @param vector_width Number of lanes to process
     */
    void calculateForcesVectorized(
        const std::vector<double>& disp_x,
        const std::vector<double>& disp_y,
        const std::vector<double>& disp_z,
        const std::vector<double>& softenings,
        const std::vector<double>& masses,
        SIMDMask& mask,
        std::vector<double>& accel_x_out,
        std::vector<double>& accel_y_out,
        std::vector<double>& accel_z_out,
        std::vector<double>& potential_out,
        size_t vector_width) const
    {
        accel_x_out.resize(vector_width);
        accel_y_out.resize(vector_width);
        accel_z_out.resize(vector_width);
        potential_out.resize(vector_width);

        // Calculate squared distances
        std::vector<double> r2(vector_width);
        for (size_t lane = 0; lane < vector_width; ++lane) {
            r2[lane] = disp_x[lane] * disp_x[lane] 
                     + disp_y[lane] * disp_y[lane]
                     + disp_z[lane] * disp_z[lane];
        }

        // Get spline coefficients with regime masking
        std::vector<double> a_coeffs, b_coeffs;
        calculateSplineForceVectorized(r2, softenings, mask, a_coeffs, b_coeffs, vector_width);

        // Assemble displacement batch
        Vector3D_Batch displacement;
        displacement.x = disp_x;
        displacement.y = disp_y;
        displacement.z = disp_z;
        
        // Compute force scale factors: force_coeff * mass per lane
        std::vector<double> force_scale(vector_width);
        for (size_t lane = 0; lane < vector_width; ++lane) {
            force_scale[lane] = b_coeffs[lane] * masses[lane];
        }
        
        // Scale displacement by force coefficients to get acceleration
        Vector3D_Batch acceleration = Vector3D_Batch::scale(displacement, force_scale);
        
        // Copy results back to output vectors
        accel_x_out = acceleration.x;
        accel_y_out = acceleration.y;
        accel_z_out = acceleration.z;
        
        // Compute potentials
        for (size_t lane = 0; lane < vector_width; ++lane) {
            potential_out[lane] = -masses[lane] * a_coeffs[lane];
        }
    }

    /**
     * Compress results based on validity mask
     * Note: With all interactions processed (including self-interactions),
     * this can be used to filter specific regime types if needed
     * 
     * Uses where-based masking instead of branching
     * 
     * @param data Input data vector
     * @param mask Validity mask
     * @param vector_width Width of data to filter
     * @return Compressed vector containing only valid entries
     */
    static std::vector<double> compressWithMask(
        const std::vector<double>& data,
        const SIMDMask& mask,
        size_t vector_width)
    {
        std::vector<double> compressed;
        
        for (size_t lane = 0; lane < vector_width; ++lane) {
            // Branchless where-based conditional: only add if valid
            if (mask.valid[lane]) {
                compressed.push_back(data[lane]);
            }
        }
        
        return compressed;
    }
};

#endif // GRAVTEST_SIMD_H
