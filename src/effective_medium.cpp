#include "emt/effective_medium.hpp"

#include <stdexcept>

namespace emt {

Tensor4 dry_compliance_contribution(const Tensor4& background_stiffness,
                                    const Tensor4& hill_tensor) {
    Tensor4 identity;
    identity.setJ();
    const Tensor4 q = background_stiffness *
                      (identity - hill_tensor * background_stiffness);
    return q.inverse();
}

Tensor4 dilute_dry_compliance(const Tensor4& background_stiffness,
                              const Tensor4& compliance_contribution,
                              double porosity) {
    if (porosity < 0.0) {
        throw std::invalid_argument("porosity must be non-negative");
    }
    Tensor4 result = compliance_contribution;
    result *= porosity;
    result += background_stiffness.inverse();
    return result;
}

Tensor4 brown_korringa_saturated_compliance(
    const Tensor4& background_stiffness,
    const Tensor4& dry_compliance,
    double porosity,
    double fluid_bulk_modulus) {
    if (porosity < 0.0) {
        throw std::invalid_argument("porosity must be non-negative");
    }
    if (fluid_bulk_modulus <= 0.0) {
        throw std::invalid_argument("fluid bulk modulus must be positive");
    }

    const Tensor4 background_compliance = background_stiffness.inverse();
    double beta0 = 0.0;
    for (int a = 1; a <= 3; ++a) {
        for (int b = 1; b <= 3; ++b) {
            beta0 += background_compliance.val(a, a, b, b);
        }
    }
    const double beta_fluid = 1.0 / fluid_bulk_modulus;

    Tensor4 saturated;
    for (int I = 1; I <= 6; ++I) {
        int i, j;
        I2ij(I, i, j);
        for (int J = I; J <= 6; ++J) {
            int k, l;
            I2ij(J, k, l);

            double sijaa = 0.0;
            double sbbkl = 0.0;
            double sccdd = 0.0;
            double sijaa0 = 0.0;
            double sbbkl0 = 0.0;
            double sccdd0 = 0.0;
            for (int a = 1; a <= 3; ++a) {
                sijaa += dry_compliance.val(i, j, a, a);
                sbbkl += dry_compliance.val(a, a, k, l);
                sijaa0 += background_compliance.val(i, j, a, a);
                sbbkl0 += background_compliance.val(a, a, k, l);
                for (int b = 1; b <= 3; ++b) {
                    sccdd += dry_compliance.val(a, a, b, b);
                    sccdd0 += background_compliance.val(a, a, b, b);
                }
            }

            const double numerator =
                (sijaa - sijaa0) * (sbbkl - sbbkl0);
            const double denominator =
                (sccdd - sccdd0) + porosity * (beta_fluid - beta0);
            double correction = 0.0;
            if (numerator != 0.0 && denominator != 0.0) {
                correction = numerator / denominator;
            }
            const double value = dry_compliance.val(I, J) - correction;
            saturated.set(I, J, value);
            if (I != J) saturated.set(J, I, value);
        }
    }
    return saturated;
}

}  // namespace emt
