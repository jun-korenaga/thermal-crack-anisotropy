#ifndef EMT_EFFECTIVE_MEDIUM_HPP
#define EMT_EFFECTIVE_MEDIUM_HPP

#include "emt/tensor.hpp"

namespace emt {

// Compliance contribution H = [C0 : (J - P : C0)]^{-1} for a dry pore.
Tensor4 dry_compliance_contribution(const Tensor4& background_stiffness,
                                    const Tensor4& hill_tensor);

// Dilute dry compliance S_dry = S0 + porosity * H.
Tensor4 dilute_dry_compliance(const Tensor4& background_stiffness,
                              const Tensor4& compliance_contribution,
                              double porosity);

// Brown-Korringa fluid substitution applied to a dry compliance tensor.
Tensor4 brown_korringa_saturated_compliance(
    const Tensor4& background_stiffness,
    const Tensor4& dry_compliance,
    double porosity,
    double fluid_bulk_modulus);

}  // namespace emt

#endif
