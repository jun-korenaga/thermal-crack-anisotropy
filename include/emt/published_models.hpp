#ifndef EMT_PUBLISHED_MODELS_HPP
#define EMT_PUBLISHED_MODELS_HPP

#include "emt/tensor.hpp"

namespace emt {

// Guo et al. (2019), equation 45, with the factor-of-four correction
// described in the manuscript. The inclination is in radians.
Tensor4 guo2019_compliance_contribution(
    const Tensor4& background_stiffness,
    double aspect_ratio,
    double inclination);

// Random-orientation compliance of Xu et al. (2020), equations 36-37,
// with the Guo et al. factor-of-four correction and the later C66 typo fixed.
Tensor4 xu2020_random_compliance_contribution(
    const Tensor4& background_stiffness,
    double aspect_ratio);

}  // namespace emt

#endif
