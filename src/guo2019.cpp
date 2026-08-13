#include "emt/published_models.hpp"

#include <cmath>
#include <stdexcept>

namespace {

struct FractureCompliances {
    double normal;
    double shear;
};

FractureCompliances corrected_compliances(const Tensor4& stiffness,
                                           double aspect_ratio) {
    if (aspect_ratio <= 0.0) {
        throw std::invalid_argument("aspect ratio must be positive");
    }

    const double c11 = stiffness.val(1, 1);
    const double c13 = stiffness.val(1, 3);
    const double c33 = stiffness.val(3, 3);
    const double c44 = stiffness.val(4, 4);
    const double c66 = stiffness.val(6, 6);
    const double c1 = std::sqrt(c11 * c33);
    const double b3 = std::sqrt(c66 / c44);
    const double b4 = std::sqrt(
        (c1 - c13) * (c1 + c13 + 2.0 * c44) / (c33 * c44));
    const double scale = 3.0 / (4.0 * M_PI * aspect_ratio);

    FractureCompliances result;
    result.normal = (8.0 * b4 / (3.0 * c1)) *
                    (scale / (1.0 - c13 * c13 / (c1 * c1)));

    // The leading 4 corrects the factor-of-four error in the published Z_T.
    result.shear = (4.0 / (3.0 * c44)) *
                   (scale /
                    (b3 + b4 -
                     2.0 * c44 * b4 / (c1 + c13 + 2.0 * c44)));
    return result;
}

}  // namespace

namespace emt {

Tensor4 guo2019_compliance_contribution(const Tensor4& stiffness,
                                        double aspect_ratio,
                                        double inclination) {
    const auto compliance = corrected_compliances(stiffness, aspect_ratio);
    const double sine = std::sin(inclination);
    const double cosine = std::cos(inclination);
    const double sine_squared = sine * sine;
    const double cosine_squared = cosine * cosine;
    const double sine_cosine = sine * cosine;

    Tensor4 result;
    result = 0.0;
    double value = compliance.shear * sine_squared;
    result.set(2, 2, value);
    result.set(6, 6, value);

    result.set(3, 3, compliance.normal * cosine_squared);
    result.set(4, 4, compliance.normal * sine_squared +
                         compliance.shear * cosine_squared);
    result.set(5, 5, compliance.shear * cosine_squared);

    value = -compliance.shear * sine_cosine;
    result.set(2, 4, value);
    result.set(4, 2, value);
    result.set(5, 6, value);
    result.set(6, 5, value);

    value = -compliance.normal * sine_cosine;
    result.set(3, 4, value);
    result.set(4, 3, value);
    return result;
}

Tensor4 xu2020_random_compliance_contribution(const Tensor4& stiffness,
                                               double aspect_ratio) {
    const auto compliance = corrected_compliances(stiffness, aspect_ratio);
    Tensor4 result;
    result = 0.0;

    result.set(1, 1, 0.25 * compliance.shear);
    result.set(2, 2, 0.25 * compliance.shear);
    result.set(3, 3, 0.50 * compliance.normal);
    result.set(4, 4, 0.25 * compliance.normal + 0.50 * compliance.shear);
    result.set(5, 5, 0.25 * compliance.normal + 0.50 * compliance.shear);
    result.set(6, 6, 0.50 * compliance.shear);
    return result;
}

}  // namespace emt
