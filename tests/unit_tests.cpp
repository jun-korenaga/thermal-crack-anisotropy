#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

#include "emt/effective_medium.hpp"
#include "emt/ellipsoid.hpp"
#include "emt/hill_tensor.hpp"
#include "emt/published_models.hpp"
#include "emt/tensor.hpp"

namespace {

void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

double relative_rms(const Tensor4& calculated, const Tensor4& reference) {
    double numerator = 0.0;
    double denominator = 0.0;
    for (int I = 1; I <= 6; ++I) {
        for (int J = 1; J <= 6; ++J) {
            const double difference = calculated.val(I, J) - reference.val(I, J);
            numerator += difference * difference;
            denominator += reference.val(I, J) * reference.val(I, J);
        }
    }
    return std::sqrt(numerator / denominator);
}

Tensor4 calculate_hill(const Tensor4& stiffness, const Ellipsoid& ellipsoid) {
    HillTensor calculation(stiffness, ellipsoid, 40, 80, 0.4, false, 20);
    Tensor4 hill;
    calculation.setTensor(hill);
    return hill;
}

}  // namespace

int main() {
    try {
        Array1d<double> constants(5);
        constants(1) = 47.0;
        constants(2) = 8.0;
        constants(3) = 5.0;
        constants(4) = 34.0;
        constants(5) = 17.0;
        const Tensor4 stiffness(constants);

        const Tensor4 hill = calculate_hill(stiffness, Ellipsoid(1.0, 1.0, 0.1));
        const Tensor4 calculated_eshelby = hill * stiffness;
        const Tensor4 exact_eshelby =
            eshelby_transverse_exact(stiffness, 1.0, 0.1);
        require(relative_rms(calculated_eshelby, exact_eshelby) < 2.0e-3,
                "adaptive integration did not reproduce the Withers benchmark");

        const Tensor4 scaled_hill =
            calculate_hill(stiffness, Ellipsoid(2.0, 2.0, 0.2));
        require(relative_rms(hill, scaled_hill) < 1.0e-12,
                "Hill tensor changed under uniform ellipsoid scaling");

        for (int I = 1; I <= 6; ++I) {
            for (int J = 1; J <= 6; ++J) {
                require(std::abs(hill.val(I, J) - hill.val(J, I)) < 2.0e-5,
                        "Hill tensor failed major-symmetry check");
            }
        }

        const Tensor4 background_compliance = stiffness.inverse();
        const Tensor4 saturated = emt::brown_korringa_saturated_compliance(
            stiffness, background_compliance, 0.0, 2.2);
        require(relative_rms(saturated, background_compliance) < 1.0e-14,
                "zero-porosity fluid substitution changed the background");

        const Tensor4 xu2020 =
            emt::xu2020_random_compliance_contribution(stiffness, 0.01);
        require(std::abs(xu2020.val(6, 6) -
                         2.0 * xu2020.val(1, 1)) < 1.0e-14,
                "Xu et al. comparison tensor was assembled incorrectly");

        std::cout << "All EMT unit tests passed.\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Test failure: " << error.what() << '\n';
        return 1;
    }
}
