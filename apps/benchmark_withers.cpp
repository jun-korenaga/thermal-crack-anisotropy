#include <algorithm>
#include <cmath>
#include <exception>
#include <iostream>

#include "cli.hpp"
#include "emt/ellipsoid.hpp"
#include "emt/hill_tensor.hpp"
#include "emt/tensor.hpp"

namespace {

double relative_rms(const Tensor4& calculated, const Tensor4& reference) {
    double squared_difference = 0.0;
    double squared_reference = 0.0;
    for (int I = 1; I <= 6; ++I) {
        for (int J = 1; J <= 6; ++J) {
            const double difference = calculated.val(I, J) - reference.val(I, J);
            squared_difference += difference * difference;
            squared_reference += reference.val(I, J) * reference.val(I, J);
        }
    }
    return std::sqrt(squared_difference / squared_reference);
}

double maximum_major_asymmetry(const Tensor4& tensor) {
    double maximum = 0.0;
    for (int I = 1; I <= 6; ++I) {
        for (int J = I + 1; J <= 6; ++J) {
            maximum = std::max(
                maximum, std::abs(tensor.val(I, J) - tensor.val(J, I)));
        }
    }
    return maximum;
}

void usage() {
    std::cout
        << "Usage: emt_withers_benchmark [options]\n"
        << "  --aspect VALUE       a3/a1 (default: 0.01)\n"
        << "  --ntheta N           initial polar divisions (default: 100)\n"
        << "  --weight VALUE       adaptive weight (default: 0.4)\n"
        << "  --max-level N        recursion cap (default: 20)\n"
        << "  --stiffness CSV      C11,C12,C13,C33,C44 in GPa\n"
        << "                       (default: 47,8,5,34,17)\n";
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (cli::has_flag(argc, argv, "--help")) {
            usage();
            return 0;
        }
        const double aspect = cli::option_double(argc, argv, "--aspect", 0.01);
        const int ntheta = cli::option_int(argc, argv, "--ntheta", 100);
        const double weight = cli::option_double(argc, argv, "--weight", 0.4);
        const int max_level = cli::option_int(argc, argv, "--max-level", 20);
        cli::require_positive(aspect, "aspect");
        cli::require_positive(ntheta, "ntheta");

        auto constants = cli::ti_parameters(cli::option(
            argc, argv, "--stiffness", "47,8,5,34,17"));
        const Tensor4 stiffness(constants);
        const Ellipsoid ellipsoid(1.0, 1.0, aspect);

        HillTensor calculation(stiffness, ellipsoid, ntheta, 2 * ntheta,
                               weight, false, max_level);
        Tensor4 hill;
        calculation.setTensor(hill);
        const Tensor4 eshelby = hill * stiffness;
        const Tensor4 exact = eshelby_transverse_exact(stiffness, 1.0, aspect);

        std::cout << "aspect,ntheta,weight,integration_points,relative_rms_error,"
                     "max_absolute_major_asymmetry,max_level_hits\n";
        std::cout << aspect << ',' << ntheta << ',' << weight << ','
                  << calculation.numIntegPoints() << ','
                  << relative_rms(eshelby, exact) << ','
                  << maximum_major_asymmetry(hill) << ','
                  << calculation.numHitMaxLevel() << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "emt_withers_benchmark: " << error.what() << '\n';
        return 2;
    }
}
