#include <cmath>
#include <cstdint>
#include <exception>
#include <iostream>
#include <random>

#include "cli.hpp"
#include "emt/effective_medium.hpp"
#include "emt/ellipsoid.hpp"
#include "emt/hill_tensor.hpp"
#include "emt/published_models.hpp"

namespace {

void usage() {
    std::cout
        << "Usage: emt_random_orientations [options]\n"
        << "  --samples N          number of crack orientations (default: 100)\n"
        << "  --seed N             deterministic random seed (default: 20260813)\n"
        << "  --output-every N     running-average interval (default: 10)\n"
        << "  --aspect VALUE       a3/a1 (default: 0.01)\n"
        << "  --ntheta N           initial polar divisions (default: 100)\n"
        << "  --weight VALUE       adaptive weight (default: 0.4)\n"
        << "  --max-level N        recursion cap (default: 20)\n"
        << "  --stiffness CSV      C11,C12,C13,C33,C44 in GPa\n"
        << "                       (default: 47,8,5,34,17)\n"
        << "Output includes adaptive running averages and the corrected\n"
        << "Xu et al. (2020) comparison tensor.\n";
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (cli::has_flag(argc, argv, "--help")) {
            usage();
            return 0;
        }
        const int samples = cli::option_int(argc, argv, "--samples", 100);
        const int output_every =
            cli::option_int(argc, argv, "--output-every", 10);
        const std::uint32_t seed = static_cast<std::uint32_t>(
            cli::option_int(argc, argv, "--seed", 20260813));
        const double aspect = cli::option_double(argc, argv, "--aspect", 0.01);
        const int ntheta = cli::option_int(argc, argv, "--ntheta", 100);
        const double weight = cli::option_double(argc, argv, "--weight", 0.4);
        const int max_level = cli::option_int(argc, argv, "--max-level", 20);
        cli::require_positive(samples, "samples");
        cli::require_positive(output_every, "output-every");
        cli::require_positive(aspect, "aspect");
        cli::require_positive(ntheta, "ntheta");

        auto constants = cli::ti_parameters(cli::option(
            argc, argv, "--stiffness", "47,8,5,34,17"));
        const Tensor4 stiffness(constants);

        std::mt19937 generator(seed);
        std::uniform_real_distribution<double> uniform(0.0, 1.0);
        Tensor4 average;
        average = 0.0;

        const Tensor4 xu2020 =
            emt::xu2020_random_compliance_contribution(stiffness, aspect);

        std::cout << "samples,seed,model";
        cli::print_tensor_header(std::cout, "H");
        std::cout << '\n';
        for (int sample = 1; sample <= samples; ++sample) {
            // For an axisymmetric crack, normals n and -n are equivalent. A
            // uniform hemisphere therefore samples azimuth uniformly and
            // cos(beta) uniformly on [0,1].
            const double alpha = 360.0 * uniform(generator);
            const double beta =
                std::acos(uniform(generator)) * 180.0 / M_PI;
            const Ellipsoid ellipsoid(1.0, 1.0, aspect, alpha, beta, 0.0);
            HillTensor calculation(stiffness, ellipsoid, ntheta, 2 * ntheta,
                                   weight, false, max_level);
            Tensor4 hill;
            calculation.setTensor(hill);
            const Tensor4 contribution =
                emt::dry_compliance_contribution(stiffness, hill);

            average *= static_cast<double>(sample - 1) / sample;
            Tensor4 increment = contribution;
            increment /= static_cast<double>(sample);
            average += increment;

            if (sample % output_every == 0 || sample == samples) {
                std::cout << sample << ',' << seed << ",adaptive";
                cli::print_tensor(std::cout, average);
                std::cout << '\n';
                std::cout << sample << ',' << seed << ",xu2020_corrected";
                cli::print_tensor(std::cout, xu2020);
                std::cout << '\n';
            }
        }
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "emt_random_orientations: " << error.what() << '\n';
        return 2;
    }
}
