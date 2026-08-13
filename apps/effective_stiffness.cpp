#include <exception>
#include <iostream>
#include <string>
#include <vector>

#include "cli.hpp"
#include "emt/effective_medium.hpp"
#include "emt/ellipsoid.hpp"
#include "emt/hill_tensor.hpp"

namespace {

void usage() {
    std::cout
        << "Usage: emt_effective_stiffness [options]\n"
        << "  --configuration NAME single, hexagonal, or horizontal\n"
        << "                       (default: single)\n"
        << "  --aspect VALUE       minor/major semi-axis (default: 0.01)\n"
        << "  --porosity VALUE     inclusion volume fraction (default: 0.001)\n"
        << "  --azimuth DEG        single/hexagonal reference azimuth (default: 0)\n"
        << "  --dip DEG            horizontal-crack rotation about y (default: 0)\n"
        << "  --fluid-bulk VALUE   fluid bulk modulus in GPa (default: 2.2)\n"
        << "  --dry                skip Brown-Korringa fluid substitution\n"
        << "  --ntheta N           initial polar divisions (default: 200)\n"
        << "  --weight VALUE       adaptive weight (default: 0.4)\n"
        << "  --max-level N        recursion cap (default: 20)\n"
        << "  --stiffness CSV      C11,C12,C13,C33,C44 in GPa\n"
        << "                       (default: 57.624,19.88,14.2,41.16,13.48)\n";
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (cli::has_flag(argc, argv, "--help")) {
            usage();
            return 0;
        }
        const std::string configuration =
            cli::option(argc, argv, "--configuration", "single");
        const double aspect = cli::option_double(argc, argv, "--aspect", 0.01);
        const double porosity =
            cli::option_double(argc, argv, "--porosity", 0.001);
        const double azimuth =
            cli::option_double(argc, argv, "--azimuth", 0.0);
        const double dip = cli::option_double(argc, argv, "--dip", 0.0);
        const double fluid_bulk =
            cli::option_double(argc, argv, "--fluid-bulk", 2.2);
        const bool dry = cli::has_flag(argc, argv, "--dry");
        const int ntheta = cli::option_int(argc, argv, "--ntheta", 200);
        const double weight = cli::option_double(argc, argv, "--weight", 0.4);
        const int max_level = cli::option_int(argc, argv, "--max-level", 20);
        cli::require_positive(aspect, "aspect");
        cli::require_positive(ntheta, "ntheta");
        if (porosity < 0.0) {
            throw std::invalid_argument("porosity must be non-negative");
        }

        auto constants = cli::ti_parameters(cli::option(
            argc, argv, "--stiffness", "57.624,19.88,14.2,41.16,13.48"));
        const Tensor4 stiffness(constants);

        std::vector<Ellipsoid> ellipsoids;
        if (configuration == "single") {
            // Vertical crack: the minor semi-axis is initially along y.
            ellipsoids.emplace_back(1.0, aspect, 1.0, azimuth, 0.0, 0.0);
        } else if (configuration == "hexagonal") {
            for (int index = 0; index < 3; ++index) {
                ellipsoids.emplace_back(1.0, aspect, 1.0,
                                        azimuth + 120.0 * index, 0.0, 0.0);
            }
        } else if (configuration == "horizontal") {
            // Horizontal at dip=0; its normal rotates toward x as dip grows.
            ellipsoids.emplace_back(1.0, 1.0, aspect, 0.0, dip, 0.0);
        } else {
            throw std::invalid_argument(
                "configuration must be single, hexagonal, or horizontal");
        }

        Tensor4 average_contribution;
        average_contribution = 0.0;
        for (const auto& ellipsoid : ellipsoids) {
            HillTensor calculation(stiffness, ellipsoid, ntheta, 2 * ntheta,
                                   weight, false, max_level);
            Tensor4 hill;
            calculation.setTensor(hill);
            average_contribution +=
                emt::dry_compliance_contribution(stiffness, hill);
        }
        average_contribution /= static_cast<double>(ellipsoids.size());

        const Tensor4 dry_compliance = emt::dilute_dry_compliance(
            stiffness, average_contribution, porosity);
        const Tensor4 effective_compliance =
            dry ? dry_compliance
                : emt::brown_korringa_saturated_compliance(
                      stiffness, dry_compliance, porosity, fluid_bulk);
        const Tensor4 effective_stiffness = effective_compliance.inverse();

        std::cout << "configuration,aspect,porosity,azimuth_deg,dip_deg,state";
        cli::print_tensor_header(std::cout, "C");
        std::cout << '\n';
        std::cout << configuration << ',' << aspect << ',' << porosity << ','
                  << azimuth << ',' << dip << ','
                  << (dry ? "dry" : "saturated");
        cli::print_tensor(std::cout, effective_stiffness);
        std::cout << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "emt_effective_stiffness: " << error.what() << '\n';
        return 2;
    }
}
