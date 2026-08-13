#include <cmath>
#include <exception>
#include <iostream>

#include "cli.hpp"
#include "emt/effective_medium.hpp"
#include "emt/ellipsoid.hpp"
#include "emt/hill_tensor.hpp"
#include "emt/published_models.hpp"

namespace {

void usage() {
    std::cout
        << "Usage: emt_guo_comparison [options]\n"
        << "  --aspect VALUE       a3/a1 (default: 0.01)\n"
        << "  --ntheta N           initial polar divisions (default: 200)\n"
        << "  --weight VALUE       adaptive weight (default: 0.4)\n"
        << "  --max-level N        recursion cap (default: 20)\n"
        << "  --angle-step DEG     inclination increment (default: 10)\n"
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
        const int ntheta = cli::option_int(argc, argv, "--ntheta", 200);
        const double weight = cli::option_double(argc, argv, "--weight", 0.4);
        const int max_level = cli::option_int(argc, argv, "--max-level", 20);
        const double angle_step =
            cli::option_double(argc, argv, "--angle-step", 10.0);
        cli::require_positive(aspect, "aspect");
        cli::require_positive(ntheta, "ntheta");
        cli::require_positive(angle_step, "angle-step");

        auto constants = cli::ti_parameters(cli::option(
            argc, argv, "--stiffness", "47,8,5,34,17"));
        const Tensor4 stiffness(constants);

        std::cout << "inclination_deg,model";
        cli::print_tensor_header(std::cout, "H");
        std::cout << '\n';
        for (double angle = 0.0; angle <= 90.0 + 1.0e-10;
             angle += angle_step) {
            const Ellipsoid ellipsoid(1.0, 1.0, aspect, 0.0, 0.0, angle);
            HillTensor calculation(stiffness, ellipsoid, ntheta, 2 * ntheta,
                                   weight, false, max_level);
            Tensor4 hill;
            calculation.setTensor(hill);
            const Tensor4 numerical =
                emt::dry_compliance_contribution(stiffness, hill);
            const Tensor4 guo = emt::guo2019_compliance_contribution(
                stiffness, aspect, angle * M_PI / 180.0);

            std::cout << angle << ",adaptive";
            cli::print_tensor(std::cout, numerical);
            std::cout << '\n';
            std::cout << angle << ",guo2019_corrected";
            cli::print_tensor(std::cout, guo);
            std::cout << '\n';
        }
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "emt_guo_comparison: " << error.what() << '\n';
        return 2;
    }
}
