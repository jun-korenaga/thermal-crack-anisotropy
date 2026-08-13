#ifndef EMT_APP_CLI_HPP
#define EMT_APP_CLI_HPP

#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "emt/tensor.hpp"

namespace cli {

inline bool has_flag(int argc, char** argv, const std::string& name) {
    for (int i = 1; i < argc; ++i) {
        if (argv[i] == name) return true;
    }
    return false;
}

inline std::string option(int argc, char** argv, const std::string& name,
                          const std::string& fallback) {
    for (int i = 1; i < argc; ++i) {
        if (argv[i] == name) {
            if (i + 1 >= argc) {
                throw std::invalid_argument("missing value after " + name);
            }
            return argv[i + 1];
        }
    }
    return fallback;
}

inline double option_double(int argc, char** argv, const std::string& name,
                            double fallback) {
    return std::stod(option(argc, argv, name, std::to_string(fallback)));
}

inline int option_int(int argc, char** argv, const std::string& name,
                      int fallback) {
    return std::stoi(option(argc, argv, name, std::to_string(fallback)));
}

inline std::vector<double> parse_csv_doubles(const std::string& value,
                                             std::size_t expected) {
    std::vector<double> result;
    std::stringstream stream(value);
    std::string item;
    while (std::getline(stream, item, ',')) result.push_back(std::stod(item));
    if (result.size() != expected) {
        throw std::invalid_argument("expected " + std::to_string(expected) +
                                    " comma-separated numbers, got: " + value);
    }
    return result;
}

inline Array1d<double> ti_parameters(const std::string& value) {
    const auto parsed = parse_csv_doubles(value, 5);
    Array1d<double> result(5);
    for (int i = 1; i <= 5; ++i) result(i) = parsed[i - 1];
    return result;
}

inline void print_tensor_header(std::ostream& out, const std::string& prefix) {
    for (int I = 1; I <= 6; ++I) {
        for (int J = 1; J <= 6; ++J) out << ',' << prefix << I << J;
    }
}

inline void print_tensor(std::ostream& out, const Tensor4& tensor) {
    out << std::setprecision(17);
    for (int I = 1; I <= 6; ++I) {
        for (int J = 1; J <= 6; ++J) out << ',' << tensor.val(I, J);
    }
}

inline void require_positive(double value, const std::string& name) {
    if (value <= 0.0) throw std::invalid_argument(name + " must be positive");
}

inline void require_positive(int value, const std::string& name) {
    if (value <= 0) throw std::invalid_argument(name + " must be positive");
}

}  // namespace cli

#endif
