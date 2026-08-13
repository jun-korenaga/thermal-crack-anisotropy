/*
 * error.cc
 */

#include <iostream>
#include <cstdlib>
#include "emt/error.hpp"

void error(const std::string& message)
{
    std::cerr << message << '\n';
    std::exit(1);
}
