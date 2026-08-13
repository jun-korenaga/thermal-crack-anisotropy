cmake_minimum_required(VERSION 3.16)

project(thermal_crack_emt
        VERSION 0.1.0
        DESCRIPTION "Hill-tensor calculations for ellipsoidal inclusions in TI media"
        LANGUAGES CXX)

include(CTest)

add_library(emt STATIC
    src/effective_medium.cpp
    src/ellipsoid.cpp
    src/error.cpp
    src/eshelby_exact.cpp
    src/guo2019.cpp
    src/hill_tensor.cpp
    src/tensor.cpp)
target_include_directories(emt PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_compile_features(emt PUBLIC cxx_std_17)

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
    target_compile_options(emt PRIVATE -Wall -Wextra -Wpedantic)
endif()

function(add_emt_app target source)
    add_executable(${target} ${source})
    target_link_libraries(${target} PRIVATE emt)
    target_include_directories(${target} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/apps)
    target_compile_features(${target} PRIVATE cxx_std_17)
endfunction()

add_emt_app(emt_withers_benchmark apps/benchmark_withers.cpp)
add_emt_app(emt_guo_comparison apps/compare_guo2019.cpp)
add_emt_app(emt_random_orientations apps/random_orientations.cpp)
add_emt_app(emt_effective_stiffness apps/effective_stiffness.cpp)

if(BUILD_TESTING)
    add_executable(emt_unit_tests tests/unit_tests.cpp)
    target_link_libraries(emt_unit_tests PRIVATE emt)
    target_compile_features(emt_unit_tests PRIVATE cxx_std_17)
    add_test(NAME unit_tests COMMAND emt_unit_tests)
    add_test(NAME withers_cli
             COMMAND emt_withers_benchmark --aspect 0.1 --ntheta 20)
    add_test(NAME guo_cli
             COMMAND emt_guo_comparison --aspect 0.1 --ntheta 20
                                         --angle-step 90)
    add_test(NAME random_cli
             COMMAND emt_random_orientations --samples 2 --output-every 1
                                              --ntheta 10 --seed 7)
    add_test(NAME effective_cli
             COMMAND emt_effective_stiffness --aspect 0.1 --porosity 0.001
                                             --ntheta 20)
endif()
