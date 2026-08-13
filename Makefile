CXX ?= c++
AR ?= ar
CPPFLAGS += -Iinclude -Iapps
CXXFLAGS ?= -O3 -std=c++17 -Wall -Wextra -Wpedantic

BUILD_DIR := build
OBJECT_DIR := $(BUILD_DIR)/obj

CORE_SOURCES := \
  src/effective_medium.cpp \
  src/ellipsoid.cpp \
  src/error.cpp \
  src/eshelby_exact.cpp \
  src/guo2019.cpp \
  src/hill_tensor.cpp \
  src/tensor.cpp
CORE_OBJECTS := $(patsubst src/%.cpp,$(OBJECT_DIR)/%.o,$(CORE_SOURCES))
LIBRARY := $(BUILD_DIR)/libemt.a

PROGRAMS := \
  $(BUILD_DIR)/emt_withers_benchmark \
  $(BUILD_DIR)/emt_guo_comparison \
  $(BUILD_DIR)/emt_random_orientations \
  $(BUILD_DIR)/emt_effective_stiffness

.PHONY: all test clean

all: $(PROGRAMS) $(BUILD_DIR)/emt_unit_tests

$(OBJECT_DIR):
	mkdir -p $@

$(OBJECT_DIR)/%.o: src/%.cpp | $(OBJECT_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

$(LIBRARY): $(CORE_OBJECTS)
	mkdir -p $(BUILD_DIR)
	$(AR) rcs $@ $^

$(BUILD_DIR)/emt_withers_benchmark: apps/benchmark_withers.cpp $(LIBRARY)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< $(LIBRARY) -o $@

$(BUILD_DIR)/emt_guo_comparison: apps/compare_guo2019.cpp $(LIBRARY)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< $(LIBRARY) -o $@

$(BUILD_DIR)/emt_random_orientations: apps/random_orientations.cpp $(LIBRARY)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< $(LIBRARY) -o $@

$(BUILD_DIR)/emt_effective_stiffness: apps/effective_stiffness.cpp $(LIBRARY)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< $(LIBRARY) -o $@

$(BUILD_DIR)/emt_unit_tests: tests/unit_tests.cpp $(LIBRARY)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< $(LIBRARY) -o $@

test: all
	$(BUILD_DIR)/emt_unit_tests
	$(BUILD_DIR)/emt_withers_benchmark --aspect 0.1 --ntheta 20 >/dev/null
	$(BUILD_DIR)/emt_guo_comparison --aspect 0.1 --ntheta 20 --angle-step 90 >/dev/null
	$(BUILD_DIR)/emt_random_orientations --samples 2 --output-every 1 --ntheta 10 --seed 7 >/dev/null
	$(BUILD_DIR)/emt_effective_stiffness --aspect 0.1 --porosity 0.001 --ntheta 20 >/dev/null

clean:
	rm -rf $(BUILD_DIR)
