CXX ?= c++
PYTHON ?= python3
LATEXMK ?= latexmk
EMT_ROOT ?= ..
CPPFLAGS += -I$(EMT_ROOT)/include -I.
CXXFLAGS ?= -O3 -std=c++17 -Wall -Wextra -Wpedantic

BUILD_DIR := build
PROGRAM := $(BUILD_DIR)/benchmark_kanaun09
EMT_LIBRARY := $(EMT_ROOT)/build/libemt.a

.PHONY: all benchmark figures note clean

all: $(PROGRAM)

$(EMT_LIBRARY):
	$(MAKE) -C $(EMT_ROOT) all

$(BUILD_DIR):
	mkdir -p $@

$(PROGRAM): benchmark_kanaun09.cc kanaun09.cc kanaun09.h $(EMT_LIBRARY) | $(BUILD_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) benchmark_kanaun09.cc kanaun09.cc \
	  $(EMT_LIBRARY) -o $@

benchmark: $(PROGRAM)
	mkdir -p results
	$(PROGRAM) -P360 -Q64 -N800 -W0.4 > results/benchmark.csv

figures:
	cd note && $(PYTHON) make_figures.py

note: figures
	cd note && $(LATEXMK) -pdf -interaction=nonstopmode -halt-on-error benchmark_note.tex

clean:
	rm -rf $(BUILD_DIR)
