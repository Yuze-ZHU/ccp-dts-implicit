CXX ?= g++
CXXFLAGS ?= -O3 -std=c++17 -Wall -Wextra -Wpedantic -Wno-unused-parameter
DEPS_DIR ?= ../../ccpDeps/install
EIGEN_DIR ?= ../../eigen-3.3.9
CPPFLAGS += -I$(DEPS_DIR)/include -I$(EIGEN_DIR)
LDFLAGS += -L$(DEPS_DIR)/lib
LDLIBS += -lsacado -lteuchoscore

TARGET = bin/ccp-dts-implicit
SOURCES = $(wildcard src/*.cpp)
OBJECTS = $(patsubst src/%.cpp,build/%.o,$(SOURCES))

all: $(TARGET)

$(TARGET): $(OBJECTS) | bin
	$(CXX) $(CXXFLAGS) $(OBJECTS) -o $@ $(LDFLAGS) $(LDLIBS)

build/%.o: src/%.cpp src/ccp-2d.hpp | build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

bin build:
	mkdir -p $@

run: $(TARGET)
	./$(TARGET) -c cases/lymberopoulos1993.case

clean:
	rm -rf build bin

.PHONY: all run clean
