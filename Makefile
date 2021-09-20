SOURCES=$(shell find . -name "*.cpp")
CXXFLAGS= -std=c++17 -Wall -O3 -ffast-math -march=native
OBJECTS=$(SOURCES:%.cpp=%.o)
TARGET=main

.PHONY: all
all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(LINK.cpp) $^ -std=c++17 $(LOADLIBES) $(LDLIBS) -o $@

.PHONY: clean
clean:
	rm -rf $(OBJECTS)