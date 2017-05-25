#

CXX       ?= g++
CFLAGS    ?= -Wall -Wno-parentheses -std=gnu++0x -O3 -fopenmp
#CFLAGS    ?= -Wall -Wno-parentheses -std=gnu++11 -O3 -fopenmp
#COPTIMIZE ?= -O3

main: main.cpp coarsening.cpp
	@echo "compiling..."
	$(CXX) $(CFLAGS) -o $@ $^

.PHONY: clean
clean:
	$(RM) *~ main
