CXX = gcc
CXXFLAGS = -fopenmp

RM = rm -f
MKDIRS = mkdir -p

TARGETS = $(basename $(wildcard omp-*.*))

all : $(TARGETS)

%:%.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

clean:
	$(RM) $(TARGETS)

