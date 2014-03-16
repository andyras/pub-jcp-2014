# name of main executable
### This line should be just about all you have to change ###
BIN = dynamix
#############################################################

# directories for headers, objects and source files
INCDIR = ./include
OBJDIR = ./obj
SRCDIR = ./src
BINDIR = ../bin

# compiler
CXX = icpc

# flags for compiler
# optimized CXXFLAGS
# Intel
#CXXFLAGS = -O0 -Wall --std=c++11 -openmp
CXXFLAGS = -O3 -fast -xHOST -no-prec-div -Wall --std=c++11 -fopenmp -mkl -no-multibyte-chars

# non-optimized CXXFLAGS
#CXXFLAGS = -Wall -g -O0 -openmp --std=c++11

# flags for linker
INCLUDES = -I$(INCDIR) -I$(MKLROOT)/include/fftw
# includes for CVODE
INCLUDES += -I$(HOME)/bin/include -I$(HOME)/bin/include/cvode -L$(HOME)/bin/lib -I/usr/include/sundials -I/usr/include/cvode -L/usr/lib
#LDFLAGS = -lsundials_cvode -lsundials_nvecserial -lgomp -lfftw3
# Intel flags for linker
LDFLAGS = -lsundials_cvode -lsundials_nvecserial -mkl

# make will look for .cpp files in $(SRCDIR)
vpath %.cpp $(SRCDIR)

# list of source files
SOURCES = $(wildcard $(SRCDIR)/*.cpp)

# object files have same name as .cpp files, but with .o extension
# TODO simplify this line
OBJECTS = $(patsubst $(SRCDIR)/%.cpp,obj/%.o,$(SOURCES))

# build the main executable; this should be listed first
$(BIN): $(OBJECTS)
	$(CXX) -o $@ $^ $(LDFLAGS) $(INCLUDES)

# automatic rule for building objects
# TODO combine this line with the next one
$(OBJDIR)/%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

# all objects depend on object directory
$(OBJECTS): | $(OBJDIR)

# make object directory if it does not exist
$(OBJDIR):
	mkdir -p $(OBJDIR)

.PHONY: clean build install uninstall $(BINDIR)
clean:
	rm -rf $(OBJDIR)
	rm -f $(BIN)

build:
	make clean
	make

install: $(BIN) | $(BINDIR)
	cp dynamix $(BINDIR)

uninstall:
	rm -rf $(BINDIR)/$(BIN)

$(BINDIR):
	mkdir -p ../bin
