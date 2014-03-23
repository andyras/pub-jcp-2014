# ---- variables for building ---- #

# the name of the executable
BIN = dynamix
# directories for headers, objects and source files
INCDIR = include
OBJDIR = obj
SRCDIR = src

# ---- [DO NOT CHANGE] base flags ---- #
LDFLAGS = -lsundials_cvode -lsundials_nvecserial
INCLUDES = -I./$(INCDIR)
CXXFLAGS = -O3 -Wall

# ---- [CHANGE] compiler-specific flags ---- #
# --- Intel compiler --- #
CXX = icpc
CXXFLAGS += -fopenmp
# -- MKL takes care of BLAS/LAPACK -- #
LDFLAGS += -mkl

# --- GCC --- #
#CXX = g++
#CXXFLAGS +=

# ---- [CHANGE] locations of libraries ---- #
# --- flags for BLAS --- #
#LDFLAGS += -lgslcblas
#INCLUDES += -I/usr/local/Cellar/gsl/1.16/include/gsl/ -L/usr/local/Cellar/gsl/1.16/lib/gsl/
# --- flags for LAPACK --- #
#LDFLAGS += -llapack
#INCLUDES += -I/usr/local/Cellar/lapack/3.4.2/include -L/usr/local/Cellar/lapack/3.4.2/lib

# make will look for .cpp files in $(SRCDIR)
vpath %.cpp $(SRCDIR)

# list of source files
SOURCES = $(wildcard $(SRCDIR)/*.cpp)

# object files have same name as .cpp files, but with .o extension
OBJECTS = $(patsubst $(SRCDIR)/%.cpp,obj/%.o,$(SOURCES))

# ---- build rules ---- #

# build the main executable; this should be listed first
$(BIN): $(OBJECTS)
	$(CXX) -o $@ $^ $(LDFLAGS) $(INCLUDES)

# automatic rule for building objects
$(OBJDIR)/%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

# all objects depend on object directory
$(OBJECTS): | $(OBJDIR)

# make object directory if it does not exist
$(OBJDIR):
	mkdir -p $(OBJDIR)

.PHONY: clean

clean:
	rm -rf $(OBJDIR)
	rm -f $(BIN)
