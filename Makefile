# ---- [DO NOT CHANGE] base flags ---- #
LDFLAGS = -lsundials_cvode -lsundials_nvecserial
INCLUDES =
CXXFLAGS = -O3 -xHOST -no-prec-div -Wall -std=c++11 -no-multibyte-chars
CPPFLAGS =

# ---- [CHANGE] compiler-specific flags ---- #
# --- Intel compiler --- #
CXX = icpc
CXXFLAGS += -fopenmp
LDFLAGS += -mkl
# --- GCC --- #
#CXX = g++
#CXXFLAGS +=
#CPPFLAGS +=
#LDFLAGS += -llapack

# ---- [CHANGE] set if SUNDIALS or LAPACK are in non-standard places ---- #
#INCLUDES += /foo/bar/includes

# ---- [DO NOT CHANGE] below this line unless you really mean it ---- #
BUILDDIR = build
BINDIR = bin
DATADIR = data
PLOTDIR = plots

# absolute path to dynamix executable
DYNAMIX = $(CURDIR)/$(BINDIR)/dynamix


#include Makevars


all: code data plots

code:
	$(MAKE) -C $(BUILDDIR)

data: code
	cd $(DATADIR)
	for datumdir in $$(find . -mindepth 1 -maxdepth 1 -type d)
	do
	  cd $${datumdir}
	  $(DYNAMIX)
	  cd ..
	done

plots: data | $(PLOTDIR)
	cd $(PLOTDIR)
	for plt in $$(ls *.plt)
	do
	  ./$${plt}
	done

pdfplots: plots
	if command -v ps2pdf &> /dev/null; then
	  cd $(PLOTDIR)
	  for eps in $$(ls *.eps)
	  do
	    ps2pdf $${eps}
	  done
	else
	  echo "ps2pdf not installed"
	fi

pngplots: plots
	if command -v convert &> /dev/null; then
	  cd $(PLOTDIR)
	  for eps in $$(ls *.eps)
	  do
	    convert -density 300 $${eps} $${eps%eps}png
	  done
	else
	  echo "convert not installed"
	fi


# here are rules for cleaning up and for making dirs if they don't exist

DIRS = $(BINDIR) $(PLOTDIR)

CLEANS = clean-code clean-data clean-plots

.PHONY: clean $(CLEANS) $(DIRS)

clean: $(CLEANS)

clean-code:
	cd build
	make clean

clean-data:
	rm -f $(DATADIR)/*/*.out

clean-plots:
	rm -f $(PLOTDIR)/*.{eps,pdf,png}

$(BINDIR):
	mkdir -p $(BINDIR)

$(PLOTDIR):
	mkdir -p $(PLOTDIR)
