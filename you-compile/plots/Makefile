# ---- set up variables ---- #

DYNAMIX = ~/git/dynamix/wavefunction/bin/dynamix

DATAPREFIX = sys
# assume that each run produces a tcprob.out file
DATADIRS = $(wildcard $(DATAPREFIX)*)
RUNFILES = $(patsubst %,%/tcprob.out,$(DATADIRS))

# assume that the names of plots match the names of the scripts
PLOTSCRIPTS = $(wildcard *plt)
PLOTS = $(patsubst %plt,%pdf,$(PLOTSCRIPTS))

# ---- define rules ---- #

.PHONY: all data plot clean

all: data plot

data: $(RUNFILES) sysA_incoherent/avg

sysA_incoherent/avg:
	@echo Making $@...
	@dir=$@; cd $${dir%/avg}; ./averaging.py $(DYNAMIX)

$(RUNFILES):
	@echo Making $@...
	@dir=$@; cd $${dir%/tcprob.out}; $(DYNAMIX)

plot: $(PLOTS)

$(PLOTS):
	./$(patsubst %pdf,%plt,$@)

clean:
	rm -f */*out
	rm -rf syA_incoherent/avg
	rm -f fig2numAna.pdf
