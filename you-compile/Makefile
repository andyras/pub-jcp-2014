# ---- set up variables ---- #

DYNAMIX = ~/git/dynamix/wavefunction/bin/dynamix

# --- names of data files to be created --- #
DATADIR = data
# assume that each run produces a tcprob.out file
DATADIRS = $(wildcard $(DATADIR)/*)
RUNFILES = $(patsubst %,%/tcprob.out,$(DATADIRS))
# data which have to be thermally averaged
AVGDATA = $(DATADIR)/sysA_incoherent/avg

# --- names of plot files to be created --- #
PLOTDIR = plots
# assume that the names of plots match the names of the scripts
PLOTSCRIPTS = $(wildcard $(PLOTDIR)/*plt)
PLOTS = $(patsubst %plt,%pdf,$(PLOTSCRIPTS))

# ---- define rules ---- #

.PHONY: all data plots clean

all: data plots

data: $(RUNFILES) $(AVGDATA)

$(AVGDATA):
	@echo Making $@...
	@dir=$@; cd $${dir%/avg}; ./averaging.py $(DYNAMIX) &> log.out

$(RUNFILES):
	@echo Making $@...
	@dir=$@; cd $${dir%/tcprob.out}; $(DYNAMIX) &> log.out

plots: $(PLOTS)

$(PLOTS):
	@echo Plotting $@...
	@cd $$(dirname $@); ../$(patsubst %pdf,%plt,$@)

clean:
	rm -f $(DATADIR)/*/*out
	rm -rf $(DATADIR)/syA_incoherent/avg
	rm -f $(PLOTDIR)/*.pdf
