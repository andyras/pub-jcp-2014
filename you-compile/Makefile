# ---- set up variables ---- #

CODEDIR = code
DYNAMIX = $(CURDIR)/$(CODEDIR)/dynamix

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

all: code data plots
	@echo Done making everything.

code: $(DYNAMIX)
	@echo Done making code.

$(DYNAMIX):
	$(MAKE) -C $(CODEDIR)

data: $(RUNFILES) $(AVGDATA)
	@echo Done making data.

$(AVGDATA):
	@echo Making $@...
	@dir=$@; cd $${dir%/avg}; ./averaging.py $(DYNAMIX) &> dynamix.log

$(RUNFILES):
	@echo Making $@...
	@dir=$@; cd $${dir%/tcprob.out}; $(DYNAMIX) &> dynamix.log

plots: $(PLOTS)
	@echo Dome making plots.

$(PLOTS):
	@echo Plotting $@...
	@cd $$(dirname $@); ../$(patsubst %pdf,%plt,$@)

clean:
	rm -f $(DATADIR)/*/*out
	rm -f $(DATADIR)/*/*log
	rm -rf $(DATADIR)/syA_incoherent/avg
	rm -f $(PLOTDIR)/*.pdf
	$(MAKE) -C $(CODEDIR) clean
