# Thin wrapper over run.sh so the common invocations are discoverable and the
# CI job is one word. Everything here is available directly from ./run.sh;
# nothing but argument spelling lives in this file.
#
# Override any setting from env/cases.conf on the command line:
#   make smoke DCS_HOME=/tmp/dcs-tools-1.2.0 THREADS=32

.PHONY: help static smoke full heavy golden list clean check

help:
	@echo "targets:"
	@echo "  static   packaging + CLI contract checks, no test data needed"
	@echo "  smoke    static plus the WGS case          (default gate)"
	@echo "  full     every case, including heavy datasets"
	@echo "  heavy    only the large/ploidy cases (30x, mouse, rice, potato)"
	@echo "  golden   re-run the full tier and install its products as baseline"
	@echo "  list     show the tiers and the available cases"
	@echo "  check    syntax-check the suite itself"
	@echo "  clean    remove work/ and reports/, keep golden/"
	@echo
	@echo "settings come from env/cases.conf and can be overridden:"
	@echo "  make smoke DCS_HOME=/path/to/dcs-tools THREADS=32"

static:
	@./run.sh --tier static

smoke:
	@./run.sh --tier smoke

full:
	@./run.sh --tier full

heavy:
	@./run.sh --tier heavy

# Deliberately not a dependency of anything: refreshing the baseline is a
# decision, not a step. A golden updated automatically on failure is how a
# regression suite quietly stops catching regressions.
golden:
	@./run.sh --tier full --update-golden

list:
	@./run.sh --list

check:
	@for f in run.sh cases/*.sh; do bash -n $$f && echo "ok  $$f"; done
	@for f in scripts/*.py scripts/lib/*.py; do \
		python3 -m py_compile $$f && echo "ok  $$f"; done

clean:
	@rm -rf work reports
	@echo "removed work/ and reports/  (golden/ kept; delete it by hand)"
