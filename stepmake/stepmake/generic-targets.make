.PHONY : all clean config default diff dist doc exe help html lib TAGS\
	 po

all:	 default
	$(LOOP)

man:
	$(LOOP)

# be careful about deletion.
clean: local-clean
	-rm -f $(outdir)/*
	$(LOOP)

ifneq ($(strip $(depth)),.)
dist:
	make -C $(depth) dist
endif

distclean: clean 
	$(LOOP)
	$(MAKE) local-distclean

maintainerclean: 
	$(LOOP)
	$(MAKE)	local-maintainerclean
	$(MAKE) local-distclean


# configure:
#
config:
	./$(depth)/configure
#


# target help:
#
generic-help:
	@echo -e "\
Makefile for $(PACKAGE_NAME) $(TOPLEVEL_VERSION)\n\
Usage: make ["VARIABLE=value"]... [TARGET]\n\
\n\
Targets:\n"

help: generic-help local-help
	@echo -e "\
  all         update everything\n\
  clean       remove all genated stuff in $(outdir)\n\
  check       run selftest\n\
  default     same as the empty target\n\
  exe         update all executables\n\
  help        this help\n\
  install     install programs and data (prefix=$(prefix))\n\
  lib         update all libraries\n\
  TAGS        genarate tagfiles\n\
\n\
Make may be invoked from any subdirectory\n\
Note that all commands recurse into SUBDIRS;\n\
prepend \`local-' to do only cwd, eg: local-clean\n\
"\

local-help:

local-dist: $(DIST_FILES) $(OUT_DIST_FILES) $(NON_ESSENTIAL_DIST_FILES)
	mkdir -p $(distdir)/$(localdir)
	$(LN) $(DIST_FILES) $(distdir)/$(localdir)

	case "$(NON_ESSENTIAL_DIST_FILES)x" in x) ;; *) \
		$(LN) $(NON_ESSENTIAL_DIST_FILES) $(distdir)/$(localdir);; \
	esac

	case "$(OUT_DIST_FILES)x" in x) ;; *) \
		mkdir -p $(distdir)/$(localdir)/$(outdir); \
		$(LN) $(OUT_DIST_FILES) $(distdir)/$(localdir)/$(outdir);; \
	esac
#	$(foreach i, $(SUBDIRS), $(MAKE) distdir=../$(distdir) localdir=$(localdir)/$(i) -C $(i) local-dist &&) true
# absolute for installed stepmake
	$(foreach i, $(SUBDIRS), $(MAKE) topdir=$(topdir) distdir=$(distdir) localdir=$(localdir)/$(notdir $(i)) -C $(i) local-dist &&) true



html: $(HTML_FILES)

TAGS:
	-if [ "$(TAGS_FILES)" != "" ]; then \
		etags $(ETAGS_FLAGS) $(TAGS_FILES) || \
		ctags $(CTAGS_FLAGS) ".h.hh.tcc.icc" $(TAGS_FILES) $(ERROR_LOG); \
	fi

	$(LOOP)

# ugh . -> $(outdir)
$(outdir)/VERSION: $(depth)/VERSION
	cp -p $< $@

$(outdir)/version.hh: $(outdir)/VERSION
	$(PYTHON) $(step-bindir)/make-version.py $< > $@

configure: configure.in aclocal.m4
	NOCONFIGURE=yes $(srcdir)/autogen.sh
	chmod +x configure

local-clean:

local-distclean:

local-maintainerclean:

install-strip:
	$(MAKE) INSTALL="$(INSTALL) -s" install

final-install:

install: local-install
	$(LOOP)
	$(MAKE) final-install

local-install:

uninstall: local-uninstall
	$(LOOP)

local-uninstall:

installextradoc:
	-$(INSTALL) -d $(prefix)/doc/$(package)
	$(foreach i, $(EXTRA_DOC_FILES),\
		cp -r $(i) $(prefix)/doc/$(package) &&) true

include $(outdir)/dummy.dep $(DEP_FILES) # expect a warning here

$(outdir)/dummy.dep:
	-mkdir -p $(outdir)
	touch $(outdir)/dummy.dep


check: local-check
	$(LOOP)

local-check:

# ugh.  ugh ugh ugh
$(depth)/$(configuration).make: $(depth)/configure
	@echo "************************************************************"
	@echo "configure changed! You should probably reconfigure manually."
	@echo "************************************************************"
	(cd $(depth); ./config.status)
	touch $@		# do something for multiple simultaneous configs.


deb:
	$(MAKE) -C $(depth)/debian
	cd $(depth) && debuild

diff:
	$(PYTHON) $(step-bindir)/package-diff.py  --outdir=$(topdir)/$(outdir) --package=$(topdir) $(makeflags)
	-ln -f $(depth)/$(outdir)/$(distname).diff.gz $(patch-dir)

release: 
	$(PYTHON) $(step-bindir)/release.py --outdir=$(topdir)/$(outdir) --package=$(topdir)

rpm-mandrake: $(depth)/$(package-icon) dist
	$(MAKE) -C $(depth)/make
	cd $(depth) && rpm -ba make/$(outdir)/lilypond.mandrake.spec

rpm: $(depth)/$(package-icon) dist
	@echo "Assuming Red Hat system"	#FIXME: check distro, then issue rpm
	$(MAKE) -C $(depth)/make
	cd $(depth) && rpm -bb make/$(outdir)/lilypond.redhat.spec
#	su -c 'rpm -tb $(depth)/$(outdir)/$(distname).tar.gz'

