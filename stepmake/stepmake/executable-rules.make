
$(EXECUTABLE): $(O_FILES) $(outdir)/version.hh $(addsuffix /$(outbase)/library.a,$(MODULE_LIBS))
	$(LD) -o $@ $(O_FILES) $(LOADLIBES) $(ALL_LDFLAGS)

%/library.a:
	$(MAKE) -C $(dir $@)
