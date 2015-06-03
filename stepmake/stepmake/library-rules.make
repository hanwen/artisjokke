
$(LIBRARY): $(O_FILES)
	$(AR) $(ARFLAGS) $@ $(O_FILES)
	$(RANLIB) $@ || $(AR) ts $@ || true



