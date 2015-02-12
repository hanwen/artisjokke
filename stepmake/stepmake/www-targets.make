
local-WWW:

ifneq ($(strip $(depth)),.)
WWW: local-WWW
	$(LOOP)

WWW-clean: local-WWW-clean
	$(LOOP)
endif

#ugh, this is (but whole make web/www/WWW is) lilypond specific
local-web:
	$(MAKE) out=www LILYPOND_BOOK_FORMAT=texi-html local-WWW

web:
	$(MAKE) out=www LILYPOND_BOOK_FORMAT=texi-html WWW

local-WWW-clean:
	rm -f $(outdir)/*

local-web-clean:
	$(MAKE) out=www local-WWW-clean

web-clean:
	$(MAKE) out=www WWW-clean

local-help: www-targets-help

www-targets-help:
	@echo -e "\
  web         update website in out-www\n\
  web-clean   clean out-www\n\
"
