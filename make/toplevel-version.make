# Toplevel_version.make

include $(depth)/VERSION
TOPLEVEL_MAJOR_VERSION:=$(MAJOR_VERSION)
TOPLEVEL_MINOR_VERSION:=$(MINOR_VERSION)
TOPLEVEL_PATCH_LEVEL:=$(PATCH_LEVEL)
# use to send patches, always empty for released version:
TOPLEVEL_MY_PATCH_LEVEL:=$(MY_PATCH_LEVEL)
