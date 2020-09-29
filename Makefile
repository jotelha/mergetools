.SILENT:

VMFILES = pkgIndex.tcl mergetools.tcl

VMVERSION = 0.1
DIR = $(PLUGINDIR)/noarch/tcl/mergetools$(VMVERSION)

bins:
win32bins:
dynlibs:
staticlibs:
win32staticlibs:

distrib:
	@echo "Copying mergetools $(VMVERSION) files to $(DIR)"
	mkdir -p $(DIR)
	cp $(VMFILES) $(DIR)
	cp README.md $(DIR)/README
