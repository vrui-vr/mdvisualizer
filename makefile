########################################################################
# Makefile for MD Visualizer.
# Copyright (c) 2019-2025 Oliver Kreylos
# 
# This file is part of the MD Visualizer (MDVisualizer).
# 
# The MD Visualizer is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation; either version 2 of the
# License, or (at your option) any later version.
# 
# The MD Visualizer is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with the MD Visualizer; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
# USA
########################################################################

# Directory containing the Vrui build system. The directory below
# matches the default Vrui installation; if Vrui's installation
# directory was changed during Vrui's installation, the directory below
# must be adapted.
VRUI_MAKEDIR = /usr/local/share/Vrui-14.1/make

# Base installation directory for the MD Visualizer. If this is set to
# the default of $(PROJECT_ROOT), the MD Visualizer does not have to be
# installed to be run. Created executables, configuration files, and
# resources will be installed in the bin, etc, and share directories
# under the given base directory, respectively.
# Important note: Do not use ~ as an abbreviation for the user's home
# directory here; use $(HOME) instead.
INSTALLDIR = $(PROJECT_ROOT)

########################################################################
# Everything below here should not have to be changed
########################################################################

# Name of the project
PROJECT_NAME = MDVisualizer
PROJECT_DISPLAYNAME = MD Visualizer

# Version number for installation subdirectories. This is used to keep
# subsequent release versions of the MD Visualizer from clobbering each
# other. The value should be identical to the major.minor version number
# found in VERSION in the root package directory.
PROJECT_MAJOR = 1
PROJECT_MINOR = 6

# Include definitions for the system environment and system-provided
# packages
include $(VRUI_MAKEDIR)/SystemDefinitions
include $(VRUI_MAKEDIR)/Packages.System
include $(VRUI_MAKEDIR)/Configuration.Vrui
include $(VRUI_MAKEDIR)/Packages.Vrui
include $(VRUI_MAKEDIR)/Configuration.Collaboration
include $(VRUI_MAKEDIR)/Packages.Collaboration

########################################################################
# Specify additional compiler and linker flags
########################################################################

CFLAGS += -Wall -pedantic

########################################################################
# List common packages used by all components of this project
# (Supported packages can be found in $(VRUI_MAKEDIR)/Packages.*)
########################################################################

PACKAGES = 

########################################################################
# Specify all final targets
########################################################################

EXECUTABLES = $(EXEDIR)/MDVisualizer

.PHONY: all
all: $(EXECUTABLES)

########################################################################
# Pseudo-target to print configuration options and configure the package
########################################################################

.PHONY: config config-invalidate
config: config-invalidate $(DEPDIR)/config

config-invalidate:
	@mkdir -p $(DEPDIR)
	@touch $(DEPDIR)/Configure-Begin

$(DEPDIR)/Configure-Begin:
	@mkdir -p $(DEPDIR)
	@echo "---- $(PROJECT_FULLDISPLAYNAME) configuration options: ----"
	@touch $(DEPDIR)/Configure-Begin

$(DEPDIR)/Configure-Package: $(DEPDIR)/Configure-Begin
	@echo "Collaborative visualization enabled"
	@touch $(DEPDIR)/Configure-Package

$(DEPDIR)/Configure-Install: $(DEPDIR)/Configure-Package
	@echo "---- $(PROJECT_FULLDISPLAYNAME) installation configuration ----"
	@echo "Installation directory: $(INSTALLDIR)"
	@echo "Executable directory: $(EXECUTABLEINSTALLDIR)"
	@touch $(DEPDIR)/Configure-Install

$(DEPDIR)/Configure-End: $(DEPDIR)/Configure-Install
	@echo "---- End of $(PROJECT_FULLDISPLAYNAME) configuration options ----"
	@touch $(DEPDIR)/Configure-End

$(DEPDIR)/config: $(DEPDIR)/Configure-End $(CONFIGFILES)
	@touch $(DEPDIR)/config

########################################################################
# Specify other actions to be performed on a `make clean'
########################################################################

.PHONY: extraclean
extraclean:

.PHONY: extrasqueakyclean
extrasqueakyclean:

# Include basic makefile
include $(VRUI_MAKEDIR)/BasicMakefile

########################################################################
# Specify build rules for executables
########################################################################

MDVISUALIZER_SOURCES = Color.cpp \
                       Atom.cpp \
                       Trace.cpp \
                       CharmmTrace.cpp \
                       LammpsTrace.cpp \
                       EspressoTrace.cpp \
                       MDVisualizer.cpp

$(MDVISUALIZER_SOURCES:%.cpp=$(OBJDIR)/%.o): | $(DEPDIR)/config

$(EXEDIR)/MDVisualizer: PACKAGES += MYCOLLABORATION2CLIENT VRUIALL
$(EXEDIR)/MDVisualizer: $(MDVISUALIZER_SOURCES:%.cpp=$(OBJDIR)/%.o)
.PHONY: MDVisualizer
MDVisualizer: $(EXEDIR)/MDVisualizer

########################################################################
# Specify installation rules for header files, libraries, executables,
# configuration files, and shared files.
########################################################################

install: $(ALL)
	@echo Installing $(PROJECT_DISPLAYNAME) in $(INSTALLDIR)...
	@install -d $(INSTALLDIR)
	@install -d $(EXECUTABLEINSTALLDIR)
	@install $(EXECUTABLES) $(EXECUTABLEINSTALLDIR)
