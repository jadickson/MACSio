# Copyright (c) 2015, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
# Written by Mark C. Miller
#
# LLNL-CODE-676051. All rights reserved.
#
# This file is part of MACSio
# 
# Please also read the LICENSE file at the top of the source code directory or
# folder hierarchy.
# 
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License (as published by the Free Software
# Foundation) version 2, dated June 1991.
# 
# This program is distributed in the hope that it will be useful, but WITHOUT 
# ANY WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the terms and conditions of the GNU General
# Public License for more details.
# 
# You should have received a copy of the GNU General Public License along with
# this program; if not, write to the Free Software Foundation, Inc., 59 Temple
# Place, Suite 330, Boston, MA 02111-1307 USA

# This floating point variable is used to order plugin objects during
# the main link for MACSio to allow dependent libraries that are common
# to multiple plugins to be placed later on the link line. Bigger
# numbers here cause them to appear later in the link line.
TYPHONIO_BUILD_ORDER = 3.0

TYPHONIO_VERSION = 1.6.0
TYPHONIO_FILE = typhonio-$(TYPHONIO_VERSION).tar.gz
TYPHONIO_URL = https://github.com/UK-MAC/typhonio

ifneq ($(TYPHONIO_HOME),)

TYPHONIO_LDFLAGS = -L$(TYPHONIO_HOME)/lib -ltyphonio
TYPHONIO_CFLAGS = -I$(TYPHONIO_HOME)/include

TYPHONIO_SOURCES = macsio_typhonio.c

# Lindstrom's ZFP compression library
ifneq ($(ZFP_HOME),)
TYPHONIO_CFLAGS += -DHAVE_ZFP -I$(ZFP_HOME)/include
TYPHONIO_LDFLAGS += -L$(ZFP_HOME)/lib -lzfp
endif

ifneq ($(SZIP_HOME),)
TYPHONIO_LDFLAGS += -L$(SZIP_HOME)/lib -lsz -Wl,-rpath,$(SZIP_HOME)/lib
TYPHONIO_CFLAGS += -DHAVE_SZIP
endif

ifneq ($(ZLIB_HOME),)
TYPHONIO_LDFLAGS += -L$(ZLIB_HOME)/lib
endif

ifneq ($(HDF5_HOME),)
TYPHONIO_LDFLAGS += -L$(HDF5_HOME)/lib -lhdf5
TYPHONIO_CFLAGS += -I$(HDF5_HOME)/include
endif

TYPHONIO_LDFLAGS += -lz -lm

PLUGIN_OBJECTS += $(TYPHONIO_SOURCES:.c=.o)
PLUGIN_LDFLAGS += $(TYPHONIO_LDFLAGS)
PLUGIN_LIST += typhonio

endif

macsio_typhonio.o: ../plugins/macsio_typhonio.c
	$(CXX) -c $(TYPHONIO_CFLAGS) $(MACSIO_CFLAGS) $(CFLAGS) ../plugins/macsio_typhonio.c

$(TYPHONIO_FILE):
	$(DLCMD) $(TYPHONIO_FILE) $(TYPHONIO_URL)

list-tpls-typhonio:
	@echo "$(TYPHONIO_FILE) ($(TYPHONIO_URL))"

download-tpls-typhonio: $(TYPHONIO_FILE)
