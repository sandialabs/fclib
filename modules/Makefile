## Copyright (2000) Sandia Corporation. Under the terms of Contract
## DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains 
## certain rights in this software.
##
## Redistribution and use in source and binary forms, with or without
## modification, are permitted provided that the following conditions are met:
##
##   * Redistributions of source code must retain the above copyright
##     notice, this list of conditions and the following disclaimer.
##
##   * Redistributions in binary form must reproduce the above copyright
##     notice, this list of conditions and the following disclaimer in the
##     documentation and/or other materials provided with the
##     distribution.
##
##   * Neither the name of Sandia nor the names of any contributors may
##     be used to endorse or promote products derived from this software
##     without specific prior written permission.
##
##   * Modified source versions must be plainly marked as such, and must
##     not be misrepresented as being the original software.
##
## THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
## AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
## IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
## ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHORS OR CONTRIBUTORS BE LIABLE FOR 
## ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
## DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
## SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
## CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT 
## LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
## OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH 
## DAMAGE.

## $Source: /home/Repositories/fcdmf/fclib/modules/Makefile,v $
## $Revision: 1.34 $
## $Date: 2006/08/30 19:20:01 $
## 
## Makefile for modules
##
include ../Makefile.include

##
## Files (in alphabetical order)
##

LIB_PRIVATEHFILES = \
	datasetP.h \
	elemdeathP.h \
	featureP.h \
	fileioP.h \
	geomP.h \
	libraryP.h \
	lsdynaioP.h \
	meshP.h \
	sequenceP.h \
	seriesP.h \
	shapeP.h \
	sierraP.h \
	storageP.h \
	subsetP.h \
	tableP.h \
	topoP.h \
	variableP.h \
	varmathP.h

# IO stuff that has external dependencies
LIB_EXO_PRIVATEHFILES = exodusioP.h

LIB_IO_PRIVATEHFILES = ${LIB_EXO_PRIVATEHFILES} 

LIB_PUBLICHFILES = \
	base.h \
	dataset.h \
	elemdeath.h \
	error.h \
	feature.h \
	fileio.h \
	geom.h \
	library.h \
	mesh.h \
        sequence.h \
	series.h \
	shape.h \
	sierra.h \
	statistics.h \
	storage.h \
	subset.h \
	threshold.h \
	topo.h \
	track.h \
	util.h \
	variable.h \
	varmath.h \
	custom.h

LIB_SUPPORTHFILES = fc.h fcP.h

LIB_SRCS = \
	base.c \
	dataset.c \
	elemdeath.c \
	feature.c \
	fileio.c \
	geom.c \
        library.c \
	lsdynaio.c \
	mesh.c \
	sequence.c \
	series.c \
	shape.c \
	sierra.c \
	statistics.c \
	storage.c \
	subset.c \
	table.c \
	threshold.c \
	topo.c \
	util.c \
	track.c \
	variable.c \
	varmath.c \
	custom.c

#IO stuff that has external dependencies
LIB_EXO_SRCS = exodusio.c
LIB_IO_SRCS = ${LIB_EXO_SRCS} 

# all files for release target
ALL_FILES = Makefile ${LIB_PRIVATEHFILES} ${LIB_IO_PRIVATEHFILES} \
	${LIB_PUBLICHFILES} ${LIB_SUPPORTHFILES} ${LIB_SRCS} ${LIB_IO_SRCS}

# files for this build
LIB_HFILES = $(LIB_PRIVATEHFILES) $(LIB_PUBLICHFILES) $(LIB_SUPPORTHFILES)
LIB_OBJS = $(LIB_SRCS:.c=.o)
ifeq (${HAVE_EXODUS},yes)
LIB_HFILES += ${LIB_EXO_PRIVATEHFILES}
LIB_OBJS += ${LIB_EXO_SRCS:.c=.o}
endif

##
## Targets (ranlib added to support OSX)
##

lib : $(FCLIB)

$(FCLIB) : $(LIB_OBJS) 
	ar -r $@ $?
	ranlib $@

clean :
	$(RM) *.o $(FCLIB) $(GCOV_FILES)

## Additional dependencies 
## Explicit depend on fc.h & fcP.h so that new lib build forces
## new program builds, so changes to them get tested
## (explicit deps too hard to maintain) (FIX? use make depends?)
$(LIB_OBJS) : $(LIB_HFILES)

#generic compiling directive

%.o : %.c
	$(CC) -c -fPIC -o $@ $(CFLAGS) $(CPPFLAGS) $<

# test changes
regtest :
	(cd ..; make regtest)
check :
	(cd ..; make check)

install : $(FCLIB)
	cp -f $(LIB_HFILES) $(INSTALL_DIR)/include
	-for i in $(LIB_HFILES); do \
	   chmod go=r $(INSTALL_DIR)/include/$$i; \
	done
	cp -f $(FCLIB) $(INSTALL_DIR)/lib
	chmod go=r $(INSTALL_DIR)/lib/$(FCLIB)
	ranlib $(INSTALL_DIR)/lib/$(FCLIB)

uninstall :
	-for i in $(LIB_HFILES); do \
	   rm -f $(INSTALL_DIR)/include/$$i; \
	done
	rm -f $(INSTALL_DIR)/lib/$(FCLIB)

release :
	-for i in $(ALL_FILES); do \
	  (cd ..; tar --append --file=$(FC_RELEASE).tar \
	    $(FC_RELEASE)/modules/$$i) \
	done

.PHONY : lib clean regtest check install uninstall release
