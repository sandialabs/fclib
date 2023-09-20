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

## $Source: /usr/local/Repositories/fcdmf/fclib/tools/Makefile,v $
## $Revision: 1.57 $
## $Date: 2007/03/07 21:31:40 $
## 
## Makefile for tools
##
include ../Makefile.include

##
## Files (in alphabetical order)
##

TOOL_SCRIPTS = \
	dyna2names.pl

TOOL_HFILES =

TOOL_SRCS = \
	analyzeSpotWelds.c \
	bb2exo.c \
	edgeStats.c \
	extents.c \
	fc2ensight.c \
	fcconvert.c \
	fcdump.c \
	fcswapAssoc.c \
	gaplines.c \
	geomExtract.c \
	loadHistory.c \
	maxStress.c \
	screwBreaks.c \
	seriesAnalyzer.c \
	sierradump.c \
	tearCompare.c \
	tears.c \
	threshBoundBox.c \
	varNormalize.c \
	varSmooth.c \
	varStats.c \
	volumeSegments.c \
	windowVars.c 
#	testrd.c \
#	testwt.c


# Wrappers - will need to be munged upon install
TOOL_WRAPPERS_SRCS = fc2exo.in fcinfo.in
TOOL_WRAPPERS = $(TOOL_WRAPPERS_SRCS:.in=)

# Stuff that has dependencies on external IO libraries
TOOL_EXO_SRCS = exodump.c
TOOL_IO_SRCS = ${TOOL_EXO_SRCS} 
 
# all files for release target
ALL_FILES = Makefile $(TOOL_SCRIPTS) $(TOOL_HFILES) $(TOOL_SRCS) \
               ${TOOL_WRAPPERS_SRCS} ${TOOL_IO_SRCS}

# files for this build
TOOL_OBJS = $(TOOL_SRCS:.c=.o)
TOOL_EXES = $(TOOL_SRCS:.c=)
ifeq (${HAVE_EXODUS},yes)
TOOL_OBJS += $(TOOL_EXO_SRCS:.c=.o)
TOOL_EXES += $(TOOL_EXO_SRCS:.c=)
endif


##
## Targets
##

first : tools

tools : $(TOOL_EXES) $(TOOL_WRAPPERS)

# tools are dependent on the fc library
# this will update the local excutables if the library is new,
# but will not automatically update the library if it needs it
$(TOOL_EXES) : ../$(FCLIB_PATH)/$(FCLIB)

clean : 
	$(RM) $(TOOL_EXES) $(TOOL_WRAPPERS) $(GCOV_FILES)

# generic compiling directive for C code
% : %.c
	$(CC) -o $@ $(CFLAGS) -I../$(FCLIB_PATH) $(CPPFLAGS) $@.c \
		-L../$(FCLIB_PATH) $(LDFLAGS) $(FCLDLIB) $(LDLIBS)

# generic compiling directive for wrappers
% : %.in
	rm -f $@
	sed -e "s@FC_VERSION@${FC_VERSION}@g" \
	    -e "s@BUILD_DATE@`date`@g" \
	    -e "s@FC_BIN@`pwd`@g" \
	    -e "s@FC_CFLAGS@-I`pwd`/../modules@g" \
	    -e "s@FC_LDFLAGS@-L`pwd`/../modules@g" \
	    -e "s@FC_LIBS@${FCLDLIB}@g" \
	    -e "s@CFLAGS@${CPPFLAGS}@g" \
	    -e "s@LDFLAGS@${LDFLAGS}@g" \
	    -e "s@LIBS@${LDLIBS}@g" \
	     $@.in > $@
	chmod ugo+x $@

install : $(TOOL_EXES)
	cp -f $(TOOL_EXES) $(TOOL_SCRIPTS) $(INSTALL_DIR)/bin
	-for i in $(TOOL_WRAPPERS); do \
	  rm -f $(INSTALL_DIR)/bin/$$i; \
	  sed -e "s@FC_VERSION@${FC_VERSION}@g" \
	      -e "s@BUILD_DATE@`date`@g" \
	      -e "s@FC_BIN@${INSTALL_DIR}/bin@g" \
	      -e "s@FC_CFLAGS@-I${INSTALL_DIR}/include@g" \
	      -e "s@FC_LDFLAGS@-L${INSTALL_DIR}/lib@g" \
	      -e "s@FC_LIBS@${FCLDLIB}@g" \
	      -e "s@CFLAGS@${CPPFLAGS}@g" \
	      -e "s@LDFLAGS@${LDFLAGS}@g" \
	      -e "s@LIBS@${LDLIBS}@g" \
	      $$i.in > $(INSTALL_DIR)/bin/$$i; \
	  chmod ugo+x $(INSTALL_DIR)/bin/$$i; \
	done
	-for i in $(TOOL_EXES) $(TOOL_SCRIPTS) $(TOOL_WRAPPERS); do \
	   chmod go=rx $(INSTALL_DIR)/bin/$$i; \
	done

uninstall : 
	-for i in $(TOOL_EXES) $(TOOL_SCRIPTS) $(TOOL_WRAPPERS); do \
	   rm -f $(INSTALL_DIR)/bin/$$i; \
	done

release :
	-for i in $(ALL_FILES); do \
	  (cd ..; tar --append --file=$(FC_RELEASE).tar \
	    $(FC_RELEASE)/tools/$$i) \
	done

.PHONY : clean install uninstall release
