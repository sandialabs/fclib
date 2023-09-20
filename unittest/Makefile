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

## $Source: /usr/local/Repositories/fcdmf/fclib/unittest/Makefile,v $
## $Revision: 1.44 $
## $Date: 2006/08/30 19:20:05 $
## 
## Makefile for unit tests (using check)
##
include ../Makefile.include

CHECK_HOME = ../util/check
CK_CPPFLAGS = ${CPPFLAGS} -I${CHECK_HOME}/src
CK_LDFLAGS = ${LDFLAGS} -L${CHECK_HOME}/src
CK_LDLIBS = ${LDLIBS} -lcheck

##
## Files (in alphabetical order)
##

CHECK_HFILES = checkall.h
CHECK_MAIN = checkall.c
CHECK_SRCS = \
	checkbase.c \
	checkcustom.c \
	checkelemdeath.c \
	checkdataset.c \
	checkfeature.c \
	checkfileio.c \
	checkgeom.c \
	checklibrary.c \
	checkmesh.c \
	checksequence.c \
	checkseries.c \
	checkshape.c \
	checkstats.c \
	checkstorage.c \
	checksubset.c \
	checktable.c \
	checkthresh.c \
	checktopo.c \
	checktrack.c \
	checkutil.c \
	checkvariable.c \
	checkvarmath.c
CHECK_SRCS2 = \
	checkfileio2.c

CHECK_OBJS = $(CHECK_SRCS:.c=.o)
CHECK_OBJS2 = $(CHECK_SRCS2:.c=.o)
CHECK_EXES = $(CHECK_MAIN:.c=)
CHECK_EXES2 = $(CHECK_EXES:=2)

ALL_FILES = Makefile $(CHECK_HFILES) $(CHECK_MAIN) $(CHECK_SRCS) $(CHECK_SRCS2)


##
## Targets
##

first : ${CHECK_EXES} ${CHECK_EXES2}

check: ${CHECK_EXES}
	./checkall

check2: datalinks ${CHECK_EXES2} 
	./checkall2

# Make links to control access to .names files
# ?FIX? This is duplicated in regtest, is there a better way?
# (I think not because regtest w/ be manipulating .names and here we don't)
datalinks:
	mkdir -p datalinks/run30-EMMI-screws
	ln -sf ../../../$(FCLIB_DATA_DIR)/lsdynaData/run30-EMMI-screws/d3plot   \
		datalinks/run30-EMMI-screws/d3plot
	ln -sf ../../../$(FCLIB_DATA_DIR)/lsdynaData/run30-EMMI-screws/d3plot01 \
		datalinks/run30-EMMI-screws/d3plot01
	ln -sf ../../../$(FCLIB_DATA_DIR)/lsdynaData/run30-EMMI-screws/d3plot02 \
		datalinks/run30-EMMI-screws/d3plot02
	ln -sf ../../../$(FCLIB_DATA_DIR)/lsdynaData/run30-EMMI-screws/d3plot03 \
		datalinks/run30-EMMI-screws/d3plot03
	ln -sf ../../../$(FCLIB_DATA_DIR)/lsdynaData/run30-EMMI-screws/d3plot04 \
		datalinks/run30-EMMI-screws/d3plot04
	ln -sf ../../../$(FCLIB_DATA_DIR)/lsdynaData/run30-EMMI-screws/d3plot05 \
		datalinks/run30-EMMI-screws/d3plot05
	ln -sf ../../../$(FCLIB_DATA_DIR)/lsdynaData/run30-EMMI-screws/d3plot06 \
		datalinks/run30-EMMI-screws/d3plot06
	ln -sf ../../../$(FCLIB_DATA_DIR)/lsdynaData/run30-EMMI-screws/d3plot07 \
		datalinks/run30-EMMI-screws/d3plot07

checkall : checkall.c $(CHECK_HFILES) $(CHECK_OBJS)
	$(CC) -o $@ $(CFLAGS) -I../$(FCLIB_PATH) $(CK_CPPFLAGS) \
		$(CHECK_MAIN) $(CHECK_OBJS) \
		-L../$(FCLIB_PATH) $(CK_LDFLAGS) $(FCLDLIB) $(CK_LDLIBS)

checkall2 : checkall.c $(CHECK_HFILES) $(CHECK_OBJS2)
	$(CC) -o $@ $(CFLAGS) -DEXTERNAL_DATA -I../$(FCLIB_PATH) $(CK_CPPFLAGS) \
		$(CHECK_MAIN) $(CHECK_OBJS2) \
		-L../$(FCLIB_PATH) $(CK_LDFLAGS) $(FCLDLIB) $(CK_LDLIBS)

# unit tests are dependent on the fc library
# this will update the local excutables if the library is new,
# but will not automatically update the library if it needs it
# (not dependent on check library since that should be static)
$(CHECK_EXES) : ../$(FCLIB_PATH)/$(FCLIB)
$(CHECK_OBJS) : ../$(FCLIB_PATH)/$(FCLIB)
$(CHECK_EXES2) : ../$(FCLIB_PATH)/$(FCLIB)
$(CHECK_OBJS) : ../$(FCLIB_PATH)/$(FCLIB)

clean : 
	$(RM) *.o $(CHECK_EXES) $(CHECK_EXES2) *.ex2 *.fcx *.bb *.names \
		$(GCOV_FILES)
	$(RM) -r datalinks

release :
	-for i in $(ALL_FILES); do \
	  (cd ..; tar --append --file=$(FC_RELEASE).tar \
	    $(FC_RELEASE)/unittest/$$i) \
	done

# generic compiling directives
%.o : %.c $(CHECK_HFILES)
	$(CC) -c -o $@ $(CFLAGS) -I../$(FCLIB_PATH) $(CK_CPPFLAGS) $<

# special compiling directives for check2 tests
checkfileio2.o : checkfileio2.c $(CHECK_HFILES)
	$(CC) -c -o $@ $(CFLAGS) -DFCLIB_DATA_DIR='"${FCLIB_DATA_DIR}"' \
                 -I../$(FCLIB_PATH) $(CK_CPPFLAGS) $<

.PHONY : first check check2 clean release
