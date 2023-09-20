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

## $Source$
## $Revision$
## $Date$
## 
## Makefile for examples
##
include ../Makefile.include

##
## Files (in alphabetical order)
##

EXAMPLE_HFILES =
EXAMPLE_SRCS = \
	elemdeath_ex1.c \
	ex1_intro.c \
	ex2_bigDataAccess.c \
	ex3_createCharactrztn.c \
	seriesDemo.c 
EXAMPLE_OBJS = $(EXAMPLE_SRCS:.c=.o)
EXAMPLE_EXES = $(EXAMPLE_SRCS:.c=)
EXAMPLE_OUTPUT = seriesDemo.out seriesDemo.gnu

ALL_FILES = Makefile $(EXAMPLE_HFILES) $(EXAMPLE_SRCS)

##
## Targets
##

first : $(EXAMPLE_EXES)

examples : $(EXAMPLE_EXES)

# examples are dependent on the fc library
# this will update the local excutables if the library is new,
# but will not automatically update the library if it needs it
$(EXAMPLE_EXES) : ../$(FCLIB_PATH)/$(FCLIB)

clean : 
	$(RM) $(EXAMPLE_EXES) $(GCOV_FILES) $(EXAMPLE_OUTPUT)

# generic compiling directive
% : %.c
	$(CC) -o $@ $(CFLAGS) -I../$(FCLIB_PATH) $(CPPFLAGS) $@.c \
		-L../$(FCLIB_PATH) $(LDFLAGS) $(FCLDLIB) $(LDLIBS)

release :
	-for i in $(ALL_FILES); do \
	  (cd ..; tar --append --file=$(FC_RELEASE).tar \
	    $(FC_RELEASE)/examples/$$i) \
	done

.PHONY : clean install uninstall release