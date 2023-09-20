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

## $Source: /home/Repositories/fcdmf/fclib/data/Makefile,v $
## $Revision: 1.20 $
## $Date: 2006/09/08 23:33:42 $
## 
## Makefile for data generation
##
include ../Makefile.include

##
## Files (by type, then in alphabetical order, mostly)
##

DATA_RAW = \
	small_tri.data \
	small_tri.elem \
	small_tri.vert \
	tiny_tri.data \
	tiny_tri.elem \
	tiny_tri.vert \
	tri_hex_maxvals.txt \
	screw.elem \
	screw.vert \
	lsdyna-names.k \
	presto-example.i \
	spotweldmodel.i \
	bbtest.bb 

# These should extract to a directory of the same name
DATA_TARS = \
	StressTests.tgz
DATA_TARS_DIRS = $(DATA_TARS:.tgz=)

# Codes that generate data
DATA_HFILES =
DATA_SRCS = \
	data_generator.c \
	dataset_generator.c \
	gaps_generator.c \
	gaussians_generator.c \
	mesh_generator.c \
	screws_generator.c \
	spotweld_generator.c
DATA_OBJS = $(DATA_SRCS:.c=.o)
DATA_EXES = $(DATA_SRCS:.c=)
DATA_SCRIPTS = makedata.sh

# All files
ALL_FILES = Makefile $(DATA_RAW) $(DATA_TARS) \
            $(DATA_HFILES) $(DATA_SRCS) $(DATA_SCRIPTS)

##
## Targets
##

data : $(DATA_EXES) $(DATA_RAW) $(DATA_TARS_DIRS) $(DATA_SCRIPTS)
	./makedata.sh
	./gaps_generator
	./gaussians_generator
	./screws_generator
	./spotweld_generator
	$(RM) data
	echo "This file is a place holder to force make dependencies" > data

# generators are dependent on the fc library
# this will update the local excutables if the library is new,
# but will not automatically update the library if it needs it
$(DATA_EXES) : ../$(FCLIB_PATH)/$(FCLIB)

# generic directive for tar'd data
% : %.tgz
	tar -xzf $<
	touch $@

# most generated files start with 'gen_' for easy cleanup
clean : 
	$(RM) $(DATA_EXES) data gen_* $(GCOV_FILES)
	$(RM) -r $(DATA_TARS_DIRS)

# generic compiling directive
% : %.c
	$(CC) -o $@ $(CFLAGS) -I../$(FCLIB_PATH) $(CPPFLAGS) $@.c \
		-L../$(FCLIB_PATH) $(LDFLAGS) $(FCLDLIB) $(LDLIBS)

release :
	-for i in $(ALL_FILES) ; do \
	  (cd ..; tar --append --file=${FC_RELEASE}.tar \
	    $(FC_RELEASE)/data/$$i) \
        done

# test changes
regtest :
	(cd ..; make regtest)
check :
	(cd ..; make check)

.PHONY : clean release regtest check