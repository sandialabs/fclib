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

## $Source: /home/Repositories/fcdmf/fclib/util/Makefile,v $
## $Revision: 1.13 $
## $Date: 2006/11/15 22:37:03 $
##
## Makefile for utilities

# 2006/11/14 - WSD - Commented this out because was having problems with
#   LDFLAGS and CPPFLAGS variables showing up when running check's configure.
#   But we don't need anything from there right now. 
#   (I tried to unexport LDFLAGS and CPPFLAGS, but that clobbered the
#    initial values in the environment)
# 2006/11/22 - WSD - DOH! We do need this, so we still have the bug but
#   it probably won't hurt anyone ...?
include ../Makefile.include

##
## Files
##

# check
# 2006/11/15 - WSD - I haven't upgraded check because I can't get the
#  newer versions (just tried 0.9.4) to pass their tests on my Mac
CHECK_TARBALL = check-0.9.0.tar.gz
CHECK_DIR = $(CHECK_TARBALL:.tar.gz=)


ALL_FILES = Makefile $(CHECK_TARBALL) 


##
## Targets
##

first : check

# 2006/11/14 - WSD - Make link before building so that if check isn't
#   building, the user can muck about and their changes won't get clobbered
check : $(CHECK_TARBALL)
	rm -rf ${CHECK_DIR} check
	tar -xzf $<
	ln -s ${CHECK_DIR} check
	(cd ${CHECK_DIR} ; ./configure; make ; make check )

clean distclean : 
	rm -rf ${CHECK_DIR} check

release :
	-for i in $(ALL_FILES); do \
	  (cd ..; tar --append --file=$(FC_RELEASE).tar \
	    $(FC_RELEASE)/util/$$i) \
	done

.PHONY : distclean release