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

## $Source: /home/Repositories/fcdmf/fclib/Makefile,v $
## $Revision: 1.190 $
## $Date: 2006/08/30 19:19:59 $
## 
## Makefile for the feature characterization library
##

## Edit 'Makefile.include' to set compiler options and paths to libraries
include Makefile.include

##
## Files
##

# doxygen
DOX_CONFIG = doxyfiles/doxygen_fc_config
DOX_SUPPORT = doxyfiles/doxy_module_defs \
              doxyfiles/doxy_page_defs
HTML_PAGES = doxyfiles/FCLib_characterizations.html \
             doxyfiles/FCLib_cvsConventions.html \
             doxyfiles/FCLib_installation.html \
             doxyfiles/FCLib_releaseNotes \
	     doxyfiles/FClibGUI.html
HTML_SUPPORT = doxyfiles/release_example_s.png \
	       doxyfiles/release_example_600.png
ALL_DOX_FILES = $(DOX_CONFIG) $(DOX_SUPPORT) $(HTML_PAGES) $(HTML_SUPPORT)

# generated files
README_FILES = INSTALL RELEASE_NOTES


##
## Targets
##

# basic = the library and it's tools
first : lib tools

# target for lazy developers
world : lib tools examples data regtest util check $(README_FILES)

lib :
	(cd modules; make lib)

tools : lib
	(cd tools; make tools)

examples : lib
	(cd examples; make examples)

data : lib
	(cd data; make data)

util :
	(cd util; make)

regtest: lib tools examples data
	(cd regtest; make regtest)

regtest2: lib tools examples data
	(cd regtest; make regtest2)

check : lib data util
	(cd unittest; make check)

check2 : lib data util
	(cd unittest; make check2)

interfaces : lib
	(cd interfaces; make)

gui : lib interfaces
	(cd gui; make)

# create plain text README files from the original html files
# (html2text src can be found in repository ascdd/ascdd-admin/util)
HTML2TEXT_CMD = html2text -style pretty -nobs -ascii
RULELINE = "=========================================================================="
INSTALL: doxyfiles/FCLib_installation.html
	rm -f INSTALL
	cat COPYRIGHT > INSTALL
	echo "" >> INSTALL
	echo ${RULELINE} >> INSTALL
	echo "" >> INSTALL
	echo "FCLib Installation Instructions" >> INSTALL
	$(HTML2TEXT_CMD) doxyfiles/FCLib_installation.html >> INSTALL
RELEASE_NOTES: doxyfiles/FCLib_releaseNotes.html
	rm -f RELEASE_NOTES
	cat COPYRIGHT > RELEASE_NOTES
	echo "" >> RELEASE_NOTES
	echo ${RULELINE} >> RELEASE_NOTES
	echo "" >> RELEASE_NOTES
	echo "FCLib Release Notes" >> RELEASE_NOTES
	$(HTML2TEXT_CMD) doxyfiles/FCLib_releaseNotes.html >> RELEASE_NOTES

# insert FClib version number in doxygen configuration file
# (Find line with PROJECT_NUMBER and replace " " with version number)
doc dox : $(README_FILES) $(DOX_CONFIG) $(DOX_SUPPORT)
	sed -e '/PROJECT_NUMBER /s/" "/"FCLib-$(FC_VERSION)"/g' \
		$(DOX_CONFIG) > doxy_temp
	rm -rf doc/html
	doxygen doxy_temp
	cp $(HTML_SUPPORT) doc/html
	rm -f doxy_temp

# documentation without code snippets for public web pages
# stuff this into a tarball
webdoc: $(README_FILES) $(DOX_CONFIG) $(DOX_SUPPORT)
	sed -e '/PROJECT_NUMBER /s/" "/"FCLib-$(FC_VERSION)"/g' \
	    -e '/SOURCE_BROWSER /s/YES/NO/g' \
	    -e '/INLINE_SOURCES /s/YES/NO/g' \
	    -e '/HTML_OUTPUT /s/html/html_no_code/g' \
		$(DOX_CONFIG) > doxy_temp
	rm -rf doc/html_no_code
	doxygen doxy_temp
	cp $(HTML_SUPPORT) doc/html_no_code
	rm -f doxy_temp
	ln -s doc/html_no_code FCLibManual
	tar -czf FCLibManual.tar.gz FCLibManual/*
	rm FCLibManual


# ranlib needed for OS X
#FIX create install directory if it doesn't exist
#FIX? keep in mind that some scripts are having paths changed w/ sed
install : lib tools
	mkdir -p $(INSTALL_DIR)/bin
	chmod go=rx $(INSTALL_DIR)/bin
	mkdir -p $(INSTALL_DIR)/include
	chmod go=rx $(INSTALL_DIR)/include
	mkdir -p $(INSTALL_DIR)/lib
	chmod go=rx $(INSTALL_DIR)/lib
	(cd modules; make install)
	(cd tools; make install)

uninstall :
	(cd modules; make uninstall)
	(cd tools; make uninstall)

clean :
	(cd data; make clean)
	(cd modules; make clean)
	(cd regtest; make clean)
	(cd tools; make clean)
	(cd examples; make clean)
	(cd unittest; make clean)

distclean : clean
	(cd util; make distclean)
	rm -rf bin include lib

reallydistclean : distclean
	rm -rf doc $(README_FILES) $(FC_RELEASE).tar.gz FCLibManual.tar.gz

# A symbolic link to the local directory is created which is used to
# essentially rename the root directory of the tarball. The files at the
# root level and for the documentation are added at creation of the tar
# file. All other files are added by 'make release' in the child directories.
release : $(README_FILES) doc
	chmod a+w Makefile.include
	rm -f $(FC_RELEASE).tar.gz
	ln -s . ${FC_RELEASE}
	tar -cvf ${FC_RELEASE}.tar \
	  $(FC_RELEASE)/Makefile \
	  $(FC_RELEASE)/Makefile.include \
	  $(FC_RELEASE)/COPYRIGHT \
          $(FC_RELEASE)/README \
	  $(FC_RELEASE)/INSTALL \
	  $(FC_RELEASE)/RELEASE_NOTES \
	  $(FC_RELEASE)/doxyfiles/doxy* \
	  $(FC_RELEASE)/doxyfiles/*.html \
	  $(FC_RELEASE)/doxyfiles/*.png \
	  $(FC_RELEASE)/doc/html
	(cd data; make release)
	(cd modules; make release)
	(cd regtest; make release)
	(cd tools; make release)
	(cd examples; make release)
	(cd unittest; make release)
	(cd util; make release)
	rm -f $(FC_RELEASE)
	gzip $(FC_RELEASE).tar


# to avoid conflicts with actual named files/directories
.PHONY : first world lib tools examples data regtest util check doc dox \
	install clean release

##
## Misc Stuff
##

# uncomment and change this to use sed to change lots of files at once
# the example changes instances of 'sad' to 'happy'
# WARNING! If you screw up the regular express, this will delete the
# contents of the the files!!! (NOTE: use '\s' instead of ' ' to match a space)
#sed :
#	-for i in $(ALL_HFILES) ${ALL_SRCS} testutils.c ; do \
#	  sed -e s/fc_sequencehandle/FC_Sequence/g $$i > sedtemp; \
#	  mv sedtemp $$i; \
#	done
