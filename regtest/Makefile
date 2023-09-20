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

## $Source: /usr/local/Repositories/fcdmf/fclib/regtest/Makefile,v $
## $Revision: 1.109 $
## $Date: 2007/05/10 19:09:37 $
## 
## Makefile for regression testing
##
include ../Makefile.include

##
## Files (in alphabetical order)
##

TEST_HFILES =
TEST_SRCS = \
        errortesthelper.c \
	printtest.c
TEST_OBJS = $(TEST_SRCS:.c=.o)
TEST_EXES = $(TEST_SRCS:.c=)

TEST_SCRIPTS = \
	errortest.pl
TEST_SUPPORT = \
	tests.out.worked \
	tests2.out.worked

ALL_FILES = Makefile $(TEST_HFILES) $(TEST_SRCS) $(TEST_SCRIPTS) \
	$(TEST_SUPPORT)

##
## Targets
##

tests : $(TEST_EXES)
$(TEST_EXES) : $(TEST_HFILES)

# Dates don't matter
DATEREGEXP='\(Mon\|Tue\|Wed\|Thu\|Fri\|Sat\|Sun\).*20[0-9][0-9]'
# A test in varnormalize has different round off on diff machines
# FIX? Is there a better way than just ignoring this line?
VARNORMALIZE='^        ave = 0.24009\(3\|4\) +/- 0.144165'
#VARNORMALIZE='         ave = 0.24009\(3\|4\) +/- 0.144165'
regtest :  tests
	rm -f tests.out
	(make runtests) >> tests.out 2>&1
	diff -I "^make" -I $(DATEREGEXP) -I $(VARNORMALIZE) \
		tests.out tests.out.worked


# additional tests that we don't want to run all of the time
regtest2 : tests
	rm -f tests2.out
	(make runtests2) >> tests2.out 2>&1
	diff -I "^make" -I $(DATEREGEXP) -I $(VARNORMALIZE) \
		tests2.out tests2.out.worked



# tools are dependent on the fc library
# this will update the local excutables if the library is new,
# but will not automatically update the library if it needs it
$(TEST_EXES) : ../$(FCLIB_PATH)/$(FCLIB)

clean : 
	$(RM)  $(TEST_EXES) *.out *.fcx *.ex2 *.g *.bb *.names \
		*.dat *.gnu *.eps *.script $(GCOV_FILES)

distclean: clean
	$(RM) d3plot d3plot??

release :
	-for i in $(ALL_FILES); do \
	  (cd ..; tar --append --file=$(FC_RELEASE).tar \
	    $(FC_RELEASE)/regtest/$$i) \
	done

# the regression tests depends on having built the tools, and will fail
# if they are not avaiable, but will not automatically update them
TOOLS_FOR_TESTS = ../tools/analyzeSpotWelds \
                  ../tools/fcdump \
	          ../tools/sierradump
# each set of commands fora test are grouped after an explanatory echo
# test are in order roughtly of dependency and complexity
runtests : $(TEST_EXES) $(TOOLS_FOR_TESTS)
	echo "Starting regression tests"
	echo "Testing fc_exitIfError() & fc_exitIfErrorPrintf() ..."
	./errortest.pl
	echo "Testing various fc_printXYZ() functions ..."
	./printtest
	echo "Testing fcdump ..."
	../tools/fcdump -A -v ../data/gen_multivar_seq.ex2
	../tools/fcdump -A -v -f "vertex ids" -s variable ../data/gen_multivar_seq.ex2
# Test the filtering on fcdump
	echo "fcdump filter: Two Vars, all steps ==> should be (0:10) and (0:10)"
	../tools/fcdump -A -s seqvar -f "/tri mesh/" ../data/gen_multivar_seq.ex2 | grep Step
	echo "fcdump filter: One Var, all steps ==> should be steps (0:10)"
	../tools/fcdump -A -s seqvar -f "/tri mesh/temperature per vertex" ../data/gen_multivar_seq.ex2 | grep Step
	echo "fcdump filter: One Var, Step [0]"
	../tools/fcdump -A -s seqvar -f "/tri mesh/temperature per vertex[0]" ../data/gen_multivar_seq.ex2 | grep Step
	echo "fcdump filter: One Var, Steps [2:3] ==> should be steps (2:3)"
	../tools/fcdump -A -s seqvar -f "/tri mesh/temperature per vertex[2:3]" ../data/gen_multivar_seq.ex2 | grep Step
	echo "fcdump filter: One Var, Last step [last] ==> should be step (10)"
	../tools/fcdump -A -s seqvar -f "/tri mesh/temperature per vertex[last]" ../data/gen_multivar_seq.ex2 | grep Step
	echo "fcdump filter: One Var, Reverse order [last:first] ==> should be steps (10:0)"
	../tools/fcdump -A -s seqvar -f "/tri mesh/temperature per vertex[last:first]" ../data/gen_multivar_seq.ex2 | grep Step
	echo "fcdump filter: One Var, All using [:last] ==> should be steps (0:10)"
	../tools/fcdump -A -s seqvar -f "/tri mesh/temperature per vertex[:last]" ../data/gen_multivar_seq.ex2 | grep Step
	echo "fcdump filter: One Var, All using [:]  ==> should be steps (0:10)"
	../tools/fcdump -A -s seqvar -f "/tri mesh/temperature per vertex[:]" ../data/gen_multivar_seq.ex2 | grep Step
	echo "fcdump filter: One Var, Last step [last] w/ all data ==> step (10)"
	../tools/fcdump -A -s seqvar -f "/tri mesh/temperature per vertex[last]" ../data/gen_multivar_seq.ex2
# Varstats
	echo "Testing varStats ..."
	../tools/varStats ../data/gen_small_tri.ex2
	../tools/varStats ../data/gen_small_tri.ex2 velocity
	../tools/varStats ../data/gen_multivar_seq.ex2
	../tools/varStats ../data/gen_multivar_seq.ex2 "temperature per vertex"
	echo "Testing edgeStats ..."
	../tools/edgeStats ../data/gen_multimesh.ex2
	../tools/edgeStats ../data/gen_multimesh.ex2 "hex mesh"
	../tools/edgeStats -d displacement ../data/gen_gaussians.ex2
	echo "Testing extents ..."
	../tools/extents ../data/gen_multimesh.ex2
	../tools/extents ../data/gen_multimesh.ex2 "hex mesh"
	../tools/extents -d "displacement per vertex" ../data/gen_small_multivar_seq.ex2
	echo "Testing sierra dump ..."
	../tools/sierradump ../data/presto-example.i
	echo "Testing volumeSegments ..."
	../tools/volumeSegments ../data/gen_multivar_seq.ex2 "temperature per vertex" ">" 100 1
	../tools/volumeSegments ../data/gen_hex_seq.ex2 "temperature per vertex" ">" 150 1
	../tools/volumeSegments -d "displacement per vertex" ../data/gen_hex_seq.ex2 "temperature per vertex" ">" 150 1
	echo "Testing analyzeSpotWelds ..."
	../tools/analyzeSpotWelds -i ../data/spotweldmodel.i -s ../data/gen_spotweldmodel.ex2 -o gen_spotweldmodel
	cat gen_spotweldmodel.dat 
	cat gen_spotweldmodel.gnu
	cat gen_spotweldmodel.script
	echo "Testing threshBoundBox ..."
	../tools/threshBoundBox  ../data/gen_hex.ex2 temperature ">" 100
	../tools/threshBoundBox ../data/gen_hex_seq.ex2 "temperature per vertex" ">" 150 
	../tools/threshBoundBox -c ../data/gen_hex_seq.ex2 "temperature per vertex" ">" 150 
	../tools/threshBoundBox -r 4 ../data/gen_hex_seq.ex2 "temperature per vertex" ">" 150 
	../tools/threshBoundBox -d displacement ../data/gen_gaussians.ex2 temperature ">" 30
	echo "Testing varNormalize ..."
	../tools/varStats ../data/gen_multivar_seq.ex2 'element ids'
	../tools/varNormalize ../data/gen_multivar_seq.ex2 'element ids' temp.ex2
	../tools/varStats temp.ex2 'element ids_nmz'
	../tools/varStats ../data/gen_multivar_seq.ex2 'temperature per element'
	../tools/varNormalize ../data/gen_multivar_seq.ex2 'temperature per element' temp.ex2
	../tools/varStats temp.ex2 'temperature per element_nmz'
	../tools/varStats ../data/gen_multivar_seq.ex2 'vertex ids'
	../tools/varNormalize ../data/gen_multivar_seq.ex2 'vertex ids' temp.ex2
	../tools/varStats temp.ex2 'vertex ids_nmz'
	../tools/varNormalize --max 2 --sparse ../data/gen_multivar_seq.ex2 'vertex ids' temp.ex2
	../tools/varStats temp.ex2 'vertex ids_nmz'
	../tools/fcdump temp.ex2
	../tools/varNormalize --permesh ../data/gen_multivar_seq.ex2 'vertex ids' temp.ex2
	../tools/varStats temp.ex2 'vertex ids_nmz'
	../tools/varNormalize --permesh -f ../data/tri_hex_maxvals.txt ../data/gen_multivar_seq.ex2 'vertex ids' temp.ex2
	../tools/varStats temp.ex2 'vertex ids_nmz'
	../tools/varStats ../data/gen_multivar_seq.ex2 'temperature per vertex'
	../tools/varNormalize ../data/gen_multivar_seq.ex2 'temperature per vertex' temp.ex2
	../tools/varStats temp.ex2 'temperature per vertex_nmz'
	../tools/varNormalize --max 2 --sparse ../data/gen_multivar_seq.ex2 'temperature per vertex' temp.ex2
	../tools/varStats temp.ex2 'temperature per vertex_nmz'
	../tools/fcdump temp.ex2
	../tools/varNormalize --permesh ../data/gen_multivar_seq.ex2 'temperature per vertex' temp.ex2
	../tools/varStats temp.ex2 'temperature per vertex_nmz'
	../tools/varNormalize --permesh -f ../data/tri_hex_maxvals.txt ../data/gen_multivar_seq.ex2 'temperature per vertex' temp.ex2
	../tools/varStats temp.ex2 'temperature per vertex_nmz'
	../tools/varNormalize --perstep ../data/gen_multivar_seq.ex2 'temperature per vertex' temp.ex2
	../tools/varStats temp.ex2 'temperature per vertex_nmz'
	../tools/varNormalize --permesh --perstep ../data/gen_multivar_seq.ex2 'temperature per vertex' temp.ex2
	../tools/varStats temp.ex2 'temperature per vertex_nmz'
	../tools/varNormalize --permesh --perstep ../data/StressTests/Tension/d3plot 'von mises stress' temp.ex2
	../tools/varStats temp.ex2 'von mises stress'
	../tools/varStats temp.ex2 'von mises stress_nmz'
	../tools/varNormalize --permesh --perstep ../data/StressTests/Tension/d3plot 'pressure' temp.ex2
	../tools/varStats temp.ex2 'pressure'
	../tools/varStats temp.ex2 'pressure_nmz'
	echo "Testing varSmooth ..."
	../tools/varStats ../data/gen_small_tri.ex2 temperature
	../tools/varSmooth ../data/gen_small_tri.ex2 temperature 1 temp.ex2
	../tools/fcdump -a temp.ex2
	../tools/varStats temp.ex2 temperature_smoothed
	../tools/varSmooth --replace ../data/gen_small_tri.ex2 temperature 1 temp.ex2
	../tools/fcdump temp.ex2
	../tools/varStats temp.ex2 temperature
	../tools/varSmooth -d displacement ../data/gen_small_tri.ex2 temperature 1 temp.ex2
	../tools/fcdump temp.ex2
	../tools/varStats temp.ex2 temperature_smoothed
	../tools/varStats ../data/gen_small_multivar_seq.ex2 "temperature per vertex"
	../tools/varSmooth ../data/gen_small_multivar_seq.ex2 "temperature per vertex" 1 temp.ex2
	../tools/fcdump temp.ex2
	../tools/varStats temp.ex2 "temperature per vertex_smoothed"
	../tools/varSmooth --replace ../data/gen_small_multivar_seq.ex2 "temperature per vertex" 1 temp.ex2
	../tools/fcdump temp.ex2
	../tools/varStats temp.ex2 "temperature per vertex"
	../tools/varSmooth -d "displacement per vertex" ../data/gen_small_multivar_seq.ex2 "temperature per vertex" 1 temp.ex2
	../tools/fcdump temp.ex2
	../tools/varStats temp.ex2 "temperature per vertex_smoothed"
	../tools/varStats ../data/gen_small_multivar_seq.ex2 "vertex ids"
	../tools/varSmooth -d "displacement per vertex" ../data/gen_small_multivar_seq.ex2 "vertex ids" 1 temp.ex2
	../tools/fcdump temp.ex2
	../tools/varStats temp.ex2 "vertex ids_smoothed"
	../tools/varSmooth --replace -d "displacement per vertex" ../data/gen_small_multivar_seq.ex2 "vertex ids" 1 temp.ex2
	../tools/fcdump temp.ex2
	../tools/varStats temp.ex2 "vertex ids"
	echo "Testing loadHistory ..."
	../tools/loadHistory -o gen_loadHistScal ../data/gen_gaussians.ex2 grid "bottom edge" temperature
	cat gen_loadHistScal.dat
	cat gen_loadHistScal.gnu
	../tools/loadHistory -o gen_loadHistVect ../data/gen_gaussians.ex2 grid "bottom edge" displacement
	cat gen_loadHistVect.dat
	cat gen_loadHistVect.gnu
	echo "Testing examples ..."
	../examples/ex1_intro
	../examples/ex2_bigDataAccess
	../examples/ex3_createCharactrztn
	echo "Testing fc2ensight ..."
	../tools/fc2ensight ../data/gen_multivar_seq.ex2
	cat gen_multivar_seq_ens/gen_multivar_seq.case
	cat gen_multivar_seq_ens/gen_multivar_seq.geo
	cat gen_multivar_seq_ens/gen_multivar_seq.temperature_per_vertex_02
	cat gen_multivar_seq_ens/gen_multivar_seq.temperature_per_element_02
	rm -rf gen_multivar_seq_ens/
	echo "Testing fcconvert ..."
	rm -f temp.ex2
	../tools/fcconvert -o temp.ex2 -t EXO ../data/gen_gaussians.ex2
	../tools/exodump temp.ex2
	../tools/exodump -a temp.ex2
	rm -f temp.ex2
	echo "Testing exodump ..."
	rm -f temp.ex2
	../tools/fcconvert -o temp.ex2 -t EXO ../data/gen_multivar_seq.ex2
	../tools/exodump temp.ex2
	../tools/exodump -a temp.ex2
	rm -f temp.ex2
	echo "Testing bb2exo ..."
	rm -rf bbtest.ex2 bbtest.names
	../tools/bb2exo ../data/bbtest.bb; ../tools/fcdump bbtest.ex2
	echo "Testing tears ..."
	rm -f *tear*.bb
	../tools/tears -n 1 damage ">" 50 ../data/gen_gaussians.ex2
	cat tear-regions.bb
	cat tear-exposeds.bb
	cat tears.fcx
	mv tears.fcx tears1.fcx
	rm -f *tear*.bb
	../tools/tears -n 2 damage ">" 50 damage "<" -30 ../data/gen_gaussians.ex2
	cat tear-regions.bb
	cat tear-exposeds.bb
	cat tears.fcx
	mv tears.fcx tears2.fcx
	rm -f *tear*.bb
	../tools/tears -d displacement -n 1 damage ">" 50 ../data/gen_gaussians.ex2
	cat tear-regions.bb
	cat tear-exposeds.bb
	cat displ-tear-exposeds.bb
	cat tears.fcx
	mv tears.fcx tears3.fcx
	rm -f *tear*.bb
	../tools/tears -w 60 -n 2 damage ">" 60 damage ">" 50 ../data/gen_gaussians.ex2
	cat tear-regions.bb
	cat tear-exposeds.bb
	cat tears.fcx
	rm -rf *tear*.bb
	../tools/tears -t 0 -n 1 damage ">" 50 ../data/gen_gaussians.ex2
	cat tear-regions.bb
	cat tear-exposeds.bb
	cat tears.fcx
	echo "Testing tearCompare ..."
	rm -f temp.out
	../tools/tears -n 1 damage ">" 1 ../data/gen_gaussians.ex2 > temp.out
	mv tears.fcx tears4.fcx
	../tools/tearCompare tears1.fcx tears3.fcx
	../tools/tearCompare tears1.fcx tears2.fcx
	../tools/tearCompare tears2.fcx tears1.fcx
	../tools/tearCompare tears1.fcx tears4.fcx
	../tools/tearCompare tears4.fcx tears1.fcx
	../tools/tearCompare -t Tear0 tears1.fcx tears2.fcx
	../tools/tearCompare -t Tear0 tears4.fcx tears1.fcx
	../tools/tearCompare -t Tear1 tears2.fcx tears1.fcx
	echo "Testing dyna parsers ..."
	rm -f *.names
	../tools/dyna2names.pl ../data/lsdyna-names.k
	cat d3plot.names
	echo "Testing screwBreak ..."
	../tools/screwBreaks -n 1 damage ">" 0 ../data/gen_screws.ex2
	echo "Testing gaplines ..."
	../tools/gaplines ../data/gen_gaps.ex2 mesh0 mesh1 displacement 0.001
	../tools/fcdump gaplines.ex2
	../tools/gaplines ../data/gen_gaps.ex2 mesh0 mesh1 displacement 0
	../tools/gaplines ../data/gen_gaps.ex2 mesh0 mesh2 displacement 0.001
	../tools/gaplines ../data/gen_gaps.ex2 mesh0 mesh2 displacement 0
	../tools/gaplines ../data/gen_gaps.ex2 mesh0 mesh3 displacement 0.001
	../tools/gaplines ../data/gen_gaps.ex2 mesh0 mesh3 displacement 0
	../tools/gaplines ../data/gen_gaps.ex2 mesh0 mesh4 displacement 0.001
	../tools/gaplines -n 1 isDead "=" 1 ../data/gen_gaps.ex2 mesh0 mesh1 displacement 0
	echo "Testing maxStress ..."
	../tools/maxStress ../data/StressTests/Compression/d3plot
	../tools/maxStress ../data/StressTests/Tension/d3plot
	echo "Testing geomExtract ..."
	../tools/geomExtract -d displacement ../data/gen_gaussians.ex2
	../tools/fcdump -a extracted.g
	echo "Testing swapAssociation ..."
	../tools/fcswapAssoc -d e2v -r _mod -o swapAssoc1.ex2 ../data/gen_multivar_seq.ex2
	../tools/fcdump swapAssoc1.ex2
	../tools/varStats swapAssoc1.ex2 "element ids_mod"
	../tools/fcswapAssoc -d v2e -r _mod -o swapAssoc2.ex2 ../data/gen_multivar_seq.ex2
	../tools/fcdump swapAssoc2.ex2
	../tools/varStats swapAssoc2.ex2 "vertex ids_mod"
	echo "Testing subsetter ..."
	../tools/subsetter -M -a ../data/gen_hex.ex2 "temperature" ">" 80 1
	../tools/fcdump subsetter_geomonly.ex2
	../tools/fcdump subsetter_subsets.ex2
	../tools/varStats subsetter_subsets.ex2 "temperature"
	../tools/varStats subsetter_subsets.ex2 "FC_ELEMENTMAP"
	../tools/varStats subsetter_subsets.ex2 "FC_VERTEXMAP"
	rm -f subsetter_geomonly.ex2
	rm -f subsetter_subsets.ex2
	../tools/subsetter -a ../data/gen_hex.ex2 "temperature" ">" 80 1
	../tools/fcdump subsetter_subsets.ex2
	rm -f subsetter_subsets.ex2
	../tools/subsetter -a ../data/gen_multivar_seq.ex2 "temperature per vertex" ">" 300 1
	../tools/fcdump subsetter_subsets.ex2
	rm -f subsetter_subsets.ex2
	echo "Testing subsetter and regionreassembler pair"
	rm -f subsetter_subsets.ex2
	../tools/subsetter -M -a ../data/gen_multivar_seq.ex2 "temperature per vertex" ">" 300 1
	../tools/fcdump subsetter_geomonly.ex2
	../tools/fcdump subsetter_subsets.ex2
	../tools/varStats subsetter_subsets.ex2 "temperature per vertex"
	../tools/varStats subsetter_subsets.ex2 "FC_ELEMENTMAP"
	../tools/varStats subsetter_subsets.ex2 "FC_VERTEXMAP"
	../tools/varStats subsetter_subsets.ex2 "element ids"
	../tools/varStats subsetter_subsets.ex2 "vertex ids"
	../tools/regionreassembler subsetter_geomonly.ex2 subsetter_subsets.ex2
	../tools/fcdump reassembled.ex2
	../tools/varStats subsetter_subsets.ex2 "temperature per vertex"
	../tools/varStats subsetter_subsets.ex2 "element ids"
	../tools/varStats subsetter_subsets.ex2 "vertex ids"
	rm -f subsetter_geomonly.ex2
	rm -f subsetter_subsets.ex2
	rm -f reassembled.ex2
	echo "Finished regresion tests"

# These are tests that depend on data in fclib-data module
# Make links to lsdyna data to control access to .names files
# (?FIX? Make a target for the data directory?)
runtests2 : $(TEST_EXES) $(TOOLS_FOR_TESTS)
	echo "Starting regression2 tests"
	echo "Setup ..."
	mkdir -p run30-EMMI-screws
	ln -sf ../../$(FCLIB_DATA_DIR)/lsdynaData/run30-EMMI-screws/d3plot run30-EMMI-screws/d3plot
	ln -sf ../../$(FCLIB_DATA_DIR)/lsdynaData/run30-EMMI-screws/d3plot01 run30-EMMI-screws/d3plot01
	ln -sf ../../$(FCLIB_DATA_DIR)/lsdynaData/run30-EMMI-screws/d3plot02 run30-EMMI-screws/d3plot02
	ln -sf ../../$(FCLIB_DATA_DIR)/lsdynaData/run30-EMMI-screws/d3plot03 run30-EMMI-screws/d3plot03
	ln -sf ../../$(FCLIB_DATA_DIR)/lsdynaData/run30-EMMI-screws/d3plot04 run30-EMMI-screws/d3plot04
	ln -sf ../../$(FCLIB_DATA_DIR)/lsdynaData/run30-EMMI-screws/d3plot05 run30-EMMI-screws/d3plot05
	ln -sf ../../$(FCLIB_DATA_DIR)/lsdynaData/run30-EMMI-screws/d3plot06 run30-EMMI-screws/d3plot06
	ln -sf ../../$(FCLIB_DATA_DIR)/lsdynaData/run30-EMMI-screws/d3plot07 run30-EMMI-screws/d3plot07
	echo "Testing exodump ..."
	../tools/exodump ../$(FCLIB_DATA_DIR)/exodusData/cancrush/can.ex2
	../tools/exodump ../$(FCLIB_DATA_DIR)/exodusData/spotweld-old/simple-assy-new-bolts.e
	echo "Testing LSDyna reader ..."
	../tools/fcdump -A -v ../$(FCLIB_DATA_DIR)/lsdynaData/weird/d3plot
	rm -f run30-EMMI-screws/d3plot.names
	../tools/dyna2names.pl 	../$(FCLIB_DATA_DIR)/lsdynaData/run30-EMMI-screws/top-emmi.k
	cat d3plot.names
	../tools/fcdump -v run30-EMMI-screws/d3plot | grep "^Mesh"
	mv d3plot.names run30-EMMI-screws
	../tools/fcdump -v run30-EMMI-screws/d3plot | grep "^Mesh"
	echo "Testing analyzeSpotWelds on real data"
	../tools/analyzeSpotWelds -o spotwelds -i ../$(FCLIB_DATA_DIR)/exodusData/spotweld-old/simple-assy-new-bolts.i -s ../$(FCLIB_DATA_DIR)/exodusData/spotweld-old/simple-assy-new-bolts.ex2
	cat spotwelds.dat
	cat spotwelds.gnu
	cat spotwelds.script
	echo "Testing tears on real data"
	rm -f *tear*.bb
	../tools/tears -d displacement run30-EMMI-screws/d3plot
	cat tear-regions.bb
	cat tear-exposeds.bb
	cat displ-tear-exposeds.bb
	rm -f *tear*.bb
	../tools/tears -ss run30-EMMI-screws/d3plot weld-top-layers foam-for-components.foam-for-components
	rm -f *tear*.bb
	../tools/tears -d displacement -m 2 run30-EMMI-screws/d3plot
	cat tear-regions.bb
	cat tear-exposeds.bb
	cat displ-tear-exposeds.bb
	rm -f *tear*.bb
	../tools/tears -w .15 -n 2 elem_death "==" 0 elem_var_8 ">" .1 	run30-EMMI-screws/d3plot weld-top-layers weld-mid-layers
	cat tear-regions.bb
	cat tear-exposeds.bb
	echo "Testing screwBreaks on real data ..."
	../tools/screwBreaks run30-EMMI-screws/d3plot
	echo "Cleaning up ..."
	rm -r run30-EMMI-screws
	echo "Finished regresion2 tests"

# generic compiling directives
%.o : %.c
	$(CC) -c -o $@ $(CFLAGS) -I../$(FCLIB_PATH) $(CPPFLAGS) $<
% : %.c
	$(CC) -o $@ $(CFLAGS) -I../$(FCLIB_PATH) $(CPPFLAGS) $@.c \
		 -L../$(FCLIB_PATH) $(LDFLAGS) $(FCLDLIB) $(LDLIBS)

.PHONY : tests regtest regtest2 runtests runtests2 clean distclean release
