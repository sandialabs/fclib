#!/bin/sh

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

## $Source: /usr/local/Repositories/fcdmf/fclib/data/makedata.sh,v $
## $Revision: 1.27 $
## $Date: 2006/08/30 19:20:00 $
##
## fcdmf script for generating test data
## most generated files are prefixed with 'gen_' for easy cleanup
## DEC-5-2002  W Koegler  Created.
## OCT-X-2007 ACG removing saf in favor of exodus and all unused cases.

# data_generator makes displacement and stress data for small & tiny
./data_generator displ -o gen_small_tri_displ.data -v small_tri.vert
./data_generator stress -o gen_small_tri_vstress.data -v small_tri.vert 
./data_generator stress -o gen_small_tri_estress.data -v small_tri.vert \
    -e small_tri.elem
./data_generator displ -o gen_tiny_tri_displ.data -v tiny_tri.vert

# little single mesh datasets with provided vert,elem & data files
./dataset_generator -m 1 -o gen_small_tri.ex2 -e 1 small_tri.elem -v 1 \
    small_tri.vert -d 1 small_tri.data -d 1 gen_small_tri_displ.data \
    -Mv -Mst -Mi
./dataset_generator -m 1 -o gen_small_tri.ex2 -e 1 small_tri.elem -v 1 \
    small_tri.vert -d 1 small_tri.data -d 1 gen_small_tri_displ.data \
    -Mv -Mst -Mi
./dataset_generator -m 1 -o gen_tiny_tri.ex2 -e 1 tiny_tri.elem -v 1 \
    tiny_tri.vert -d 1 tiny_tri.data -Mst

# little dataset with variety of seq vars and stuff
./dataset_generator -m 2 -o gen_small_multivar_seq.ex2 \
    -v 1 small_tri.vert -e 1 small_tri.elem -d 1 small_tri.data \
    -d 1 gen_small_tri_displ.data \
    -v 2 tiny_tri.vert  -e 2 tiny_tri.elem  -d 2 tiny_tri.data \
    -d 2 gen_tiny_tri_displ.data \
    -Mt -Mi

# mesh_generator makes vertex & element files for each celltype (8)
./mesh_generator

# data_generator makes dipslacement
./data_generator displ -o gen_point_line_tri_quad_displ.data \
    -v gen_point_line_tri_quad.vert
./data_generator displ -o gen_tet_prism_hex_displ.data \
    -v gen_tet_prism_hex.vert
./data_generator displ -o gen_pyramid_displ.data -v gen_pyramid.vert
# data_generator makes temp data with 2 gaussian peaks
# vertex associated data is generated for all of the meshes
./data_generator temp2D -o gen_point_line_tri_quad.data \
    -v gen_point_line_tri_quad.vert
./data_generator temp3D -o gen_tet_prism_hex.data -v gen_tet_prism_hex.vert
./data_generator temp3D -o gen_pyramid.data -v gen_pyramid.vert
# element associated data for a few meshes, appended with '_perElem.data'
./data_generator temp2D -o gen_tri_perElem.data \
    -v gen_point_line_tri_quad.vert -e gen_tri.elem
./data_generator temp3D -o gen_hex_perElem.data \
    -v gen_tet_prism_hex.vert -e gen_hex.elem

# single mesh datasets using generated meshes and generated data and
# default subsets
./dataset_generator -m 1 -o gen_point.ex2 -v 1 gen_point_line_tri_quad.vert \
    -e 1 gen_point.elem -d 1 gen_point_line_tri_quad.data
./dataset_generator -m 1 -o gen_line.ex2 -v 1 gen_point_line_tri_quad.vert \
    -e 1 gen_line.elem -d 1 gen_point_line_tri_quad.data
./dataset_generator -m 1 -o gen_tri.ex2 -v 1 gen_point_line_tri_quad.vert \
    -e 1 gen_tri.elem -d 1 gen_point_line_tri_quad.data
./dataset_generator -m 1 -o gen_quad.ex2 -v 1 gen_point_line_tri_quad.vert \
    -e 1 gen_quad.elem -d 1 gen_point_line_tri_quad.data
./dataset_generator -m 1 -o gen_tet.ex2 -v 1 gen_tet_prism_hex.vert \
    -e 1 gen_tet.elem -d 1 gen_tet_prism_hex.data
./dataset_generator -m 1 -o gen_prism.ex2 -v 1 gen_tet_prism_hex.vert \
    -e 1 gen_prism.elem -d 1 gen_tet_prism_hex.data
./dataset_generator -m 1 -o gen_pyramid.ex2 -v 1 gen_pyramid.vert \
    -e 1 gen_pyramid.elem -d 1 gen_pyramid.data
./dataset_generator -m 1 -o gen_hex.ex2 -v 1 gen_tet_prism_hex.vert \
    -e 1 gen_hex.elem -d 1 gen_tet_prism_hex.data -Msub


# multiple mesh dataset with full variety of meshes and default subsets
./dataset_generator -m 8 -o gen_multimesh.ex2 \
    -v 1 gen_point_line_tri_quad.vert   -e 1 gen_point.elem \
	-d 1 gen_point_line_tri_quad.data \
    -v 2 gen_point_line_tri_quad.vert   -e 2 gen_line.elem \
	-d 2 gen_point_line_tri_quad.data \
    -v 3 gen_point_line_tri_quad.vert   -e 3 gen_tri.elem \
	-d 3 gen_point_line_tri_quad.data \
    -v 4 gen_point_line_tri_quad.vert   -e 4 gen_quad.elem \
	-d 4 gen_point_line_tri_quad.data \
    -v 5 gen_tet_prism_hex.vert         -e 5 gen_tet.elem \
	-d 5 gen_tet_prism_hex.data \
    -v 6 gen_tet_prism_hex.vert         -e 6 gen_prism.elem \
	-d 6 gen_tet_prism_hex.data \
    -v 7 gen_pyramid.vert               -e 7 gen_pyramid.elem \
	-d 7 gen_pyramid.data \
    -v 8 gen_tet_prism_hex.vert         -e 8 gen_hex.elem \
	-d 8 gen_tet_prism_hex.data \
    -Msub

# single mesh datasets with time series and subsets
./dataset_generator -m 1 -o gen_hex_seq.ex2 \
    -v 1 gen_tet_prism_hex.vert -e 1 gen_hex.elem \
    -d 1 gen_tet_prism_hex_displ.data \
    -d 1 gen_tet_prism_hex.data -Mt -Msub

# double mesh data set with variety of seq vars (2 on 1st mesh, 1 on 2nd)
#  and some non-seq vars
./dataset_generator -m 2 -o gen_multivar_seq.ex2 \
    -v 1 gen_point_line_tri_quad.vert -e 1 gen_tri.elem \
    -d 1 gen_point_line_tri_quad.data -d 1 gen_tri_perElem.data \
    -v 2 gen_tet_prism_hex.vert -e 2 gen_hex.elem  \
    -d 2 gen_tet_prism_hex.data -Mt -Mi -Msub
