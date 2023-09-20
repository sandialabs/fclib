/*
 * Copyright (2000) Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
 * certain rights in this software.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 *   * Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *
 *   * Redistributions in binary form must reproduce the above copyright
 *     notice, this list of conditions and the following disclaimer in the
 *     documentation and/or other materials provided with the
 *     distribution.
 *
 *   * Neither the name of Sandia nor the names of any contributors may
 *     be used to endorse or promote products derived from this software
 *     without specific prior written permission.
 *
 *   * Modified source versions must be plainly marked as such, and must
 *     not be misrepresented as being the original software.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHORS OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
 * DAMAGE.
 */

/** 
 * \file mesh_generator.c
 * \brief Make vertex and element files for variety of meshes.
 *
 * $Source: /usr/local/Repositories/fcdmf/fclib/data/mesh_generator.c,v $
 * $Revision: 1.15 $
 * $Date: 2006/09/19 00:57:58 $
 *
 * \description
 *
 *     Quick program to make all sorts of meshes
 *     points, lines, tris & quads are on a gaussian surface
 *     hex, tet, prism & pryamid are in blocks (boring)
 *
 *     The range is from (0,0,0) to (A,B,C)
 *
 *     1st make the gaussian surface w/ quad mesh then divide 
 *         up for tris, lines & points
 *     2nd make the hex block and then divide up for tets, prisms & pyramids
 *
 *     Quad & hex elements are ordered with the x dimension increasing fastes
 *     and z slowest. i.e.
 *
 * \modifications
 *    - 5/3/02 created by Wendy Koegler
 *    - 5/21/04 WSK, renamed struct mesh to typedef MeshInfo to help
 *      clear up doxygen documentation.
 *    - 8/24/04 WSK, Changed to write coords with more precision
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

typedef struct {
  int numVert;
  int numElem;
  int numVperE;
  char elemType[1028];
  char name[1028];
  double* coords;
  int* conns;
} MeshInfo;

double* make_coords(MeshInfo *mesh);  // makes & returns pointer to coords
int* make_conns(MeshInfo *mesh); // makes & returns pointer to conns

void make_mesh(MeshInfo *mesh);  // initializes mesh
void free_mesh(MeshInfo *mesh);  // deletes coords & conns

int write_vertex_file(char* vertFileName, MeshInfo mesh);
int write_element_file(char* elemFileName, MeshInfo mesh);

int main(void) {
  int M, N, P;    // number of elements in x, y & z directions
  int i, j, k;    // indices for x, y & z directions
  int idx;      // composite index, e.g. i + j*M + k*N*M
  int numVert2d, numElem2d;   // number verts, number elements for surfaces
  int numVert3d, numElem3d;   // number verts, number elements for volumes
  int* conn_p;            // pointer to connectivity list
  int* conn_p_prev;       // another pointer to previously made conn list
  double A, B, C;         // Range of block: (0, 0, 0) to (A, B, C)
  double sigma;           // parameter for gaussian surface
  double x, y;
  MeshInfo quad;
  MeshInfo tri;
  MeshInfo line;
  MeshInfo point;
  MeshInfo hex;          
  MeshInfo tet;
  MeshInfo prism;
  MeshInfo pyramid;

  // ****************************************************************
  // setup
  // ****************************************************************
  // set the range
  A = 10;
  B = 8;
  C = 5;
  sigma = 1./7.;

  // set the number of elements in each dimension, for surfaces & for vols
  M = 11;
  N = 12;
  P = 6;
  numVert2d = (M+1)*(N+1);
  numElem2d = M*N;
  numVert3d = (M+1)*(N+1)*(P+1);
  numElem3d = M*N*P;

  // ****************************************************************
  // make QUAD vert & element file
  // ****************************************************************
  // basics for QUAD
  make_mesh(&quad);
  quad.numVert = numVert2d;
  quad.numElem = numElem2d;
  quad.numVperE = 4;
  strcpy(quad.elemType, "QUAD");
  strcpy(quad.name, "quad mesh");

  // make coordinates
  make_coords(&quad);
  for (i = 0; i < M+1; i++) 
    for (j = 0; j < N+1; j++) { 
      idx = j*(M+1) + i;
      // FIX? I don't know how portable this is
      x = quad.coords[idx*3] = A*i/M;
      y = quad.coords[idx*3+1] = B*j/N;
      quad.coords[idx*3+2] = C * exp(-0.5 * pow( (x/A - 0.5)/sigma, 2) ) * 
        exp(-0.5 * pow( (y/B - 0.5)/sigma, 2) );
    }
  
  // make connections
  conn_p = make_conns(&quad);
  for (j = 0; j < N; j++) {
    for (i = 0; i < M; i++) {
      conn_p[0] = j*(M+1) + i;
      conn_p[1] = j*(M+1) + i+1;
      conn_p[2] = (j+1)*(M+1) + i+1;
      conn_p[3] = (j+1)*(M+1) + i;
      conn_p += 4;
    }
  }

  // write vert & connection files
  write_vertex_file("gen_point_line_tri_quad.vert", quad);
  write_element_file("gen_quad.elem", quad);

  // ***************************************************************
  // make TRI element file (uses QUAD vert file)
  // ***************************************************************
  // basics for TRI
  make_mesh(&tri);
  tri.numVert = quad.numVert;
  tri.numElem = quad.numElem*2;  // divide each quad into two tris
  tri.numVperE = 3;
  strcpy(tri.elemType, "TRI");
  strcpy(tri.name, "tri mesh");

  // make connections
  conn_p = make_conns(&tri);
  conn_p_prev = quad.conns;
  for (i = 0; i < quad.numElem; i++) { // for each quad make 2 tris
    // alternate where the diagonal goes
    if (i % 2) {
      // 1st triangle
      conn_p[0] = conn_p_prev[0];
      conn_p[1] = conn_p_prev[1];
      conn_p[2] = conn_p_prev[2];
      // 2nd triangle
      conn_p[3] = conn_p_prev[0];
      conn_p[4] = conn_p_prev[2];
      conn_p[5] = conn_p_prev[3];
    }
    else {
      // 1st triangle
      conn_p[0] = conn_p_prev[0];
      conn_p[1] = conn_p_prev[1];
      conn_p[2] = conn_p_prev[3];
      // 2nd triangle
      conn_p[3] = conn_p_prev[1];
      conn_p[4] = conn_p_prev[2];
      conn_p[5] = conn_p_prev[3];
    }
    // advance pointers
    conn_p += 2*3;   // advance 2 triangles
    conn_p_prev += 4;  // advance 1 quad
  }

  // write vert & connection files
  write_element_file("gen_tri.elem", tri);

  // ****************************************************************
  // make LINE element file (uses QUAD vert file)
  // ****************************************************************
  // basics for LINE
  make_mesh(&line);
  line.numVert = quad.numVert;
  line.numElem = quad.numElem*3;  // divide each quad into three lines
  line.numVperE = 2;
  strcpy(line.elemType, "LINE");
  strcpy(line.name, "line mesh");

  // make connections
  conn_p = make_conns(&line);
  conn_p_prev = quad.conns;
  for (i = 0; i < quad.numElem; i++) { // for each quad make 3 lines
    // 1st line
    conn_p[0] = conn_p_prev[0];
    conn_p[1] = conn_p_prev[1];
    // 2nd line
    conn_p[2] = conn_p_prev[0];
    conn_p[3] = conn_p_prev[2];
    // 3rd line
    conn_p[4] = conn_p_prev[3];
    conn_p[5] = conn_p_prev[0];
    // advance pointers
    conn_p += 3*2;   // advance 3 lines
    conn_p_prev += 4;  // advance 1 quad
  }

  // write vert & connection files
  write_element_file("gen_line.elem", line);

  // ****************************************************************
  // make POINT element file (uses QUAD vert file)
  // ****************************************************************
  // basics for POINT
  make_mesh(&point);
  point.numVert = quad.numVert;
  point.numElem = quad.numVert;
  point.numVperE = 1;
  strcpy(point.elemType, "POINT");
  strcpy(point.name, "point mesh");

  // make connections
  conn_p = make_conns(&point);
  for (i = 0; i < point.numElem; i++) // each vert is it's own element
    conn_p[i] = i;

  // write vert & connection files
  write_element_file("gen_point.elem", point);

  // ****************************************************************
  // make Hex mesh
  // ****************************************************************
  // basics for HEX
  make_mesh(&hex);
  hex.numVert = numVert3d;
  hex.numElem = numElem3d;
  hex.numVperE = 8;
  strcpy(hex.elemType, "HEX");
  strcpy(hex.name, "hex mesh");

  // make coordinates
  make_coords(&hex);
  for (i = 0; i < M+1; i++) 
    for (j = 0; j < N+1; j++) 
      for (k = 0; k < P+1; k++) {
        idx = (k*(N+1) + j)*(M+1) + i;
        hex.coords[idx*3] = A*i/M;
        hex.coords[idx*3+1] = B*j/N;
        hex.coords[idx*3+2] = C*k/P;
    }

  // make connections
  conn_p = make_conns(&hex);
  for (k = 0; k < P; k++)
    for (j = 0; j < N; j++)
      for (i = 0; i < M; i++) {
        conn_p[0] = (k*(N+1) + j) * (M+1) + i;
        conn_p[1] = (k*(N+1) + j) * (M+1) + i+1;
        conn_p[2] = (k*(N+1) + (j+1)) * (M+1) + i+1;
        conn_p[3] = (k*(N+1) + (j+1)) * (M+1) + i;
        conn_p[4] = ((k+1)*(N+1) + j) * (M+1) + i;
        conn_p[5] = ((k+1)*(N+1) + j) * (M+1) + i+1;
        conn_p[6] = ((k+1)*(N+1) + (j+1)) * (M+1) + i+1;
        conn_p[7] = ((k+1)*(N+1) + (j+1)) * (M+1) + i;
        conn_p += 8;
      }

  // write vert & connection files
  write_vertex_file("gen_tet_prism_hex.vert", hex);
  write_element_file("gen_hex.elem", hex);

  // ****************************************************************
  // make TET element file (has same vert file as hex)
  // ****************************************************************
  // There are a number of different decompositions of cubes to tets
  // I chose the one way to get 5 tets from each cube
  // basics for TET
  make_mesh(&tet);
  tet.numVert = hex.numVert;
  tet.numElem = hex.numElem*5;   // 5 tets per hex
  tet.numVperE = 4;
  strcpy(tet.elemType, "TET");
  strcpy(tet.name, "tet mesh");

  // make connections
  conn_p = make_conns(&tet);
  conn_p_prev = hex.conns;
  for (i = 0; i < hex.numElem; i++) { // for each hex make 5 tets
    // tet 1
    conn_p[0] = conn_p_prev[0];
    conn_p[1] = conn_p_prev[1];
    conn_p[2] = conn_p_prev[2];
    conn_p[3] = conn_p_prev[5];
    // tet 2
    conn_p[4] = conn_p_prev[0];
    conn_p[5] = conn_p_prev[2];
    conn_p[6] = conn_p_prev[3];
    conn_p[7] = conn_p_prev[7];
    // tet 3
    conn_p[8] = conn_p_prev[0];
    conn_p[9] = conn_p_prev[5];
    conn_p[10] = conn_p_prev[7];
    conn_p[11] = conn_p_prev[4];
    // tet 4
    conn_p[12] = conn_p_prev[2];
    conn_p[13] = conn_p_prev[7];
    conn_p[14] = conn_p_prev[5];
    conn_p[15] = conn_p_prev[6];
    // tet 5
    conn_p[16] = conn_p_prev[0];
    conn_p[17] = conn_p_prev[5];
    conn_p[18] = conn_p_prev[2];
    conn_p[19] = conn_p_prev[7];
    // advance pointers
    conn_p += 5*4;   // advance 5 tets
    conn_p_prev += 8;  // advance 1 hex
  }

  // print connection file
  write_element_file("gen_tet.elem", tet);
  
  // ****************************************************************
  // make PRISM element file (has same vert file as hex)
  // ****************************************************************
  // basics for PRISM
  make_mesh(&prism);
  prism.numVert = hex.numVert;
  prism.numElem = hex.numElem*2;   // divide hex's in half
  prism.numVperE = 6;
  strcpy(prism.elemType, "PRISM");
  strcpy(prism.name, "prism mesh");

  // make connections
  conn_p = make_conns(&prism);
  conn_p_prev = hex.conns;
  for (i = 0; i < hex.numElem; i++) { // for each hex make 2 prism
    // base of first prism
    conn_p[0] = conn_p_prev[0];
    conn_p[1] = conn_p_prev[1];
    conn_p[2] = conn_p_prev[3];
    // top of first prism
    conn_p[3] = conn_p_prev[4];
    conn_p[4] = conn_p_prev[5];
    conn_p[5] = conn_p_prev[7];
    // base of second prism
    conn_p[6] = conn_p_prev[1];
    conn_p[7] = conn_p_prev[2];
    conn_p[8] = conn_p_prev[3];
    // top of second prism
    conn_p[9] = conn_p_prev[5];
    conn_p[10] = conn_p_prev[6];
    conn_p[11] = conn_p_prev[7];
    // advance pointers
    conn_p += 2*6;  // advance 2 prisms
    conn_p_prev += 8;  // advance 1 hex
  }

  // print connection file
  write_element_file("gen_prism.elem", prism);
  
  // ****************************************************************
  // make PYRAMID vert & element file (almost same vert file as hex)
  // ****************************************************************
  // make by adding another vert to center of each hex
  // basics for PYRAMID
  make_mesh(&pyramid);
  pyramid.numVert = hex.numVert + hex.numElem;
  pyramid.numElem = hex.numElem*6;   // a new pyramid for each face of a hex
  pyramid.numVperE = 5;
  strcpy(pyramid.elemType, "PYRAMID");
  strcpy(pyramid.name, "pyramid mesh");

  // make coords--copy hex's and then add coords that are centers of hex's
  make_coords(&pyramid);
  for (i = 0; i < hex.numVert*3; i++)
    pyramid.coords[i] = hex.coords[i];
  for (i = 0; i < M; i++) 
    for (j = 0; j < N; j++) 
      for (k = 0; k < P; k++) {
        idx = hex.numVert + (k*N + j)*M + i; // start at hex.numVert
        pyramid.coords[idx*3] = A*(i+0.5)/M;
        pyramid.coords[idx*3+1] = B*(j+0.5)/N;
        pyramid.coords[idx*3+2] = C*(k+0.5)/P;
    }
    
  // make connections
  conn_p = make_conns(&pyramid);
  conn_p_prev = hex.conns;
  for (i = 0; i < hex.numElem; i++) { // for each hex make 6 pyramids
    // 1st pyramid
    conn_p[0] = conn_p_prev[0];
    conn_p[1] = conn_p_prev[1];
    conn_p[2] = conn_p_prev[2];
    conn_p[3] = conn_p_prev[3];
    conn_p[4] = hex.numVert + i;
    // 2nd pyramid
    conn_p[5] = conn_p_prev[0];
    conn_p[6] = conn_p_prev[4];
    conn_p[7] = conn_p_prev[5];
    conn_p[8] = conn_p_prev[1];
    conn_p[9] = hex.numVert + i;
    // 3rd pyramid
    conn_p[10] = conn_p_prev[1];
    conn_p[11] = conn_p_prev[5];
    conn_p[12] = conn_p_prev[6];
    conn_p[13] = conn_p_prev[2];
    conn_p[14] = hex.numVert + i;
    // 4th pyramid
    conn_p[15] = conn_p_prev[2];
    conn_p[16] = conn_p_prev[6];
    conn_p[17] = conn_p_prev[7];
    conn_p[18] = conn_p_prev[3];
    conn_p[19] = hex.numVert + i;
    // 5th pyramid
    conn_p[20] = conn_p_prev[0];
    conn_p[21] = conn_p_prev[3];
    conn_p[22] = conn_p_prev[7];
    conn_p[23] = conn_p_prev[4];
    conn_p[24] = hex.numVert + i;
    // 6th pyramid
    conn_p[25] = conn_p_prev[4];
    conn_p[26] = conn_p_prev[7];
    conn_p[27] = conn_p_prev[6];
    conn_p[28] = conn_p_prev[5];
    conn_p[29] = hex.numVert + i;
    // advance pointers
    conn_p += 6*5;  // advance 6 pyramids
    conn_p_prev += 8;  // advance 1 hex
  }

  // print connection file
  write_vertex_file("gen_pyramid.vert", pyramid);
  write_element_file("gen_pyramid.elem", pyramid);
  

  // ****************************************************************
  // cleanup & exit
  // ****************************************************************
  // cleanup
  free_mesh(&quad);
  free_mesh(&tri);
  free_mesh(&line);
  free_mesh(&point);
  free_mesh(&hex);
  free_mesh(&tet);
  free_mesh(&prism);
  free_mesh(&pyramid);
  
  return 0;
}

int write_vertex_file(char* vertFileName, MeshInfo mesh) {
  int i;
  FILE *outfile;

  // open coordinate file
  outfile = fopen(vertFileName, "w");
  if (outfile == NULL) {
    fprintf(stderr, "cannot open '%s' for writing\n", vertFileName);
    exit (-1);
  }

  // write number of verts
  fprintf(outfile, "%d\n", mesh.numVert);

  // write coordinates
  for (i = 0; i < mesh.numVert; i++) {
    fprintf(outfile, "%.15g %.15g %.15g\n", mesh.coords[i*3], 
	    mesh.coords[i*3+1], mesh.coords[i*3+2]);
  }

  // close file
  fclose(outfile);

  return 0;
}

int write_element_file(char* elemFileName, MeshInfo mesh) {
  int i, j;
  FILE *outfile;

  // open connections file
  outfile = fopen(elemFileName, "w");
  if (outfile == NULL) {
    fprintf(stderr, "cannot open '%s' for writing\n", elemFileName);
    return -1;
  }

  // write vert & element count
  fprintf(outfile, "%d %d\n", mesh.numVert, mesh.numElem);

  // write element type
  fprintf(outfile, "%s\n", mesh.elemType);

  // write mesh name
  fprintf(outfile, "%s\n", mesh.name);

  // write element vert ids
  for (i = 0; i < mesh.numElem; i++) {
    for (j = 0; j < mesh.numVperE-1; j++) 
      fprintf(outfile, "%d ", mesh.conns[i*mesh.numVperE + j]);
    fprintf(outfile, "%d\n", mesh.conns[i*mesh.numVperE + mesh.numVperE - 1]);
  }

  // close file
  fclose(outfile);

  return 0;
}

double* make_coords(MeshInfo *mesh) {
  mesh->coords = (double*)malloc(mesh->numVert*3*sizeof(double));
  return mesh->coords;
}

int* make_conns(MeshInfo *mesh) {
  mesh->conns = (int*)malloc(mesh->numElem*mesh->numVperE*sizeof(int));
  return mesh->conns;
}

void make_mesh(MeshInfo *mesh) {
  // set to 0
  mesh->numVert = 0;
  mesh->numElem = 0;
  mesh->numVperE = 0;
  strcpy(mesh->elemType, "");
  strcpy(mesh->name, "");
  mesh->coords = NULL;
  mesh->conns = NULL;
}
void free_mesh(MeshInfo *mesh) {
  // set to 0
  mesh->numVert = 0;
  mesh->numElem = 0;
  mesh->numVperE = 0;
  strcpy(mesh->elemType, "");
  strcpy(mesh->name, "");

  // free dynamically allocated
  if(mesh->coords)
    free(mesh->coords);
  mesh->coords = NULL;
  if(mesh->conns)
    free(mesh->conns);
  mesh->conns = NULL;
}







