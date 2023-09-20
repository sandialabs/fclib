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
 * \file screws_generator.c
 * \brief A simple program to create a dataset suitable for screwBreak testing.
 *
 * $Source: /usr/local/Repositories/fcdmf/fclib/data/screws_generator.c,v $
 * $Revision: 1.4 $
 * $Date: 2006/11/07 23:49:01 $
 *
 * \description
 *
 *   Copy the screw mesh and each have different damage fields so can test for
 *     - breaking through shaft
 *     - breaking by end erosion
 *     - breaking by shaft side erosion
 *
 *   The first mesh will have 1 screw and will have a tear. The 2nd mesh 
 *   will have 2 screws and one will erode at the bottom, and the other
 *   will erod at the side.
 * 
 * \modifications:
 *    - 07/18/2006 Created, WSD.
 */
#include <string.h>
#include "fc.h"
#include "fcP.h"

int main(void) {
  FC_ReturnCode rc;
  int i, j, k;
  int numDim;
  FC_Dataset dataset;
  FC_Sequence sequence;
  int numMesh = 2, numStep = 11, temp_numStep;
  int numVerts[2], numElems[2], numVertPerElem;
  FC_ElementType elemType;
  FC_Mesh meshes[2];
  char* file_name_exo = "../data/gen_screws.ex2";
  char* mesh_names[2] = { "screw-tear", "screws-erode" };
  FC_Variable *damages[2];
  char* var_name = "damage";
  double* time_coords;
  int *conns, *newConns;
  double *coords, *newCoords;
  double *data;
  FC_Coords lowers, uppers;
  double spacing;
  int breakStepID;
  // the tear IDs  will be a cross section through shaft and parrallel to the
  // end, therefore the IDs are equal to endIDs + some multiple of numEnd > 1.
  int numTear = 140;
  int tearIDs[140] = { 
    1971, 1972, 1973, 1974, 1975, 1976, 1977, 1978, 1979, 1980,
    1981, 1982, 1983, 1984, 1985, 1986, 1987, 1988, 1989, 1990,
    1991, 1992, 1993, 1994, 1995, 1996, 1997, 1998, 1999, 2000,
    2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010,
    2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020,
    2021, 2022, 2023, 2024, 2025, 2026, 2027, 2028, 2029, 2030,
    2031, 2032, 2033, 2034, 2035, 2036, 2037, 2038, 2039, 2040,
    2041, 2042, 2043, 2044, 2045, 2046, 2047, 2048, 2049, 2050,
    2051, 2052, 2053, 2054, 2055, 2056, 2057, 2058, 2059, 2060,
    2061, 2062, 2063, 2064, 2065, 2066, 2067, 2068, 2069, 2070,
    2071, 2072, 2073, 2074, 2075, 2076, 2077, 2078, 2079, 2080,
    2081, 2082, 2083, 2084, 2085, 2086, 2087, 2088, 2089, 2090,
    2091, 2092, 2093, 2094, 2095, 2096, 2097, 2098, 2099, 2100,
    2101, 2102, 2103, 2104, 2105, 2106, 2107, 2108, 2109, 2110 };
  int numEnd = 140;
  int endIDs[140] = { 
    0, 2251, 2252, 2253, 2254, 2255, 2256, 2257, 2258, 2259,
    2260, 2261, 2262, 2263, 2264, 2265, 2266, 2267, 2268, 2269,
    2270, 2271, 2272, 2273, 2274, 2275, 2276, 2277, 2278, 2279,
    2280, 2281, 2282, 2283, 2284, 2285, 2286, 2287, 2288, 2289,
    2290, 2291, 2292, 2293, 2294, 2295, 2296, 2297, 2298, 2299,
    2300, 2301, 2302, 2303, 2304, 2305, 2306, 2307, 2308, 2309,
    2310, 2311, 2312, 2313, 2314, 2315, 2316, 2317, 2318, 2319,
    2320, 2321, 2322, 2323, 2324, 2325, 2326, 2327, 2328, 2329,
    2330, 2331, 2332, 2333, 2334, 2335, 2336, 2337, 2338, 2339,
    2340, 2341, 2342, 2343, 2344, 2345, 2346, 2347, 2348, 2349,
    2350, 2351, 2352, 2353, 2354, 2355, 2356, 2357, 2358, 2359,
    2360, 2361, 2362, 2363, 2364, 2365, 2366, 2367, 2368, 2369,
    2370, 2371, 2372, 2373, 2374, 2375, 2376, 2377, 2378, 2379,
    2380, 2381, 2382, 2383, 2384, 2385, 2386, 2387, 2388, 2389 };
  int numSide = 228;
  int sideIDs[228] = { 
    1551, 1552, 1553, 1554, 1555, 1556, 1557, 1558, 1559, 1560,
    1561, 1562, 1563, 1564, 1565, 1566, 1567, 1568, 1569, 1570,
    1571, 1572, 1573, 1574, 1575, 1576, 1577, 1578, 1579, 1580,
    1581, 1582, 1583, 1584, 1585, 1586, 1587, 1588, 1691, 1692,
    1693, 1694, 1695, 1696, 1697, 1698, 1699, 1700, 1701, 1702,
    1703, 1704, 1705, 1706, 1707, 1708, 1709, 1710, 1711, 1712,
    1713, 1714, 1715, 1716, 1717, 1718, 1719, 1720, 1721, 1722,
    1723, 1724, 1725, 1726, 1727, 1728, 1831, 1832, 1833, 1834,
    1835, 1836, 1837, 1838, 1839, 1840, 1841, 1842, 1843, 1844,
    1845, 1846, 1847, 1848, 1849, 1850, 1851, 1852, 1853, 1854,
    1855, 1856, 1857, 1858, 1859, 1860, 1861, 1862, 1863, 1864,
    1865, 1866, 1867, 1868, 1971, 1972, 1973, 1974, 1975, 1976,
    1977, 1978, 1979, 1980, 1981, 1982, 1983, 1984, 1985, 1986,
    1987, 1988, 1989, 1990, 1991, 1992, 1993, 1994, 1995, 1996,
    1997, 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006,
    2007, 2008, 2111, 2112, 2113, 2114, 2115, 2116, 2117, 2118,
    2119, 2120, 2121, 2122, 2123, 2124, 2125, 2126, 2127, 2128,
    2129, 2130, 2131, 2132, 2133, 2134, 2135, 2136, 2137, 2138,
    2139, 2140, 2141, 2142, 2143, 2144, 2145, 2146, 2147, 2148,
    2251, 2252, 2253, 2254, 2255, 2256, 2257, 2258, 2259, 2260,
    2261, 2262, 2263, 2264, 2265, 2266, 2267, 2268, 2269, 2270,
    2271, 2272, 2273, 2274, 2275, 2276, 2277, 2278, 2279, 2280,
    2281, 2282, 2283, 2284, 2285, 2286, 2287, 2288 };

  // Init library & create dataset
  rc = fc_setLibraryVerbosity(FC_WARNING_MESSAGES);
  fc_exitIfError(rc);
  rc = fc_initLibrary();
  fc_exitIfError(rc);
  rc = fc_createDataset("Breaking screw examples", &dataset);
  fc_exitIfErrorPrintf(rc, "Failed to create dataset '%s'", file_name_exo);

  // create a sequence
  time_coords = malloc(numStep*sizeof(double));
  for (i = 0; i < numStep; i++)
    time_coords[i] = i + .1*i;
  fc_createSequence(dataset, "time", &sequence);
  fc_setSequenceCoordsPtr(sequence, numStep, FC_DT_DOUBLE, 
			  (void*)time_coords);

  // create the first screw from files
  rc = _fc_makeMeshFromFiles(dataset, "screw.vert", "screw.elem", &meshes[0]);
  fc_exitIfErrorPrintf(rc, "Failed to create screw mesh from fc files");
  fc_changeMeshName(meshes[0], mesh_names[0]);
  fc_getMeshInfo(meshes[0], NULL, &numDim, &numVerts[0], &numElems[0], &elemType);
  fc_getMeshCoordsPtr(meshes[0], &coords);
  fc_getMeshElementConnsPtr(meshes[0], &conns);

  // Make the next two screws from first, space over 1 1/2 times
  rc = fc_getMeshBoundingBox(meshes[0], &numDim, &lowers, &uppers);
  numVertPerElem = fc_getElementTypeNumVertex(elemType);
  fc_exitIfError(rc);
  spacing = uppers[0] - lowers[0];
  spacing *= 1.5;
  numVerts[1] = 2*numVerts[0];
  numElems[1] = 2*numElems[0];
  newCoords = (double*)malloc(numVerts[1]*numDim*sizeof(double));
  newConns = (int*)malloc(numElems[1]*numVertPerElem*sizeof(int));
  if (!newCoords || !newConns)
    fc_exitIfError(rc);
  for (i = 0; i < numVerts[0]; i++) {
    newCoords[i*numDim] = coords[i*numDim] + spacing;
    newCoords[i*numDim+1] = coords[i*numDim+1];
    newCoords[i*numDim+2] = coords[i*numDim+2];
  }
  for (i = 0; i < numVerts[0]; i++) {
    newCoords[numVerts[0]*numDim + i*numDim] = coords[i*numDim] + 2*spacing;
    newCoords[numVerts[0]*numDim + i*numDim+1] = coords[i*numDim+1];
    newCoords[numVerts[0]*numDim + i*numDim+2] = coords[i*numDim+2];
  }
  memcpy(newConns, conns, numElems[0]*numVertPerElem*sizeof(int));
  for (i = 0; i < numElems[0]; i++) {
    for (j = 0; j < numVertPerElem; j++)
      newConns[numElems[0]*numVertPerElem + i*numVertPerElem+j] =
	conns[i*numVertPerElem+j] + numVerts[0];
  }
  rc = fc_createMesh(dataset, mesh_names[1], &meshes[1]);
  fc_exitIfError(rc);
  rc = fc_setMeshCoordsPtr(meshes[1], numDim, numVerts[1], newCoords);
  fc_exitIfError(rc);
  rc = fc_setMeshElementConnsPtr(meshes[1], elemType, numElems[1], newConns);
  fc_exitIfError(rc);

  // temp stuff to get our ids for dead elements
  // figured out the end because it had smallest # of elems
  // figured out side by checking in ensight
  // figured out tear by subtracking end elements off twice 
  //   and taking the new end
  /*
  {
    int numShape;
    FC_Shape *shapes;
    FC_Subset subset0, subset1, subset2;
    int sideID, numSide;

    // erode & side
    fc_getMeshShapes(meshes[0], 45, 0, &numShape, &shapes);
    for (i = 0; i < 2; i++) {
      int numMember, *memberIDs;
      fc_getSubsetMembersAsArray(shapes[0].elems[i], &numMember, &memberIDs);
      printf("mesh side %d:\n", i);
      printf("int numMember = %d;\n", numMember);
      printf("int memberIDs[%d] = { ", numMember);
      for (j = 0; j < numMember; j++) {
	printf("%d, ", memberIDs[j]);
	if ((j+1)%10 == 0)
	  printf("\n");
      }
      printf(" };\n");
    }
    
    // tear
    fc_createSubset(meshes[0], "whole", FC_AT_ELEMENT, &subset0);
    for (i = 0; i < numElem; i++)
      fc_addMemberToSubset(subset0, i);
    fc_createSubsetIntersection(subset0, "XOR", shapes[0].elems[0], "temp", 
				&subset1);
    fc_getSubsetShapes(subset1, 45, 0, &numShape, &shapes);
    fc_createSubsetIntersection(subset1, "XOR", shapes[0].elems[0], "temp",
				&subset2);
    fc_getSubsetShapes(subset2, 45, 0, &numShape, &shapes);
    for (i = 0; i < 1; i++) {
      int numMember, *memberIDs;
      fc_getSubsetMembersAsArray(shapes[0].elems[i], &numMember, &memberIDs);
      printf("subset side %d:\n", i);
      printf("int numMember = %d;\n", numMember);
      printf("int memberIDs[%d] = { ", numMember);
      for (j = 0; j < numMember; j++) {
	printf("%d, ", memberIDs[j]);
	if ((j+1)%10 == 0)
	  printf("\n");
      }
      printf(" };\n");
    }
  }
  */

  // create the damage seq vars
  breakStepID = 8;
  for (i = 0; i < numMesh; i++) {
    rc = fc_createSeqVariable(meshes[i], sequence, var_name, &temp_numStep,
			 &damages[i]);
    fc_exitIfError(rc);
    for (j = 0; j < numStep; j++) {
      double fraction;
      data = (double*)calloc(numElems[i], sizeof(double)); // init to zero
      //switch (i) {
      //case 0: num = numTear; ids = tearIDs; break;
      //case 1: num = numEnd; ids = endIDs;   break;
      //case 2: num = numSide; ids = sideIDs; break;
      // }
      if (j == 0) 
	fraction = 0;
      else if (j < breakStepID)
	fraction = ((double)j)/numStep;
      else
	fraction = 1;
      if (i == 0) {
	for (k = 0; k < fraction*numTear; k++)
	  data[tearIDs[k]] = 1;
      }
      else {
	// stagger the breaking time of the two screw mesh
	for (k = 0; k < fraction*numEnd; k++) 
	  data[endIDs[k]] = 1;
	if (j+1 == breakStepID)  // make the 2nd screw break a timestep early
	  fraction = 1;
	for (k = 0; k < fraction*numSide; k++) 
	  data[sideIDs[k]+numElems[0]] = 1;
      }
      rc = fc_setVariableDataPtr(damages[i][j], numElems[i], 1, FC_AT_ELEMENT,
				 FC_MT_SCALAR, FC_DT_DOUBLE, (void*)data);
      fc_exitIfErrorPrintf(rc, "failed to set variable data");
    }
  }

  // Write the dataset (have to make a copy first)
  fc_writeDataset(dataset, file_name_exo, FC_FT_EXODUS);

  // Final library
  fc_finalLibrary();

  exit(0);  
}
