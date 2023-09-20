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
 * \file seriesAnalyzer.c
 * \brief reads in series and does some analyses on them.
 *
 * $Source: /usr/local/Repositories/fcdmf/fclib/tools/seriesAnalyzer.c,v $
 * $Revision: 1.28 $
 * $Date: 2006/08/30 19:20:04 $
 *
 * \description
 *   Given an inputfile containing the datafilenames it
 *   reads in series from file and does some analyses on them.
 *   Wired into the code right now is the option of doing the calcs
 *   on the data itself, or on its normal form.
 *   MinMax and integral on each. Euclidean distance 
 *   and dimensionless area between all pairs.
 *
 *  commented out - prints file of normal form of data set and its derivative - see
 *   note about this in todo and the normalwriteoutfunction for warning.
 *  (der isnt quite right but it jsut means its off by some const)
 *
 *  Only works on 1 data point 1 component data.
 *
 *  reads in from inputfile in format:
 *   numfiles
 *   datafilename1
 *   datafilename2
 *   datafilenameN
 *
 *  data files in format:
 *   #dataname seqname numseqpts varname  (e.g., can1 displacement 20 Force)
 *   time1 val1 
 *   time2 val2 
 *   time3 valN 
 *
 *  writes results to "<filename>.out" in format
 *   //#regularized intersecting sequences
 *   //analyses results - min, max, integral, euc, and klein for
 *   //each pair ordered by closest klein
 *   //two blank lines at the end of a seq set then start the next one
 *
 * \modifications
 *   - 04/19/05 ACG Created.
 *   - 07/21/05 ACG removing gnu option for now
 *   - 07/27/05 ACG trimmed mods list so only hot current stuff listed
 *   - 07/27/05 ACG parameter passed for derivative is wired into the code.
 *   - 07/27/05 ACG having it print out normal forms of data and its derivative for plotting purposes
 *
 * \todo  
 *    - fix normal form of derivative to account for fist and last pts being zero
 *    - is there a better way to get the derivative tolerance parameter
 *    - change wired in raw data vs normal form option.
 *    - relax assumptions on single input data point
 *    - more FT (skip over and continue)
 *      on missing files or bad files or badnames in inputdeck
 */

#include <fc.h>
#include <math.h>
#include <string.h>

//globals
FC_Mesh glob_mesh;
FC_Dataset glob_dataset;
int glob_meshDataPoint;
//must be single datapoint single component
int numComponent = 1;
int glob_meshDim = 1;


//reading in data
FILE *fp, *fpout, *fpder,*fpgnu, *lfp;
int maxchar = 80;

typedef enum{
  VMAX = 0,
  VMIN = 1,
  VINT = 2,    
  VEUC = 3,
  VKLEIN = 4,
} stattype;


typedef struct{
  char* name; //will be responsible for freeing this memory
  //  char* refname; //will be responsible for freeing this memory
  double max;
  double min;
  int maxindex;
  int minindex;
  double integral;
  double eucdist;
  double kleindist;
} Stats;



static void statsInit(Stats *x){
  x->name = NULL;
  //  x.refname = NULL;
  x->min = 0.;
  x->max = 0.;
  x->minindex = 0;
  x->maxindex = 0;
  x->integral = 0.;
  x->eucdist = 0;
  x->kleindist = 0;
}


static void setupFiles(char* dirin, char* namein, char* basedir,
                       char* filein, char* fileout, char* fileder, char* filegnu){
  strcpy(basedir,dirin);
  strcat(basedir,"/");

  strcpy(filein, basedir);
  strcat(filein, namein);

  strcpy(fileout, filein);
  strcat(fileout, ".out");

  strcpy(fileder, filein);
  strcat(fileder, ".der");

  strcpy(filegnu, filein);
  strcat(filegnu, ".gnu");

  printf("inputfile <%s> outputfile <%s>\n",filein,fileout);
  fflush(NULL);

  fp = fopen(filein,"r");
  if (fp == NULL){
    printf("failed to open input file\n");
    fflush(NULL);
    exit (-1);
  }
  fpout = fopen(fileout,"w");
  if (fpout == NULL){
    printf("failed to open output file\n");
    fflush(NULL);
    exit (-1);
  }
  fprintf(fpout, "inputfile <%s> outputfile <%s>\n",filein,fileout);
  fflush(NULL);

  fpder = fopen(fileder,"w");
  if (fpder == NULL){
    printf("failed to open derivative file\n");
    fflush(NULL);
    exit (-1);
  }
  /*
  fpgnu = fopen(filegnu,"w");      
  if (fpgnu == NULL){
    printf("failed to open gnu output file\n");
    fflush(NULL);
    exit (-1);
  }
  */
}


static void setupDataset(int k){
  //parameter doesnt matter

  FC_ReturnCode rc;
  double* glob_meshcoords;
  int i =k;

  //make a data set to play in
  rc = fc_createDataset("glob dataset", &glob_dataset);
  fc_exitIfErrorPrintf(rc, "abort: failed to create dataset");

  //make a mesh
  rc = fc_createMesh(glob_dataset," glob mesh", &glob_mesh);
  fc_exitIfErrorPrintf(rc, "abort: failed to create mesh");
  //make mesh coords - doesnt matter what these values are
  //hand these off to the mesh
  glob_meshcoords = (double*) malloc(glob_meshDataPoint*sizeof(double));
  for (i = 0; i < glob_meshDataPoint; i++){
    glob_meshcoords[i] = (double)i;
  }
  rc = fc_setMeshCoords(glob_mesh, glob_meshDim, glob_meshDataPoint,
                        glob_meshcoords);
  fc_exitIfErrorPrintf(rc, "abort: failed to set vertex coords");

  free (glob_meshcoords);
}


static int readOneData(int id, int* numSeqStep,
                     FC_Sequence* seq, FC_Variable** seqVar){
  //read in data. exits entirely if there is problem

  FC_ReturnCode rc;
  int i,k;
  int ii = id;
  char dataname[maxchar];
  char name[maxchar];
  char seqname[maxchar];
  char varname[maxchar];
  double* timecoords;
  double* data;

  //   #dataname seqname numseqpts varname  (e.g., can1 displacement 20 Force)

  //read in data
  fscanf(lfp,"#%s %s %d %s\n",dataname,seqname,&numSeqStep[ii],
         varname);
  printf("dataname = %s seqname = %s numStep = %d varname = %s\n",
         dataname,seqname,numSeqStep[ii],varname);
  fprintf(fpout,"dataname = %s seqname = %s numStep = %d varname = %s\n",
          dataname,seqname,numSeqStep[ii],varname);
  fflush(NULL);

  //the memory for these gets handed off
  timecoords = (double*)malloc(numSeqStep[ii]*sizeof(double));
  data = (double*)malloc(numSeqStep[ii]*sizeof(double));
  
  //now read in its data
  for (i = 0; i < numSeqStep[ii]; i++){
    if (fscanf(lfp,"%lg %lg\n",&timecoords[i],&data[i]) !=2){
      printf("error in input data format");
      fprintf(fpout, "error in input data format");
      fflush(NULL);
      exit (-1);
    }
  }

  // Create the sequence
  strcpy(name,dataname);
  strcat(name,"_");
  strcat(name,seqname);
  rc = fc_createSequence(glob_dataset, name, &seq[ii]);
  fc_exitIfErrorPrintf(rc, "failed to create sequence");
  rc = fc_setSequenceCoordsPtr(seq[ii], numSeqStep[ii],
                               FC_DT_DOUBLE, (void*)timecoords);
  fc_exitIfErrorPrintf(rc, "failed to set time values");
  
  //create seq var
  rc = fc_createSeqVariable(glob_mesh,seq[ii],dataname,&k,&seqVar[ii]);
  fc_exitIfErrorPrintf(rc, "failed to create sv1");
  
  for (i = 0; i < numSeqStep[ii]; i++){
    rc = fc_setVariableData(seqVar[ii][i],glob_meshDataPoint,
                            numComponent,
                            FC_AT_VERTEX,FC_MT_SCALAR,
                            FC_DT_DOUBLE,
                            &data[i*glob_meshDataPoint*numComponent]);
    fc_exitIfErrorPrintf(rc, "failed to set data for sv");
  }
  free (data);

  return 1;
}


static int readData(char *basedir, int numSeq, int* numSeqStep,
                    FC_Sequence* seq, FC_Variable** seqVar){
  int i = 0;
  char lfilein[maxchar];

  while(fscanf(fp,"%s\n",lfilein) && i < numSeq){
    char fname[maxchar];
    strcpy(fname,basedir);
    strcat(fname,lfilein);

    lfp = fopen(fname,"r");
    if (lfp == NULL){
      fc_printfErrorMessage("failed to open input file %s", lfilein);
      i++;
      return 0; //make this more robust to missing files later
      //      continue;
    }
    //now read in sets of data
    if (readOneData(i,numSeqStep,seq,seqVar) == 0){
      fc_printfErrorMessage("failed to read input file %s", lfilein);
      return 0;
    }
    fclose(lfp); 
    i++;
  }
  
  fclose(fp);
  return 1;
}


static int makeIntersecting(double tmin, double tmax, int numSeq,
                             int* numSeqStep, FC_Sequence* seq, FC_Variable** seqVar,
                             int*  numStepIntersect, FC_Sequence* seqintersect,
                             FC_Variable** seqVarIntersect){
  //make comparison seqs over teh range
  FC_ReturnCode rc;
  int i, ii;


  printf("Making Intersecting SeqVars.....\n");
  fflush(NULL);

  //default
  *numStepIntersect = 0;

  //first get intersecting sequence
  if (FC_DBL_EQUIV(tmin,0.0) && FC_DBL_EQUIV(tmax,0.0)){
    //then do it for the whole range
    rc = fc_createIntersectingRegularSequence(numSeq,seq,0,
                                              "seqintersect", seqintersect);
    if (rc != FC_SUCCESS){
      fprintf(fpout,"Can't get intersecting range over whole range - exiting\n");
      fflush(NULL);
      fc_exitIfErrorPrintf(rc, "Can't get intersecting range over whole range - exiting\n");
      return 0; //this wont matter
    }

  } else {
    //do it for the bounds
    rc = fc_createIntersectingRangeRegularSequence(numSeq,seq,tmin,tmax,
                                                   "seqintersect", seqintersect);
    if (rc != FC_SUCCESS){
      fprintf(fpout,"Can't get intersecting range over this range\n");
      printf("Can't get intersecting range over this range\n");
      fflush(NULL);
      return 0;
    }
  }

  rc = fc_getSequenceNumStep(*seqintersect,numStepIntersect);
  if (rc != FC_SUCCESS){
    fprintf(fpout,"failed to get numsteps on new seq\n");
    printf("failed to get numsteps on new seq\n");
    fflush(NULL);
    fc_deleteSequence(*seqintersect);
    return 0;
  }
  
  //now put them all on the same seq
  for (ii = 0; ii < numSeq; ii++){
    char *name1, *name2;
    
    //now get linear interpolation 
    rc = fc_getVariableName(seqVar[ii][0],&name1);
    name2 = (char*)malloc((strlen(name1))+5*sizeof(char));
    strcpy(name2,name1);
    strcat(name2,"_new");
    rc =  fc_linearInterpolation (numSeqStep[ii], seqVar[ii], *seqintersect,
                                  name1,&seqVarIntersect[ii]);

    free(name1);
    free(name2);
    if (rc != FC_SUCCESS){
      fprintf(fpout,"failed to do linear interpolation\n");
      printf("failed to do linear interpolation\n");
      fflush(NULL);
      //clean up all the ones weve made so far
      for ( i = 0; i <= ii; i++){
        rc = fc_deleteSeqVariable(*numStepIntersect,seqVarIntersect[i]);
      }
      return 0;
    }
  }

  return 1;
}


static int makeNormal(int numStepIntersect, 
                       int numSeq, FC_Variable **seqVarIntersect, 
                      FC_Variable **seqVarIntersectNormal){
  FC_ReturnCode rc;
  int i,ii;


  printf("Making NormalForms of Intersecting SeqVars.....\n");
  fflush(NULL);
  for (ii = 0; ii < numSeq; ii++){
    char *name1,*name2;

    rc = fc_getVariableName(seqVarIntersect[ii][0],&name1);
    name2 = (char*)malloc((strlen(name1))+5*sizeof(char));
    strcpy(name2,name1);
    strcat(name2,"_N");

    rc = fc_normalForm(numStepIntersect,seqVarIntersect[ii],name2,&seqVarIntersectNormal[ii]);
    if (rc != FC_SUCCESS){
      fprintf(fpout,"#error:cant get normalform for this var\n");
      fflush(NULL);
      //clean up all the ones weve made so far
      for ( i = 0; i <= ii; i++){
        rc = fc_deleteSeqVariable(numStepIntersect,seqVarIntersectNormal[i]);
      }
      free(name1);
      free(name2);
      return 0;
    }


    free(name1);
    free(name2);
  }
  return 1;
}

static void doSingles(int numStepIntersect,int numSeq, FC_Variable** seqVarIntersect,
                      Stats** svstats){
  //do stats that are on individual sequence vars
  FC_ReturnCode rc;
  int ii;

  printf("Performing Single Seq Analyses\n");
  fflush(NULL);
  for (ii = 0; ii < numSeq; ii++){
    FC_Variable minvar, maxvar, minindexvar, maxindexvar, integralvar;
    double *mindata, *maxdata, *integraldata;
    int *minindexdata, *maxindexdata;
    char* vname;
    
    rc = fc_getVariableName(seqVarIntersect[ii][0],&vname);
    rc =  fc_getSeqVariableSeriesMinMax (numStepIntersect, seqVarIntersect[ii],
					&minvar, &maxvar, &minindexvar,
					&maxindexvar); 
    fc_exitIfErrorPrintf(rc, "failed to get minmax for sv");
    rc = fc_getVariableDataPtr(minindexvar,(void**)&minindexdata);
    fc_exitIfErrorPrintf(rc, "failed to get min indexvariable data");
    rc = fc_getVariableDataPtr(minvar,(void**)&mindata);
    fc_exitIfErrorPrintf(rc, "failed to get min variable data");
    rc = fc_getVariableDataPtr(maxindexvar,(void**)&maxindexdata);
    fc_exitIfErrorPrintf(rc, "failed to get max indexvariable data");
    rc = fc_getVariableDataPtr(maxvar,(void**)&maxdata);
    fc_exitIfErrorPrintf(rc, "failed to get max variable data");
    
    rc = fc_integral_TR(numStepIntersect, seqVarIntersect[ii],
                                    "intvar",&integralvar);
    fc_exitIfErrorPrintf(rc, "failed to get integral for sv");
    rc = fc_getVariableDataPtr(integralvar,(void**)&integraldata);
    fc_exitIfErrorPrintf(rc, "failed to get integral data");


    //fill in stats
    (*svstats)[ii].name = (char*)malloc((strlen(vname)+1)*sizeof(char));
    strcpy((*svstats)[ii].name,vname);
    //know only 1DP 1component data
    (*svstats)[ii].min = mindata[0];
    (*svstats)[ii].max = maxdata[0];
    (*svstats)[ii].minindex = minindexdata[0];
    (*svstats)[ii].maxindex = maxindexdata[0];
    (*svstats)[ii].integral = integraldata[0];
    

    //clean up
    free(vname);
    fc_deleteVariable(minvar);
    fc_deleteVariable(maxvar);
    fc_deleteVariable(minindexvar);
    fc_deleteVariable(maxindexvar);
    fc_deleteVariable(integralvar);
  }

}


static void doDoubles(int numStepIntersect,FC_Variable *seqVarref, FC_Variable *seqVarcomp,
                      Stats *statscomp){

  FC_ReturnCode rc;
  double *data1, *data2;
  FC_Variable euc_var, area_var;

  printf("Performing Comparison Seq Analyses\n");
  fflush(NULL);
  //    svstats[ii].refname = (char*)malloc((strlen(name1))*sizeof(char));
  //    strcpy(svstats[ii].refname,name1);
  
  //first do euclidean distance
  rc =  fc_euclideanDistanceBetweenCurves(numStepIntersect, seqVarref,
                                          seqVarcomp,
                                        "euc_dist", &euc_var);
  fc_exitIfErrorPrintf(rc, "failed to calc euclidean dist");
  rc = fc_getVariableDataPtr(euc_var,(void**)&data1); //know this is a double
  fc_exitIfErrorPrintf(rc, "failed to get euclidean dist data");
        
  //now do area between curves
  rc =  fc_dimensionlessAreaBetweenCurves (numStepIntersect,
                                               seqVarref,seqVarcomp,
                                               "euc_dist", &area_var);
  fc_exitIfErrorPrintf(rc, "failed to calc dimensionless area");
  rc = fc_getVariableDataPtr(area_var,(void**)&data2); //know this is a double
  fc_exitIfErrorPrintf(rc, "failed to get euclidean dist data");
        
  //single component, single data point
  statscomp->eucdist = data1[0];
  statscomp->kleindist = data2[0];
        
  //clean up
  fc_deleteVariable(euc_var);
  fc_deleteVariable(area_var);
}


static int origWriteout(int numSeq, int* numSeqStep,
                         FC_Sequence* seq, FC_Variable **seqVar){
  //writeout of orig data for gnuplot. gnuplot file is based on filein, not fileout!
  //if there is an error in there then it jus aborts on the writeout, btu the program continues.

  char *name;
  double *dataa, *coords;
  int i, ii;
  FC_ReturnCode rc;
  int ret = 1;
  
  fprintf(fpout,"\n\n#Orig sequences\n");
  fflush(NULL);
  for (ii = 0; ii < numSeq; ii++){
    rc = fc_getVariableName(seqVar[ii][0],&name);
    if (rc != FC_SUCCESS){
      fprintf(fpout,"#error:cant get name for this var\n");
      fflush(NULL);
      ret = 0;
      continue;
    }
    fprintf(fpout,"#%10s\n ",name);
    fflush(NULL);
    //    fprintf(fpgnu,"%s \"%s\" index %d title \"%s\" with lines\n",
    //      (ii == 0 ? "plot": "replot"),filein,ii,name);
    
    free(name);
    
    rc = fc_getSequenceCoordsPtr(seq[ii],(void**)&coords);
    if (rc != FC_SUCCESS){
      fprintf(fpout,"#error:cant get coords for this var\n");
      fflush(NULL);
      ret = 0;
      continue;
    }
    
    for (i = 0; i < numSeqStep[ii]; i++){
      rc = fc_getVariableDataPtr(seqVar[ii][i],(void**)&dataa);
      if (rc != FC_SUCCESS){
        fprintf(fpout,"#error:cant get variable data\n");
        fflush(NULL);
        ret = 0;
        break;
      }
      //know its one data point one component
      fprintf(fpout,"%10g %10g\n",coords[i],dataa[0]);
      fflush(NULL);
    }
    fprintf(fpout,"\n\n");
    fflush(NULL);
  }

  return ret;
}


static void writeoutCompare(int numStepIntersect, FC_Sequence seqintersect, 
                                 int numSeq, FC_Variable **seqVarIntersect){
  FC_ReturnCode rc;
  char *name;
  double *dataa, *coords;
  int i,ii;


  rc = fc_getSequenceCoordsPtr(seqintersect,(void**)&coords);
  fc_exitIfErrorPrintf(rc, "failed to get seq coords");      
  fprintf(fpout,"\n\n#Time Range: %10g - %10g\n",coords[0],coords[numStepIntersect-1]);
  printf("\n\n#Time Range: %10g - %10g\n",coords[0],coords[numStepIntersect-1]);
  fprintf(fpout,"#comparison sequences\n#%9s ","time");
  fflush(NULL);
  //      printf("#comparison sequences\n#%9s ","time");
  for (ii = 0; ii < numSeq; ii++){
    rc = fc_getVariableName(seqVarIntersect[ii][0],&name);
    fc_exitIfErrorPrintf(rc, "failed to get var name");
    fprintf(fpout,"%10s ",name);
    //  printf("%10s ",name);
    fflush(NULL);
    free(name);
  }
  fprintf(fpout,"\n");
  //      printf("\n");
  fflush(NULL);
  
  for (i = 0; i < numStepIntersect; i++){
    fprintf(fpout,"%10g ", coords[i]);
    //  printf("%10g ", coords[i]);
    fflush(NULL);
    for (ii = 0; ii < numSeq; ii++){
      rc = fc_getVariableDataPtr(seqVarIntersect[ii][i],(void**)&dataa);
      fc_exitIfErrorPrintf(rc, "failed to get variable data");
      //know its one data point one component
      fprintf(fpout,"%10g ",dataa[0]);
      //          printf("%10g ",dataa[0]);
      fflush(NULL);
    }
    fprintf(fpout,"\n");
    //  printf("\n");
    fflush(NULL);
  }
  fprintf(fpout,"\n\n");
  //      printf("\n\n");
  fflush(NULL);
}


/*
static void writeoutNormal(int numStepIntersect, FC_Sequence seqintersect, 
                                 int numSeq, FC_Variable **seqVarIntersect){
  //writeing out to derivate file normal forms of orig seq and derivate. 
  //recall that derivative is given as zero for the first two vals and last two vals

  FC_ReturnCode rc;
  int i,ii;
  for (ii = 0; ii < numSeq; ii++){
    FC_Variable *derivativevar, *normalseqvar, *normaldervar;//, *tempdervar;
      double *derivativedata, *seqvardata, *coords;


      //picking huge tolerance cuase we know its regular
      rc = fc_firstDerivative_REA(numStepIntersect, seqVarIntersect[ii], FLT_EPSILON,
                                  "dervar",&derivativevar);
      fc_exitIfErrorPrintf(rc, "failed to get derivative for sv");

      //first and last 2 dp of derivate var are zero, so we want to shift and unshift the data for the normal form
      //havent doe this yet

      rc = fc_getSequenceCoordsPtr(seqintersect,(void**)&coords);
      if (rc != FC_SUCCESS){
        fprintf(fpout,"#error:cant get coords for this var\n");
        fflush(NULL);
        continue;
      }

      rc = fc_normalForm(numStepIntersect,seqVarIntersect[ii],"normalintersect",&normalseqvar);
      if (rc != FC_SUCCESS){
        fprintf(fpout,"#error:cant get normalform for this var\n");
        fflush(NULL);
        continue;
      }


      rc = fc_normalForm(numStepIntersect,derivativevar,"normalder",&normaldervar);
      if (rc != FC_SUCCESS){
        fprintf(fpout,"#error:cant get normalform for this var\n");
        fflush(NULL);
        continue;
      }

      
      fprintf(fpder,"#%d\n",ii);
      fflush(NULL);
      for (i = 0; i < numStepIntersect; i++){
        rc = fc_getVariableDataPtr(normaldervar[i],(void**)&derivativedata);
        if (rc != FC_SUCCESS){
          fprintf(fpout,"#error:cant get variable data\n");
          fflush(NULL);
          break;
        }

        rc = fc_getVariableDataPtr(normalseqvar[i],(void**)&seqvardata);
        if (rc != FC_SUCCESS){
          fprintf(fpout,"#error:cant get variable data\n");
          fflush(NULL);
          break;
        }
        //know its one data point one component
        fprintf(fpder,"%10g %10g %10g\n",coords[i],seqvardata[0],derivativedata[0]);
        fflush(NULL);
      }
      fprintf(fpder,"\n\n");
      fflush(NULL);

      rc = fc_deleteSeqVariable(numStepIntersect,derivativevar);
      rc = fc_deleteSeqVariable(numStepIntersect,normaldervar);
      rc = fc_deleteSeqVariable(numStepIntersect,normalseqvar);
      free(derivativevar);
      free(normaldervar);
      free(normalseqvar);
  }
  
  fclose(fpder);
}

*/


static void writeoutChart(int kk, char* name1, FC_Sequence seqintersect, int numSeq, Stats* svstats){
  //now output chart of quantities - write them in order of increasing klein dist
  //to the ref one

  FC_ReturnCode rc;
  int i, j, tmp;
  int order[numSeq];
  double* coords;

  //get order
  for (i = 0; i < numSeq; i++){
    order[i] = i;
  }
  
  //  printf("sorting....");
  for (i=0; i < numSeq-1 ; i++) {
    for (j=0; j < numSeq-1-i; j++){
      if (fabs(svstats[order[j+1]].kleindist)
          < fabs(svstats[order[j]].kleindist)) {  /* compare the two neighbors */
        tmp = order[j];         /* swap a[j] and a[j+1]      */
        order[j] = order[j+1];
        order[j+1] = tmp;
      }
    }
  }
  //  printf("done\n");


  rc = fc_getSequenceCoordsPtr(seqintersect,(void**)&coords);
  fc_exitIfErrorPrintf(rc, "failed to get seq coords");      
  
 
  //writeout
  fprintf(fpout,"#%14s: %10s %10s %10s %10s %10s %10s %10s\n",
          "NAME", "MIN", "MINTIME", "MAX", "MAXTIME", "INTEGRAL",
          "EUCDIST", "KLEINDIST");
  printf("#%14s: %10s %10s %10s %10s %10s %10s %10s\n",
         "NAME", "MIN", "MINTIME", "MAX", "MAXTIME", "INTEGRAL",
         "EUCDIST", "KLEINDIST");
  
  fprintf(fpout,"REF= %10s: %10g %10g %10g %10g %10g\n",
          svstats[kk].name, svstats[kk].min, coords[svstats[kk].minindex],
          svstats[kk].max, coords[svstats[kk].maxindex],svstats[kk].integral);
  printf("REF= %10s: %10g %10g %10g %10g %10g\n",
         svstats[kk].name, svstats[kk].min, coords[svstats[kk].minindex],
         svstats[kk].max, coords[svstats[kk].maxindex],svstats[kk].integral);
  fflush(NULL);
  for (i = 0; i <numSeq; i++){
    int ii = order[i];
    if (strcmp(svstats[ii].name,name1)!=0){//same as if ii DNE kk
      fprintf(fpout,"%15s: %10g %10g %10g %10g %10g %10g %10g\n", 
              svstats[ii].name, svstats[ii].min, coords[svstats[ii].minindex],
              svstats[ii].max, coords[svstats[ii].maxindex], svstats[ii].integral,
              svstats[ii].eucdist, svstats[ii].kleindist);
      printf("%15s: %10g %10g %10g %10g %10g %10g %10g\n",
             svstats[ii].name, svstats[ii].min, coords[svstats[ii].minindex],
             svstats[ii].max, coords[svstats[ii].maxindex], svstats[ii].integral,
             svstats[ii].eucdist, svstats[ii].kleindist);
      fflush(NULL);
    }
  }
  fprintf(fpout,"\n");
  printf("\n");
  fflush(NULL);
}


//------------------------------------------MAIN---------------------------------//
int main(int argc, char** argv) {
  FC_ReturnCode rc;
  int i, ii, kk;
  int ret;

  //seqvars
  int numSeq;
  FC_Variable **seqVar;
  int *numSeqStep;
  FC_Sequence *seq;

  //looping
  int validbounds = 1;
  double tmin = 0.;
  double tmax = 0.;
  //  int iteration = 0;

  char filein[maxchar],fileout[maxchar],fileder[maxchar],filegnu[maxchar], basedir[maxchar];

  //get filename
  if (argc != 3) {
    printf("usage: seriesAnalyzer inputDir inputFile\n");
    fflush(NULL);
    exit(-1);
  }

  //-----------------------------FILES---------------------------------------//
  setupFiles(argv[1], argv[2], basedir, filein, fileout,fileder,filegnu);

  //------------------------- FC SETUP--------------------------------------//
  rc = fc_setLibraryVerbosity(FC_ERROR_MESSAGES);
  rc = fc_initLibrary();
  fc_exitIfErrorPrintf(rc, "failed to initialize fc library");

  
  glob_meshDataPoint =1;

  setupDataset(3); //set up dataset, mesh etc. these vals dont matter

  //------------------------ORIG DATA--------------------------------------//
  //get number of datafiles

  fscanf(fp,"%d\n",&numSeq);
  if (numSeq < 1){
    fc_printfErrorMessage("invalid numseq");
    return 0;
  }
    
  //allocate out here so not passing *** 
  seq = (FC_Sequence*)malloc(numSeq*sizeof(FC_Sequence));
  numSeqStep = (int*)malloc(numSeq*sizeof(int));
  seqVar = (FC_Variable**)malloc(numSeq*sizeof(FC_Variable*));

  i = readData(basedir,numSeq,numSeqStep,seq,seqVar);
  if (i == 0){
    fc_printfErrorMessage("prob reading data - exiting" );
    return 0;
  }

  //--------------------WRITE OUT ORIG DATA  --------------------------//
  //dont care if there is an error in writeout
  ret = origWriteout(numSeq,numSeqStep,seq,seqVar);
  //07/21/05 ACG removing gnu option for now
  //  ret = origWriteout(filein,numSeq,numSeqStep,seq,seqVar);
  //  fclose(fpgnu);


  //---------------------------------ANALYSES LOOP-------------------------------//
  //now loop thru - analyses are repeatable over different ranges
  while (validbounds){
  //analyses
    int numStepIntersect;
    FC_Sequence seqintersect;
    FC_Variable **seqVarIntersect, **seqVarIntersectNormal;
    FC_Variable **seqVarCompare;
    Stats* svstats;

    seqVarIntersect = (FC_Variable**)malloc(numSeq*sizeof(FC_Variable*));
    seqVarIntersectNormal = (FC_Variable**)malloc(numSeq*sizeof(FC_Variable*));
    seqVarCompare = (FC_Variable**)malloc(numSeq*sizeof(FC_Variable*));
    svstats = (Stats*)malloc(numSeq*sizeof(Stats));
    for (i = 0; i < numSeq; i++){
      statsInit(&svstats[i]);
    }

    //----------------PUT THEM ON SAME SEQ FOR COMPARISONS---------------------//
    ret = makeIntersecting(tmin, tmax, numSeq, numSeqStep,seq, seqVar, &numStepIntersect,
                     &seqintersect, seqVarIntersect);

    if (ret == 0){
      printf("ERROR: cant process these bounds. Continue with new bounds.\n");
      fprintf(fpout,"#ERROR: cant process these bounds.\n");
      fflush(NULL);
    } else {
      //makes normal forms of the intersecting seqvars
      ret = makeNormal(numStepIntersect,numSeq,seqVarIntersect,
                       seqVarIntersectNormal);
      if (ret == 0){
        printf("ERROR: cant process these bounds. Continue with new bounds.\n");
        fprintf(fpout,"#ERROR: cant process these bounds.\n");
        fflush(NULL);
        exit(-1);
      } else {
        
        //decide which set we are going to use
        for (i = 0; i < numSeq; i++){
          seqVarCompare[i] = seqVarIntersectNormal[i];
        }

        //---------------WRITEOUT INTERSECTING SEQS---------------------//
        writeoutCompare(numStepIntersect, seqintersect, numSeq, seqVarCompare);

        //---------------------------------ANALYSES----------------------------------//

        //get min and max and integral on the intersecting range:
        doSingles(numStepIntersect,numSeq,seqVarCompare,&svstats); //sinlges

        /*
        //FIX  - pass normal in now that im calculating it
        if (iteration == 0){
          //only do this the first time. it closes fpder
          writeoutNormal(numStepIntersect,seqintersect, numSeq, seqVarIntersect);
          iteration++;
        }
        */
        
        //now do euclidean distance and dimensionless area
        //take each one as the reference one in turn:
        for (kk = 0; kk < numSeq; kk++){
          char* name1;
          rc = fc_getVariableName(seqVarCompare[kk][0],&name1);
          fc_exitIfErrorPrintf(rc, "failed to get var name");
          
          for (ii = 0; ii < numSeq; ii++){
            if (ii == kk) continue;
            doDoubles(numStepIntersect,seqVarCompare[kk],seqVarCompare[ii],
                      &(svstats[ii])); //comparisons
          }
        
          //----------------------WRITEOUT CHART----------------------------------//
          writeoutChart(kk, name1, seqintersect, numSeq, svstats);
        
          //clean up for this iteration of ref variable
          free(name1);
        }
      }


      //clean up for all iterations
      for (i = 0; i< numSeq; i++){
        fc_deleteSeqVariable(numStepIntersect,seqVarIntersect[i]);
        fc_deleteSeqVariable(numStepIntersect,seqVarIntersectNormal[i]);
        free(seqVarIntersect[i]);
        free(seqVarIntersectNormal[i]);
        free (svstats[i].name);
        //      free (svstats[i].refname);
      }
      //inner one is just a pointer
      //do we need to delete the seq var ?
      free(seqVarCompare);
      
      fc_deleteSequence(seqintersect);
    }
    free (seqVarIntersect);
    free (seqVarIntersectNormal);
    free (svstats);

    //----------------------GET NEW BOUNDS-----------------------------------------------//
    validbounds = 0;
    while (!validbounds){
      printf("Enter new Tmin and Tmax. Use 0 0 to do calc over entire range. Use -1 -1 to exit.\n");
      //      printf("You can load %s in gnuplot to see the full range of data curves.\n",filegnu);
      fflush(NULL);
      scanf("%lf %lf",&tmin,&tmax);
      if (tmin < tmax || (FC_DBL_EQUIV(tmin,0.0) && FC_DBL_EQUIV(tmax,0.0))){
        validbounds = 1;
      } else {
        if (FC_DBL_EQUIV(tmin,-1.0) && FC_DBL_EQUIV(tmax,-1.0)){
          printf("Exiting\n");
          fflush(NULL);
          break; //this should break out of while
        } else {
          printf("Invalid bounds\n");
          fflush(NULL);
        }
      }
    }
    if (!validbounds){
      //exiting
      break;
    }
  }

  //final clean up
  fclose(fpout);

  for (i = 0; i< numSeq; i++){
    fc_deleteSeqVariable(numSeqStep[i],seqVar[i]);
    free(seqVar[i]);
    fc_deleteSequence(seq[i]);
  }

  free (numSeqStep);
  free (seq);
  free (seqVar);

  fc_deleteMesh(glob_mesh);
  fc_deleteDataset(glob_dataset);

  fc_finalLibrary();
  exit(0);

}
