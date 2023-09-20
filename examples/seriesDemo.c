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
 * \file seriesDemo.c
 * \brief demo on how to do shiftAndScale of SeqVars
 *        and some analyses.
 *
 * $Source: /usr/local/Repositories/fcdmf/fclib/examples/seriesDemo.c,v $
 * $Revision: 1.9 $
 * $Date: 2006/10/19 03:14:49 $
 *
 * \description
 *
 *   Demo of doing shift and scale on some seqvars and
 *   some analyses. Writes gnuplotable results into
 *   seriesDemo.out and gnuplot file seriesDemo.gnu.
 *   Run output file by starting gnuplot and loading
 *   seriesDemo.gnu.
 *
 *   Suggested usage:
 *     1) run the demo and save STDOUT to a file: ./seriesDemo > output
 *     2) start gnuplot: gnuplot 
 *     3) in gnuplot: load "seriesDemo.gnu"
 *     4) read the instruction as they come up on screen to
 *        see what the demo is about
 *     5) then look through the seriesDemo.c code and STDOUT file
 *         
 *
 * \modifications
 *   - 06/14/05 ACG Created.
 *   - 01/20/06 ACG changed data to use sin wave, added eucdist
 *     analysis as well, made output more verbose for gnuplot,
 *     and enhanced the output to STDOUT
 */

#include <fc.h>
#include <math.h>
#include <string.h>

//int main(int argc, char** argv) {
int main(void){
  FC_ReturnCode rc;
  FC_Dataset dataset;
  
  FC_Sequence *seqs, *compseqs, seqintersect, ss_seq;
  int numSeqs = 2;
  int numStepIntersect;
  FC_Mesh mesh;
  FC_Variable **seqvars, **compseqvars, **seqvarsIntersect, *ss_seqvar;
  int numSeqVars = numSeqs;

  int numStep = 9;
  //seq 1 is seq 0 shifted by 5 and with twice the frequency
  double seqvals0[9] = {0,1,2,3,4,5,6,7,8};
  double seqvals1[9] = {5.0,5.5,6.0,6.5,7.0,7.5,8.0,8.5,9.0};
  double datavals[9];
  void* seqcoords;

  FC_Variable areavar,eucvar;
  double *areadata, *eucdata;
  double results[3][2];
  int i,j,casex;

  FILE *fpout, *fpgnu; 
  char *fileout = "seriesDemo.out";
  char *filegnu = "seriesDemo.gnu";


  //get datavals for a sin wave
  for (i = 0; i < numStep; i++){
    datavals[i] = sin(seqvals0[i]*2*FC_PI/(double)(numStep-1))+3;
  }

  rc = fc_setLibraryVerbosity(FC_WARNING_MESSAGES);
  if (rc != FC_SUCCESS)
    exit(rc);
  rc = fc_initLibrary();
  if (rc != FC_SUCCESS)
    exit(rc);

  // make a dataset to play in
  rc = fc_createDataset("temp dataset", &dataset);
  if (rc != FC_SUCCESS){
    fc_printfErrorMessage("abort: failed to create dataset");
    exit(FC_ERROR);
  }
  rc = fc_createMesh(dataset,"mesh",&mesh);
  if (rc != FC_SUCCESS){
    fc_printfErrorMessage("abort: failed to create mesh");
    exit(FC_ERROR);
  }

  //now make original sequences and seqvars
  //seq 1 is seq 0 shifted and scaled, so the tests cases will
  // unshift and unscale it and see what happens.
  seqs = (FC_Sequence*)malloc(numSeqs*sizeof(FC_Sequence));
  compseqs = (FC_Sequence*)malloc(numSeqs*sizeof(FC_Sequence));
  seqvars = (FC_Variable**)malloc(numSeqVars*sizeof(FC_Variable*));
  compseqvars = (FC_Variable**)malloc(numSeqVars*sizeof(FC_Variable*));
  seqvarsIntersect = (FC_Variable**)malloc(numSeqVars*sizeof(FC_Variable*));

  for (i = 0; i < numSeqs; i++){
    char *name,*seqname;
    double* scoords;
    name = (i == 0 ? "sv0":"sv1");
    seqname = (i == 0 ? "s0":"s1");
    scoords = (i == 0 ? seqvals0:seqvals1);
    rc = fc_createSequence(dataset, seqname, &(seqs[i]));
    if (rc != FC_SUCCESS){
      fc_printfErrorMessage("abort: failed to create sequence");
      exit(FC_ERROR);
    }
    rc = fc_setSequenceCoords(seqs[i], numStep, FC_DT_DOUBLE,
				 (void*)scoords);
    if (rc != FC_SUCCESS){
      fc_printfErrorMessage("abort: failed to set sequence coords");
      exit(FC_ERROR);
    }
    rc = fc_createSeqVariable(mesh,seqs[i],name,&numStep,&(seqvars[i]));
    if (rc != FC_SUCCESS){
      fc_printfErrorMessage("abort: failed to create seqvar");
      exit(FC_ERROR);
    }
    for (j = 0; j < numStep; j++){
      rc = fc_setVariableData(seqvars[i][j],1,1,FC_AT_WHOLE_MESH,
			      FC_MT_SCALAR,FC_DT_DOUBLE,
			      &datavals[j]);
      if (rc != FC_SUCCESS){
	fc_printfErrorMessage("abort: failed to set seqvar coords");
	exit(FC_ERROR);
      }
    }
  }

  fpout = fopen(fileout,"w");

  //print orig sequences & write to file
  printf("ORIGINAL SEQS & VARS:\n");
  fprintf(fpout,"# ORIGINAL SEQS & VARS:\n");
  for (i = 0; i < numSeqVars; i++){
    rc = fc_printSequence(seqs[i],"ORIG SEQ",1);
    rc = fc_printSeqVariable(numStep,seqvars[i],"ORIG SEQVAR",1);
    rc = fc_getSequenceCoordsPtr(seqs[i],&seqcoords);
    fc_exitIfErrorPrintf(rc, "failed to get seq coords");
    fprintf(fpout,"#\tsv%d: (time val)\n",i);
    for (j = 0; j < numStep; j++){
      double *dataa;
      rc = fc_getVariableDataPtr(seqvars[i][j],(void**)&dataa);
      fc_exitIfErrorPrintf(rc, "failed to get seq var vals");
      fprintf(fpout,"%10g %10g\n",((double*)seqcoords)[j],dataa[0]);
    }
    printf("\n");
    fprintf(fpout,"\n\n");
  }


  // examine 3 cases: 1) direct diff 2) shift and diff 3) shift and scale and diff
  // (case 3 will be identical seqvars)
  for (casex = 0; casex <3; casex++){
    printf("\n\n------------------------------------------------------\n");
    printf("CASE %d:\n",casex);
    fprintf(fpout,"#CASE %d:\n",casex);

    //first get the starting seq vars for the case:
    compseqvars[0] = seqvars[0];
    compseqs[0] = seqs[0];
    switch (casex){
    case 0:
      //do nothing - direct diff
      printf("This case uses the original seqvars.\n");
      compseqvars[1] = seqvars[1];
      compseqs[1] = seqs[1];
      break;
    default:
      {
	double shift,scale;
	shift = -5.0;
	scale =  (casex == 2? 2.0 : 1.0);
	printf("This case shifts the second seqvar by %g and scales it by %g.\n",
	       shift, scale);
	// shift and scale seqvar 1:
	rc = fc_shiftAndScaleSequence(seqs[1],dataset,shift,scale,"ss_seq",&ss_seq);
	fc_exitIfErrorPrintf(rc, "failed to get shifted and scaled sequence");
	//now get new seq var on the new seq
	rc = fc_copySeqVariable (numStep,seqvars[1],mesh,ss_seq, "ss_seqvar",
				 &ss_seqvar);
	fc_exitIfErrorPrintf(rc, "failed to get shifted and scaled sequence var");
	compseqvars[1] = ss_seqvar;
	compseqs[1] = ss_seq;
	break;
      }
    }

    //display seq vars used for this case & write to file
    printf("\n\ninitial seqs & seqvars for this case:\n");
    fprintf(fpout,"#initial seqs & seqvars for this case:\n");
    for (i = 0; i < numSeqVars; i++){
      rc = fc_printSequence(compseqs[i],"ORIG SEQ",1);
      rc = fc_printSeqVariable(numStep,compseqvars[i],"ORIG SEQVAR",1);
      rc = fc_getSequenceCoordsPtr(compseqs[i],&seqcoords);
      fc_exitIfErrorPrintf(rc, "failed to get seq coords");
      fprintf(fpout,"#\tsv%d: (time val)\n",i);
      for (j = 0; j < numStep; j++){
	double *dataa;
	rc = fc_getVariableDataPtr(compseqvars[i][j],(void**)&dataa);
	fc_exitIfErrorPrintf(rc, "failed to get seq var vals");
	fprintf(fpout,"%10g %10g\n",((double*)seqcoords)[j],dataa[0]);
      }
      printf("\n");
      fprintf(fpout,"\n\n");
    }
    
    //now need to put them on the same seq in order to do the analyses
    //frist get intersecting sequence
    rc = fc_createIntersectingRegularSequence(numSeqs,compseqs,0,
					      "seqintersect", &seqintersect);
    fc_exitIfErrorPrintf(rc, "failed to get intersecting sequence");
    //    printf("\n\nseq after intersection:\n"); 
    //    rc = fc_printSequence(seqintersect,"seqintersect",1);
    rc = fc_getSequenceNumStep(seqintersect,&numStepIntersect);
    fc_exitIfErrorPrintf(rc, "failed to get numsteps on new sequence");
    
    //then put them all on the same seq and display them
    printf("\n\norig seqvars for this case interpolated onto the same sequence:\n");
    for (i = 0; i < numSeqVars; i++){
      char *name1,*name2;

      //now do linear interpolation 
      name1 = NULL;
      rc = fc_getVariableName(compseqvars[i][0],&name1);
      fc_exitIfErrorPrintf(rc, "failed to get name of comp seqvar");
      name2 = malloc(sizeof(char)*strlen(name1)+20);
      strcpy(name2,name1);
      strcat(name2,"_newcase");
      rc = fc_linearInterpolation(numStep, compseqvars[i], seqintersect,
				  name1,&seqvarsIntersect[i]);
      fc_exitIfErrorPrintf(rc, "failed to create linear interpolation of seqvar");
      rc = fc_printSequence(seqintersect,"COMMON SEQUENCE",1);
      fc_printSeqVariable(numStepIntersect,compseqvars[i],
			  (i == 0 ? "INTERPOLATED REF SEQVAR" : "INTERPOLATED COMP SEQVAR"),1);
      free(name1);
      free(name2);
    }


    //writing them out to output file as well in gnuplottable format
    fprintf(fpout,"#comp seqs & seqvars for this case:\n");
    rc = fc_getSequenceCoordsPtr(seqintersect,&seqcoords);
    fc_exitIfErrorPrintf(rc, "failed to get seq coords");
    for (j = 0; j < numStepIntersect; j++){
      double *dataa, *datab;
      rc = fc_getVariableDataPtr(seqvarsIntersect[0][j],(void**)&dataa);
      fc_exitIfErrorPrintf(rc, "failed to get seq var vals");
      rc = fc_getVariableDataPtr(seqvarsIntersect[1][j],(void**)&datab);
      fc_exitIfErrorPrintf(rc, "failed to get seq var vals");
      fprintf(fpout,"%10g %10g %10g\n",((double*)seqcoords)[j],dataa[0],datab[0]);
    }
    fprintf(fpout,"\n\n");

    //now do dimensionless area and eucdistance
    rc =  fc_dimensionlessAreaBetweenCurves(numStepIntersect,
						 seqvarsIntersect[0],
						 seqvarsIntersect[1],
						 "area", &areavar);
    fc_exitIfErrorPrintf(rc, "failed to calc dimensionless area");
    rc = fc_getVariableDataPtr(areavar,(void**)&areadata); //know this is a double
    fc_exitIfErrorPrintf(rc, "failed to get area data");
    rc =  fc_euclideanDistanceBetweenCurves(numStepIntersect,
						 seqvarsIntersect[0],
						 seqvarsIntersect[1],
						 "euc", &eucvar);
    fc_exitIfErrorPrintf(rc, "failed to calc euc distance");
    rc = fc_getVariableDataPtr(eucvar,(void**)&eucdata); //know this is a double
    fc_exitIfErrorPrintf(rc, "failed to get area data");
    rc=fc_getSequenceCoordsPtr(seqintersect,&seqcoords);
    fc_exitIfErrorPrintf(rc, "failed to get seq coords");
    //area data and eucdata have only 1 data point since the original seq vars 
    //only had 1 data point per timestep
    results[casex][0] = areadata[0]; //for gnuplot file    
    results[casex][1] = eucdata[0]; //for gnuplot file    
    printf("\n\ncase %d: DP = %d range = [%4g,%4g] dimensionless area  = %10g eucdistance = %10g\n\n",
	   casex, 
	   numStepIntersect,
	   ((double*)seqcoords)[0],
	   ((double*)seqcoords)[numStepIntersect-1],
	   areadata[0],eucdata[0]);
    fprintf(fpout,"#RESULT: case %d: DP = %d range = [%4g,%4g] dimensionlessarea = %10g eucdistance = %10g\n\n\n",
	    casex,numStepIntersect,
	   ((double*)seqcoords)[0],
	   ((double*)seqcoords)[numStepIntersect-1],
	   areadata[0],eucdata[0]);

    //clean up
    rc = fc_deleteVariable(areavar);
    fc_exitIfErrorPrintf(rc, "failed to clean up area var");
    for (i = 0; i< numSeqVars; i++){
      rc = fc_deleteSeqVariable(numStepIntersect,seqvarsIntersect[i]);
      fc_exitIfErrorPrintf(rc, "failed to clean up seqvars intersect");
      free(seqvarsIntersect[i]);
    }

    rc = fc_deleteSequence(seqintersect);
    fc_exitIfErrorPrintf(rc, "failed to clean up seq intersect");

    //clean up if we used a shiftandscale
    if (casex !=0){
      rc = fc_deleteSeqVariable(numStep,ss_seqvar);
      fc_exitIfErrorPrintf(rc, "failed to clean up ss seqvar");
      free(ss_seqvar);
      rc = fc_deleteSequence(ss_seq);
      fc_exitIfErrorPrintf(rc, "failed to clean up ss seq");
    }
  }

  {
    //writing out gnuplot file
    fpgnu = fopen(filegnu,"w");

    //first print orig seqs
    fprintf(fpgnu,"set title \"ORIG SEQVARS\"\n");
    fprintf(fpgnu,"plot \"%s\" index 0 title \"%s\" with lines\n",
	    fileout,"ORIG SV0");
    fprintf(fpgnu,"replot \"%s\" index 1 title \"%s\" with lines\n",
	      fileout,"ORIG SV1");
    fprintf(fpgnu,"pause 0 \"This example demonstrates series analyses. Analyses for\" \n");
    fprintf(fpgnu,"pause 0 \"3 different pairs of seq vars are considered.\" \n");
    fprintf(fpgnu,"pause 0 \"\" \n");
    fprintf(fpgnu,"pause 0 \"Two sequence variables are shown in the plot in red and green.\" \n");
    fprintf(fpgnu,"pause 0 \"For each case we use the red seqvar and a variation on the green\" \n");
    fprintf(fpgnu,"pause 0 \"seqvar, in order to see how the variation alters the results\" \n");
    fprintf(fpgnu,"pause 0 \"of the analyses. The altered version of the green curve will be\" \n");
    fprintf(fpgnu,"pause 0 \"shown in blue.\" \n");
    fprintf(fpgnu,"pause 0 \"\" \n");
    fprintf(fpgnu,"pause 0 \"We want to do the analyses using the red and blue seqvars.\" \n");
    fprintf(fpgnu,"pause 0 \"However, in order to do the analyses, the two seq vars \" \n");
    fprintf(fpgnu,"pause 0 \"must be on the same sequence. Therefore, we find the\" \n");
    fprintf(fpgnu,"pause 0 \"overlapping time region of the two seqvars and interpolate\" \n");
    fprintf(fpgnu,"pause 0 \"the seqvars to make new seq vars over this region. The\" \n");
    fprintf(fpgnu,"pause 0 \"resulting seqvars that are actually used for the analyses will\" \n");
    fprintf(fpgnu,"pause 0 \"be shown in magenta and teal, for the red and green seqvars,\" \n");
    fprintf(fpgnu,"pause 0 \"respectively.\" \n"); 
    fprintf(fpgnu,"pause 0 \"\" \n");
    fprintf(fpgnu,"pause 0 \"Results of the dimensionlessArea and euclideanDistance \" \n");
    fprintf(fpgnu,"pause 0 \"analyses for each case will be shown in the title of the plot.\" \n");
    fprintf(fpgnu,"pause 0 \"Note that the dimensionlessArea result is dependent on the\" \n");
    fprintf(fpgnu,"pause 0 \"max value of the reference function over that range.\" \n");
    fprintf(fpgnu,"pause 0 \"Also note that the euclidean distance is dependent upon.\" \n");
    fprintf(fpgnu,"pause 0 \"the number of datapoints in the range. Therefore, these values\" \n");
    fprintf(fpgnu,"pause 0 \"should be considered with caution.\" \n");
    fprintf(fpgnu,"pause 0 \"\" \n");
    fprintf(fpgnu,"pause -1 \"Hit return to see the cases...\" \n");

    for (i = 0; i < 3; i++){
      char* casestring;
      fprintf(fpgnu,"pause 0\"CASE %d:\" \n",i);
      fprintf(fpgnu,"pause 0 \"\" \n");

      fprintf(fpgnu,"set title \"CASE %d\" \n",i);
      fprintf(fpgnu,"plot \"%s\" index 0 title \"%s\" with lines\n",
	      fileout,"ORIG SV0");
      fprintf(fpgnu,"replot \"%s\" index 1 title \"%s\" with lines\n",
	    fileout,"ORIG_SV1");

      fprintf(fpgnu,"pause 0 \"First we see the original red and green seqvars.\" \n");
      fprintf(fpgnu,"pause -1 \"Hit return to show the variation on the green seqvar...\"\n");
      fprintf(fpgnu,"pause 0 \"\" \n");

      switch(i){
      case 0:
	casestring = "THISCASE_SV1 = ORIG_SV1";
	break;
      case 1:
	casestring = "THISCASE_SV1 = ORIG_SV1 shifted by -5";
	break;
      case 2:
	casestring = "THISCASE_SV1 = ORIG_SV1 shifted by -5 and scaled by 2";
	break;
      default:
	casestring = "";
	break;
      }

      fprintf(fpgnu,"replot \"%s\" index %d title \"%s\" with lines\n",
	    fileout,3+3*i,casestring);

      fprintf(fpgnu,"pause 0 \"The variation on the green seqvar is now shown in blue.\"\n");
      fprintf(fpgnu,"pause -1 \"Hit return to put vars on intersecting seq...\"\n");
      fprintf(fpgnu,"pause 0 \"\" \n");

      fprintf(fpgnu,"set title \"CASE %d - DIM AREA = %10g EUC DIST = %10g\" \n",i,results[i][0],results[i][1]);
      fprintf(fpgnu,"replot \"%s\" index %d using 1:2 title \"%s\" with lines\n",
	    fileout,4+3*i,"CASE_INTSV0");
      fprintf(fpgnu,"replot \"%s\" index %d using 1:3 title \"%s\" with lines\n",
	      fileout,4+3*i,"CASE_INT1SV1");

      fprintf(fpgnu,"pause 0 \"The red and blue curves are interpolated over\"\n");
      fprintf(fpgnu,"pause 0 \"the overlapping time region. The resulting curves\"\n");
      fprintf(fpgnu,"pause 0 \"are shown in magneta (red) and teal (blue).\"\n");
      fprintf(fpgnu,"pause 0 \"\" \n");

      fprintf(fpgnu,"pause 0 \"Values for the analyses based on these curves are\"\n");
      fprintf(fpgnu,"pause 0 \"given in the plot title.\"\n");
      fprintf(fpgnu,"pause 0 \"\" \n");
      if (i !=2 ){
	fprintf(fpgnu,"pause -1 \"Hit return to move to case %d...\"\n",i+1);
	fprintf(fpgnu,"pause 0 \"\" \n");
      }
    }
    fprintf(fpgnu,"pause 0 \"Finished.\" \n");
  }


  printf("FINAL CLEANUP\n");
  //clean up

  fclose(fpout);
  free(seqvarsIntersect);

  for (i = 0; i< numSeqVars; i++){
    rc = fc_deleteSeqVariable(numStep,seqvars[i]);
    fc_exitIfErrorPrintf(rc, "failed to clean up seqvar");
    rc = fc_deleteSequence(seqs[i]);
    fc_exitIfErrorPrintf(rc, "failed to clean up seq");
    free(seqvars[i]);
  }

  free(seqvars);
  free(seqs);
  free(compseqvars); // just pointers
  free(compseqs); // just pointers

  rc = fc_deleteMesh(mesh);
  fc_exitIfErrorPrintf(rc, "failed to clean up mesh");
  rc = fc_deleteDataset(dataset);
  fc_exitIfErrorPrintf(rc, "failed to clean up dataset");

  rc = fc_finalLibrary();
  fc_exitIfErrorPrintf(rc, "failed to finalize library");
  exit(0);
}
