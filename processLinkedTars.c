#include "log.h"
#include "format.h"
#include "linestream.h"
#include "transcriptAssemblyConfig.h"
#include <stdio.h>



typedef struct {
  char* chromosome;
  char strand;
  Array starts;
  Array ends;
  int start;
  int end;
  int transcriptLength;
  int numExons;
  double averageNumSupportingReads;
  double averageTarExpression;
  int numRecursions;
} Transcript;



typedef struct {
  Transcript *transcript;
  double averageNumSupportingReads;
  double averageTarExpression;
} TranscriptClusterEntry;



static void printTranscriptClusterEntries (FILE *fp1, FILE*fp2, Array transcriptClusterEntries)
{
  int i,j;
  int stop;
  static int transcriptNumber = 1;
  static int isoformNumber = 1;
  TranscriptClusterEntry *currTCE;
    
  stop = MIN (arrayMax (transcriptClusterEntries),MAX_NUM_ISOFORMS);
  for (i = 0; i < stop; i++) {
    currTCE = arrp (transcriptClusterEntries,i,TranscriptClusterEntry);
    fprintf (fp1,"%d\tInterval_%d_%.1f_%.1f\n",
             isoformNumber,
             transcriptNumber,
             currTCE->transcript->averageNumSupportingReads,
             currTCE->transcript->averageTarExpression);
    fprintf (fp2,"Interval_%d_%.1f_%.1f\t%s\t%c\t%d\t%d\t%d\t",
             transcriptNumber,
             currTCE->transcript->averageNumSupportingReads,
             currTCE->transcript->averageTarExpression,
             currTCE->transcript->chromosome,
             currTCE->transcript->strand,
             arru (currTCE->transcript->starts,0,int),
             arru (currTCE->transcript->ends,arrayMax (currTCE->transcript->ends) - 1,int),
             arrayMax (currTCE->transcript->starts));
    for (j = 0; j < arrayMax (currTCE->transcript->starts); j++) {
      fprintf (fp2,"%d%s",arru (currTCE->transcript->starts,j,int),j < arrayMax (currTCE->transcript->starts) - 1 ? "," : "\t");
    }
    for (j = 0; j < arrayMax (currTCE->transcript->ends); j++) {
      fprintf (fp2,"%d%s",arru (currTCE->transcript->ends,j,int),j < arrayMax (currTCE->transcript->ends) - 1 ? "," : "\n");
    }
    transcriptNumber++;
  }
  isoformNumber++;
}



static void clearTranscripts (Array transcripts) 
{
  Transcript *currTranscript;
  int i;

  for (i = 0; i < arrayMax (transcripts); i++) {
    currTranscript = arrp (transcripts,i,Transcript);
    hlr_free (currTranscript->chromosome);
    arrayDestroy (currTranscript->starts);
    arrayDestroy (currTranscript->ends);
  }
  arrayClear (transcripts);
}



static int sortTranscripts (Transcript *a, Transcript *b) 
{
  int diff;

  diff = a->strand - b->strand;
  if (diff != 0) {
    return diff;
  } 
  diff = a->start - b->start;
  if (diff != 0) {
    return diff;
  }
  return b->end - a->end;
}



static int sortTranscriptClusterEntries (TranscriptClusterEntry *a, TranscriptClusterEntry *b)
{
  double diff;

  diff = b->transcript->averageNumSupportingReads - a->transcript->averageNumSupportingReads;
  if (diff < 0) {
    return -1;
  }
  if (diff > 0) {
    return 1;
  }
  diff = b->transcript->averageTarExpression - a->transcript->averageTarExpression;
  if (diff < 0) {
    return -1;
  }
  if (diff > 0) {
    return 1;
  }
  return 0;
}



static void processTranscripts (FILE *fp1, FILE*fp2, Array transcripts)
{
  Transcript *currTranscript,*nextTranscript;
  int i,j;
  TranscriptClusterEntry *currTCE;
  static Array transcriptClusterEntries = NULL;

  if (transcriptClusterEntries == NULL) {
    transcriptClusterEntries = arrayCreate (1000,TranscriptClusterEntry);
  }
  arraySort (transcripts,(ARRAYORDERF)sortTranscripts);
  i = 0;
  while (i < arrayMax (transcripts)) {
    currTranscript = arrp (transcripts,i,Transcript);
    arrayClear (transcriptClusterEntries);
    currTCE = arrayp (transcriptClusterEntries,arrayMax (transcriptClusterEntries),TranscriptClusterEntry);
    currTCE->transcript = currTranscript;
    j = i + 1;
    while (j < arrayMax (transcripts)) {
      nextTranscript = arrp (transcripts,j,Transcript);
      if (currTranscript->strand == nextTranscript->strand &&
          currTranscript->start <= nextTranscript->start && 
          currTranscript->end >= nextTranscript->end) { 
        currTCE = arrayp (transcriptClusterEntries,arrayMax (transcriptClusterEntries),TranscriptClusterEntry);
        currTCE->transcript = nextTranscript;
      }
      else {
        break;
      }
      j++;
    }
    arraySort (transcriptClusterEntries,(ARRAYORDERF)sortTranscriptClusterEntries);
    printTranscriptClusterEntries (fp1,fp2,transcriptClusterEntries);
    i = j;
  }
  clearTranscripts (transcripts); 
}



int main (int argc, char *argv[]) 
{
  LineStream ls;
  char *line;
  Texta tokens;
  int i;
  int start,end;
  Array transcripts;
  Transcript *currTranscript;
  Stringa buffer;
  FILE *fp1,*fp2;

  if (argc != 2) {
    usage ("%s <prefix>",argv[0]);
  }
  buffer = stringCreate (100);
  stringPrintf (buffer,"%s_isoformMap.txt",argv[1]);
  fp1 = fopen (string (buffer),"w");
  stringPrintf (buffer,"%s_transcriptCoordinates.interval",argv[1]);
  fp2 = fopen (string (buffer),"w");
  if (fp1 == NULL || fp2 == NULL) {
    die ("Unable to open files.");
  }
  transcripts = arrayCreate (10000,Transcript);
  ls = ls_createFromFile ("-");
  while (line = ls_nextLine (ls)) {
    if  (strEqual (line,"###")) {
      processTranscripts (fp1,fp2,transcripts);
      continue;
    }
    currTranscript = arrayp (transcripts,arrayMax (transcripts),Transcript);
    tokens = textFieldtok (line," ");
    currTranscript->numRecursions = atoi (textItem (tokens,0));
    currTranscript->averageNumSupportingReads = atof (textItem (tokens,1));
    currTranscript->averageTarExpression = atof (textItem (tokens,2));
    currTranscript->chromosome = hlr_strdup (textItem (tokens,3));
    currTranscript->strand = textItem (tokens,4)[0];
    currTranscript->starts = arrayCreate (100,int);
    currTranscript->ends = arrayCreate (100,int);
    currTranscript->transcriptLength = 0;
    for (i = 6; i < arrayMax (tokens); i = i + 2) {
      start = atoi (textItem (tokens,i - 1));
      end = atoi (textItem (tokens,i));
      currTranscript->transcriptLength += (end - start + 1);
      array (currTranscript->starts,arrayMax (currTranscript->starts),int) = start;
      array (currTranscript->ends,arrayMax (currTranscript->ends),int) = end;
    }
    currTranscript->start = arru (currTranscript->starts,0,int);
    currTranscript->end = arru (currTranscript->ends,arrayMax (currTranscript->ends) - 1,int);
    textDestroy (tokens);
    if (arrayMax (currTranscript->starts) != arrayMax (currTranscript->ends)) {
      die ("Expected same number of starts (n = %d) and ends (n = %d)!",arrayMax (currTranscript->starts),arrayMax (currTranscript->ends));
    }
    currTranscript->numExons = arrayMax (currTranscript->starts);
    currTranscript->averageNumSupportingReads = currTranscript->averageNumSupportingReads / currTranscript->numExons;
    currTranscript->averageTarExpression = currTranscript->averageTarExpression / currTranscript->numExons;
  }
  processTranscripts (fp1,fp2,transcripts);
  ls_destroy (ls);
  stringDestroy (buffer);
  fclose (fp1);
  fclose (fp2);
  return 0;
}
