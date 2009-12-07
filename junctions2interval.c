#include "log.h"
#include "format.h"
#include "mrf.h"
#include "transcriptAssemblyConfig.h"
#include <stdio.h>



typedef struct {
  Stringa buffer;
  int junctionStart;
  int junctionEnd;
  int junctionNumber;
} Entry;



static int sortEntries (Entry *a, Entry *b)
{
  int diff;

  diff = a->junctionStart - b->junctionStart;
  if (diff != 0) {
    return diff;
  }
  return a->junctionEnd - b->junctionEnd;
}


static void processRead (Array entries, MrfRead *currMrfRead) 
{
  static int junctionNumber = 1;
  Entry *currEntry;
  MrfBlock *currMrfBlock1,*currMrfBlock2;
  
  if (arrayMax (currMrfRead->blocks) != 2) {
    return;
  }
  currMrfBlock1 = arrp (currMrfRead->blocks,0,MrfBlock);
  currMrfBlock2 = arrp (currMrfRead->blocks,1,MrfBlock);
  if (!strEqual (currMrfBlock1->targetName,currMrfBlock2->targetName) ||
      currMrfBlock1->strand != currMrfBlock2->strand) {
    die ("Alignment blocks must have same targetName and strand");
  }
  currEntry = arrayp (entries,arrayMax (entries),Entry);
  currEntry->buffer = stringCreate (200);
  currEntry->junctionStart = currMrfBlock1->targetEnd;
  currEntry->junctionEnd = currMrfBlock2->targetStart;
  currEntry->junctionNumber = junctionNumber;
  stringAppendf (currEntry->buffer,"SpliceJunction_%d\t%s\t%c\t%d\t%d\t2\t%d,%d\t%d,%d",
                 junctionNumber,currMrfBlock1->targetName,currMrfBlock1->strand,currMrfBlock1->targetStart,currMrfBlock2->targetEnd,
                 currMrfBlock1->targetStart,currMrfBlock2->targetStart,currMrfBlock1->targetEnd,currMrfBlock2->targetEnd);
  junctionNumber++;
}



int main (int argc, char *argv[]) 
{
  MrfEntry *currMrfEntry;
  Array entries;
  Entry *currEntry,*nextEntry;
  int countNumSupportingReads;
  Stringa buffer;
  FILE *fp1,*fp2;
  int i,j;
  
  if (argc != 2) {
    usage ("%s <prefix>",argv[0]);
  }
  buffer = stringCreate (100);
  stringPrintf (buffer,"%s_spliceJunctions.interval",argv[1]);
  fp1 = fopen (string (buffer),"w");
  stringPrintf (buffer,"%s_numSupportingReads.txt",argv[1]);
  fp2 = fopen (string (buffer),"w");
  if (fp1 == NULL || fp2 == NULL) {
    die ("Unable to open files!");
  }
  entries = arrayCreate (100000,Entry);
  mrf_init ("-");
  while (currMrfEntry = mrf_nextEntry ()) {
    processRead (entries,&currMrfEntry->read1);
    if (currMrfEntry->isPairedEnd) {
      processRead (entries,&currMrfEntry->read2);
    }  
  }
  mrf_deInit ();
  arraySort (entries,(ARRAYORDERF)sortEntries);
  i = 0;
  while (i < arrayMax (entries)) {
    currEntry = arrp (entries,i,Entry);
    countNumSupportingReads = 1;
    j = i + 1;
    while (j < arrayMax (entries)) {
      nextEntry = arrp (entries,j,Entry);
      if (currEntry->junctionStart == nextEntry->junctionStart && 
          currEntry->junctionEnd == nextEntry->junctionEnd) {
        countNumSupportingReads++;
      }
      else {
        break;
      }
      j++;
    }
    i = j;
    if (countNumSupportingReads >= MIN_NUM_SUPPORTING_READS_FOR_SPLICE_JUNCTION) {
      fprintf (fp1,"%s\n",string (currEntry->buffer));
      fprintf (fp2,"SpliceJunction_%d\t%d\n",currEntry->junctionNumber,countNumSupportingReads);
    }
  }
  fclose (fp1);
  fclose (fp2);
  stringDestroy (buffer);
  return 0;
}
