#include "log.h"
#include "format.h"
#include "linestream.h"
#include "intervalFind.h"



typedef struct {
  int id;
  char *transcriptName;
} Isoform;



static Array readIsoforms (char* fileName) 
{
  LineStream ls;
  char* line;
  char* pos;
  Array isoforms;
  Isoform *currIsoform;

  isoforms = arrayCreate (30000,Isoform);
  ls = ls_createFromFile (fileName);
  while (line = ls_nextLine (ls)) {
    if (line[0] == '\0') {
      continue;
    }
    if (pos = strchr (line,'\t')) {
      *pos = '\0';
      currIsoform = arrayp (isoforms,arrayMax (isoforms),Isoform);
      currIsoform->id = atoi (line);
      currIsoform->transcriptName = hlr_strdup (pos + 1);
    }
  }
  ls_destroy (ls);
  return isoforms;
}



static int sortIsoforms (Isoform *a, Isoform *b)
{
  return a->id - b->id;
}



int main (int argc, char *argv[]) 
{
  Array isoforms;
  Array transcripts;
  Isoform *currIsoform,*nextIsoform;
  int i,j;
  int count;
  FILE *fp1,*fp2;
  Stringa buffer;
  Interval *currTranscript;
  SubInterval *currExon;
  int transcriptLength;

  if (argc != 4) {
    usage ("%s <prefix> <isoformMap.txt> <transcriptCoordinates.interval>",argv[0]);
  }
  isoforms = readIsoforms (argv[2]); 
  arraySort (isoforms,(ARRAYORDERF)sortIsoforms);
  buffer = stringCreate (100);
  stringPrintf (buffer,"%s_numIsoformsPerGene.txt",argv[1]);
  if (!(fp1 = fopen (string (buffer),"w"))) {
    die ("Unable to open file: %s",string (buffer));
  }
  i = 0;
  while (i < arrayMax (isoforms)) {
    currIsoform = arrp (isoforms,i,Isoform);
    count = 1;
    j = i + 1;
    while (j < arrayMax (isoforms)) {
      nextIsoform = arrp (isoforms,j,Isoform);
      if (currIsoform->id != nextIsoform->id) {
        break;
      }
      count++;
      j++;
    }
    i = j;
    fprintf (fp1,"%d\n",count);
  }
  fclose (fp1);
  stringPrintf (buffer,"%s_numExonsPerTranscript.txt",argv[1]);
  if (!(fp1 = fopen (string (buffer),"w"))) {
    die ("Unable to open file: %s",string (buffer));
  }
  stringPrintf (buffer,"%s_transcriptLengths.txt",argv[1]);
  if (!(fp2 = fopen (string (buffer),"w"))) {
    die ("Unable to open file: %s",string (buffer));
  }
  intervalFind_addIntervalsToSearchSpace (argv[3],0);
  transcripts = intervalFind_getAllIntervals ();
  for (i = 0; i < arrayMax (transcripts); i++) {
    currTranscript = arrp (transcripts,i,Interval);
    transcriptLength = 0;
    for (j = 0; j < arrayMax (currTranscript->subIntervals); j++) {
      currExon = arrp (currTranscript->subIntervals,j,SubInterval);
      transcriptLength += (currExon->end - currExon->start + 1);
    }
    fprintf (fp1,"%d\n",arrayMax (currTranscript->subIntervals));
    fprintf (fp2,"%d\n",transcriptLength);
  }
  fclose (fp1);
  fclose (fp2);
  stringDestroy (buffer);
  return 0;
}
