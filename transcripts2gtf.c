#include "log.h"
#include "format.h"
#include "linestream.h"
#include "intervalFind.h"
#include "transcriptAssemblyConfig.h"



typedef struct {
  int id;
  char* transcriptName;
} KnownIsoform;



typedef struct {
  char* transcriptName;
  double rpkm;
} Isoquant;



static Array readKnownIsoforms (char* fileName) 
{
  LineStream ls;
  char* line;
  char* pos;
  Array knownIsoforms;
  KnownIsoform *currKnownIsoform;

  knownIsoforms = arrayCreate (30000,KnownIsoform);
  ls = ls_createFromFile (fileName);
  while (line = ls_nextLine (ls)) {
    if (line[0] == '\0') {
      continue;
    }
    if (pos = strchr (line,'\t')) {
      *pos = '\0';
      currKnownIsoform = arrayp (knownIsoforms,arrayMax (knownIsoforms),KnownIsoform);
      currKnownIsoform->id = atoi (line);
      currKnownIsoform->transcriptName = hlr_strdup (pos + 1);
    }
  }
  ls_destroy (ls);
  return knownIsoforms;
}



static int sortIsoquants (Isoquant *a, Isoquant *b)
{
  return strcmp (a->transcriptName,b->transcriptName);
}



static Array readIsoquants (char *fileName) 
{
  Isoquant *currIsoquant;
  LineStream ls;
  char* line;
  Array isoquants;
  Texta tokens;

  isoquants = arrayCreate (30000,Isoquant);
  ls = ls_createFromFile (fileName);
  while (line = ls_nextLine (ls)) {
    if (line[0] == '\0') {
      continue;
    }
    tokens = textFieldtok (line,"\t");
    currIsoquant = arrayp (isoquants,arrayMax (isoquants),Isoquant);
    currIsoquant->transcriptName = hlr_strdup (textItem (tokens,3));
    currIsoquant->rpkm = atof (textItem (tokens,5));
    textDestroy (tokens);
  }
  ls_destroy (ls);
  return isoquants;
}



static int sortKnownIsoformsById (KnownIsoform *a, KnownIsoform *b) 
{
  return a->id - b->id;
}



static int sortIntervalsByName (Interval *a, Interval *b) 
{
  return strcmp (a->name,b->name);
}



static Interval* getTranscriptFromName (Array transcripts, char *name) 
{
  Interval *currTranscript,testTranscript;
  int index;

  testTranscript.name = hlr_strdup (name);
  if (arrayFind (transcripts,&testTranscript,&index,(ARRAYORDERF)sortIntervalsByName)) {
    currTranscript = arrp (transcripts,index,Interval);
  }
  else {
    die ("Expected to find transcript: %s",name);
  }
  hlr_free (testTranscript.name);
  return currTranscript;
}



static double getRpkm (Array isoquants, char *transcriptName) 
{
  Isoquant testIsoquant,*currIsoquant;
  int index;

  testIsoquant.transcriptName = hlr_strdup (transcriptName);
  if (arrayFind (isoquants,&testIsoquant,&index,(ARRAYORDERF)sortIsoquants)) {
    currIsoquant = arrp (isoquants,index,Isoquant);
  }
  else {
    die ("Expected to find RPKM for %s!",transcriptName);
  }
  hlr_free (testIsoquant.transcriptName);
  return currIsoquant->rpkm;
}



static void printOutput (Array transcripts, Array isoquants, Texta transcriptNames) 
{
  int i,j;
  Interval *currTranscript;
  SubInterval *currExon;
  static int geneId = 1;
  double rpkm;
  int count;

  count = 0;
  for (i = 0; i < arrayMax (transcriptNames); i++) {
    rpkm = getRpkm (isoquants,textItem (transcriptNames,i));
    if (rpkm < MIN_RPKM) {
      continue;
    }
    count++;
    currTranscript = getTranscriptFromName (transcripts,textItem (transcriptNames,i));
    printf ("%s\tprediction\ttranscript\t%d\t%d\t.\t%c\t.\tgene_id \"gene_%d\"; transcript_id \"transcript_%d_%d\"; RPKM \"%.2f\";\n",
            currTranscript->chromosome,currTranscript->start,currTranscript->end,currTranscript->strand,geneId,geneId,count,rpkm);
    for (j = 0; j < arrayMax (currTranscript->subIntervals); j++) {
      currExon = arrp (currTranscript->subIntervals,j,SubInterval);
       printf ("%s\tprediction\texon\t%d\t%d\t.\t%c\t.\tgene_id \"gene_%d\"; transcript_id \"transcript_%d_%d\";\n",
               currTranscript->chromosome,currExon->start,currExon->end,currTranscript->strand,geneId,geneId,count);
    }
  }
  if (count > 0) {
    geneId++;
  }
}



int main (int argc, char *argv[]) 
{
  Array knownIsoforms;
  Array transcripts;
  KnownIsoform *currKnownIsoform,*nextKnownIsoform;
  Texta transcriptNames;
  int i,j;
  Array isoquants;

  if (argc != 4) {
    usage ("%s <prefix_isoformMap.txt> <prefix_transcriptCoordinates.interval> <prefix.isoquants>",argv[0]);
  }
  knownIsoforms = readKnownIsoforms (argv[1]); 
  arraySort (knownIsoforms,(ARRAYORDERF)sortKnownIsoformsById);
  transcripts = intervalFind_parseFile (argv[2]);
  arraySort (transcripts,(ARRAYORDERF)sortIntervalsByName);
  isoquants = readIsoquants (argv[3]);
  arraySort (isoquants,(ARRAYORDERF)sortIsoquants);
  transcriptNames = textCreate (100);
  i = 0;
  while (i < arrayMax (knownIsoforms)) {
    currKnownIsoform = arrp (knownIsoforms,i,KnownIsoform);
    textClear (transcriptNames);
    textAdd (transcriptNames,currKnownIsoform->transcriptName);
    j = i + 1;
    while (j < arrayMax (knownIsoforms)) {
      nextKnownIsoform = arrp (knownIsoforms,j,KnownIsoform);
      if (currKnownIsoform->id == nextKnownIsoform->id) {
        textAdd (transcriptNames,nextKnownIsoform->transcriptName);
      }
      else {
        break;
      }
      j++;
    }
    i = j;
    printOutput (transcripts,isoquants,transcriptNames);
  } 
  return 0;
}
