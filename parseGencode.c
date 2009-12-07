#include "log.h"
#include "format.h"
#include "linestream.h"
#include <stdio.h>



// <seqname> <source> <feature> <start> <end> <score> <strand> <frame> [attributes] [comments] 



typedef struct {
  char* chromosome;
  char* source;
  char* feature;
  int start;
  int end;
  char* score;
  char* strand;
  int frame;
  char* group;
  char* geneId;
  char* transcriptId;
  char* transcriptType;
  char* geneType;
  char* geneName;
  char* transcriptName;
} GtfEntry;



static int sortGtfEntries (GtfEntry *a, GtfEntry *b)
{
  int diff;

  diff = strcmp (a->geneId,b->geneId);
  if (diff != 0) {
    return diff;
  }
  return strcmp (a->transcriptId,b->transcriptId);
} 



static void writeAnnotation (FILE *fp, GtfEntry *currGtfEntry, Array items)
{
  int i;
  static Array starts = NULL; 
  static Array ends = NULL;
  GtfEntry *thisGtfEntry;

  if (starts == NULL) {
    starts = arrayCreate (100,int);
  }
  else {
    arrayClear (starts);
  }
  if (ends == NULL) {
    ends = arrayCreate (100,int);
  }
  else {
    arrayClear (ends);
  }
  for (i = 0; i < arrayMax (items); i++) {
    thisGtfEntry = arru (items,i,GtfEntry*);
    if (strEqual (thisGtfEntry->feature,"CDS")) {
      array (starts,arrayMax (starts),int) = thisGtfEntry->start - 1; // convert to zero-based intervals
      array (ends,arrayMax (ends),int) = thisGtfEntry->end;
    }
  }
  if (arrayMax (starts) != arrayMax (ends)) {
    die ("Unequal number of starts and ends");
  }
  if (strEqual (currGtfEntry->chromosome,"chrMT")) {
    currGtfEntry->chromosome[4] = '\0';
  }
  arraySort (starts,(ARRAYORDERF)arrayIntcmp);
  arraySort (ends,(ARRAYORDERF)arrayIntcmp);
  if (strEqual (currGtfEntry->strand,"+")) {
    arru (ends,arrayMax (ends) - 1,int) = arru (ends,arrayMax (ends) - 1,int) + 3;
  }
  else if (strEqual (currGtfEntry->strand,"-")) {
    arru (starts,0,int) = arru (starts,0,int) - 3;
  }
  fprintf (fp,"%s|%s|%s|%s\t%s\t%s\t",
          currGtfEntry->geneId,currGtfEntry->transcriptId,currGtfEntry->geneName,currGtfEntry->transcriptName,
          currGtfEntry->chromosome,currGtfEntry->strand);
  fprintf (fp,"%d\t%d\t%d\t",arru (starts,0,int),arru (ends,arrayMax (ends) - 1,int),arrayMax (starts));
  for (i = 0; i < arrayMax (starts); i++) {
    fprintf (fp,"%d%s",arru (starts,i,int),i < arrayMax (starts) - 1 ? "," : "\t");
  }
  for (i = 0; i < arrayMax (ends); i++) {
    fprintf (fp,"%d%s",arru (ends,i,int),i < arrayMax (ends) - 1 ? "," : "\n");
  }
}



static char* getAttribute (Texta tokens, char *attributeName)
{
  int i;
  char *pos1,*pos2;

  i = 0; 
  while (i < arrayMax (tokens)) {
    if (strstr (textItem (tokens,i),attributeName)) {
      pos1 = strchr (textItem (tokens,i),'"');
      pos2 = strrchr (textItem (tokens,i),'"');
      if (pos1 == NULL || pos2 == NULL) {
        die ("Unexpected token: %s",textItem (tokens,i));
      } 
      *pos2 = '\0';
      return (pos1 + 1);
    }
    i++;
  }
  die ("Expected to find attribute: %s",attributeName);
  return NULL;
}



static int isValid (Array items) 
{
  int i;
  GtfEntry *currGtfEntry;
  int hasCds,hasStart;
  
  hasCds = 0;
  hasStart = 0;
  for (i = 0; i < arrayMax (items); i++) {
    currGtfEntry = arru (items,i,GtfEntry*);
    if (strEqual (currGtfEntry->feature,"CDS")) {
      hasCds = 1;
    }
    if (strEqual (currGtfEntry->feature,"start_codon")) {
      hasStart = 1;
    }
  }
  if (hasCds == 1 && hasStart == 1) {
    return 1;
  }
  return 0;
}



int main (int argc, char *argv[]) 
{
  LineStream ls;
  char* line;
  WordIter w;
  GtfEntry *currGtfEntry,*nextGtfEntry,*prevGtfEntry;
  Array gtfEntries;
  int i,j;
  Texta tokens;
  Array items;
  int isoformId;
  FILE *fp1,*fp2;
  Stringa buffer;
  
  if (argc != 2) {
    usage ("%s <prefix>",argv[0]);
  }
  buffer = stringCreate (100);
  stringPrintf (buffer,"%s_isoformMap.txt",argv[1]);
  fp1 = fopen (string (buffer),"w");
  stringPrintf (buffer,"%s_transcriptCoordinates.interval",argv[1]);
  fp2 = fopen (string (buffer),"w");
  if (fp1 == NULL || fp2 == NULL) {
    die ("Unable to open files");
  }
  gtfEntries = arrayCreate (10000,GtfEntry);
  ls = ls_createFromFile ("-");
  while (line = ls_nextLine (ls)) {
    w = wordIterCreate (line,"\t",0);
    currGtfEntry = arrayp (gtfEntries,arrayMax (gtfEntries),GtfEntry);
    currGtfEntry->chromosome = hlr_strdup (wordNext (w));
    currGtfEntry->source = hlr_strdup (wordNext (w));
    currGtfEntry->feature = hlr_strdup (wordNext (w));
    currGtfEntry->start = atoi (wordNext (w));
    currGtfEntry->end = atoi (wordNext (w));
    currGtfEntry->score = hlr_strdup (wordNext (w));
    currGtfEntry->strand = hlr_strdup (wordNext (w));
    currGtfEntry->frame = atoi (wordNext (w));
    currGtfEntry->group = hlr_strdup (wordNext (w));
    tokens = textStrtokP(currGtfEntry->group,";");
    currGtfEntry->geneId = hlr_strdup (getAttribute (tokens,"gene_id"));
    currGtfEntry->transcriptId = hlr_strdup (getAttribute (tokens,"transcript_id"));
    currGtfEntry->transcriptType = hlr_strdup (getAttribute (tokens,"transcript_type"));
    currGtfEntry->geneType = hlr_strdup (getAttribute (tokens,"gene_type"));
    currGtfEntry->geneName = hlr_strdup (getAttribute (tokens,"gene_name"));
    currGtfEntry->transcriptName = hlr_strdup (getAttribute (tokens,"transcript_name"));
    textDestroy (tokens);
    wordIterDestroy (w);
  }
  ls_destroy (ls);
  arraySort (gtfEntries,(ARRAYORDERF)sortGtfEntries);
  items = arrayCreate (100000,GtfEntry*);
  i = 0;
  isoformId = 1;
  while (i < arrayMax (gtfEntries)) {
    currGtfEntry = arrp (gtfEntries,i,GtfEntry);
    if (i == 0) {
      prevGtfEntry = currGtfEntry;
    }
    arrayClear (items);
    array (items,arrayMax (items),GtfEntry*) = currGtfEntry;
    j = i + 1;
    while (j < arrayMax (gtfEntries)) {
      nextGtfEntry = arrp (gtfEntries,j,GtfEntry);
      if (strEqual (currGtfEntry->geneId,nextGtfEntry->geneId) &&
          strEqual (currGtfEntry->transcriptId,nextGtfEntry->transcriptId)) {
        array (items,arrayMax (items),GtfEntry*) = nextGtfEntry;
      }
      else { 
        break;
      }
      j++;
    }
    i = j;
    if (isValid (items)) {
      writeAnnotation (fp2,currGtfEntry,items);
      if (!strEqual (prevGtfEntry->geneId,currGtfEntry->geneId)) {
        isoformId++;
        prevGtfEntry = currGtfEntry;
      }
      fprintf (fp1,"%d\t%s|%s|%s|%s\n",isoformId,currGtfEntry->geneId,currGtfEntry->transcriptId,currGtfEntry->geneName,currGtfEntry->transcriptName);
    }
  }
  fclose (fp1);
  fclose (fp2);
  return 0;
}
