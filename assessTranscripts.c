#include "log.h"
#include "format.h"
#include "linestream.h"
#include "intervalFind.h"



#define TYPE_TRANSCRIBED 0
#define TYPE_NON_TRANSCRIBED 1



static int getType (Interval *currTranscript, int position) 
{
  int i;
  SubInterval *currExon;

  if (position < currTranscript->start || currTranscript->end < position) {
    return TYPE_NON_TRANSCRIBED;
  }
  for (i = 0; i < arrayMax (currTranscript->subIntervals); i++) {
    currExon = arrp (currTranscript->subIntervals,i,SubInterval);
    if (currExon->start <= position && position <= currExon->end) {
      return TYPE_TRANSCRIBED;
    }
  }
  return TYPE_NON_TRANSCRIBED;
}



static int sortTranscriptByName (Interval *a, Interval *b) 
{
  return strcmp (a->name,b->name);
}



int main (int argc, char *argv[]) 
{
  Array annotationTranscripts,predictionTranscripts;
  int i,j,k,l;
  Interval *currPT,*currAT;
  int matrix[2][2];
  int start,end;
  double sensitivity,specificity,bestSensitivity,bestSpecificity;
   
  if (argc != 3) {
    usage ("%s <transcriptCoordinatesAnnotation.interval> <transcriptCoordinatesPrediction.interval>",argv[0]);
  } 
  intervalFind_addIntervalsToSearchSpace (argv[1],0);
  predictionTranscripts = intervalFind_parseFile (argv[2]);
  arraySort (predictionTranscripts,(ARRAYORDERF)sortTranscriptByName);
  for (i = 0; i < arrayMax (predictionTranscripts); i++) {
    currPT = arrp (predictionTranscripts,i,Interval);
    annotationTranscripts = intervalFind_getOverlappingIntervals (currPT->chromosome,currPT->start,currPT->end);
    bestSensitivity = 0;
    bestSpecificity = 0;
    for (j = 0; j < arrayMax (annotationTranscripts); j++) {
      currAT = arru (annotationTranscripts,j,Interval*);
      if (currPT->strand != currAT->strand) {
        continue;
      }
      for (k = 0; k < 2; k++) {
        for (l = 0; l < 2; l++) {
          matrix [k][l] = 0;
        }
      }
      start = MIN (currPT->start,currAT->start);
      end = MAX (currPT->end,currAT->end);
      for (k = start; k <= end; k++) {
        matrix[getType (currPT,k)][getType (currAT,k)]++;
      }
      sensitivity = 1.0 * matrix[0][0] / (matrix[0][0] + matrix[1][0]);
      specificity = 1.0 * matrix[1][1] / (matrix[1][1] + matrix[0][1]);
      if ((bestSensitivity * bestSpecificity) < (sensitivity * specificity)) {
        bestSensitivity = sensitivity;
        bestSpecificity = specificity;
      }
    }
    printf ("%s\t%.4f\t%.4f\t%.4f\n",currPT->name,bestSpecificity,1-bestSpecificity,bestSensitivity);
  }
  return 0;
}
