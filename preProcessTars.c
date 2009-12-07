#include "format.h"
#include "log.h"
#include "linestream.h"
#include "intervalFind.h"
#include "mrfUtil.h"



#define TYPE_JUNCTION_START 1
#define TYPE_JUNCTION_END 2



typedef struct {
  int position;
  int type;
} Boundary;



static void getJunctionCoordinates (Interval *currJunction, int *leftStart, int *leftEnd, int *rightStart, int *rightEnd)
{
  SubInterval *leftPart;
  SubInterval *rightPart;
   
  leftPart = arrp (currJunction->subIntervals,0,SubInterval);
  rightPart = arrp (currJunction->subIntervals,1,SubInterval);
  *leftStart = leftPart->start;
  *leftEnd = leftPart->end;
  *rightStart = rightPart->start;
  *rightEnd = rightPart->end;
}



static int sortBoundariesByPosition (Boundary *a, Boundary *b) 
{
  return a->position - b->position;
}



int main (int argc, char *argv[]) 
{
  Array tars;
  Tar *currTar;
  int i,j,k;
  Array intervals;
  Interval *currJunction;
  int currLeftStart,currLeftEnd,currRightStart,currRightEnd;
  Array boundaries;
  Boundary *currBoundary,*prevBoundary;
  Array tarCoordinates;

  if (argc != 2) {
    usage ("%s <junctions.interval>",argv[0]);
  }
  intervalFind_addIntervalsToSearchSpace (argv[1],0);
  tars = readTarsFromBedFile ("-");
  boundaries = arrayCreate (1000,Boundary);
  tarCoordinates = arrayCreate (1000,int);
  for (i = 0; i < arrayMax (tars); i++) {
    currTar = arrp (tars,i,Tar);
    arrayClear (boundaries);
    intervals = intervalFind_getOverlappingIntervals (currTar->targetName,currTar->start,currTar->end);
    for (j = 0; j < arrayMax (intervals); j++) {
      currJunction = arru (intervals,j,Interval*);
      getJunctionCoordinates (currJunction,&currLeftStart,&currLeftEnd,&currRightStart,&currRightEnd);
      if (currTar->start <= currLeftEnd && currLeftEnd <= currTar->end) {
        currBoundary = arrayp (boundaries,arrayMax (boundaries),Boundary);
        currBoundary->position = currLeftEnd;
        currBoundary->type = TYPE_JUNCTION_START;
      }
      if (currTar->start <= currRightStart && currRightStart <= currTar->end) {
        currBoundary = arrayp (boundaries,arrayMax (boundaries),Boundary);
        currBoundary->position = currRightStart;
        currBoundary->type = TYPE_JUNCTION_END;
      }
    }
    arraySort (boundaries,(ARRAYORDERF)sortBoundariesByPosition);
    arrayClear (tarCoordinates);
    array (tarCoordinates,arrayMax (tarCoordinates),int) = currTar->start;
    for (k = 1; k < arrayMax (boundaries); k++) {
      prevBoundary = arrp (boundaries,k - 1,Boundary);
      currBoundary = arrp (boundaries,k,Boundary);
      if (prevBoundary->type == TYPE_JUNCTION_START && currBoundary->type == TYPE_JUNCTION_END) {
        array (tarCoordinates,arrayMax (tarCoordinates),int) = prevBoundary->position;
        array (tarCoordinates,arrayMax (tarCoordinates),int) = currBoundary->position;
      }
    }
    array (tarCoordinates,arrayMax (tarCoordinates),int) = currTar->end;
    if ((arrayMax (tarCoordinates)) % 2 != 0) {
      die ("Expected an even number of TAR coordinates");
    }
    for (k = 0; k < arrayMax (tarCoordinates) - 1; k = k + 2) {
      printf ("%s\t%d\t%d\n",currTar->targetName,arru (tarCoordinates,k,int),arru (tarCoordinates,k + 1,int));
    }
  }
  return 0;
}
