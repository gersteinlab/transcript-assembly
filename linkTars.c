#include "format.h"
#include "log.h"
#include "linestream.h"
#include "intervalFind.h"
#include "numUtil.h"
#include "transcriptAssemblyConfig.h"



typedef struct {
  char *junctionId;
  int numSupportingReads;
} JunctionSupport;



typedef struct {
  char *tarId;
  double expressionValue;
} TarExpression;



typedef struct {
  Interval *tar;
  double tarExpression;
  Interval *junction;
  int numSupportingReads;
} Link;



typedef struct {
  Interval *junction;
  int flag;
} Flag;



static int sortLinksByTar (Link *a, Link *b)
{
  return a->tar - b->tar;
}



static int sortLinksByJunction (Link *a, Link *b)
{
  return a->junction - b->junction;
}



static int sortFlagsByJunction (Flag *a, Flag *b)
{
  return a->junction - b->junction;
}



static int sortLinks (Link **a, Link **b)
{
  return (*a)->junction - (*b)->junction;
}



static int sortJunctionSupportsByJunctionId (JunctionSupport *a, JunctionSupport *b) 
{
  return strcmp (a->junctionId,b->junctionId);
}



static int sortTarExpressionsByTarId (TarExpression *a, TarExpression *b) 
{
  return strcmp (a->tarId,b->tarId);
}



static Array readJunctionSupport (char *fileName) 
{
  LineStream ls;
  char *line,*pos;
  Array junctionSupports;
  JunctionSupport *currJS;

  junctionSupports = arrayCreate (100000,JunctionSupport);
  ls = ls_createFromFile (fileName);
  while (line = ls_nextLine (ls)) {
    pos = strchr (line,'\t');
    *pos = '\0';
    currJS = arrayp (junctionSupports,arrayMax (junctionSupports),JunctionSupport);
    currJS->junctionId = hlr_strdup (line);
    currJS->numSupportingReads = atoi (pos + 1);
  }
  ls_destroy (ls);
  return junctionSupports;
}



static int getNumSupportingReads (Array junctionSupports, char *junctionId)
{
  JunctionSupport testJS;
  int index;

  testJS.junctionId = hlr_strdup (junctionId);
  if (!arrayFind (junctionSupports,&testJS,&index,(ARRAYORDERF)sortJunctionSupportsByJunctionId)) { 
    die ("Unable to find junction support for %s",junctionId);
  }
  hlr_free (testJS.junctionId);
  return arrp (junctionSupports,index,JunctionSupport)->numSupportingReads;
}



static Array readTarExpression (char *fileName) 
{
  LineStream ls;
  char *line,*pos;
  Array tarExpressions;
  TarExpression *currTE;
  
  tarExpressions = arrayCreate (100000,TarExpression);
  ls = ls_createFromFile (fileName);
  while (line = ls_nextLine (ls)) {
    pos = strchr (line,'\t');
    *pos = '\0';
    currTE = arrayp (tarExpressions,arrayMax (tarExpressions),TarExpression);
    currTE->tarId = hlr_strdup (line);
    currTE->expressionValue = atof (pos + 1);
  }
  ls_destroy (ls);
  return tarExpressions;
}



static double getTarExpression (Array tarExpressions, char *tarId)
{
  TarExpression testTE;
  int index;

  testTE.tarId = hlr_strdup (tarId);
  if (!arrayFind (tarExpressions,&testTE,&index,(ARRAYORDERF)sortTarExpressionsByTarId)) { 
    die ("Unable to find tar expression for %s",tarId);
  }
  hlr_free (testTE.tarId);
  return arrp (tarExpressions,index,TarExpression)->expressionValue;
}



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



static int isLeftOverlap (int tarStart, int tarEnd, Interval *currJunction) 
{
  SubInterval *leftPart;
 
  leftPart = arrp (currJunction->subIntervals,0,SubInterval);
  if (rangeIntersection (leftPart->start,leftPart->end,tarStart,tarEnd) > 0) {
    return 1;
  }
  return 0;
}



static int isRightOverlap (int tarStart, int tarEnd, Interval *currJunction) 
{
  SubInterval *rightPart;
 
  rightPart = arrp (currJunction->subIntervals,1,SubInterval);
  if (rangeIntersection (rightPart->start,rightPart->end,tarStart,tarEnd) > 0) {
    return 1;
  }
  return 0;
}



static void setFlag (Array flags, Interval *currJunction)
{
  Flag testFlag,*currFlag;
  int index;
  testFlag.junction = currJunction;
  if (!arrayFind (flags,&testFlag,&index,(ARRAYORDERF)sortFlagsByJunction)) {
    die ("Unable to find flag for junction!");
  }
  currFlag = arrp (flags,index,Flag);
  currFlag->flag = 1;
}



static int getFlag (Array flags,  Interval *currJunction)
{
  Flag testFlag,*currFlag;
  int index;

  testFlag.junction = currJunction;
  if (!arrayFind (flags,&testFlag,&index,(ARRAYORDERF)sortFlagsByJunction)) {
    die ("Unable to find flag for junction!");
  }
  currFlag = arrp (flags,index,Flag);
  return currFlag->flag;
}



static void getLinksForTar (Array links, Interval *currTar, Array linksSortedByTar)
{
  Link testLink,*thisLink;
  int index;
  int i;

  testLink.tar = currTar;
  if (arrayFind (linksSortedByTar,&testLink,&index,(ARRAYORDERF)sortLinksByTar)) {
    i = index;
    while (i >= 0) {
      thisLink = arrp (linksSortedByTar,i,Link);
      if (thisLink->tar == currTar) { 
        array (links,arrayMax (links),Link*) = thisLink;
      }
      else {
        break;
      }
      i--;
    }
    i = index + 1;
    while (i < arrayMax (linksSortedByTar)) {
      thisLink = arrp (linksSortedByTar,i,Link);
      if (thisLink->tar == currTar) {
        array (links,arrayMax (links),Link*) = thisLink;
      }
      else {
        break;
      }
      i++;
    }
  }
  else {
    die ("Expected to find tar!");
  }
}



static void getTarsForJunction (Array tars, Interval *currJunction, Array linksSortedByJunction)
{
  Link testLink,*thisLink;
  int index;
  int i;

  testLink.junction = currJunction;
  if (arrayFind (linksSortedByJunction,&testLink,&index,(ARRAYORDERF)sortLinksByJunction)) {
    i = index;
    while (i >= 0) {
      thisLink = arrp (linksSortedByJunction,i,Link);
      if (thisLink->junction == currJunction) {
        array (tars,arrayMax (tars),Interval*) = thisLink->tar;
      }
      else {
        break;
      }
      i--;
    }
    i = index + 1;
    while (i < arrayMax (linksSortedByJunction)) {
      thisLink = arrp (linksSortedByJunction,i,Link);
      if (thisLink->junction == currJunction) {
        array (tars,arrayMax (tars),Interval*) = thisLink->tar;
      }
      else {
        break;
      }
      i++;
    }
  }
  else {
    die ("Expected to find junction!");
  }
}



static void linkTars (Link *currLink, Array linksSortedByTar, Array linksSortedByJunction, Array flags, int level, char *str, 
                      int numSupportingReads, double tarExpression) 
{
  int i;
  Array links;
  Stringa local;
  int currLeftStart,currLeftEnd,currRightStart,currRightEnd;
  int newLeftStart,newLeftEnd,newRightStart,newRightEnd;
  Array tars;
  Link *newLink;
  int moreTarsToLink;
  Interval *currTar;
  int rightTarEnd;

  if (level > MAX_NUM_RECURSIONS) {
    return;
  }
  local = stringCreate (100);
  stringCat (local,str);
  setFlag (flags,currLink->junction);
  getJunctionCoordinates (currLink->junction,&currLeftStart,&currLeftEnd,&currRightStart,&currRightEnd);
  tars = arrayCreate (100,Interval*);
  links = arrayCreate (100,Link*);
  getTarsForJunction (tars,currLink->junction,linksSortedByJunction);
  rightTarEnd = currLink->tar->end;
  for (i = 0; i < arrayMax (tars); i++) {
    currTar = arru (tars,i,Interval*);
    if (rightTarEnd < currTar->end) {
      rightTarEnd = currTar->end;
    }
    getLinksForTar (links,currTar,linksSortedByTar);
  }
  arraySort (links,(ARRAYORDERF)sortLinks);
  arrayUniq (links,NULL,(ARRAYORDERF)sortLinks);
  if (level == 0 ) {
    if (currLink->tar->start < currLeftEnd) {
      stringAppendf (local,"%d ",currLink->tar->start);
    }
    else {
      stringAppendf (local,"%d ",currLeftStart);
    }
  }  
  stringAppendf (local,"%d %d ",currLeftEnd,currRightStart);
  moreTarsToLink = 0;
  for (i = 0; i < arrayMax (links); i++) {
    newLink = arru (links,i,Link*);
    getJunctionCoordinates (newLink->junction,&newLeftStart,&newLeftEnd,&newRightStart,&newRightEnd);
    if (currRightStart < newLeftEnd && currLink->junction->strand == newLink->junction->strand) {
      linkTars (newLink,linksSortedByTar,linksSortedByJunction,flags,level + 1,string (local),
                numSupportingReads + newLink->numSupportingReads,tarExpression + newLink->tarExpression);
      moreTarsToLink = 1;
    }
  }
  if (moreTarsToLink == 0) {
    if (rightTarEnd > currRightStart) {
      printf ("%d %d %0.1f %s %c %s%d\n",level + 1,numSupportingReads,tarExpression,currLink->tar->chromosome,currLink->junction->strand,string (local),rightTarEnd);
    }
    else {
      printf ("%d %d %0.1f %s %c %s%d\n",level + 1,numSupportingReads,tarExpression,currLink->tar->chromosome,currLink->junction->strand,string (local),currRightEnd);
    }
  }
  arrayDestroy (links);
  arrayDestroy (tars);
  stringDestroy (local);
}



int main (int argc, char *argv[]) 
{
  Array junctionSupports;
  Array tarExpressions;
  Array tars;
  Interval *currTar;
  int i,j;
  Array intervals;
  Interval *currJunction;
  Array links;
  Link *currLink,*nextLink;
  Array linksSortedByTar,linksSortedByJunction;
  Flag *currFlag;
  Array flags;
  char *str;
  static char* prevChromosome = NULL;
  
  if (argc != 5) {
    usage ("%s <prefix_spliceJunctions.interval> <prefix_numSupportingReads.txt> <tars.interval> <tars.expression>",argv[0]);
  }
  links = arrayCreate (10000,Link);
  intervalFind_addIntervalsToSearchSpace (argv[1],0);
  junctionSupports = readJunctionSupport (argv[2]); 
  arraySort (junctionSupports,(ARRAYORDERF)sortJunctionSupportsByJunctionId);
  tars = intervalFind_parseFile (argv[3]);
  tarExpressions = readTarExpression (argv[4]);
  arraySort (tarExpressions,(ARRAYORDERF)sortTarExpressionsByTarId);
  for (i = 0; i < arrayMax (tars); i++) {
    currTar = arrp (tars,i,Interval);
    intervals = intervalFind_getOverlappingIntervals (currTar->chromosome,currTar->start,currTar->end);
    for (j = 0; j < arrayMax (intervals); j++) {
      currJunction = arru (intervals,j,Interval*);
      if (isLeftOverlap (currTar->start,currTar->end,currJunction) || 
          isRightOverlap (currTar->start,currTar->end,currJunction)) {
        currLink = arrayp (links,arrayMax (links),Link);
        currLink->junction = currJunction;
        currLink->numSupportingReads = getNumSupportingReads (junctionSupports,currJunction->name);
        currLink->tar = currTar;
        currLink->tarExpression = getTarExpression (tarExpressions,currTar->name);
      }
    }
  }
  linksSortedByJunction = arrayCopy (links);
  arraySort (linksSortedByJunction,(ARRAYORDERF)sortLinksByJunction);
  flags = arrayCreate (100000,Flag);
  i = 0;
  while (i < arrayMax (linksSortedByJunction)) {
    currLink = arrp (linksSortedByJunction,i,Link);
    currFlag = arrayp (flags,arrayMax (flags),Flag);
    currFlag->junction = currLink->junction;
    currFlag->flag = 0;
    j = i + 1;
    while (j < arrayMax (linksSortedByJunction)) {
      nextLink = arrp (linksSortedByJunction,j,Link);
      if (currLink->junction != nextLink->junction) {
        break;
      }
      j++;
    }
    i = j;
  }
  arraySort (flags,(ARRAYORDERF)sortFlagsByJunction);
  linksSortedByTar = arrayCopy (links);
  arraySort (linksSortedByTar,(ARRAYORDERF)sortLinksByTar);
  str = hlr_strdup ("");
  for (i = 0; i < arrayMax (links); i++) {
    currLink = arrp (links,i,Link);
    if (i == 0) {
      strReplace (&prevChromosome,currLink->tar->chromosome);
    }
    if (getFlag (flags,currLink->junction) == 1) {
      continue;
    }
    if (!strEqual (prevChromosome,currLink->tar->chromosome)) {
      puts ("###");
    }
    linkTars (currLink,linksSortedByTar,linksSortedByJunction,flags,0,str,currLink->numSupportingReads,currLink->tarExpression);
    strReplace (&prevChromosome,currLink->tar->chromosome);
  }
  hlr_free (str);
  return 0;
}
