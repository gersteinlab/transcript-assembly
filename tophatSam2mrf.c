#include "format.h"
#include "log.h"
#include "linestream.h"
#include "mrf.h"



int main (int argc, char *argv[]) 
{
  LineStream ls;
  char *line;
  Texta tokens;
  char *chromosome = NULL;
  char *cigar = NULL;
  char strand;
  int position;
  int intronic;
  int offset;
  int match1,match2;
  int i;

  puts (MRF_COLUMN_NAME_BLOCKS);
  ls = ls_createFromFile ("-");
  while (line = ls_nextLine (ls)) {
    strand = '.';
    tokens = textFieldtok (line,"\t");
    strReplace (&chromosome,textItem (tokens,2));
    position = atoi (textItem (tokens,3));
    strReplace (&cigar,textItem (tokens,5));
    if (arrayMax (tokens) == 13) {
      strand = textItem (tokens,12)[5];
    }
    textDestroy (tokens);
    tokens = textFieldtokP (cigar,"MN");
    if (arrayMax (tokens) == 2) {
      match1 = atoi (textItem (tokens,0));
      printf ("%s:%c:%d:%d:%d:%d\n",
              chromosome,strand,position,position + match1 - 1,1,match1);
    }
    else {
      offset = position;
      for (i = 0; i < arrayMax (tokens) - 2 - 1; i = i + 2) {
        match1 = atoi (textItem (tokens,i));
        intronic = atoi (textItem (tokens,i + 1));
        match2 = atoi (textItem (tokens,i + 2));
        printf ("%s:%c:%d:%d:%d:%d,%s:%c:%d:%d:%d:%d\n",
                chromosome,strand,offset,offset + match1 - 1,1,match1,
                chromosome,strand,offset + match1 + intronic,offset + match1 + intronic + match2 - 1,match1 + 1,match1 + match2);
        offset +=  (match1 + intronic);
      }
    }
    textDestroy (tokens);
  }
  ls_destroy (ls);
  return 0;
}
