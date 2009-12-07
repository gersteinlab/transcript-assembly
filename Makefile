# required input from environment:
# $(BIOINFOCONFDIR) and $(BICOSN)

# ---------------- setup symbols needed in the rest of the file --------

include $(BIOINFOCONFDIR)/biosdefs.make
.SUFFIXES:
SHELL = /bin/sh



ROOTFLAGS = -D_REENTRANT -pthread -m64
ROOTLIBS  = -L$(ROOTSYS)/lib -lCore -lCint -lRIO -lNet -lHist -lGraf \
            -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix \
            -lPhysics -lMathCore -lThread -lfreetype -pthread -lm -ldl
ROOTINC = -I$(ROOTSYS)/include



MRFDIR = /home1/lh372/SOFT/mrf



# ----------------------- entry points --------------

PROGRAMS = junctions2interval tophatSam2mrf preProcessTars linkTars processLinkedTars parseGencode transcriptStatistics assessTranscripts transcripts2gtf 


all: $(PROGRAMS)

clean:
	rm -rf $(PROGRAMS)



junctions2interval: junctions2interval.c $(MRFDIR)/mrf.o $(BIOSLIB)
	-@/bin/rm -f junctions2interval
	$(CC) $(CFLAGSO) $(BIOSINC) -I$(MRFDIR) junctions2interval.c $(MRFDIR)/mrf.o -o junctions2interval $(BIOSLNK)

tophatSam2mrf: tophatSam2mrf.c $(MRFDIR)/mrf.h $(BIOSLIB)
	-@/bin/rm -f tophatSam2mrf
	$(CC) $(CFLAGSO) $(BIOSINC) -I$(MRFDIR) tophatSam2mrf.c -o tophatSam2mrf $(BIOSLNK)

preProcessTars: preProcessTars.c $(MRFDIR)/mrfUtil.o $(BIOSLIB)
	-@/bin/rm -f preProcessTars
	$(CC) $(CFLAGSO) $(BIOSINC) -I$(MRFDIR) preProcessTars.c $(MRFDIR)/mrfUtil.o -o preProcessTars $(BIOSLNK) -lm

linkTars: linkTars.c transcriptAssemblyConfig.h $(BIOSLIB)
	-@/bin/rm -f linkTars
	$(CC) $(CFLAGSO) $(BIOSINC) linkTars.c -o linkTars $(BIOSLNK) -lm

processLinkedTars: processLinkedTars.c transcriptAssemblyConfig.h $(BIOSLIB)
	-@/bin/rm -f processLinkedTars
	$(CC) $(CFLAGSO) $(BIOSINC) processLinkedTars.c -o processLinkedTars $(BIOSLNK) -lm

parseGencode: parseGencode.c $(BIOSLIB)
	-@/bin/rm -f parseGencode
	$(CC) $(CFLAGSO) $(BIOSINC) parseGencode.c -o parseGencode $(BIOSLNK)

transcriptStatistics: transcriptStatistics.c $(BIOSLIB)
	-@/bin/rm -f transcriptStatistics
	$(CC) $(CFLAGSO) $(BIOSINC) transcriptStatistics.c -o transcriptStatistics $(BIOSLNK) -lm

assessTranscripts: assessTranscripts.c $(BIOSLIB)
	-@/bin/rm -f assessTranscripts
	$(CC) $(CFLAGSO) $(BIOSINC) assessTranscripts.c -o assessTranscripts $(BIOSLNK) -lm

transcripts2gtf: transcripts2gtf.c transcriptAssemblyConfig.h $(BIOSLIB)
	-@/bin/rm -f transcripts2gtf
	$(CC) $(CFLAGSO) $(BIOSINC) transcripts2gtf.c -o transcripts2gtf $(BIOSLNK) -lm

