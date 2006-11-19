
if DEBUG
D_CXXFLAGS = -Wall -g -D__GNDEBUG__ 
else
endif
AM_CXXFLAGS = -DCOMMAND_LINE $(D_CXXFLAGS) 

bin_PROGRAMS = sgEvolver scoreAlignment2 mfaToMultiFiles maf2xmfa
sgEvolver_SOURCES = \
sgEvolver.cpp Mutator.cpp PhyloTree.cpp Alignment.cpp \
gnGapSource.cpp gnStringSource.cpp Alignment.h \
gnStringSource.h  PhyloTree.h GisSplayTree.h \
GappedIntervalSequence.h  gnGapSource.h  Mutator.h \
twister.c twister.h

LIBRARY_CL = $(DEPS_LIBS) 

sgEvolver_LDFLAGS = @STATIC_FLAG@
sgEvolver_LDADD = $(LIBRARY_CL)
EXTRA_sgEvolver_SOURCES = getopt.c getopt.h getopt1.c
sgEvolver_DEPENDENCIES = @GETOPT_LONG_SYSTEM@

scoreAlignment2_SOURCES = scoreAlignment2.cpp
scoreAlignment2_LDFLAGS = @STATIC_FLAG@
scoreAlignment2_LDADD = $(LIBRARY_CL)

mfaToMultiFiles_SOURCES = mfaToMultiFiles.cpp
mfaToMultiFiles_LDFLAGS = @STATIC_FLAG@
mfaToMultiFiles_LDADD = $(LIBRARY_CL)

maf2xmfa_SOURCES = maf2xmfa.cpp
maf2xmfa_LDFLAGS = @STATIC_FLAG@
maf2xmfa_LDADD = $(LIBRARY_CL)

INCLUDES = @DEPS_CFLAGS@ 

bin_SCRIPTS = simujobrun.pl simujobgen.pl scoregen.pl rgradientlegend.R rgradientplot.R
data_DATA = simujobparams.pm

