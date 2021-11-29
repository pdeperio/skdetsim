C Changes from standard GEANT include:
C KWBANK increased to 3000000 words long.  Apparently gcalor needs the space
*
*
#ifndef CERNLIB_GEANT321_GCBANK_INC
#define CERNLIB_GEANT321_GCBANK_INC
*
*
* gcbank.inc
*
#if !defined(CERNLIB_INCLUDE)
#if defined(CERNLIB_TYPE)
#include "geant321/gtbank.inc"
#endif
      PARAMETER (KWBANK=3000000,KWWORK=5200)
      COMMON/GCBANK/NZEBRA,GVERSN,ZVERSN,IXSTOR,IXDIV,IXCONS,FENDQ(16)
     +             ,LMAIN,LR1,WS(KWBANK)
      DIMENSION IQ(2),Q(2),LQ(8000),IWS(2)
      EQUIVALENCE (Q(1),IQ(1),LQ(9)),(LQ(1),LMAIN),(IWS(1),WS(1))
      EQUIVALENCE (JCG,JGSTAT)
#include "geant321/gclink.inc"
#endif
*
* gcbank.inc
*
#if defined(CERNLIB_INCLUDE)
*
#include "geant321/incbeg.inc"
     +   gcbank
#include "geant321/incend.inc"
#endif

#endif
