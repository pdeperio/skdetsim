*
* $Id: gskpho.F,v 1.2 2012/02/02 15:11:19 shimpeit Exp $
*
* $Log: gskpho.F,v $
* Revision 1.2  2012/02/02 15:11:19  shimpeit
* Make sure array indices do not go out of bound.
*
* Revision 1.1  2012/01/27 02:08:57  shimpeit
* Scattered/reflected photons are tagged.
* The information is stored in ISTORY variable in the GCRACK common block:
* ISTORY=(# of refl)*1000+(# of scat)
*
* Revision 1.1.1.1  1995/10/24 10:21:43  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 03/07/94  19.07.28  by  S.Giani
*-- Author :
      SUBROUTINE GSKPHO (IGK)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *                                                                *
C.    *   Stores in stack JSTAK either the IGKth Cherenkov photon of   *
C.    *   /GCKIN2/, or the NPHOT tracks when IGK is 0.                 *
C.    *                                                                *
C.    *   Called by : 'User'                                           *
C.    *   Authors   : F.Carminati                                      *
C.    *                                                                *
C.    ******************************************************************
C.
#include "gcbank.h"
#include "geant321/gcjloc.inc"
#include "geant321/gcking.inc"
#include "geant321/gcnum.inc"
#include "geant321/gcunit.inc"
#include "geant321/gctrak.inc"
#include "geant321/gcmzfo.inc"
#include "geant321/gcstak.inc"
#include "geant321/gckine.inc"

#include "skdetsim_parm.h"

      INTEGER*4 IGK,FLGCPHOT
*
      PARAMETER (NWPOLA=3)
C.
C.    ------------------------------------------------------------------
C.
*
* *** Make sure the request is reasonable
      IF (NGPHOT.LE.0) GO TO 999
      FLGCPHOT=0
      IF (IGK.EQ.0) THEN
         N1 = 1
         N2 = NGPHOT
      ELSE
         IF ((OUTPUT_SCATTABLE.eq.1).and.(IGK.LT.0)) THEN ! photon is indirect
            IF (IGK.LE.-20000000) THEN ! reflected photon
               FLGCPHOT = -IGK - 20000000 + 1000
            ELSE IF ((IGK.GT.-20000000).and.(IGK.LE.-10000000)) THEN ! scattered photon
               FLGCPHOT = -IGK - 10000000 + 1
            ELSE ! revived photon
               FLGCPHOT = -IGK
            ENDIF
            N1 = 1
            N2 = NGPHOT
         ELSE
            IF (IGK.LT.0.OR.IGK.GT.NGPHOT) THEN
               WRITE(CHMAIL,10000) IGK, NGPHOT
10000 FORMAT(' **** GSKPHO: ',I6,'th photon requested, ',
     + I6,' in stack')
               CALL GMAIL(0,0)
               GO TO 999
            ENDIF
            N1 = IGK
            N2 = IGK
         ENDIF
      ENDIF
*
* *** Store photons in stack JSTAK
*
* *** Check that the bank is there and big enough
      IF (JSTAK.EQ.0) THEN
         NDBOOK = NTSTKP*NWSTAK +3
         NDPUSH = NTSTKS*NWSTAK
         CALL MZBOOK (IXCONS,JSTAK,JSTAK,1,'STAK', 1,1,NDBOOK, IOSTAK,3)
         IQ(JSTAK+2) = NTSTKP
      ENDIF
      IF(IQ(JSTAK-2).EQ.0) THEN
         CALL MZPUSH(IXCONS,JSTAK,1,0,'I')
      ENDIF
      IF(LQ(JSTAK-1).EQ.0) THEN
         CALL MZBOOK(IXCONS,JPOLAR,JSTAK,-1,'POLA',0,0,
     +               NTSTKP*NWPOLA,3,-1)
      ENDIF
      MISSNG=IQ(JSTAK+1)+N2-N1+1-IQ(JSTAK+2)
      IF (MISSNG.GT.0) THEN
         CALL MZPUSH (IXCONS, JSTAK, 0, MISSNG*NWSTAK+NDPUSH, 'I')
         IQ(JSTAK+2) = IQ(JSTAK+2) +NTSTKS+MISSNG
      ENDIF
      MSSPOL = IQ(JSTAK+1)+N2-N1+1-IQ(JPOLAR-1)/3
      IF(MSSPOL.GT.0) THEN
         CALL MZPUSH (IXCONS,JPOLAR,0,(NTSTKS+MSSPOL)*NWPOLA, 'I')
      ENDIF
*
* *** Now we can go
      DO 20  NN=N1,N2
         JST = JSTAK  +IQ(JSTAK+1)*NWSTAK +3
         JPO = JPOLAR +IQ(JSTAK+1)*NWPOLA
         IQ(JSTAK+1) = IQ(JSTAK+1) +1
         IF (IQ(JSTAK+3).EQ.0) IQ(JSTAK+3) = IQ(JSTAK+1)
         IF (IQ(JSTAK+1).GT.NSTMAX) NSTMAX = IQ(JSTAK+1)
*
         IQ(JST+1) = -ITRA
         IQ(JST+2) = 50

         if (OUTPUT_SCATTABLE.eq.1) then
            IQ(JST+3) = FLGCPHOT
         else
            IQ(JST+3) = 0
         endif

         DO 10 I = 1,3
            Q(JST+3+I) = XPHOT(I  ,NN)
            Q(JST+6+I) = XPHOT(I+3,NN)*XPHOT(7,NN)
            Q(JPO+  I) = XPHOT(I+7,NN)
   10    CONTINUE
         Q(JST+10) = XPHOT(11,NN)
         Q(JST+11) = SAFETY
         Q(JST+12) = UPWGHT
*
         NALIVE = NALIVE +1
   20 CONTINUE
*                                                             END GSKPHO
  999 END
