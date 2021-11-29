CDECK  ID>, GDRELA.
*CMZ :  3.21/02 29/03/94  15.41.21  by  S.Giani
*-- Author :
      SUBROUTINE GDRELA
C.
C.    ******************************************************************
C.    *                                                                *
C.    *  Initialise ionisation energy loss by filling proton DE/DX     *
C.    *  tables for each material.                                     *
C.    *                                                                *
C.    *  For chemical mixtures,compounds & molecules the approximation *
C.    *  is made that                                                  *
C.    *                                                                *
C.    *          DE/DX = W(1)*DE/DX(1)+W(2)*DE/DX(2)...+W(N)*DE/DX(N)  *
C.    *  with,                                                         *
C.    *          DE/DX(i) appropriate to the i'th constituent.         *
C.    *                                                                *
C.    *  For mixtures    W(i) = fractional wght of each element.       *
C.    *  For molecules   W(i) = No. atoms*atomic wght/molecular wght.  *
C.    *                                                                *
C.    *    ==>Called by : GPHYSI                                       *
C.    *       Authors    R.Brun, G.Patrick  *********                  *
C.    *                                                                *
C.    ******************************************************************
C.
C  2001.Apr.24 modifyed by Y.Itow to consider correct density effect for H2O
C  2001.JUN.14 modified by H.Kim to consider correct density effect for H2O
C
      IMPLICIT NONE
#include "gcbank.h"
#include "geant321/gcmulo.inc"
#include "geant321/gcjloc.inc"
#include "geant321/gctrak.inc"
#include "geant321/gconsp.inc"
*
CYI Newly added paramteres to calc T1L1
      real T1L,T2L,EMPROT,T1LIM,T2LIM
      PARAMETER (T1L=0.00001,T2L=0.002)
      PARAMETER (EMPROT=0.9382723)
CYI
#include "geant321/gcmate.inc"
CYI
CYI   ======> new variable for apdetsim 
      real AATMP(2),ZZTMP(2)
      integer iwater
#include "gumate.h"
CYI
CYI
C Local defines for implicit none:
      integer nlmat,jel1,jel2,jel3,ichan,l
      real t,t1,dedx,dedx1,hmass,aa,zz,wght,dedxc,dedxc1
C
C.
C.    ------------------------------------------------------------------
C.
C            Number of constituents(ie. element,mixture or compound).
C
      NLMAT = Q(JMA+11)
      NLMAT  = IABS(NLMAT)
      IF (NLMAT.EQ.0) GO TO 999
      IF(Z.LT.1.) GO TO 999
C
      JEL1=LQ(JMA-1)
      JEL2=LQ(JMA-2)
      JEL3=LQ(JMA-3)
C
      ICHAN=IEKBIN
      T      = ELOW(ICHAN)
      T1     = 10.**(EKBIN(1)+(ICHAN-0.5)/GEKA)
CYI  ====> this medium is water?
CYI
         if(NLMAT.eq.2)then
           Do L=1,NLMAT
             AATMP(L) = Q(JMIXT+L)
             ZZTMP(L) = Q(JMIXT+NLMAT+L)
           enddo
           if((AATMP(1).eq.AWTR1.and.AATMP(2).eq.AWTR2).or.
     $(AATMP(2).eq.AWTR1.and.AATMP(1).eq.AWTR2))then
              iwater = 1
           else
              iwater = 0
           endif
CHK    
         ELSE
            IWATER =0

         ENDIF

CYI
CYI
C
C==========>  A,Ionisation losses for particles other than electrons
C            Simple element.
C
      IF (NLMAT.EQ.1) THEN
         CALL GDRELP(A,Z,DENS,T,DEDX)
         IF(DEDX.LT.0.)DEDX=0.
*
* *** auxiliary integration point for Range tables
         IF(ICHAN.NE.NEK1) THEN
            CALL GDRELP(A,Z,DENS,T1,DEDX1)
            IF(DEDX1.LT.0.)DEDX1=0.
         ENDIF
      ELSE
C
C            Mixture/compound : Loop over chemical constituents.
C

         DEDX = 0.
         DEDX1 = 0.
         HMASS = real(PMASS)
         T1LIM=HMASS*T1L/EMPROT
         T2LIM=HMASS*T2L/EMPROT
         if(iwater.eq.1.and.T.gt.T1LIM.and.T.gt.T2LIM)then
CYI This is water. And T is higher enough to use Bethe-Broch.
CYI Use special dedx routine with correct density effect
CYI Otherwise get into usual loop.
            CALL GDRELP_H2O(A,Z,DENS,T,DEDX)
            if(DEDX.LT.0)DEDX=0.
* *** auxiliary integration point for Range tables
            IF(ICHAN.NE.NEK1) THEN
               CALL GDRELP_H2O(A,Z,DENS,T1,DEDX1)
               IF(DEDX1.LT.0.)DEDX1=0.
            ENDIF
         else

         DO 10 L=1,NLMAT
            AA = Q(JMIXT+L)
            ZZ = Q(JMIXT+NLMAT+L)
            WGHT = Q(JMIXT+2*NLMAT+L)
            CALL GDRELP(AA,ZZ,DENS*WGHT,T,DEDXC)
            IF(DEDXC.LT.0.)DEDXC=0.
            DEDX = DEDX + WGHT*DEDXC
*
* *** auxiliary integration point for Range tables
            IF(ICHAN.NE.NEK1) THEN
               CALL GDRELP(AA,ZZ,DENS*WGHT,T1,DEDXC1)
               IF(DEDXC1.LT.0.)DEDXC1=0.
               DEDX1 = DEDX1 + WGHT*DEDXC1
            ENDIF
   10    CONTINUE
         endif
      ENDIF
C
      Q(JEL3+ICHAN)=Q(JEL3+ICHAN)+DEDX*DENS
      IF(ICHAN.NE.NEK1) THEN
         WS(NEKBIN*3+ICHAN) = WS(NEKBIN*3+ICHAN)+DEDX1*DENS
      ENDIF
C
C===========>  B, Ionisation losses for muons
C

C     Simple element
      IF(NLMAT.EQ.1) THEN
         CALL GDRELM(A,Z,DENS,T,DEDX)
         IF(DEDX.LT.0.)DEDX=0.
*
* *** auxiliary integration point for Range tables
         IF(ICHAN.NE.NEK1) THEN
            CALL GDRELM(A,Z,DENS,T1,DEDX1)
            IF(DEDX1.LT.0.)DEDX1=0.
         ENDIF
      ELSE
C
C     Mixture/compound
C
CHK  FOR MUON
         DEDX = 0.
         DEDX1 = 0.
         HMASS = real(EMMU)
         T1LIM=HMASS*T1L/EMPROT
         T2LIM=HMASS*T2L/EMPROT
CHK 
         if(iwater.eq.1.and.T.gt.T1LIM.and.T.gt.T2LIM)then 
CYI This is water. Use special dedx routine with correct density effect
            CALL GDRELM_H2O(A,Z,DENS,T,DEDX)
            if(DEDX.LT.0)DEDX=0.
* *** auxiliary integration point for Range tables
            IF(ICHAN.NE.NEK1) THEN
               CALL GDRELM_H2O(A,Z,DENS,T1,DEDX1)
               IF(DEDX1.LT.0.)DEDX1=0.
            ENDIF
         else
         DO 20 L=1,NLMAT
            AA=Q(JMIXT+L)
            ZZ=Q(JMIXT+NLMAT+L)
            WGHT=Q(JMIXT+2*NLMAT+L)
            CALL GDRELM(AA,ZZ,DENS*WGHT,T,DEDXC)
            IF(DEDXC.LT.0.)DEDXC=0.
            DEDX=DEDX+WGHT*DEDXC
*
* *** auxiliary integration point for Range tables
            IF(ICHAN.NE.NEK1) THEN
               CALL GDRELM(AA,ZZ,DENS*WGHT,T1,DEDXC1)
               IF(DEDXC1.LT.0.)DEDXC1=0.
               DEDX1 = DEDX1 + WGHT*DEDXC1
            ENDIF
   20    CONTINUE
         endif
      ENDIF
C
      Q(JEL2+ICHAN)=Q(JEL2+ICHAN)+DEDX*DENS
      IF(ICHAN.NE.NEK1) THEN
         WS(NEKBIN*2+ICHAN) = WS(NEKBIN*2+ICHAN)+DEDX1*DENS
      ENDIF
C
C===========>  C, Ionisation losses for electrons/positrons
C
CYI ==> Still same treatment of density effect for electron/positrons so far
CYI 
CHK ==> NEW DEDX FOR ELECTRON/POSITRON

      IF(IWATER.EQ.1) THEN
      CALL GDRELE_H2O(T,-1.,JMA,DEDX)
      Q(JEL1+ICHAN)=Q(JEL1+ICHAN)+DEDX
      CALL GDRELE_H2O(T,+1.,JMA,DEDX)
      Q(JEL1+ICHAN+NEK1)=Q(JEL1+ICHAN+NEK1)+DEDX
*
* *** auxiliary integration point for Range tables
      IF(ICHAN.NE.NEK1) THEN
         CALL GDRELE_H2O(T1,-1.,JMA,DEDX1)
         WS(ICHAN)=WS(ICHAN)+DEDX1
         CALL GDRELE_H2O(T1,+1.,JMA,DEDX1)
         WS(NEKBIN+ICHAN)=WS(NEKBIN+ICHAN)+DEDX1
      ENDIF
CHK
      ELSE
      CALL GDRELE(T,-1.,JMA,DEDX)
      Q(JEL1+ICHAN)=Q(JEL1+ICHAN)+DEDX
      CALL GDRELE(T,+1.,JMA,DEDX)
      Q(JEL1+ICHAN+NEK1)=Q(JEL1+ICHAN+NEK1)+DEDX
*
* *** auxiliary integration point for Range tables
      IF(ICHAN.NE.NEK1) THEN
         CALL GDRELE(T1,-1.,JMA,DEDX1)
         WS(ICHAN)=WS(ICHAN)+DEDX1
         CALL GDRELE(T1,+1.,JMA,DEDX1)
         WS(NEKBIN+ICHAN)=WS(NEKBIN+ICHAN)+DEDX1
      ENDIF
      ENDIF
C
  999 END


