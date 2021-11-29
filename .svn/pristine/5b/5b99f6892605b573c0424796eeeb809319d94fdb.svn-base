CDECK  ID>, GTMUON.
*CMZ :  3.21/02 29/03/94  15.41.24  by  S.Giani
*-- Author :
      SUBROUTINE GTMUON
C.
C.    ******************************************************************
C.    *                                                                *
C.    *   Muon track. Computes step size and propagates particle       *
C.    *    through step.                                               *
C.    *                                                                *
C.    *   ==>Called by : GTRACK                                        *
C.    *      Authors     R.Brun, F.Bruyant, M.Maire ********           *
C.    *                                                                *
C.    ******************************************************************
C.
      implicit none

#include "gcbank.h"
#include "geant321/gccuts.inc"
#include "geant321/gcjloc.inc"
#include "geant321/gckine.inc"
#include "geant321/gcmate.inc"
#include "geant321/gcmate.inc"
#include "geant321/gcmulo.inc"
#include "geant321/gconsp.inc"
#include "geant321/gcphys.inc"
#include "geant321/gcstak.inc"
#include "geant321/gctmed.inc"
#include "geant321/gctrak.inc"
#include "geant321/gcking.inc"
#include "geant321/gcjump.inc"
#include "geant321/gcvolu.inc"
#include "geant321/gcvdma.inc"

#include "apcard.h"

#include "skdetsim_parm.h"
C
C

      integer iproc,ist,jst,i,ier,iek1,iek2,inear,isame,i1,ikcut
      real gekrt1,smulof,vectmp,corr,dems,gekint,rndm
      save ikcut

#include "epsmac.h"

      REAL GAM1,GAM2,DIRGAM(3)
      INTEGER ngam,ii
      REAL PGAM(3,2),EGAM(2)

      REAL MOMIN(4)

c    decay of mu+-
      real p_eminus(3),e_eminus
      real p_eplus(3),e_eplus

#ifdef WithG4
C    Inclusive neutron production of mu-
C    capture on 16O
      REAL*8 secVec(9,250), pmuminus(3)
      INTEGER Nscnd,iscnd 
#endif


C.
C.    ------------------------------------------------------------------
*
* *** Particle below energy threshold ? short circuit
*
      IF (GEKIN.LE.CUTMUO) GO TO 100

C     Save input momentum for storage in secondary stack later
      if (IAPMORETRUTH.ne.0) then
         do i=1,4
            MOMIN(i) = VECT(i+3)
         end do
      endif
*
* *** Update local pointers if medium has changed (only first time)
      IF (IUPD.EQ.0) THEN
         IUPD  = 1
         JLOSS = LQ(JMA-2)
         JBREM = LQ(JMA-9)
         JPAIR = LQ(JMA-10)
         JDRAY = LQ(JMA-11)
         JMUNU = LQ(JMA-14)
         JRANG = LQ(JMA-16)
         JCOEF = LQ(JMA-18)
         JMULOF= LQ(JTM-2)
         IF(IMCKOV.EQ.1) THEN
            JTCKOV = LQ(JTM-3)
            JABSCO = LQ(JTCKOV-1)
            JEFFIC = LQ(JTCKOV-2)
            JINDEX = LQ(JTCKOV-3)
            JCURIN = LQ(JTCKOV-4)
            NPCKOV = Q(JTCKOV+1)
         ENDIF
         OMCMOL= Q(JPROB+21)
         CHCMOL= Q(JPROB+25)
         IKCUT = Q(JMULOF+NEK1+1)
         STOPC = Q(JMULOF+NEK1+2)
         IF(ISTRA.GT.0) THEN
            JTSTRA = LQ(JMA-19)
            JTSTCO = LQ(JTSTRA-1)
            JTSTEN = LQ(JTSTRA-2)
         ENDIF
      ENDIF
*
* *** Compute current step size
*
      STEP   = STEMAX
      IPROC  = 103
      GEKRT1 = 1. -GEKRAT
      IEK1   = IEKBIN+NEK1
      IEK2   = IEKBIN+2*NEK1
*
*  **   Step limitation due to bremsstrahlung ?
*
      IF (IBREM.GT.0) THEN
         STEPBR = GEKRT1*Q(JBREM+IEK2) +GEKRAT*Q(JBREM+IEK2+1)
         SBREM  = STEPBR*ZINTBR
         IF (SBREM.LT.STEP) THEN
            STEP  = SBREM
            IPROC = 9
         ENDIF
      ENDIF
*
*  **   Step limitation due to pair production ?
*
      IF (IPAIR.GT.0) THEN
         STEPPA = GEKRT1*Q(JPAIR+IEK1) +GEKRAT*Q(JPAIR+IEK1+1)
         SPAIR  = STEPPA*ZINTPA
         IF (SPAIR.LT.STEP) THEN
            STEP  = SPAIR
            IPROC = 6
         ENDIF
      ENDIF
*
*  **   Step limitation due to decay ?
*
      IF (IDCAY.NE.0) THEN
         SDCAY = SUMLIF*VECT(7)/AMASS
         IF (SDCAY.LT.STEP) THEN
            STEP  = SDCAY
            IPROC = 5
         ENDIF
      ENDIF
*
*  **   Step limitation due to delta-ray ?
*
      IF (IDRAY.GT.0) THEN
         STEPDR = GEKRT1*Q(JDRAY+IEK2) +GEKRAT*Q(JDRAY+IEK2+1)
         SDRAY  = STEPDR*ZINTDR
         IF (SDRAY.LT.STEP) THEN
            STEP  = SDRAY
            IPROC = 10
         ENDIF
      ENDIF
*
*  **   Step limitation due to nuclear interaction ?
*
      IF (IMUNU.GT.0) THEN
         IF(GEKIN.GE.5.)THEN
            STEPMU = GEKRT1*Q(JMUNU+IEKBIN) +GEKRAT*Q(JMUNU+IEKBIN+1)
            SMUNU = STEPMU*ZINTMU
            IF (SMUNU.LT.STEP) THEN
               STEP  = SMUNU
               IPROC = 21
            ENDIF
         ELSE
            STEPMU = BIG
         ENDIF
      ENDIF
*
      IF (STEP.LE.0.) THEN
         STEP = 0.
         GO TO 90
      ENDIF
*
*  **   Step limitation due to energy-loss,multiple scattering
*             or magnetic field ?
*
      IF (JMULOF.NE.0) THEN
         SMULOF  = GEKRT1*Q(JMULOF+IEKBIN) +GEKRAT*Q(JMULOF+IEKBIN+1)
         IF (SMULOF.LT.STEP) THEN
            STEP  = SMULOF
            IPROC = 0
         ENDIF
      ENDIF
*
*  **   Step limitation due to geometry ?
*
      IF (STEP.GE.0.95*SAFETY) THEN
         CALL GTNEXT
         IF (IGNEXT.NE.0) THEN
            STEP  = SNEXT + PREC
            IPROC = 0
         ENDIF
*
*        Update SAFETY in stack companions, if any
         IF (IQ(JSTAK+3).NE.0) THEN
            DO 10 IST = IQ(JSTAK+3),IQ(JSTAK+1)
               JST    = JSTAK + 3 + (IST-1)*NWSTAK
               Q(JST+11) = SAFETY
   10       CONTINUE
            IQ(JSTAK+3) = 0
         ENDIF
      ELSE
         IQ(JSTAK+3) = 0
      ENDIF
*
* *** Linear transport when no field or very short step
*
      IF (IFIELD.EQ.0.OR.STEP.LE.PREC) THEN
*
         IF (IGNEXT.NE.0) THEN
            DO 20 I = 1,3
               VECTMP  = VECT(I) +STEP*VECT(I+3)
               IF(VECTMP.EQ.VECT(I)) THEN
*
* *** Correct for machine precision
*
                  IF(VECT(I+3).NE.0.) THEN
                     VECTMP =
     +               VECT(I)+ABS(VECT(I))*SIGN(1.,VECT(I+3))*EPSMAC
                     IF(NMEC.GT.0) THEN
                        IF(LMEC(NMEC).EQ.104) NMEC=NMEC-1
                     ENDIF
                     NMEC=NMEC+1
                     LMEC(NMEC)=104
                  ENDIF
               ENDIF
               VECT(I) = VECTMP
   20       CONTINUE
            INWVOL = 2
            NMEC = NMEC +1
            LMEC(NMEC) = 1
         ELSE
            DO 30 I = 1,3
               VECT(I)  = VECT(I) +STEP*VECT(I+3)
   30       CONTINUE
         ENDIF
      ELSE
*
* ***   otherwise, swim particle in magnetic field
*
         if(mycoun.gt.1.and.nfmany.gt.0.and.step.ge.safety)then
           nlevel=manyle(nfmany)
           do 99 i=1,nlevel
             names(i)=manyna(nfmany,i)
             number(i)=manynu(nfmany,i)
 99        continue
           call glvolu(nlevel,names,number,ier)
           if(ier.ne.0)print *,'Fatal error in GLVOLU'
           ingoto=0
         endif
         NMEC = NMEC +1
         LMEC(NMEC) = 4
*
   40    CALL JUMPT4(JUSWIM, CHARGE, STEP, VECT, VOUT)
*
*  ** When near to boundary, take proper action (cut-step,crossing...)
*
         IF(STEP.GE.SAFETY)THEN
            INEAR = 0
            IF (IGNEXT.NE.0) THEN
               DO 50 I = 1,3
                  VNEXT(I+3) = VECT(I+3)
                  VNEXT(I) = VECT(I) +SNEXT*VECT(I+3)
   50          CONTINUE
               DO 60 I = 1,3
                  IF (ABS(VOUT(I)-VNEXT(I)).GT.EPSIL) GO TO 70
   60          CONTINUE
               INEAR = 1
            ENDIF
*
   70       CALL GINVOL (VOUT, ISAME)
            IF (ISAME.EQ.0)THEN
               IF ((INEAR.NE.0).OR.(STEP.LT.EPSIL)) THEN
                  INWVOL = 2
                  NMEC = NMEC +1
                  LMEC(NMEC) = 1
               ELSE
*              Cut step
                  STEP = 0.5*STEP
                  IF (LMEC(NMEC).NE.24) THEN
                     NMEC = NMEC +1
                     LMEC(NMEC) = 24
                  ENDIF
                  GO TO 40
               ENDIF
            ENDIF
         ENDIF
*
         DO 80 I = 1,6
            VECT(I) = VOUT(I)
   80    CONTINUE
*
      ENDIF
*
* *** Correct the step due to multiple scattering
      IF (IMULL.NE.0) THEN
         STMULS = STEP
         CORR=0.0001*(STEP/RADL)*(GETOT/(VECT(7)*VECT(7)))**2
         IF (CORR.GT.0.25) CORR = 0.25
         STEP  = (1.+CORR)*STEP
      ENDIF
*
      SLENG = SLENG + STEP
*
* *** Generate Cherenkov photons if required
*
      IF(IMCKOV.EQ.1) THEN
         CALL GGCKOV
         IF(NGPHOT.NE.0) THEN
            NMEC=NMEC+1
            LMEC(NMEC)=105
         ENDIF
      ENDIF
*
* *** apply energy loss : find the kinetic energy corresponding
*      to the new stopping range = stopmx - step
*
      IF (ILOSL.NE.0) THEN
         NMEC = NMEC +1
         LMEC(NMEC) = 3
         IF(GEKRAT.LT.0.7) THEN
            I1 = MAX(IEKBIN-1,1)
         ELSE
            I1 = MIN(IEKBIN,NEKBIN-1)
         ENDIF
         I1 = 3*(I1-1)+1
         XCOEF1 = Q(JCOEF+I1)
         XCOEF2 = Q(JCOEF+I1+1)
         XCOEF3 = Q(JCOEF+I1+2)
         IF(XCOEF1.NE.0.) THEN
            STOPMX = -XCOEF2+SIGN(ONE,XCOEF1)*SQRT(XCOEF2**2 - (XCOEF3-
     +      GEKIN/XCOEF1))
         ELSE
            STOPMX = - (XCOEF3-GEKIN)/XCOEF2
         ENDIF
         STOPRG = STOPMX - STEP
         IF (STOPRG.LT.STOPC) THEN
            STEP = STOPMX - STOPC
            GO TO 100
         ENDIF
*
         IF(XCOEF1.NE.0.) THEN
            DEMEAN=GEKIN-XCOEF1*(XCOEF3+STOPRG*(2.*XCOEF2+STOPRG))
         ELSE
            DEMEAN=GEKIN-XCOEF2*STOPRG-XCOEF3
         ENDIF
         IF(DEMEAN.LE.5.*GEKIN*EPSMAC) THEN
            DEMEAN=(GEKRT1*Q(JLOSS+IEKBIN)+GEKRAT*Q(JLOSS+IEKBIN+1))
     +      *STEP
         ENDIF
         IF (ILOSS.EQ.4.OR.IEKBIN.LE.IKCUT+1) THEN
            DESTEP = DEMEAN
         ELSE
            DEMS = DEMEAN
            CALL GFLUCT(DEMS,DESTEP)
         ENDIF
         IF (DESTEP.LT.0.) DESTEP = 0.
         GEKINT = GEKIN -DESTEP
         IF (GEKINT.LE.(1.01*CUTMUO)) GO TO 100
         DESTEL = DESTEP
         GEKIN  = GEKINT
         GETOT  = GEKIN +AMASS
         VECT(7)= SQRT((GETOT+AMASS)*GEKIN)
         CALL GEKBIN
      ENDIF
*
* *** Apply multiple scattering.
*
      IF (IMULL.NE.0) THEN
         NMEC   = NMEC +1
         LMEC(NMEC) = 2
         CALL GMULTS
      ENDIF
*
* *** Update time of flight
*
      SUMLIF = SUMLIF -STEP*AMASS/VECT(7)
      TOFG = TOFG +STEP*GETOT/(VECT(7)*CLIGHT)
      IF (TOFG.GE.TOFMAX) THEN
         ISTOP = 4
         NMEC  = NMEC +1
         LMEC(NMEC) = 22
         GO TO 999
      ENDIF
*
* *** Update interaction probabilities
*
      IF (IBREM.GT.0) ZINTBR = ZINTBR -STEP/STEPBR
      IF (IPAIR.GT.0) ZINTPA = ZINTPA -STEP/STEPPA
      IF (IDRAY.GT.0) ZINTDR = ZINTDR -STEP/STEPDR
      IF (IMUNU.GT.0) ZINTMU = ZINTMU -STEP/STEPMU
*
* ***   otherwise, apply the selected process if any
*
   90 IF (IPROC.EQ.0) GO TO 999
      NMEC = NMEC +1
      LMEC(NMEC) = IPROC
*
*  **   Bremsstrahlung ?
*
      IF (IPROC.EQ.9) THEN
         CALL GBREMM

         if (IAPMORETRUTH.ne.0) goto 998
*
*  **   Pair production ?
*
      ELSE IF (IPROC.EQ.6) THEN
         CALL GPAIRM
*
*  **   Decay ?
*
      ELSE IF (IPROC.EQ.5) THEN
         ISTOP = 1
         CALL GDECAY
*
*  **   Delta-ray ?
*
      ELSE IF (IPROC.EQ.10) THEN
         CALL GDRAY
*
*  **   Nuclear interaction ?
*
      ELSE IF (IPROC.EQ.21) THEN
         CALL GMUNU
         if (IAPMORETRUTH.ne.0) goto 998
      ENDIF
      
      GO TO 999

C     Stop current track and put it on the secondaries stack
C     so we may track what happened
 998  if (IAPMORETRUTH.ne.0) then
         NGKINE = NGKINE+1
         do i=1,3
            GKIN(i,NGKINE) = VECT(i+3)*VECT(7)
            GPOS(i,NGKINE) = VECT(i)
         end do
         GKIN(4,NGKINE) = GETOT
         GKIN(5,NGKINE) = IPART
         TOFD(NGKINE) = 0.
         ISTOP = 1
      
C     Restore parent track momentum
         do i=1,4
            VECT(i+3) = MOMIN(i)
         end do

         GO TO 9999
      endif
*
* *** Special treatment for overstopped tracks
*
  100 DESTEP = GEKIN
      DESTEL = DESTEP
      GEKIN  = 0.
      GETOT  = AMASS
      VECT(7)= 0.
      INWVOL = 0
      NMEC   = NMEC +1
      LMEC(NMEC) = 30

#ifdef WithG4
C pmuminus : added by RA 2018 Sep.
C to handle neutrons from mu- capture on 16O
      pmuminus(1) = VECT(4)*VECT(7)
      pmuminus(2) = VECT(5)*VECT(7)
      pmuminus(3) = VECT(6)*VECT(7)
#endif

      IF (IDCAY.EQ.0) THEN
         ISTOP = 2
      ELSE
         NMEC   = NMEC +1
         TOFG   = TOFG +SUMLIF/CLIGHT
         SUMLIF = 0.
         IF (TOFG.GE.TOFMAX) THEN
            ISTOP = 4
            LMEC(NMEC) = 22
            GO TO 999
         ENDIF
         LMEC(NMEC) = 5
         ISTOP = 1
         call grndm(rndm,1)
c

         if (ipart .eq. 5) then
c           CALL GDECAY
           call skdcy_muplus(p_eplus,e_eplus)
            NGKINE = NGKINE+1
               GKIN(1,NGKINE) = p_eplus(1)
               GKIN(2,NGKINE) = p_eplus(2)
               GKIN(3,NGKINE) = p_eplus(3)
               GKIN(4,NGKINE) = e_eplus
               GKIN(5,NGKINE) = 2  ! positron
               TOFD(NGKINE) = 0.
               GPOS(1,NGKINE)=VECT(1)
               GPOS(2,NGKINE)=VECT(2)
               GPOS(3,NGKINE)=VECT(3)            

         elseif (rndm.lt.0.816 .and. IPART.eq.6) then

           call skdcy_muminus(p_eminus,e_eminus)
            NGKINE = NGKINE+1
               GKIN(1,NGKINE) = p_eminus(1)
               GKIN(2,NGKINE) = p_eminus(2)
               GKIN(3,NGKINE) = p_eminus(3)
               GKIN(4,NGKINE) = e_eminus
               GKIN(5,NGKINE) = 3  ! electron
               TOFD(NGKINE) = 0.
               GPOS(1,NGKINE)=VECT(1)
               GPOS(2,NGKINE)=VECT(2)
               GPOS(3,NGKINE)=VECT(3)            
         else

#ifdef WithG4
C   If you want to consider correlation between gammas and
C   neutrons, then the following implementation is not correct.
C   This is valid only for ** INCLUSIVE ** analyses.
C   e.g. Can not do **gamma - neutron coincidence**
C   ------ RA 2018 Sep.
          if (USEG4.eq.1) then
            CALL g4_mumicap(pmuminus(1), pmuminus(2),
     +                       pmuminus(3), Nscnd, secVec)
            IF ( Nscnd.NE.0 ) THEN
               DO iscnd=1, Nscnd
                   NGKINE=NGKINE+1
                   GKIN(1,NGKINE) =secVec(6,iscnd)
                   GKIN(2,NGKINE) =secVec(7,iscnd)
                   GKIN(3,NGKINE) =secVec(8,iscnd)
                   GKIN(4,NGKINE) =secVec(9,iscnd)
                   GKIN(5,NGKINE) =secVec(1,iscnd)
                   GPOS(1,NGKINE) =VECT(1)
                   GPOS(2,NGKINE) =VECT(2)
                   GPOS(3,NGKINE) =VECT(3)
                   TOFD(NGKINE)   = 0
               END DO
            ENDIF
          endif
#endif

           CALL SKGAMMA(ngam,pgam,egam)
c           write(18,*)ngam,pgam,egam
            if(ngam.eq.0)goto 358
            do ii=1,ngam
               NGKINE = NGKINE+1
               GKIN(1,NGKINE) = pgam(1,ii)
               GKIN(2,NGKINE) = pgam(2,ii)
               GKIN(3,NGKINE) = pgam(3,ii)
               GKIN(4,NGKINE) = egam(ii)
               GKIN(5,NGKINE) = 1
               TOFD(NGKINE) = 0.
               GPOS(1,NGKINE)=VECT(1)
               GPOS(2,NGKINE)=VECT(2)
               GPOS(3,NGKINE)=VECT(3)
            end do



  358       continue
         end if
c
      ENDIF
      
  999 IF(NGPHOT.GT.0) THEN
         IF(ITCKOV.EQ.2.AND.ISTOP.EQ.0) THEN
*
*  The muon has produced Cerenkov photons and it is still alive
*  we put it in the stack and we let the photons to be tracked
            NGKINE = NGKINE+1
            GKIN(1,NGKINE) = VECT(4)*VECT(7)
            GKIN(2,NGKINE) = VECT(5)*VECT(7)
            GKIN(3,NGKINE) = VECT(6)*VECT(7)
            GKIN(4,NGKINE) = GETOT
            GKIN(5,NGKINE) = IPART
            TOFD(NGKINE) = 0.
            ISTOP = 1
c----put position as well
            GPOS(1,NGKINE)=VECT(1)
            GPOS(2,NGKINE)=VECT(2)
            GPOS(3,NGKINE)=VECT(3)
         ENDIF
      ENDIF

 9999 END
