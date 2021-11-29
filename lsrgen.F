CDECK  ID>, GGCKOV.
      SUBROUTINE LSRGEN
C.
C.    ******************************************************************
C     *   This routine generates photons to simulate laser beam.       
C     *   It is based on ggckov and called by gtnino.
C     *
C.    ******************************************************************
C.    *                                                                *
C.    *   This routine is called for each tracking step of a charged   *
C.    *   particle in a radiator. A Poisson-distributed number of      *
C.    *   photons is generated according to the Cherenkov formula,     *
C.    *   distributed evenly along the track segment and uniformly     *
C.    *   azimuth w.r.t. the particle direction. The parameters are    *
C.    *   then transformed into the Master Reference System, and they  *
C.    *   are added to the particle stack.                             *
C.    *                                                                *
C.    *   ==>Called by : GTMUON, GTHADR, GTELEC                        *
C.    *      Authors     R.Jones, F.Carminati ********                 *
C.    *                                                                *
C.     *    2006.05.10 Corrected a bug for OD segmentation F.Dufour
C.    *                                                                *
C.    *                                                                *
C.    *                                                                *
C.    ******************************************************************
C.
      IMPLICIT NONE

CMM.. copied from ggckov 03/03/12 14:50
CMM.. copied partly from ggckov B. Yang  Mar. 4 2009
#include "gcbank.h"
#include "geant321/gcjloc.inc"
#include "geant321/gctmed.inc"
#include "geant321/gcunit.inc"
#include "geant321/gctrak.inc"
#include "geant321/gckine.inc"
#include "geant321/gcking.inc"
#include "geant321/gconsp.inc"
#include "skdetsim_parm.h"
#include "guflag.h"
#ifdef UPMU
#include "gumed_u.h"
#else
#include "gumed.h"
#endif
#include "aptrkseg.h"
#include "skparm.h"
#include "skdetsim_test.h"
#include "skhead.h"

CC for laser
#include "laser.h"
      integer nsave
      REAL phdir(3),spdir(3)
      INTEGER NPHOTON

CC to reduce the the number of photons
#include "guref.h"
      REAL PRABSL0
      REAL REDUC

C Local variables
*

      REAL RPHOT(2)
      LOGICAL ROTATE
*
      real costh,sinth,cosph,sinph,costmx,sinmx2,cost,sint2,sint,THPHO
      integer j
      real ds,phi
      real sinp,cosp

C Functions called
      real chersg,refsg,pmtqesg,cheracsg,acrnsg
*   ------------------------------------------------------------------
*

CMM      write(*,*) 'LSRGEN is called $$$ Nphoto = ',nphoto
CMM       write(*,*) 'NUMED = ',numed

      NGPHOT = 0
CC for laser
      NPHOTON = 0
      nsave=0
*F.Dufour bug correction: add IMDWOE
      IF(IFANTI.EQ.0.AND.(NUMED.EQ.IMDWTO.OR.NUMED.EQ.IMDWOE)) GOTO 999

*
* *** WATER and ACRYLIC MEDIUM Cherenkov Photon Generation
**F.Dufour bug correction: add IMDWOE

      IF(NUMED.EQ.IMDWTI.OR.NUMED.EQ.IMDWTO.
     & OR.NUMED.EQ.IMDWOE.OR.NUMED.EQ.IMDWTS.OR.NUMED.EQ.IMDACR) THEN


*
* ***  Set up rotation to Particle frame
*
         CALL GFANG(VECT(4),COSTH,SINTH,COSPH,SINPH,ROTATE)
C
C+ Exit without tossing photons if requested
C
         if (TRACK_PHOTON.eq.NOT_TRACK_PHOTON) return
C
         DO 40 J=1,NPHOTO
*
* *** Sample the wave length of the photon
 20         CALL GRNDM(RPHOT, 2)
*
* *** Q.E.
ccc         IF(RNDM().GT.PMTQESG(BLAMDA)) GOTO 40
C         print *, 'rnd,qe,lamda=',rphot(2),PMTQESG(BLAMDA),blamda
C-- Algorithm of reducing the number of photons is improved and can be applied also LE mode. Byeongsu Yang 2008/04/23
            IF ( (TRACK_PHOTON.eq.TRACK_QE_OVER_ABS_PHOTON) .or.
     & (TRACK_PHOTON.eq.TRACK_QE_OVER_ABS_SAVE_MEMORY) )   then
                  REDUC = PMTQESG(BLAMDA)*REDUCFCT/PRABSL0(BLAMDA)
                  IF( (numed.eq.2 .or. numed.eq.8 .or. numed.eq.7) .and. RPHOT(2) .GE. REDUC ) GOTO 40
            ELSEIF(TRACK_PHOTON.eq.TRACK_QE_PHOTON) then
C            IF(TRACK_PHOTON.eq.TRACK_QE_PHOTON) then
C-------------------------------------------------------------------
C increasing # of photons by 20% for the SK-II MC tuning
C                          Shoei NAKAYAMA (15-APR-2003)
C-from here---------------------------------------------------------
                  REDUC = PMTQESG(BLAMDA)*REDUCFCT/0.75
                  IF( (numed.eq.2 .or. numed.eq.8 .or. numed.eq.7 ) .and. RPHOT(2) .GT. REDUC ) GOTO 40

C               if(SK_GEOMETRY.eq.SK_I
C     &              .and.REFLECTION_TUNE.eq.AP_TUNE) then
C                  IF(numed.eq.2
C     &                 .and.RPHOT(2).GT.PMTQESG(BLAMDA)/(1.0-0.25))
C     &                 GOTO 40
C                  IF(numed.eq.8
C     &                 .and.RPHOT(2).GT.PMTQESG(BLAMDA)/(1.0-0.25))
C     &                 GOTO 40
C                  IF(numed.eq.7
C     &                 .and.RPHOT(2).GT.PMTQESG(BLAMDA)/(1.0-0.25))
C     &                 GOTO 40
C
C               else if(SK_GEOMETRY.eq.SK_II
C     &              .and.REFLECTION_TUNE.eq.AP_TUNE) then
C                  IF(numed.eq.2
C     &                 .and.RPHOT(2).GT.PMTQESG(BLAMDA)*1.2/(1.0-0.25))
C     &                 GOTO 40
C                  IF(numed.eq.8
C     &                 .and.RPHOT(2).GT.PMTQESG(BLAMDA)*1.2/(1.0-0.25))
C     &                 GOTO 40
C                  IF(numed.eq.7
C     &                 .and.RPHOT(2).GT.PMTQESG(BLAMDA)*1.2/(1.0-0.25))
C     &                 GOTO 40

C-- This case added so SK-III would use SK-II data.  Kirk Bays 2006/04/18
C     An amount of photon increase is changed to 1.25 (by G.Mitsuka, Oct 2007)
C               else if(SK_GEOMETRY.eq.SK_III
C     &              .and.REFLECTION_TUNE.eq.AP_TUNE) then
C                  IF(numed.eq.2
C     &                 .and.RPHOT(2).GT.PMTQESG(BLAMDA)*1.25/(1.0-0.25))
C     &                 GOTO 40
C                  IF(numed.eq.8
C     &                 .and.RPHOT(2).GT.PMTQESG(BLAMDA)*1.25/(1.0-0.25))
C     &                 GOTO 40
C                  IF(numed.eq.7
C     &                 .and.RPHOT(2).GT.PMTQESG(BLAMDA)*1.25/(1.0-0.25))
C     &                 GOTO 40
C
C               endif
C-to here-----------------------------------------------------------
            ENDIF

CMM... SAVE TO BANK EVERY MXPHOT EVENTS
            IF(NGPHOT.EQ.MXPHOT) THEN
                CALL GSKPHO(0)
                NGPHOT = 0
                nsave=nsave+1
CMM                write(*,*) 'LSRGEN:: Saved !;',nsave
            ENDIF

CMM...Test
C            IF(NGPHOT.GT.MXPHOT) THEN
C                CALL GSKPHO(0)
C                NGPHOT = 1
C            ENDIF 


            NGPHOT=NGPHOT+1
            NPHOTON=NPHOTON+1
CC            IF(NGPHOT.GT.MXPHOT) THEN
CC               WRITE(6,*) 'Warning in GGCKOV; Overflow in the photon stack'
CC               GOTO 999
CC            ENDIF

* *** Find position, direction, etc
            CALL GRNDM(RPHOT, 2)
            IF(IGNEXT.NE.0) THEN
               DS=(STEP-PREC)*RPHOT(1)+PREC
            ELSE
               DS = STEP*RPHOT(1)
            ENDIF
            XPHOT(1,NGPHOT) = VECT(1)-VECT(4)*DS
            XPHOT(2,NGPHOT) = VECT(2)-VECT(5)*DS
            XPHOT(3,NGPHOT) = VECT(3)-VECT(6)*DS
            XPHOT(11,NGPHOT)= TOFG+(STEP-DS)*GETOT/(VECT(7)*CLIGHT)

CMM         THPHO = RNGAUSK(0,SIGBM)
CYOO         CALL sankaku(0.,SIGBM,THPHO)  ! SIGBM set in ledetsim.F
CYOO         CALL sankaku(THPHO)           ! YOO think it's irrational! but
            CALL rngausyoo(0.,SIGBM,THPHO) ! YOO 07feb2000
CYOO          CALL RNGAUSK2(0.,SIGBM,THPHO)
CYOO         CALL RNGAUS(0.,SIGBM,THPHO)        ! YOO changed 06FEB2000

            SINT = SQRT(SINT2)
            PHI  = TWOPI*RPHOT(2)
            SINT = SIN(THPHO)
            COST = COS(THPHO)
            PHI  = TWOPI*RPHOT(2)
            SINP = SIN(PHI)
            COSP = COS(PHI)

CMM..Get polarization
            phdir(1) = SINT*COSP
            phdir(2) = SINT*SINP
            phdir(3) = COST
   
C Get random direction
C            call rndir(phdir)
   
            call rncone(phdir,90.,spdir)
   
   
            XPHOT(4,NGPHOT) = SINT*COSP
            XPHOT(5,NGPHOT) = SINT*SINP
            XPHOT(6,NGPHOT) = COST
C Get random direction
C            XPHOT(4,NGPHOT) = phdir(1)
C            XPHOT(5,NGPHOT) = phdir(2)
C            XPHOT(6,NGPHOT) = phdir(3)
            XPHOT(7,NGPHOT) = 1237.8E-09 / BLAMDA
C            XPHOT(8,NGPHOT) = COST*COSP
C            XPHOT(9,NGPHOT) = COST*SINP
C            XPHOT(10,NGPHOT) = -SINT
            XPHOT(8,NGPHOT) = spdir(1)
            XPHOT(9,NGPHOT) = spdir(2)
            XPHOT(10,NGPHOT) =spdir(3)
*      
            IF(ROTATE) THEN
               CALL GDROT(XPHOT(8,NGPHOT),COSTH,SINTH,COSPH,SINPH)
               CALL GDROT(XPHOT(4,NGPHOT),COSTH,SINTH,COSPH,SINPH)
            ENDIF
   40    CONTINUE
      ENDIF

*      write(*,*) 'NGPHOT = ',NGPHOT
*      write(*,*) 'NPHOTO = ',NPHOTO
*      write(*,*) 'nsave = ',nsave
*      write(*,*) 'MAXPHOTO = ',MXPHOT
      write(*,*) 'No. of GENERATED PHOTONS = ',NPHOTON
      
       IF(NGPHOT.GT.0) CALL GSKPHO(0)
C
C+ Standard Diagnostic Test Histograms
C
      if (OUTPUT_HBOOK .eq. OUT_TEST) then
         TEST_gen_phot=TEST_Gen_Phot+NPHOTON
      endif
       NGPHOT=0
 999   RETURN
      END
