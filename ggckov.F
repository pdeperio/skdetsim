CDECK  ID>, GGCKOV.
      SUBROUTINE GGCKOV
C.
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
C.    *      Authors     R.Jones, F.Carminati ********     
C.    *      
C.     *    2006.05.10 Corrected a bug for OD segmentation F.Dufour
C.    * 
C.    *
C.    *                                                                *
C.    ******************************************************************
C.
      IMPLICIT NONE

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

CC to reduce the the number of photons
#include "guref.h"
      REAL PRABSL0
      REAL REDUC

C Local variables

*
      REAL RPHOT(2)
      LOGICAL ROTATE
*
      real beta,chern,pho
      real costh,sinth,cosph,sinph,costmx,sinmx2,cost,sint2,sint
      real blammn,blamda,freq
      integer nphoto,j
      real ds,phi
      real sinp,cosp

C Functions called
      real chersg,refsg,pmtqesg,cheracsg,acrnsg

*   ------------------------------------------------------------------
*
      NGPHOT = 0
*F.Dufour bug correction: add IMDWOE     
      IF(IFANTI.EQ.0.AND.(NUMED.EQ.IMDWTO.OR.NUMED.EQ.IMDWOE)) GOTO 999
      IF(CHARGE.EQ.0.) GOTO 999

*
* *** WATER MEDIUM Cherenkov Photon Generation
**F.Dufour bug correction: add IMDWOE 

      IF(NUMED.EQ.IMDWTI.OR.NUMED.EQ.IMDWTO.
     &            OR.NUMED.EQ.IMDWOE.OR.NUMED.EQ.IMDWTS) THEN

         CALL GRNDM(RPHOT,1)
         BETA = VECT(7)/GETOT
         CHERN  = CHERSG(BETA)
         NPHOTO = IFIX(CHERN*STEP)
         PHO    = CHERN*STEP - FLOAT(NPHOTO)
         IF(PHO.GE.RPHOT(1)) NPHOTO=NPHOTO+1
         IF(NPHOTO.LE.0) GOTO 999
*
* ***  Set up rotation to Particle frame
*
         CALL GFANG(VECT(4),COSTH,SINTH,COSPH,SINPH,ROTATE)
*
* ***  Distribute the photons in origin, direction, momemtum
         BLAMMN = 300.
         COSTMX = 1./(BETA*REFSG(BLAMMN))
C Typo on following line corrected by D.Casper, 28-7-98 (COSMX->COSTMX)
         SINMX2 = (1.-COSTMX)*(1.+COSTMX)
C
C+ Store track segment if requested
C
         if (OUTPUT_HBOOK.eq.OUT_SEG .and. sinmx2.gt.0) call seg(itra, vect(4), tofg, 
     &        nphoto, asin(sqrt(sinmx2))*180.0/3.1415927)
C
C+ Exit without tossing photons if requested
C
         if (TRACK_PHOTON.eq.NOT_TRACK_PHOTON) return
C

         DO 40 J=1,NPHOTO
*
* *** Sample the wave length of the photon
 20         CALL GRNDM(RPHOT, 2)
            FREQ   = (4.283+RPHOT(1)*5.710)*1.E14
            BLAMDA = 3.E17/FREQ
         
*
* *** Find the density function value corresponding to the
* *** wave length sampled
            COST   = 1./(BETA*REFSG(BLAMDA))
            IF(COST.GT.1.) GOTO 999
            SINT2  = (1.-COST)*(1.+COST)
*
* *** Perform hit-and-miss
ccc         IF(rndm()*SINMX2.GT.SINT2) GO TO 20
*
* *** Q.E. 
ccc         IF(RNDM().GT.PMTQESG(BLAMDA)) GOTO 40
C         print *, 'rnd,qe,lamda=',rphot(2),PMTQESG(BLAMDA),blamda
C-- Algorithm of reducing the number of photons is improved and can be applied also LE mode. Byeongsu Yang 2008/04/23
            IF ( (TRACK_PHOTON.eq.TRACK_QE_OVER_ABS_PHOTON) .or.
     & (TRACK_PHOTON.eq.TRACK_QE_OVER_ABS_SAVE_MEMORY) )   then
                  REDUC = PMTQESG(BLAMDA)*REDUCFCT/PRABSL0(BLAMDA)
                  IF( (numed.eq.2 .or. numed.eq.8) .and. RPHOT(2) .GE. REDUC ) GOTO 40
            ELSEIF(TRACK_PHOTON.eq.TRACK_QE_PHOTON) then
C            IF(TRACK_PHOTON.eq.TRACK_QE_PHOTON) then
C-------------------------------------------------------------------
C increasing # of photons by 20% for the SK-II MC tuning
C                          Shoei NAKAYAMA (15-APR-2003)
C-from here---------------------------------------------------------
                  REDUC = PMTQESG(BLAMDA)*REDUCFCT/0.75
                  IF( (numed.eq.2 .or. numed.eq.8) .and. RPHOT(2) .GT. REDUC ) GOTO 40

C               if(SK_GEOMETRY.eq.SK_I
C     &              .and.REFLECTION_TUNE.eq.AP_TUNE) then
C                  IF(numed.eq.2
C     &                 .and.RPHOT(2).GT.PMTQESG(BLAMDA)/(1.0-0.25)) 
C     &                 GOTO 40
C                  IF(numed.eq.8
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
C               
C               endif
C-to acrylic--------------------------------------------------------
            ENDIF

            NGPHOT=NGPHOT+1
            IF(NGPHOT.GT.MXPHOT) THEN
               WRITE(6,*) 'Warning in GGCKOV; Overflow in the photon stack'
               GOTO 999
            ENDIF
*
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
            
            SINT = SQRT(SINT2)
            PHI  = TWOPI*RPHOT(2)
            SINP = SIN(PHI)
            COSP = COS(PHI)
            XPHOT(4,NGPHOT) = SINT*COSP
            XPHOT(5,NGPHOT) = SINT*SINP
            XPHOT(6,NGPHOT) = COST
            XPHOT(7,NGPHOT) = 1237.8E-09 / BLAMDA
            XPHOT(8,NGPHOT) = COST*COSP
            XPHOT(9,NGPHOT) = COST*SINP
            XPHOT(10,NGPHOT) = -SINT
c            write(6,*)'ggckov',ngphot,xphot(7,ngphot),blamda,numed
*
            IF(ROTATE) THEN
               CALL GDROT(XPHOT(8,NGPHOT),COSTH,SINTH,COSPH,SINPH)
               CALL GDROT(XPHOT(4,NGPHOT),COSTH,SINTH,COSPH,SINPH)
            ENDIF
            if (OUTPUT_SCATTABLE.eq.1) then
               call recordphotoninipos(XPHOT(1,NGPHOT),XPHOT(2,NGPHOT),XPHOT(3,NGPHOT))
            endif
c            print *,'Y-direction cosine = ',XPHOT(5,NGPHOT)
c            if (XPHOT(5,NGPHOT).lt.0.5) then
c               XPHOT(4,NGPHOT) = 0.
c               XPHOT(5,NGPHOT) = 1.
c               XPHOT(6,NGPHOT) = 0.
c            endif
 40      CONTINUE
*
* ***  ACRYLIC MEDUIM Cherenkov Photon Generation
*
      ELSE IF(NUMED.EQ.IMDACR) THEN
         CALL GRNDM(RPHOT,1)
         BETA = VECT(7)/GETOT
         CHERN  = CHERACSG(BETA)
         NPHOTO = IFIX(CHERN*STEP)
         PHO    = CHERN*STEP - FLOAT(NPHOTO)
         IF(PHO.GE.RPHOT(1)) NPHOTO=NPHOTO+1
         IF(NPHOTO.LE.0) GOTO 999
*
* ***  Set up rotation to Particle frame
*
         CALL GFANG(VECT(4),COSTH,SINTH,COSPH,SINPH,ROTATE)
*
* ***  Distribute the photons in origin, direction, momemtum
         BLAMMN = 300.
         COSTMX = 1./(BETA*ACRNSG(BLAMMN))
         SINMX2 = (1.-COSTMX)*(1.+COSTMX)
C
C+ Store track segment if requested
C
         if (OUTPUT_HBOOK.eq.OUT_SEG .and. sinmx2.gt.0) call seg(itra, vect(4), tofg, 
     &        nphoto, asin(sqrt(sinmx2))*180.0/3.1415927)
C
C+ Exit without tossing photons if requested
C
         if (TRACK_PHOTON.eq.NOT_TRACK_PHOTON) return
C

         DO 140 J=1,NPHOTO
*
* *** Sample the wave length of the photon
 120        CALL GRNDM(RPHOT, 2)
            FREQ   = (4.283+RPHOT(1)*5.710)*1.E14
            BLAMDA = 3.E17/FREQ
         
*
* *** Find the density function value corresponding to the
* *** wave length sampled
            COST   = 1./(BETA*ACRNSG(BLAMDA))
            IF(COST.GT.1.) GOTO 999
            SINT2  = (1.-COST)*(1.+COST)
*
* *** Perform hit-and-miss
ccc         IF(rndm()*SINMX2.GT.SINT2) GO TO 120
*
* *** Q.E. 
ccc         IF(RNDM().GT.PMTQESG(BLAMDA)) GOTO 140
C         print *, 'rnd,qe,lamda=',rphot(2),PMTQESG(BLAMDA),blamda
C-- Algorithm of reducing the number of photons is improved and can be applied also LE mode. Byeongsu Yang 2008/04/23
            IF( (TRACK_PHOTON.eq.TRACK_QE_OVER_ABS_PHOTON) .or.
     & (TRACK_PHOTON.eq.TRACK_QE_OVER_ABS_SAVE_MEMORY) )   then
                  REDUC = PMTQESG(BLAMDA)*REDUCFCT/PRABSL0(BLAMDA)
                  IF( numed.eq.7 .and. RPHOT(2) .GE. REDUC ) GOTO 140
            ELSEIF(TRACK_PHOTON.eq.TRACK_QE_PHOTON) then
C            IF(TRACK_PHOTON.eq.TRACK_QE_PHOTON) then
C-------------------------------------------------------------------
C increasing # of photons by 20% for the SK-II MC tuning
C                          Shoei NAKAYAMA (15-APR-2003)
C-from water--------------------------------------------------------
                  REDUC = PMTQESG(BLAMDA)*REDUCFCT/0.75
                  IF( numed.eq.7 .and. RPHOT(2) .GT. REDUC ) GOTO 140

C               if(SK_GEOMETRY.eq.SK_I
C     &              .and.REFLECTION_TUNE.eq.AP_TUNE) then
C                  IF(numed.eq.7.and.RPHOT(2).GT.PMTQESG(BLAMDA)/(1.0-0.25)) 
C     &                 GOTO 140                     
C               else if(SK_GEOMETRY.eq.SK_II
C     &              .and.REFLECTION_TUNE.eq.AP_TUNE) then
C                  IF(numed.eq.7.and.RPHOT(2).GT.PMTQESG(BLAMDA)*1.2/(1.0-0.25)) 
C     &                 GOTO 140                     
C     An amount of photon increase is changed to 1.25 (by G.Mitsuka, Oct 2007)
C               else if(SK_GEOMETRY.eq.SK_III
C     &              .and.REFLECTION_TUNE.eq.AP_TUNE) then
C                  IF(numed.eq.7.and.RPHOT(2).GT.PMTQESG(BLAMDA)*1.25/(1.0-0.25)) 
C     &                 GOTO 140                     
C 
C               endif
C-to here-----------------------------------------------------------
            ENDIF

            NGPHOT=NGPHOT+1
            IF(NGPHOT.GT.MXPHOT) THEN
               WRITE(6,*) 'Warning in GGCKOV; Overflow in the photon stack'
               GOTO 999
            ENDIF
*
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
            
            SINT = SQRT(SINT2)
            PHI  = TWOPI*RPHOT(2)
            SINP = SIN(PHI)
            COSP = COS(PHI)
            XPHOT(4,NGPHOT) = SINT*COSP
            XPHOT(5,NGPHOT) = SINT*SINP
            XPHOT(6,NGPHOT) = COST
            XPHOT(7,NGPHOT) = 1237.8E-09 / BLAMDA
            XPHOT(8,NGPHOT) = COST*COSP
            XPHOT(9,NGPHOT) = COST*SINP
            XPHOT(10,NGPHOT) = -SINT
c            write(6,*)'ggckov',ngphot,blamda, numed
*
            IF(ROTATE) THEN
               CALL GDROT(XPHOT(8,NGPHOT),COSTH,SINTH,COSPH,SINPH)
               CALL GDROT(XPHOT(4,NGPHOT),COSTH,SINTH,COSPH,SINPH)
            ENDIF

            if (OUTPUT_SCATTABLE.eq.1) then
               call recordphotoninipos(XPHOT(1,NGPHOT),XPHOT(2,NGPHOT),XPHOT(3,NGPHOT))
            endif

 140     CONTINUE
                 
      ENDIF
C
C+ Standard Diagnostic Test Histograms
C
      if (OUTPUT_HBOOK .eq. OUT_TEST) then
         TEST_gen_phot=TEST_Gen_Phot+NGPHOT
      endif

 999  RETURN
      END



