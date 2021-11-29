************************************************************************
*     -------------------------------------
      SUBROUTINE DSANTI
*     -------------------------------------
*
*     (Purpose)
*        Calculate TQ of the ANTI PMT
*
*     (Input)
*
*     (Output)
*        'TQA' bank
*
*     (Creation Date and Author)
*        1994.11.02 ; First version by Koshio,Y
*        1998.02.16 : added calib. constants from tqa_13 by Kibayashi,A
*        2001. 2. 9 : RNDM -> GRNDM
*        2002. 7.15 : cleanup and tuning Kearns, Walter, Wang
*        2008. 1.17 : remove refences to geometry, rely solely on tube 
*                   : R. Wendell (rvw)
*
************************************************************************
      IMPLICIT NONE

#include "skparm.h"
#include "sktq.h"
#include "dsparma.h"
#include "dsparm.h"
#include "skmc.h"
#include "dsoffset.h"
#include "skdetsim_parm.h"
#include "guflag.h"
#include "skhead.h"
#include "geopma.h"
#include "geotnk.h"

      INTEGER IBUF(MAXRNG*2+1)
      REAL RBUF(MAXRNG*2+1)
      EQUIVALENCE (IBUF(1),RBUF(1))
      REAL TLIST(MAXPMA), QLIST(MAXPMA)
      INTEGER ICFLAG(MAXPMA), NUM(MAXPMA)
      INTEGER IOLD1, IOLD2
      INTEGER I, IERR, NW, NDATA 
      INTEGER NRING, IRING, ICABLE, NHIT, NWRD, IHTMX, IHITS
      REAL RDUMMY, RADRANG, QCNVSK_OLD, AOFFST
      REAL TNEW, RDARK, DARKT, PE1
      REAL GAIN(MAXPMA), SIGAIN
      REAL SATPARAM(MAXPMA)
      REAL EFF_FAC(MAXPMA)
      REAL PED(MAXPMA), SIGPED
      REAL SIGSPE
! made into parameters in the card file, see dscard.F, dsparma.h
! 20080213 -rvw
!      REAL GAIN_PLATEAU
!      REAL GAINFACT
      REAL GAINPLT(MAXPMA), PORTION
      REAL SMEAR, SIGSMEAR
      LOGICAL SET_GAIN, EXP_RSP, PLTEXP_RSP, GAUS_RSP
      LOGICAL SET_SAT, SET_EFF
!      DATA GAIN_PLATEAU/0.001/
!      DATA GAIN_PLATEAU/0.6/
      DATA SIGAIN, SIGPED, SIGSMEAR/0.25,0.12,0.2/
      DATA SIGSPE /0.43/
!      DATA GAINFACT/0.49/
!      DATA GAINFACT/0.9/
      DATA SET_GAIN/.TRUE./
      DATA SET_SAT/.TRUE./
      DATA SET_EFF/.TRUE./
      DATA EXP_RSP, PLTEXP_RSP, GAUS_RSP/.FALSE.,.TRUE., .FALSE./

      integer nhit_anti,index_anti(MAXPM),icabaz_tmp(MAXPMA),itaskz_tmp(MAXPMA)
      real taskz_tmp(MAXPMA),qaskz_tmp(MAXPMA)
      

      external GRNDM

      ! prevent OD simulation when requested -rvw
      if( ifanti.eq.0 ) return;
C
C     Set random gain and pedestal for each tube. We want the same
C     random tube gains so that certain gain checking histograms make
C     sense (the physics works if the tubes are scrambled for each
C     event, but gain distributions appear too narrow due to
C     averaging). Until such time as we make/measure a gain table, we
C     will set the random gains by storing the random number seed,
C     generating random numbers from a fixed seed, and then restoring
C     the old saved seed. We store/restore the seed because we do not
C     want to override the Monte Carlo use of a different seed each time
C     it runs; that is we trap our use of a fixed seed. -EK
C

**
      call vzero(itaskz_tmp,MAXPMA)
      call vzero(icabaz_tmp,MAXPMA)
      call vzero(taskz_tmp, MAXPMA)
      call vzero(qaskz_tmp, MAXPMA)
      call vzero(index_anti, MAXPMA)
      nhit_anti = 0
      
**      
      IF (SET_GAIN) THEN
         CALL GRNDMQ(IOLD1,IOLD2,0,'G')     ! store old seed
         CALL GRNDMQ(20030412,92612,1,'S')  ! initialize new seed
         CALL RNORMX(GAIN,MAXPMA,GRNDM)     ! generate random gains
         CALL RNORMX(PED,MAXPMA,GRNDM)      ! generate random pedestals
         DO I = 1, MAXPMA
            PED(I) = PED(I)*SIGPED
            IF (EXP_RSP.OR.PLTEXP_RSP) THEN
               GAIN(I) = 1.0 + GAIN(I)*SIGAIN
               GAINPLT(I) = GAIN_PLATEAU
            ENDIF
            IF (GAIN(I) .LT. 0.0) GAIN(I) = 0.0
         END DO
         CALL GRNDMQ(IOLD1,IOLD2,1,'S')     ! restore old seed
         SET_GAIN = .FALSE.
      END IF


      IF (SET_SAT) THEN
         CALL GRNDMQ(IOLD1,IOLD2,0,'G')     ! store old seed
         CALL GRNDMQ(20030442,921612,1,'S')  ! initialize new seed
         DO I = 1, MAXPMA
            CALL RNORMX(RDUMMY,1,GRNDM)

! 20080113 rvw
!           IF (SK_GEOMETRY.EQ.SK_I.OR.OD_TUBE_TYPE(I).EQ.2) THEN
            IF (OD_TUBE_TYPE(I).EQ.2) THEN

* Old tube
               SATPARAM(I) = ADSAPO(2)+ADSAPO(3)*RDUMMY

            ELSE
* New tube
               SATPARAM(I) = ADSAPN(2)+ADSAPN(3)*RDUMMY

            ENDIF
         ENDDO


         CALL GRNDMQ(IOLD1,IOLD2,1,'S')     ! restore old seed
         SET_SAT = .FALSE.
      END IF

      IF (SET_EFF) THEN
         CALL GRNDMQ(IOLD1,IOLD2,0,'G')     ! store old seed
         CALL GRNDMQ(20050442,921612,1,'S')  ! initialize new seed
         DO I = 1, MAXPMA

c Apply extra efficiency factor for old tubes.  Note this is only
c   correct if old tubes have smaller efficiency, apmtefr<1
c  (if not, should modify sgpmta.F)


            IF (OD_TUBE_TYPE(I).EQ.2) THEN
               EFF_FAC(I) = APMTEFR
            ELSE
               EFF_FAC(I) = 1.
            ENDIF

c Now apply efficiency factor for top/bottom/wall

            IF (XYZPMA(3,I).GE.ZPINTK) THEN
               EFF_FAC(I) = EFF_FAC(I)*APMTEFT
            ELSE IF (XYZPMA(3,I).LE.-ZPINTK) THEN
               EFF_FAC(I) = EFF_FAC(I)*APMTEFB
            ELSE
               EFF_FAC(I) = EFF_FAC(I)*APMTEFW
            ENDIF

         ENDDO
         CALL GRNDMQ(IOLD1,IOLD2,1,'S')     ! restore old seed
         SET_EFF = .FALSE.

      ENDIF

C
C     Initialize arrays
C
      CALL VZERO(TLIST,MAXPMA)  ! list of time for each photon hit
      CALL VZERO(QLIST,MAXPMA)  ! list of charge of each photon hit
      CALL VZERO(ICFLAG,MAXPMA) ! hit flag for each PMT
C
C     Read MCRAWTQA bank
C
      CALL KZGET2('MCRAWTQA',0,0,NW,IBUF)
      NRING=IBUF(1)
      DO I=1,NRING
         NUM(I) = IBUF(2*I+1)
      END DO
C
C-J.G. set conversion to 3.5ns per pe, then reset qcnvsk at end
C      qcnvsk = 3.5 
C-A.K. gain measured is 6.67 ns/pe or 13.34 Tcounts/pe
C
      qcnvsk_old = qcnvsk
      qcnvsk = 13.34
C
C     Loop over M.C. particles
C
      DO IRING=1,NRING
         DO I=1,NUM(IRING)

            CALL KZGET2('MCRAWTQA',IRING,I,NW,IBUF)

            ICABLE = IBUF(1)

            AOFFST=COMOFST
            RBUF(2) = RBUF(2)*1.e9 + AOFFST

            IF(ICFLAG(ICABLE).EQ.0)THEN
               ICFLAG(ICABLE)=1
               TLIST(ICABLE) =RBUF(2)*TCNVSK
            ELSE
               TNEW = RBUF(2)*TCNVSK
               IF(TNEW.LE.TLIST(ICABLE)) TLIST(ICABLE)=TNEW
            ENDIF
            QLIST(ICABLE) = QLIST(ICABLE)+FLOAT(IBUF(3))*QCNVSK

         END DO
      END DO
C
C     Create the TQA bank to hold the data
C     (If bank is already there, delete it before creating a new one)
      CALL KZBLOC('TQA',IERR)
      IF (IERR.NE.0) CALL KZBDEL('TQA')
      CALL KZBCR0('TQA',IERR)

      NHIT = 0
      NWRD = 0

      RDARK = DARKDSA(1)*FULWIDA*1.e-9
C
C     Loop over all tubes, add up hits for each tube
C
      DO 400 I=1,MAXPMA
C
C     For SK_I use exp+plt dark noise description for
C     all tubes (default values).  For SK_II use exp+plt description for
C     old tubes and gaussian description for new tubes. --JC
C     
C     If SK_II geometry reset logical parameters. --JC
C
!            write(*,*) 'OD_TUBE_TYPE', OD_TUBE_TYPE(I)
! 20080113 rvw, SK-I defaults to old tubes
!        IF (SK_GEOMETRY.EQ.SK_II) THEN
            IF (OD_TUBE_TYPE(I).EQ.1) THEN     ! new tubes
               EXP_RSP = .FALSE.
               PLTEXP_RSP = .FALSE.
               GAUS_RSP = .TRUE.
            ENDIF
            IF (OD_TUBE_TYPE(I).EQ.2) THEN     ! old tubes
               EXP_RSP = .FALSE.
               PLTEXP_RSP = .TRUE.
               GAUS_RSP = .FALSE.
            ENDIF
            IF (OD_TUBE_TYPE(I).EQ.0) THEN
               EXP_RSP = .FALSE.
               PLTEXP_RSP = .TRUE.
               GAUS_RSP = .FALSE.
            ENDIF
!        ENDIF
C
         CALL GRNDM(RDUMMY,1)
         IF (RDUMMY .LT. RDARK) THEN
C
C           Make dark hits:
C
            CALL GRNDM(RDUMMY,1)
            DARKT=RDUMMY*FULWIDA-DSANTT0+AOFFST
            IF (ICFLAG(I).EQ.0) THEN
C
C              This tube was not hit by Cherenkov light;
C              allow dark hit anywhere in time window.
C
               ICFLAG(I)=1
               QLIST(I)=QCNVSK
               TLIST(I)=DARKT*TCNVSK
            ELSE 
C
C              This tube was hit by Cherenkov light;
C              only allow dark hit preceding the trigger time of 1000
C
               IF ((TLIST(I)/TCNVSK-(DARKT-AOFFST)).le.1000) THEN
                  QLIST(I)=QCNVSK
                  TLIST(I)=DARKT*TCNVSK
               ENDIF
            ENDIF
         ENDIF

C
C        Tube has no light or dark hit, skip it
C
         IF (ICFLAG(I).EQ.0) GOTO 400
C
C        Prepare to add up hits for each tube
C
         IBUF(1) = I+2**15
         IHTMX = NINT(QLIST(I)/QCNVSK)
         QLIST(I)=0.

C
C        Loop over hits for this tube
C
         DO 401 IHITS=1,IHTMX

c Apply extra efficiency factor to this hit
            CALL GRNDM(RDUMMY,1)
            if (rdummy.gt.EFF_FAC(I)) goto 401

            IF(EXP_RSP) THEN
C              Draw single PE from exponential distribution
C              (Yes, that is a good description for the OD tubes -EK)
               CALL GRNDM(RDUMMY,1)
               PE1 =  - GAIN(I)*LOG(RDUMMY)
            ELSEIF(PLTEXP_RSP) THEN
C              Draw single PE from plateau&exponential distribution.
C              This is even better than a pure exponential distribution.
C              The joint point of the flat distribution at lower end
C              and the exponential distribution is an ajustable parameter,
C              GAIN_PLATEAU.
C                                       --- BU SuperK group 9/12/2002
               CALL GRNDM(RDUMMY,1)
               PORTION=GAINPLT(I)/(GAIN(I)+GAINPLT(I))
               IF (RDUMMY.LT.PORTION) THEN
                  CALL GRNDM(RDUMMY,1)
                  PE1 = GAINPLT(I)*RDUMMY
               ELSE
                  CALL GRNDM(RDUMMY,1)
                  PE1 = GAINPLT(I)-GAINFACT*GAIN(I)*LOG(RDUMMY)
               ENDIF               
            ELSEIF(GAUS_RSP) THEN
C     Draw single PE from Gaussian distribution.  -- Jodi Cooley
               CALL RNORMX(RDUMMY,1,GRNDM)
               PE1 = GAIN(I) + RDUMMY*SIGSPE
               
            ENDIF
            QLIST(I)=QLIST(I)+PE1
 401     CONTINUE
C
C        Apply PMT threshold 
C
c         write(6,*) 'NPMAFAC ADSCTHR',NPMAFAC,ADSCTHR

! 20080113 rvw
!        IF (SK_GEOMETRY.EQ.SK_II.AND.OD_TUBE_TYPE(I).EQ.1) THEN
         IF (OD_TUBE_TYPE(I).EQ.1) THEN
            IF (QLIST(I) .LE. ADSCTHR*ADSTHRF) THEN
               QLIST(I) = 0.
               TLIST(I) = 0.
               ICFLAG(I) = 0
               GO TO 400
            ENDIF
         ELSE
            IF (QLIST(I) .LE. ADSCTHR) THEN
               QLIST(I) = 0.
               TLIST(I) = 0.
               ICFLAG(I) = 0
               GO TO 400
            ENDIF
         ENDIF
C     
C        Final noise smearing after threshold for old tubes
C
CCS.03         CALL RNORML(RDUMMY,1)

! 20080113 rvw
!        IF (SK_GEOMETRY.EQ.SK_I.OR.OD_TUBE_TYPE(I).EQ.2) THEN
         IF (OD_TUBE_TYPE(I).EQ.2) THEN

            CALL RNORMX(RDUMMY,1,GRNDM)
            SMEAR = 1. + RDUMMY*SIGSMEAR
            QLIST(I) = QLIST(I)*SMEAR
         endif
C
C        Pedestal smearing
C
         QLIST(I) = QLIST(I) + PED(I)
C
C        MAX PE correction
C
C         CALL RDMOUT(IOSEED)
C         ISEED=I
C         CALL RDMIN(ISEED*147)

c New saturation curve 

! 20080113 rvw
!        IF (SK_GEOMETRY.EQ.SK_I.OR.OD_TUBE_TYPE(I).EQ.2) THEN
         IF (OD_TUBE_TYPE(I).EQ.2) THEN

c Old tubes
            QLIST(I) = ADSAPO(1)*QLIST(I)/
     &           (1+SATPARAM(I)*QLIST(I))
            
            
         ELSE

c New tubes
            QLIST(I) = ADSAPN(1)*QLIST(I)/
     &           (1+SATPARAM(I)*QLIST(I))
            

         ENDIF


C         CALL RDMIN(IOSEED)
C         write(*,*) "CORRECTED ANTI Charge(",I,")=",QLIST(I)

         QLIST(I)=QLIST(I)*QCNVSK
C         IBUF(2) = INT(TLIST(I)+0.5)*2**16 + INT(QLIST(I)+0.5)
         IF ((TLIST(I).gt.32767.).or.(TLIST(I).lt.-32767.)) THEN
           write(*,*) "TIME is OUT OF BOUNDS"
           GO TO 400
         ENDIF
C
C        Suppress negative values (need to understand real data behavior- ETK)
C
         IF (QLIST(I) .LT. 0.0) GO TO 400

         NHIT = NHIT + 1

C         write(*,*) "Cable #",I,":Q=",QLIST(I),"/T=",TLIST(I)
         IF (TLIST(I).ge.0.) THEN
           IBUF(2) = INT(TLIST(I)+0.5)*2**16 + INT(QLIST(I)+0.5)
         ELSE
C           IBUF(2) = INT((REAL(2**16)+TLIST(I))+0.5)*2**16
C     $              +INT(QLIST(I)+0.5)
           IBUF(2) = (INT(TLIST(I)+0.5))*2**16
     $              -INT(QLIST(I)+0.5)
         ENDIF

         NDATA = 2

         CALL KZREP0('TQA',1,'I',NDATA,IBUF)
         NWRD = NWRD+2

         nhit_anti = nhit_anti + 1
         ITASKZ_tmp(nhit_anti) =  int(TLIST(I)/tdc_unit_sk4) + 96485
         ICABAZ_tmp(nhit_anti) =  I

 400  CONTINUE

      IBUF(1) = NHIT
C      write(*,*) "---------------------NHIT=",NHIT
      IBUF(2) = NWRD
      NDATA = 2
      CALL KZREP0('TQA',0,'I',NDATA,IBUF)

C-J.G.
      QCNVSK = QCNVSK_OLD
      nhitaz = nhit

C ---
      RETURN
      END
