************************************************************************
*     -------------------------------------
      SUBROUTINE DSTRIG_MULTI
*     -------------------------------------
*
*     (Purpose)
*        Consider trigger condition etc.
*
*     (Input)
*        'DATATMP' bank        
*
*     (Output)
*        'TQ' bank
*
*     (Creation Date and Author)
*       1992.05.21 ; First version by K.S.Hirata
*
************************************************************************
      IMPLICIT NONE

#include "skhead.h"
#include "skdetsim_parm.h"
#include "skparm.h"
#include "dsparm.h"
#include "skmc.h"
#include "dsdata.h"
#include "guflag.h"
#include "dsoffset.h"

      INTEGER IBUF(MAXPM*2), ICTQ(3,MXPHDS)
      REAL RBUF(MAXPM*2), RCTQ(3,MXPHDS)
      EQUIVALENCE (IBUF,RBUF)
      EQUIVALENCE (ICTQ,RCTQ)
      INTEGER ICFLAG(MAXPM)
      CHARACTER*32 CNAME
      INTEGER NOTRIG,JPMT,IPMT,ISEGM,NSSEG,NDATA,ICNT,IMAX,IPO,J
      INTEGER NHIT,NWRD,IERR,I,ITCNT,IQCNT,ISSEG,ICSV
      REAL TMAX, TMIN, GMIN, GMAX
      REAL OFFSET,DELTATQ,OFFSETTQ
      COMMON /NOTRIG/NOTRIG
      INTEGER NHITR,MHIT,NWMAST,LSEGM,LDAT,NENT
      REAL Q,T
C-Y.H.-Support TQREAL BANK
      INTEGER*4 ITQBANK(3,MAXPM)
      REAL*4    RTQBANK(3,MAXPM)
      EQUIVALENCE(ITQBANK(1,1),RTQBANK(1,1))
      INTEGER*4 ITQRHD(3)
      REAL*4    RTQRHD(3)
      EQUIVALENCE (ITQRHD(1),RTQRHD(1))
#include <zbs.h>

      REAL SATMEAN,SATSIGMA,SATMEAN1,SATSIGMA1,SATMEAN2,SATSIGMA2
      DATA SATMEAN1/247.4/, SATSIGMA1/7.662/
      DATA SATMEAN2/220.9/, SATSIGMA2/14.16/

      PARAMETER(DELTATQ = -300.)

C --- Sort TDC's counts

      CNAME = 'DATATMP'
      ISEGM = 1
      CALL KZNSSG('DATATMP',ISEGM,NSSEG)

      IF(NSSEG.LE.0) THEN
         CALL KZBDEL('DATATMP')
         RETURN
      ENDIF

      DO 10 ISSEG=1,NSSEG
         
         CALL KZGET2('DATATMP',ISEGM,ISSEG,NDATA,IBUF)

         ICTQ(1,ISSEG)=IBUF(1)
         RCTQ(2,ISSEG)=RBUF(2)
         RCTQ(3,ISSEG)=RBUF(3)

 10   CONTINUE   

      CALL KZBDEL('DATATMP')
      ICNT = NSSEG
      IF(ICNT.GT.0) CALL SORTR(RCTQ,3,ICNT,2)

C --- Consider trigger condition

      IF(IFTRIG.EQ.1) THEN


         IMAX=ICNT-INT(TRIGDS)+1
         
         DO 30 IPO=1,IMAX
            
            CALL VZERO(ICFLAG,MAXPM)
            
            TMIN = RCTQ(2,IPO)
            TMAX = TMIN+TWINDS
            JPMT = 0

            DO 40 J=IPO,ICNT

               IF(RCTQ(2,J).GT.TMAX) THEN
                  IPMT = J-1
                  GOTO 50
               ENDIF
               IF(RCTQ(3,J).GE.DTHRSK) THEN
                 IF(ICFLAG(ICTQ(1,J)).EQ.0) JPMT= JPMT+1
                 ICFLAG(ICTQ(1,J))=ICFLAG(ICTQ(1,J))+1
               ENDIF
               
 40         CONTINUE
            IPMT=ICNT

 50         IF(JPMT.GE.INT(TRIGDS+0.5)) GOTO 60

 30      CONTINUE   

C --- If there is no trigger,
     
C         PRINT*,'NO TRIGGER !'

         GOTO 70

C --- If trigger is accepted,

 60      CONTINUE
C         PRINT*,'TRIGGER !'
         NOTRIG = NOTRIG + 1

         TRGDS = RCTQ(2,IPMT)

      ENDIF

 70   CONTINUE

C --- Apply the trigger condition
      GMAX = TRGDS + GLTRDS
      GMIN = GMAX - GATEDS
      OFFSET = OFFSET0 + (DTRGOFF - TRGDS)
      OFFSETTQ = OFFSET + DELTATQ

CC for OD
      COMOFST=OFFSET
C
C     If bank is already present delete it before creating new one
C     
      CNAME = 'TQ'
      CALL KZBLOC(CNAME,IERR)
      IF (IERR.NE.0) CALL KZBDEL(CNAME)
      CALL KZBCR0(CNAME,IERR)

      IF(ICNT.GT.0) CALL SORTI(RCTQ,3,ICNT,1)

C --  set saturation point
C     saturation charge of SK-II is changed mainly due to new pc2pe
C     changed by Shiozawa and data given by Nakayama
      IF (SK_GEOMETRY.EQ.SK_I) THEN
         SATMEAN = SATMEAN1
         SATSIGMA = SATSIGMA1
      ELSE IF (SK_GEOMETRY.EQ.SK_II) THEN
         SATMEAN = SATMEAN2
         SATSIGMA = SATSIGMA2
      ELSE
         SATMEAN = SATMEAN2
         SATSIGMA = SATSIGMA2
      ENDIF

      NHIT = 0
      NHITR = 0
      NWRD = 0
      ICSV = -1

      DO 110 I=1,ICNT
C         IF(IFTRIG.EQ.1.AND.
C     &        (RCTQ(2,I).LT.GMIN  .OR. RCTQ(2,I).GT.GMAX)) GOTO 110
         IF(RCTQ(2,I).LT.GMIN  .OR. RCTQ(2,I).GT.GMAX) GOTO 110

         ITCNT = 4095 - INT(TCNVSK*(RCTQ(2,I)+OFFSETTQ)+0.5)
         IQCNT = INT(RCTQ(3,I)*QCNVSK+0.5) + 1000


         IF(ITCNT.GE.4096.OR.ITCNT.LE.0) GOTO 110
         IF(IQCNT.GE.4096) IQCNT=4096
         IF(IQCNT.LE.0)    IQCNT=0

         NHIT = NHIT+1

C-    -Support TQREAL BANK
c-------------consider overflow
         Q = RCTQ(3,I)
         T = RCTQ(2,I)+OFFSET
         
         if (Q.ge.230.) then
            Call rngausk(satmean,satsigma,Q)
         Endif

         IF(ICTQ(1,I) .EQ. ICSV) THEN
             
            NWRD = NWRD+1
            IBUF(NWRD) = 2**29 + ITCNT*2**16 + 2**14 + 2**13 + IQCNT

            MHIT = MHIT + 1
            IF(MHIT.EQ.2) THEN
              WRITE(6,*) 'Warning A/B hit cable = ',ICTQ(1,I)
            ELSEIF(MHIT.GT.2) THEN
              IF(MOD(MHIT,2).EQ.1) THEN
                WRITE(6,*) 'Warning: double hit Cable = ',ICTQ(1,I),' A'
              ELSE
                WRITE(6,*) 'Warning: double hit Cable = ',ICTQ(1,I),' B'
              ENDIF
            ENDIF

         ELSE IF(ICTQ(1,I) .EQ. (ICSV+1)) THEN

            NWRD = NWRD+1
            IBUF(NWRD) = 2**30 + ITCNT*2**16 + IQCNT

            MHIT = 1
            NHITR = NHITR+1
            ITQBANK(1,NHITR) = ICTQ(1,I)
            RTQBANK(2,NHITR) = T
            RTQBANK(3,NHITR) = Q

         ELSE

            NWRD = NWRD+1
            IBUF(NWRD) = ICTQ(1,I)
            NWRD = NWRD+1
            IBUF(NWRD) = ITCNT*2**16 + IQCNT

            MHIT = 1
            NHITR = NHITR+1
            ITQBANK(1,NHITR) = ICTQ(1,I)
            RTQBANK(2,NHITR) = T
            RTQBANK(3,NHITR) = Q

         END IF


         ICSV = ICTQ(1,I)

 110  CONTINUE

      IF(NWRD.GT.MAXPM*2) THEN
         WRITE(6,*)'ERROR in DSTRIG'
         STOP
      ENDIF


      CALL KZREP0('TQ',1,'I',NWRD,IBUF)

      IBUF(1) = NHIT
      IBUF(2) = NWRD
      NDATA = 2
      CALL KZREP0('TQ',0,'I',NDATA,IBUF)

C
C     If bank is already present, delete it before creating new one
C     
      CALL KZBLOC('TQREAL',IERR)
      IF (IERR.NE.0) CALL KZBDEL('TQREAL')
      NWMAST = 3
      LSEGM  = 3
      CALL KZBCR1('TQREAL',NWMAST,'1I 1F 1I',NHIT,LSEGM,'1I 2F',IERR)
      IF (IERR.NE.0) THEN
         WRITE(0,*) 'DSTRIG_MULTI: ERROR IN CREATING TQREAL BANK'
         RETURN
      ENDIF

      ITQRHD(1)=NHITR
      RTQRHD(2)=2.055
      ITQRHD(3)=1
      ISEGM = 0
      NDATA = 3
      CALL KZREP1('TQREAL',ISEGM,NDATA,ITQRHD)
      
      IF (NHITR.gt.0) THEN
         CALL KZLDAT('TQREAL',LDAT)
         NENT=3*NHITR
         CALL LBFCOPY(ITQBANK,IZBS(LDAT+3+1),NENT)
      ENDIF

C ---
      RETURN
      END

