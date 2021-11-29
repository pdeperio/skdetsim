************************************************************************
*     -------------------------------------
      SUBROUTINE DSTRIG
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

#include "skparm.h"
#include "dsparm.h"
#include "skmc.h"
#include "dsdata.h"
#include "guflag.h"

      DIMENSION IBUF(MAXPM*2),RBUF(MAXPM*2)
      EQUIVALENCE (IBUF(1),RBUF(1))
      DIMENSION ICTQ(3,MXPHDS),RCTQ(3,MXPHDS)
      EQUIVALENCE (ICTQ,RCTQ)
      DIMENSION ICFLAG(MAXPM)

      character*32 cname

      COMMON /NOTRIG/NOTRIG

C --- Sort TDC's counts

      CNAME = 'DATATMP'
      ISEGM = 1
      CALL KZNSSG('DATATMP',ISEGM,NSSEG)

      if(nsseg.le.0) then
         CALL KZBDEL('DATATMP')
         return
      endif

      DO 10 ISSEG=1,NSSEG
         
         CALL KZGET2('DATATMP',ISEGM,ISSEG,NDATA,IBUF)

         ICTQ(1,ISSEG)=IBUF(1)
         RCTQ(2,ISSEG)=RBUF(2)
         RCTQ(3,ISSEG)=RBUF(3)

 10   CONTINUE   

      CALL KZBDEL('DATATMP')
      ICNT = NSSEG

C --- Consider trigger condition

      IF(IFTRIG.EQ.1) THEN

         IF(ICNT.GT.0) CALL SORTI(RCTQ,3,ICNT,2)

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
               IF(RCTQ(3,J).GE.DTHRSK
     &              .AND. ICFLAG(ICTQ(1,J)).EQ.0) JPMT= JPMT+1
               IF(RCTQ(3,J).GE.DTHRSK) 
     &              ICFLAG(ICTQ(1,J))=ICFLAG(ICTQ(1,J))+1
               
 40         CONTINUE
            IPMT=ICNT

 50         IF(JPMT.GE.INT(TRIGDS+0.5)) GOTO 60

 30      CONTINUE   

C --- If there is no trigger,
     
C         PRINT*,'NO TRIGGER !'
         GMIN = 999999.
         GMAX = -1.

         GOTO 70

C --- If trigger is accepted,

 60      CONTINUE
C         PRINT*,'TRIGGER !'
         NOTRIG = NOTRIG + 1

         GMAX = RCTQ(2,IPMT) + GLTRDS
         GMIN = GMAX - GATEDS

      ENDIF

C --- Apply the trigger condition

 70   CONTINUE

      CNAME = 'TQ'
      CALL KZBCR0('TQ',IERR)

      NHIT = 0
      NWRD = 0
      IF(ICNT.GT.0) CALL SORTI(RCTQ,3,ICNT,1)

      DO 110 I=1,ICNT

         IF(IFTRIG.EQ.1.AND.
     &        (RCTQ(2,I).LT.GMIN  .OR. RCTQ(2,I).GT.GMAX)) GOTO 110

ccc         ITCNT = 1000+INT(TCNVSK*(GMAX - RCTQ(2,I))+0.5)
ccc         ITCNT = INT(RCTQ(2,I)*TCNVSK+0.5)
ccc         IQCNT = INT(RCTQ(3,I)*QCNVSK+0.5)
         ITCNT = 4095 - INT(TCNVSK*RCTQ(2,I)+0.5)
         IQCNT = INT(RCTQ(3,I)*QCNVSK+0.5) + 1000

c         IF(ITCNT.GE.4096) ITCNT=4096
c         IF(ITCNT.LE.0)    ITCNT=0
         IF(ITCNT.GE.4096.OR.ITCNT.LE.0) GOTO 110
         IF(IQCNT.GE.4096) IQCNT=4096
         IF(IQCNT.LE.0)    IQCNT=0

         NHIT = NHIT+1

         IF(ICTQ(1,I) .EQ. ICSV) THEN
             
            NWRD = NWRD+1
            IBUF(NWRD) = 2**29 + ITCNT*2**16 
     &                 + 2**14 + 2**13 + IQCNT

         ELSE IF(ICTQ(1,I) .EQ. (ICSV+1)) THEN

            NWRD = NWRD+1
            IBUF(NWRD) = 2**30 + ITCNT*2**16 + IQCNT

         ELSE

            NWRD = NWRD+1
            IBUF(NWRD) = ICTQ(1,I)
            NWRD = NWRD+1
            IBUF(NWRD) = ITCNT*2**16 + IQCNT

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

C ---
      RETURN
      END

