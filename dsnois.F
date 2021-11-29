************************************************************************
*     -------------------------------------
      SUBROUTINE DSNOIS
*     -------------------------------------
*
*     (Purpose)
*        Dark noise of PMT's
*        Channel conversion of T & Q
*
*     (Creation Date and Author)
*       1992.05.28 ; First version by K.S.Hirata
*       1994.02.19 ; Added for dark noise study by K.S.Hirata
*       1995.02.03 ; Modified by Koshio
*
************************************************************************

#include "skparm.h"
#include "dsparm.h"
#include "skmc.h"
#include "dsdata.h"
#include "skmisch.h"
#include "dsoffset.h"

      CHARACTER*32 CNAME
      DIMENSION JBUF(3),SBUF(3),KBUF(MAXRNG*2+3)
      EQUIVALENCE(JBUF(1),SBUF(1))
C moved to dsinip Sep. 17 2008 and changed to OFFSET0
Cc      PARAMETER (OFFSET = 1000.)
C      PARAMETER (OFFSET = 540.)
Cc      PARAMETER (OFFSET = 840.)

CC for OD
      COMOFST=OFFSET0+300.

C --- Define bank name

      CNAME = 'MCPHOTON'
      CALL KZBLOC('MCPHOTON',LBIR)
      IF(LBIR.EQ.0) CNAME = 'MCRAWTQ'

C --- Read data in 'MCPHOTON' or 'MCRAWTQ' bank 
C                             and create 'MCTQ' bank

      CALL KZGET2(CNAME,0,0,NDATA,KBUF)

      NRNGDS = KBUF(1)
      NPHADS = 0
      DO 1 IRING=1,NRNGDS
         IPARDS(IRING) = KBUF(2*IRING)
         NPHRDS(IRING) = KBUF(2*IRING+1)
         NPHADS = NPHADS+NPHRDS(IRING)
 1    CONTINUE

C --- Read random number seed stored in bank and
C     use it again.

      IF(CNAME.EQ.'MCPHOTON') THEN
         ISEED1 = KBUF(2+2*NRNGDS)
         ISEED2 = KBUF(3+2*NRNGDS)
c         CALL RDMIN(ISEED)
         call grndmq(ISEED1,ISEED2,1,'S')
      ENDIF

      NWMAST = (NRNGDS+1)*2+1
      NSEGM = NRNGDS+1
      LSEGM = 3
      CALL KZBCR2('MCTQ',NWMAST,' *I ',
     &     NSEGM,LSEGM,' 1I 1F 1I',IERR)

      ICNT = 0
          
      DO 10 IRING=1,NRNGDS

         DO 20 IPH=1,NPHRDS(IRING)

            CALL KZGET2(CNAME,IRING,IPH,NDATA,JBUF)

            ICNT = ICNT+1

            IF(ICNT.GT.MXPHDS) THEN
               WRITE(6,*) 'Please use HIGH flag'
               STOP
            ENDIF

            SBUF(2)  = SBUF(2) + OFFSET0

            IF(CNAME.EQ.'MCPHOTON') THEN
               JBUF(3)  = 1
            ELSE IF(CNAME.NE.'MCRAWTQ') THEN
               WRITE(6,*) 'ERROR in DSNOIS ; 1'
               STOP
            ENDIF

            IF(ICNT.EQ.1) THEN
               TMIN = SBUF(2)
               TMAX = SBUF(2)
            ELSE
               IF(SBUF(2).LE.TMIN) TMIN = SBUF(2)
               IF(SBUF(2).GE.TMAX) TMAX = SBUF(2)
            ENDIF
               
            NDATAN = 3
            CALL KZREP2('MCTQ',IRING,IPH,NDATAN,JBUF)

 20      CONTINUE

 10   CONTINUE   

      IF(ICNT.NE.NPHADS) THEN
         WRITE(6,*) 'ERROR in DSNOIS ; 2'
         STOP
      ENDIF

      IF(ICNT.GE.1) THEN
c         TMIN = TMIN - (GATEDS - GLTRDS)
         TMIN = TMIN - DEADDS
         TMAX = TMAX + GLTRDS
      ELSE
         TMIN = 0.
         TMAX = GATEDS
c         IF(IMWATR.NE.-1)
c     &        WRITE(6,*),'Dark current ONLY; Dark rate is ',darkds,' Hz'
      END IF   

C --- Add dark noise of PMT's as (nrngds+1)-th ring

      EXPD = DARKDS*(TMAX-TMIN)*(1.E-09)

      IPH=0
       
      ISEGM = NRNGDS+1
      NDATA = 3

      DO 30 IPMT=1,MAXPM

         call grndm(dummy,1)
         IF(DUMMY.GE.EXPD) GOTO 30
         if(imis(ipmt).ne.0) goto 30
         IPH=IPH+1
          
         IF(ICNT+IPH.GT.MXPHDS) THEN
            WRITE(6,*) 'Please use HIGH flag'
            STOP
         ENDIF

         JBUF(1) = IPMT
         call grndm(dummy,1)
         SBUF(2) = DUMMY*(TMAX-TMIN)+TMIN
         JBUF(3) = 1

         CALL KZREP2('MCTQ',ISEGM,IPH,NDATA,JBUF)
          
 30   CONTINUE   

      ISEGM = 0
      ISSEG = 1
      NDATA = 2*(NRNGDS+1)+1
      
      KBUF(2*NRNGDS+2) = 0
      KBUF(2*NRNGDS+3) = IPH

      CALL KZREP2('MCTQ',ISEGM,ISSEG,NDATA,KBUF)

      NPHRDS(NRNGDS+1)=IPH

C ---

      RETURN
      END
