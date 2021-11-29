************************************************************************
*     -------------------------------------
      SUBROUTINE DSNOIS_MULTI
*     -------------------------------------
*
*     (Purpose)
*        Dark noise of PMT's
*        Channel conversion of T & Q
*
*
*     (Creation Date and Author)
*
*    2008.07.03 ; First version by 
*                 Jun Kameda based on dsnois.F
*
*                 Variable gate width is supported.
*                 Multi dark hit in single PMT is supported.
*
************************************************************************
      implicit none

#include "skhead.h"
#include "skparm.h"
#include "skdetsim_parm.h"
#include "dsparm.h"
#include "skmc.h"
#include "dsdata.h"
#include "skmisch.h"
#include "geant321/gcflag.inc"

      CHARACTER*32 CNAME
      integer JBUF(3), KBUF(MAXRNG*2+3)
      real sbuf(3)
      EQUIVALENCE(JBUF(1),SBUF(1))
      real dummy,tmin,tmax,expd
      integer ierr,ierr_poi, ndark,lbir,iring,iseed1,iseed2,i,isseg
      integer nwmast, nsegm, lsegm, icnt, iph, ndata, ndatan, isegm, ipmt


C --- Define bank name

      CNAME = 'MCPHOTON'
      CALL KZBLOC('MCPHOTON',LBIR)
      IF(LBIR.EQ.0) then
        WRITE(6,*) 'ERROR in DSNOIS_multi: No MCPHOTON bank. Stop'
        STOP 1
      ENDIF

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

C
C --- Read random number seed stored in bank and
C     use it again.

       if ((TRACK_PHOTON .eq. TRACK_QE_OVER_ABS_PHOTON) .or. 
     & ((TRACK_PHOTON .eq. TRACK_ALL_PHOTON).and.(SK_GEOMETRY.ge.SK_IV))
     &)then
         ISEED1 = KBUF(2+2*(NRNGDS+1))
         ISEED2 = KBUF(3+2*(NRNGDS+1))
       else
         ISEED1 = KBUF(2+2*NRNGDS)
         ISEED2 = KBUF(3+2*NRNGDS)
       endif

c         CALL RDMIN(ISEED)
         call grndmq(ISEED1,ISEED2,1,'S')

c      NWMAST = (NRNGDS+1)*2+1 
c      NSEGM = NRNGDS+1
c      LSEGM = 3
c      CALL KZBCR2('MCTQ',NWMAST,' *I ',
c     &     NSEGM,LSEGM,' 1I 1F 1I',IERR)


      ICNT = 0
          
      DO 10 IRING=1,NRNGDS

         DO 20 IPH=1,NPHRDS(IRING)

            CALL KZGET2(CNAME,IRING,IPH,NDATA,JBUF)

            ICNT = ICNT+1

            IF(ICNT.GT.MXPHDS) THEN
               write(6,*) '******************************************'
               WRITE(6,*) 'DSNOIS_MULTI: No. of photon exceed MXPHDS.'
               WRITE(6,*) 'Please Use Other tracking mode. Or,'
               WRITE(6,*) 'Please extend MXPHDS by hand.'
               WRITE(6,*) 'Stop.'
               write(6,*) '******************************************'
               STOP
            ENDIF

c            SBUF(2)  = SBUF(2) + OFFSET

c               JBUF(3)  = 1

            IF(ICNT.EQ.1) THEN
               TMIN = SBUF(2)
               TMAX = SBUF(2)
            ELSE
               IF(SBUF(2).LE.TMIN) THEN
                 TMIN = SBUF(2)
               ELSEIF(SBUF(2).GE.TMAX) THEN
                 TMAX = SBUF(2)
               ENDIF
            ENDIF
               
c            NDATAN = 3
c            CALL KZREP2('MCTQ',IRING,IPH,NDATAN,JBUF)

 20      CONTINUE

 10   CONTINUE   

      IF(ICNT.NE.NPHADS) THEN
         WRITE(6,*) 'ERROR in DSNOIS_MULTI ; 2'
         STOP
      ENDIF


c set time interval in which the dark hit will be generated.


      if (SK_GEOMETRY .ge. SK_IV) then
c
c    Dark hit will be generated from (N-2)-th block to (N+3)-th block.
c    LE events (-5usec,35usec) will be contained within it.
c    We should enlarge for T2K events(-500 usec, 500usec).
c    Tmin, Tmax is measured by GEANT time. 
c

      tmin = -trg_offset - twid_block*2
#ifdef NEUTRON
      tmax = twid_block
#else
      tmax = (twid_block - trg_offset) + twid_block*3
#endif

      elseif (
     & (SK_GEOMETRY .ge. SK_I) .and. (SK_GEOMETRY .le. SK_III) 
     & ) then
      
        IF(ICNT.GE.1) THEN
           TMIN = TMIN - DEADDS - GATEDS
           TMAX = TMAX + DEADDS + GATEDS
        ELSE
           TMIN = -DEADDS - GATEDS
           TMAX =  DEADDS + GATEDS
        END IF   

      else
         write(6,*) 'ERROR: DSNOISE_MULTI: SK_GEOMETRY is irregal !!'
      endif

         
C --- Add dark noise of PMT's as (nrngds+1)-th ring

      IPH=0
c
c      EXPD = DARKDS*(TMAX-TMIN)*(1.E-09)
       
      ISEGM = NRNGDS+1
c      NDATA = 3
      NDATA = 2

      DO 30 IPMT=1,MAXPM

c     skip missing PMTs

      if(imis(ipmt).ne.0) goto 30

c     calculate expected # of dark hit within the window.
      if (DOPMTDARKRT.eq.0) then
         EXPD = DARKDS*(TMAX-TMIN)*(1.E-09)
      else
         EXPD = DARKRT(IPMT)*(TMAX-TMIN)*(1.E-09)
      endif
 
c      call poissn(expd,ndark,ierr_poi)
      call gpoiss(expd,ndark,1)

      if (ndark .gt. 0) then
         do 100 i = 1,ndark

         IPH=IPH+1
          
         IF(ICNT+IPH.GT.MXPHDS) THEN
               write(6,*) '******************************************'
               WRITE(6,*) 'DSNOIS_MULTI: No. of photon exceed MXPHDS in '
               WRITE(6,*) 'dark hit generation. '
               write(6,*) 'Please Use Other tracking mode. Or,'
               WRITE(6,*) 'Please extend MXPHDS by hand.'
               WRITE(6,*) 'Stop.'
               write(6,*) '******************************************'
               STOP
         ENDIF

c         JBUF(1) = IPMT
c         call grndm(dummy,1)
c         SBUF(2) = DUMMY*(TMAX-TMIN)+TMIN
c         JBUF(3) = 1
c
c         CALL KZREP2('MCTQ',ISEGM,IPH,NDATA,JBUF)

         JBUF(1) = IPMT
         call grndm(dummy,1)
         SBUF(2) = DUMMY*(TMAX-TMIN)+TMIN
c         JBUF(3) = 1

         CALL KZREP2('MCPHOTON',ISEGM,IPH,NDATA,JBUF)

 100     continue

      endif
          
 30   CONTINUE   

      ISEGM = 0
      ISSEG = 1
      NDATA = 2*(NRNGDS+1)+1
      
      KBUF(2*NRNGDS+2) = 0
      KBUF(2*NRNGDS+3) = IPH

c      CALL KZREP2('MCTQ',ISEGM,ISSEG,NDATA,KBUF)
      CALL KZREP2('MCPHOTON',ISEGM,ISSEG,NDATA,KBUF)

      NPHRDS(NRNGDS+1)=IPH

C ---

      RETURN
      END
