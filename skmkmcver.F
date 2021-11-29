C
      subroutine skmkmcver
C
C     make detsim version bank
C 

c     !!!!! IMPORTANT NOTICE (S.Mine on Aug.4, 2011) !!!!
c     "apdslibver" used by skgetwater/water.aveX is hard-coded here
c     this file and water.aveX should be consistently updated
c     please contact to me before modifying the version here


      implicit none

#include "skhead.h"
      integer lenocc
      external lenocc
      integer NW, NSEGM, LSEGM, ISEG, IERR, IDATA(10)
      Real RDATA(10)
      CHARACTER*32 CNAME, REVTEXT, REVTEXT2
      equivalence(idata(1), rdata(1))


C---- For backward compatibility
C     APMCVER
C     
C     -- SET MONTE CARLO VERSION NUMBER
C     
C     RDATA(1) ; apdetsim version
C     
C
      IF(SK_GEOMETRY.eq.1) then
         RDATA(1)=5.2
      ELSEIF(SK_GEOMETRY.eq.2) then
         RDATA(1)=6.1
      ELSEIF(SK_GEOMETRY.eq.3) then
         RDATA(1)=7.2
      ELSEIF(SK_GEOMETRY.eq.4) then
         RDATA(1)=8.3
csm      ELSEIF(SK_GEOMETRY.eq.5) then ! just copy. Need update J.Kameda 2019/10/29
csm         RDATA(1)=8.3         
      ELSEIF(SK_GEOMETRY.eq.5) then ! add by S.Mine (Jul.9, 2021)
         RDATA(1)=8.4
      ELSE
         print*, 'SK_GEOMETRY is not set correctly in skmkmcver.F'
         stop
      ENDIF
      RDATA(2)=APDETSIM_VERSION
      RDATA(3)=0.
      RDATA(4)=0.
      RDATA(5)=0.
      RDATA(6)=0.
C     
C     -- MAKE 'MCHEAD ' BANK AND WRITE
C     
      CNAME='APMCVER'
      ISEG=0
      NW=8

      CALL KZBCR1(CNAME,NW,'-I',0,0,' ',IERR)       
      IF (IERR.NE.0) GOTO 1000
      CALL KZREP1(CNAME,ISEG,NW,RDATA)      
C----


C     SKMCVER
C     
C     -- SET MONTE CARLO VERSION NUMBER
C     
C     -- MAKE 'SKMCVER ' BANK AND WRITE
C
      CNAME='SKMCVER'

      NW=0
      NSEGM=3
      LSEGM=8
      CALL KZBCR1(CNAME,NW,' ',NSEGM,LSEGM,'8F',IERR)
      IF (IERR.NE.0) GOTO 1000
C     
C     segment #0
C     RDATA(1) ; apdetsim version
C     
C
      IF(SK_GEOMETRY.eq.1) then
         RDATA(1)=5.2
      ELSEIF(SK_GEOMETRY.eq.2) then
         RDATA(1)=6.1
      ELSEIF(SK_GEOMETRY.eq.3) then
         RDATA(1)=7.2
      ELSEIF(SK_GEOMETRY.eq.4) then
         RDATA(1)=8.3
csm      ELSEIF(SK_GEOMETRY.eq.5) then ! just copy from SK-IV
csm         RDATA(1)=8.3         
      ELSEIF(SK_GEOMETRY.eq.5) then ! add by S.Mine (Jul.12, 2021)
         RDATA(1)=8.4
      ELSE
         print*, 'SK_GEOMETRY is not set correctly in skmkmcver.F'
         stop
      ENDIF
      RDATA(2)=APDETSIM_VERSION
      RDATA(3)=0.
      RDATA(4)=0.
      RDATA(5)=0.
      RDATA(6)=0.
      RDATA(7)=0.
      RDATA(8)=0.

      ISEG=1
      CALL KZREP1(CNAME,ISEG,LSEGM,IDATA)
C
C     segment #1
C     RDATA(1) ; skdetsim version
C     
C
      RDATA(1)=SKDETSIM_VERSION
C
C     Write current CVS version of this file (skmkmcver.F)
C     into the SKMCVER bank. The text string beginning
C     with $Revision is automatically filled in during cvs update.
C     There is no planned use for this now, but it may be
C     helpful in future forensics of forgotten MC files. 
C     The nice feature is that it is fully automatic. ETK (11/2007)
C
      REVTEXT = "$Revision: 1.10 $"
      REVTEXT2 = REVTEXT(11:LENOCC(REVTEXT)-1)
      READ(REVTEXT2,*) RDATA(2)
C
      RDATA(3)=0.
      RDATA(4)=0.
      RDATA(5)=0.
      RDATA(6)=0.
      RDATA(7)=0.
      RDATA(8)=0.

      ISEG=2
      CALL KZREP1(CNAME,ISEG,LSEGM,IDATA)

      RETURN
C     
C     ++ ERROR MESSAGE
C     
 1000 WRITE(6,900)
  900 FORMAT(' *** ERROR *** ( BANK CREATION IN SKMKMCVER )')
      RETURN
      END
