***
*** Modified 9/02 Alec Habig
*** Updated to also store Saji's MUVTX information
***
*** Added 7/00 Alec Habig
*** If the particle is a penetrating particle - mu or charged pi
*** then store information about where and why it finished tracking
C
      subroutine SKMKAPMCSTOP(VERTEX,TIME,INTCODE,PARTID,TRACK,REPLACE)
      IMPLICIT NONE
#include "gt2pd.h"
C
      REAL    VERTEX(3)
      REAL    TIME
      INTEGER INTCODE
      INTEGER PARTID
      INTEGER TRACK
      LOGICAL REPLACE
C
      INTEGER    LSEGM
      PARAMETER (LSEGM  = 7)
      INTEGER    NWMAST
      PARAMETER (NWMAST = 1)
C
C     Make penetrating particle stopping info bank
C     Add the arguments to a new segment
C 
      INTEGER IDATA(LSEGM)
      REAL    RDATA(LSEGM)
      EQUIVALENCE(IDATA(1),RDATA(1))
      CHARACTER*32 CNAME
      INTEGER IERR
      INTEGER NUMVERTS /0/
C NUMVERTS needs to be static
      SAVE    NUMVERTS


      CNAME='APMCSTOP'
      

C If this is the first call, then create the master bank.  
C
      IF (NUMVERTS.EQ.0) THEN
         CALL KZBCR1(CNAME,NWMAST,'1I',0,LSEGM,'4F 3I ',IERR)      
         IF (IERR.NE.0) GOTO 1000
      ENDIF
C
C Put the number of segments present into the master bank data field
C
      IF (.NOT.REPLACE) THEN
         NUMVERTS = NUMVERTS + 1
         CALL KZREP1(CNAME,0,NWMAST,NUMVERTS)      
      ENDIF
C     
C     -- Fill the new SEGMENT
C
C        IND+ 1 : X                - vertex of stopping point
C           + 2 : Y
C           + 3 : Z
C           + 4 : TIME             - of particle stopping
C           + 5 : INTERACTION CODE - why did the particle stop?
C           + 6 : PDG PARTICLE CODE
C           + 7 : TRACK NUMBER
C     
      RDATA( 1)=VERTEX(1)
      RDATA( 2)=VERTEX(2)
      RDATA( 3)=VERTEX(3)
      RDATA( 4)=TIME*1.0e9
      IDATA( 5)=INTCODE
      IDATA( 6)=IPDG(PARTID)
      IDATA( 7)=TRACK
      CALL KZREP1(CNAME,NUMVERTS,LSEGM,IDATA(1))         

      RETURN
C     
C     ++ ERROR MESSAGE
C     
 1000 WRITE(6,900)
  900 FORMAT(' *** ERROR *** ( BANK CREATION IN SKMKAPMCSTOP )')
      RETURN

      entry skapmcclr
      numverts = 0
      call clear_stop_pointers
      return

      END


      subroutine check_for_stops
      implicit none

*** Added 7/00 Alec Habig
*** If the particle is a penetrating particle - mu or charged pi
*** then store information about where and why it finished tracking
#include "geant321/gckine.inc"
#include "geant321/gctrak.inc"

      INTEGER mecs
      real*4  vertex(3)
      integer last_itra /0/
      save    last_itra 
      integer last_ipart /0/
      save    last_ipart 
      logical replace_last

C Pions are hard, as they can make more of the same particle in hadronic
C interactions.  I have not figured out how to tell which exact particle
C is being tracked, so the logic of following pion chains gets munged.
      if ((ipart.ge.5).and.(ipart.le.9).and. 
     &     (ipart.ne.7)) then
*** We don't care about uninteresting reasons for stopping the tracking,
*** like secondary production.  We want things that:
***   nmec  = 0  (exited SK)
***   lmec = 1  (left volume) - always first.  The rest are last
***   lmec = 5  (decay)
***   lmec = 12 (hadronic interaction)
***   lmec = 13 (hadronic elastic coherent scatter)
***   lmec = 16 (nuclear absorption)
***   lmec = 19 (hadronic elastic incoherent scatter)
***   lmec = 20 (hadronic inelastic scatter)
***   lmec = 21 (muon-nuclear interaction)
***   lmec = 30 (stop)   - not always last
         if (nmec.ne.0) then
            mecs = lmec(nmec)
            if (lmec(1).eq.1) mecs = 1
         else
            mecs = 1
         endif
         if ((mecs.eq.1).or.
     &       (mecs.eq.5).or.
     &       (mecs.eq.12).or.
     &       (mecs.eq.13).or.
     &       (mecs.eq.16).or.
     &       (mecs.eq.19).or.
     &       (mecs.eq.20).or.
     &       (mecs.eq.21).or.
     &       (mecs.eq.30)) then
C Check for the same particle coming along as the last track.
C This happens, for example, if the pion didn't get stopped by the 
C hadronic interaction which got us to this point before
            if ((itra.eq.last_itra).and.(ipart.eq.last_ipart)) then
               replace_last = .TRUE.
            else
               replace_last = .FALSE.
            endif 
            last_itra   = itra
            last_ipart  = ipart
C            print *,vect(1),vect(2),vect(3),tofg
c            print *,itra,istop,ipart
c            print *,nmec,(lmec(i),i=1,nmec)
c            if (REPLACE_LAST)  print *,'Supercedes last entry'
c            print *
            vertex(1) = vect(1)
            vertex(2) = vect(2)
            vertex(3) = vect(3)
            call skmkapmcstop(vertex,tofg,mecs,ipart,itra,replace_last)
         endif
      endif

      return

      entry clear_stop_pointers

      last_itra   = 0
      last_ipart  = 0

      return
      end


C Below is stuff for Saji's MUVTX bank

C Raise MUVTX bank containing info from muvtx common
      subroutine mkmubnk
      implicit none

#include "muvtx.h"
      
      real rbuff(7)
      integer ibuff(2)
      integer ierr,i

      call kzbcr1('MUVTX',2,'1I',nvtx,7,'7F',ierr)

      ibuff(1)=idmu
      ibuff(2)=nvtx
      call kzrep1('MUVTX',0,2,ibuff)

      do i=1,nvtx
        rbuff(1)=mupos(i,1)
        rbuff(2)=mupos(i,2)
        rbuff(3)=mupos(i,3)
        rbuff(4)=mudir(i,1)
        rbuff(5)=mudir(i,2)
        rbuff(6)=mudir(i,3)
        rbuff(7)=mumom(i)

        call kzrep1('MUVTX',i,7,rbuff)
      enddo

      return
      end


C Clear the muvtx common
      subroutine clrmuver
      implicit none

#include "muvtx.h"
      integer i,j

      idmu=0
      first_id_mu=0
      last_medium = 0
      do i=1,nvtx
        mumom(i)=0.
        do j=1,3
          mupos(i,j)=0.
          mudir(i,j)=0.
        enddo
      enddo

      return
      end
