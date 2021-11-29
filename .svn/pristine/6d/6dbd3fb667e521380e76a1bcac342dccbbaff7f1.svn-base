C moved to its own file Jan 2003
C Why?  CERN type definitions only can be included once per file.
C Fortran, on the other hand (unlike C), needs the types declared once per
C routine.  So one is forced to split routines into one per file.

      subroutine save_muon_path
      implicit none

C Alec, Sept. 2002
C Save muon kinematics if a muon enters or leaves ID or OD
C Note that we check for volume change via medium type.  This is
C not so good, because lots of sub-volumes in the detector have
C the same medium type.  So, we check to see if the "water" changed
C flavor from OD to Dead Space or ID.  Thus, we check only upon entering 
C the new volume, since we can remember what the last medium was, but 
C checking ahead to the new medium upon exit is a pain the the ass.

#include "geant321/gctrak.inc"
#include "geant321/gckine.inc"
#ifdef UPMU
#include "gumed_u.h"
#else
#include "gumed.h"
#endif
#include "geant321/gctmed.inc"
#include "muvtx.h"

      integer iii

C Entering ID?
      if ((numed.eq.imdwti).and.(idmu.eq.0.or.idmu.eq.1)
     &     .and.(numed.ne.last_medium)) then
         idmu=idmu+10
c         print *,'*** ID entrance point ',
c     &           ipart,itra,(vect(iii),iii=1,7),idmu
         do iii=1,3
            mupos(1,iii)=vect(iii)
            mudir(1,iii)=vect(iii+3)
         enddo
         mumom(1)=vect(7)
         last_medium = numed
C Entering OD?  Can only start here, if it's already been in ID it won't
C get tagged as "entering" OD
      else if ((numed.eq.imdwto).and.(idmu.eq.0)
     &        .and.(numed.ne.last_medium)) then
         idmu=idmu+1
c         print *,'*** OD entrance point ',
c     &           ipart,itra,(vect(iii),iii=1,7),idmu
         do iii=1,3
            mupos(4,iii)=vect(iii)
            mudir(4,iii)=vect(iii+3)
         enddo
         mumom(4)=vect(7)
         last_medium = numed
C Exiting ID into dead space?
      else if ((numed.eq.imdwtm).and.(idmu.eq.10.or.idmu.eq.11)
     &        .and.(numed.ne.last_medium)) then
         idmu=idmu+100
c         print*,'*** ID exit point ',
c     &          ipart,itra,(vect(iii),iii=1,7),idmu
         do iii=1,3
            mupos(2,iii)=vect(iii)
            mudir(2,iii)=vect(iii+3)
         enddo
         mumom(2)=vect(7)
         last_medium = numed
C Exiting OD into master volume?
C Might never have "entered" OD, if the muon started 
C life in the ID.  Note this would happen twice, once upon exit of muon
C out the other side of detector, if we didn't trap at idmu.lt.1000
C However, by requiring idmu.gt.1, we require that the particle has
C entered ID.  Yes, this fails for OD only muons, but we don't care 
C about those. 
      else if ((numed.eq.imdmst).and.
     &        ((idmu.gt.1).and.(idmu.lt.1000))
     &        .and.(numed.ne.last_medium)) then
         idmu=idmu+1000
c         print *,'*** OD exit point ',
c     &            ipart,itra,(vect(iii),iii=1,7),idmu
         do iii=1,3
            mupos(5,iii)=vect(iii)
            mudir(5,iii)=vect(iii+3)
         enddo
         mumom(5)=vect(7)
         last_medium = numed
      end if                    ! of which volume was being exited

      return
      end
