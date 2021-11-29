      subroutine skphonuc
c*****************************************************************
c     Generate final state of photonuclear interaction.
c     Just kill the photon for now, but provide interface for more
c     complicated final states
c*****************************************************************
c   6/2012 H. A. Tanaka
c*****************************************************************
      implicit none
#include "geant321/gcking.inc"
#include "geant321/gconsp.inc"
#include "geant321/gttrak.inc"
#include "geant321/gctrak.inc"

      integer j
*
*   Kill photon
*
      istop = 1

* Code stolen from sgpintneut.F absorption
* add proton at rest to the stack

C     Generate placeholder proton to store in secondaries stack
      NGKINE=1

      DO J=1,3
         GKIN(J,NGKINE)=0
         GPOS(J,NGKINE)=VECT(J)
      end do
      GKIN(5,NGKINE) = 14
      TOFD(NGKINE) = 0.
      IFLGK(NGKINE) = -1

      return
      end
