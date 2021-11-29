      subroutine seg(ipart, dir, time, nphot, thopen)
*
* Record track segments from particle ipart, in direction dir(3),
* generating nphot cerenkov photons at the opening angle thopen
*
      implicit none
      integer ipart, nphot, ihier
      real dir(*), thopen, time
#include "aptrkseg.h"
      integer treet(0:255), treep(0:255)
      common /tree/ treet, treep
      integer invtree(0:15,0:15)
      common /inv_tree/ invtree

      real theta, phi
      real vmod
C
C+ Internal segment buffer
C
      integer ipix(100*mxseg), npho(100*mxseg), ippix(100*mxseg)
      real thpix(100*mxseg), phpix(100*mxseg)
      real aopen(100*mxseg), tpix(100*mxseg)
      integer npix
      integer i, nscan, ihmax
      integer ipout(mxseg)

      logical initialized
      data initialized/.false./

C
C+  Entry point
C
      theta = acos(dir(3)/vmod(dir,3))
      if (dir(1).eq.0 .and. dir(2).eq.0) then
         phi = 0.
      else
         phi = atan2(dir(2),dir(1))
      endif
      do while (phi.lt.0.)
         phi = phi + 8. * atan(1.)
      enddo
      phi = mod(phi,8.*atan(1.))
      call ang2pix_hier(ihier, igloo, theta, phi)
      if (npix.ge.100*mxseg) then ! Overflow
         write (*,*) 'seg: maximum number of segments exceeded'
         return
      endif
      npix = npix + 1
      ipix(npix) = ihier
      npho(npix) = nphot
      ippix(npix) = ipart
      aopen(npix) = thopen
      tpix(npix) = time
      return

      entry packseg

      nseg = 0
      nphtot = 0
      do i = 1, npix
         nphtot = nphtot + npho(i)
         nscan = 1
         do while (nscan.le.nseg .and. ipout(nscan) .ne. ipix(i))
            nscan = nscan + 1
         enddo
         if (nscan.le.nseg) then
            nphseg(nscan) = nphseg(nscan) + npho(i)
            if (tseg(nscan).gt.tpix(i)) tseg(nscan) = tpix(i)
         else
            nseg = nseg + 1
            call pix2ang_hier(ipix(i), igloo, thseg(nseg), phseg(nseg))
            nphseg(nseg) = npho(i)
            ipout(nseg) = ipix(i)
            tseg(nseg) = tpix(i)
         endif
      enddo
      return

      entry iniseg

      npix = 0
      ihmax = 0
      if (.not.initialized) then
        call init_tree
        call init_inv_tree
        initialized = .TRUE.
      endif

      return
      end
