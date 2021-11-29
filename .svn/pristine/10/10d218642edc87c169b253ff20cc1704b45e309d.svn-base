      real function skphonucxs(egam)
C**********************************************************************
C   Given gamma energy in GeV, returns 16O photonuclear cross section
C   in mb.
C         H. A. Tanaka June 2012
C**********************************************************************
      implicit none

#include "apcard.h"

      integer i, ibin
      real egam, offset
      logical first

*     Cross section file 
      character*200 pnxsecfile
      data          pnxsecfile/'o16_phonuc_smoothxsec.dat'/
      character*128 fname

*     Numbe of points
      integer       nxsecpts
      parameter(nxsecpts = 44)

*     data points
      real          pne(nxsecpts), pnxs(nxsecpts)

      data first / .true. /

c   If first call, read in xsec table.

      if(first)then
         first = .false.

         call findconsts(fname,'Photonuclear XS table',pnxsecfile,0)

         write(*,*) ' *****SKPHONUCXS -- Reading file ', 
     +        pnxsecfile
         write(*,*) ' *****SKPHONUCXS -- Scale factor', 
     +        pnscale

         open (unit=1, file=fname, status='OLD')

         do i = 1, nxsecpts
            read (1,*) pne(i), pnxs(i)
c   Convert: E to GeV
            pne(i) = pne(i)/1000.

         enddo

         close (1)
      endif

*  Linearly interpolate: Set to 0 if outside the range
      ibin = 0

      if     ( egam .lt. pne(1)        )then
         skphonucxs = 0.0
      else if( egam .ge. pne(nxsecpts) )then
         skphonucxs = 0.0
      else 
         do i = 2, nxsecpts
            if(egam .lt. pne(i))then
               ibin = i
               goto 10
            endif
         enddo

 10      continue
         
         offset     = ( egam - pne(ibin-1) ) / ( pne(ibin) - pne(ibin-1) )
         skphonucxs = (1.0-offset) * pnxs(ibin-1) + offset * pnxs(ibin)
         skphonucxs = skphonucxs * pnscale

      endif

      end


      real function skphonucstep( egam )
c******************************************************
c   Compute photonuclear interaction mean free path (cm)
c  for photon of energ egam (GeV).
c   On first call, prints parameter and xsec summary.
c                  H.A. Tanaka 06/2012
c******************************************************
      implicit none

#include "gconsp.inc"
#include "gumate.h"
#include "gttmed.inc"
#include "gctmed.inc"

      integer i, IDUM1, IDUM2, IDUM3, numat
      real dens, e, egam, rate, RDUM1, RDUM2, RDUM3, RDUM4, 
     *     RDUM5, RDUM6, RDUM7, skphonucxs, xsec
      logical first
      CHARACTER*20 CDUM

      data first / .true. /
      save

ccc      write(6,*)'numed, avo, big =', numed,avo,big 

*   Find tracking medium number.
      CALL GFTMED( NUMED, CDUM, NUMAT, IDUM1, IDUM2, RDUM1, RDUM2, 
     *                        RDUM3, RDUM4, RDUM5, RDUM6, 
     *                        RDUM7, IDUM3 ) ! Find material of this med.

*   If we are indeed in water, get the density and proceed
      if(numat .eq. IMTWTR)then

         CALL GFMATE(numat, cdum, RDUM1, RDUM2, dens,
     +               RDUM3, RDUM4, RDUM5, IDUM1 )

C   Get photonuclear cross section in mb)
C  
         if(first)then    ! Write table to logfile
            first = .false.

            xsec = skphonucxs(egam)

            write(*,*)
            write(*,*) ' *****SKPHONUCSTEP --',
     &            ' Table of photonuclear cross section'
            write(*,*) ' EGamma (MeV)     xsec (mb)'
 10         format(3x,f6.1,3x,f10.3)

            do i=0,150,5
               e = float(i)
               write(*,10) e, skphonucxs(0.001*e)
            enddo

         endif

* Now get the cross section and calculate rate, stepsize

         xsec = skphonucxs(egam)
         rate = avo*dens*xsec/16./1000.
*                              /1000 because avo assumes barn, not mb

         if(rate .gt. 0.)then
            skphonucstep = 1./rate
         else
            skphonucstep = BIG
         endif
      else
         skphonucstep = BIG
      endif
ccc      write(6,*)'Egam, xsec, uphotonucStep = ', Egam,xsec,uphotonucStep

      return
      end
