************************************************************************
*     -------------------------------------
      subroutine sgpope_read
*     -------------------------------------
*
*     (Purpose)
*       To read in values of ICAM measurement for long wavelength
*         water absorption
*
************************************************************************

      implicit none

#include "pope_inc.h"

      integer u,i,maxarray,isave/0/
      character*200 fname
      save isave

      if(isave.eq.0) then
         i = 0
         u = 94
         maxarray = 140
ccyt         fname = '/home/sklowe/skrep/const/popefry_data.dat'
         call findconsts(fname,'popefry_data.dat','popefry_data.dat',2)
         open(u, FILE=fname)

         do i=1,maxarray
            read(u,*) wlength(i), coeff(i)
            coeff(i) = 100.*coeff(i)
         enddo

         close(u)

         isave = 1
      endif

      return
      end
