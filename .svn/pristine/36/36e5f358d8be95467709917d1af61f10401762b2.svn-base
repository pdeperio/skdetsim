      function popefrysg(lamda)

      implicit none

#include "pope_inc.h"

      integer i, j
      real lamda, term1, term2, popefrysg

      if (lamda.lt.380) then
         popefrysg = 0.
         goto 200
      endif

      do i = 1,139
         if (lamda.ge.wlength(i).and.lamda.lt.wlength(i+1)) then
            j = i
            goto 100
         endif
      enddo

 100  continue

      term1 = (lamda - wlength(j+1))/(wlength(j)-wlength(j+1))
      term2 = (lamda - wlength(j))/(wlength(j+1)-wlength(j))
      popefrysg = coeff(j)*(term1) + coeff(j+1)*term2

 200  continue

      return
      end
