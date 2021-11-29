      subroutine skdcy_muminus(p_e,e_eminus)

      implicit none

#include "weakparm.h"

      integer z
      real prob, height, pp, p_e(3), e_eminus
      parameter(z=8)

 334     continue

         Call grndm(prob,1)
         pp = (78.-e_mass)*prob + e_mass
         e_eminus = sqrt(pp**2+e_mass**2)*0.001

         call muminus(z,pp,height)
         Call grndm(prob,1)

         if(prob .gt. height) goto 334

         CALL RNDIR(p_e)
         p_e(1) = p_e(1)*pp*0.001
         p_e(2) = p_e(2)*pp*0.001
         p_e(3) = p_e(3)*pp*0.001

      return
      end
