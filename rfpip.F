      SUBROUTINE rfpip

#include "geant321/gctrak.inc"
#include "geant321/gcvolu.inc"
#include "geant321/gcking.inc"

#include "linpara.h"

      REAL ABC(3)

      ISTOP=1

*  Calculate the incident angle

      rr=(vect(1)-sorx)**2+(vect(2)-sory)**2
      abc(1)=0.
      abc(2)=0.
      abc(3)=-1.

      COST1 = -(ABC(1)*VECT(4) + ABC(2)*VECT(5) + ABC(3)*VECT(6))
      if(cost1.lt.0.)write(*,*)'error in rfpip cost1=',cost1

      call refstn(cost1,abc)
      
      return
      end
c
c
