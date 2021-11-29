      SUBROUTINE rfpip2

#include "geant321/gctrak.inc"
#include "geant321/gcvolu.inc"
#include "geant321/gcking.inc"

#include "linpara.h"

      REAL ABC(3)

      ISTOP=1

*  Calculate the incident angle

      aa=abs(rpip2o-sqrt((vect(1)-sorx)**2+(vect(2)-sory)**2))
      bb=abs(vect(3)-sorz-2.*pip1l)
      if(aa.lt.bb)then
         rr=(vect(1)-sorx)**2+(vect(2)-sory)**2
         abc(1)=(vect(1)-sorx)/sqrt(rr)
         abc(2)=(vect(2)-sory)/sqrt(rr)
         abc(3)=0.
         if(aa.gt.1.e-2)write(*,*) 'rfpip2 aa=',aa
      else
         abc(1)=0.
         abc(2)=0.
         abc(3)=-1.
         if(bb.gt.1.e-2)write(*,*) 'rfpip2 bb=',bb         
      endif

      COST1 = -(ABC(1)*VECT(4) + ABC(2)*VECT(5) + ABC(3)*VECT(6))
      if(cost1.lt.0.)write(*,*)'error in rfpip2 cost1=',cost1

      call refstn(cost1,abc)
      
      return
      end
c
c
c
