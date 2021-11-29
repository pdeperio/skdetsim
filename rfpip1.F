      SUBROUTINE rfpip1

#include "geant321/gctrak.inc"
#include "geant321/gcvolu.inc"
#include "geant321/gcking.inc"

#include "linpara.h"

      REAL ABC(3)

      ISTOP=1

*  Calculate the incident angle

      aa=abs(rpip1o-sqrt((vect(1)-sorx)**2+(vect(2)-sory)**2))
      bb=abs(vect(3)-sorz)
      cc=abs(rpip1i-sqrt((vect(1)-sorx)**2+(vect(2)-sory)**2))
      if(aa.lt.bb.and.aa.lt.cc)then
         rr=(vect(1)-sorx)**2+(vect(2)-sory)**2
         abc(1)=(vect(1)-sorx)/sqrt(rr)
         abc(2)=(vect(2)-sory)/sqrt(rr)
         abc(3)=0.
         if(aa.gt.1.e-2)write(*,*) 'rfpip1 aa=',aa
      else if(bb.lt.cc)then
         abc(1)=0.
         abc(2)=0.
         abc(3)=-1.
         if(bb.gt.1.e-2)write(*,*) 'rfpip1 bb=',bb         
      else 
         rr=(vect(1)-sorx)**2+(vect(2)-sory)**2
         abc(1)=-(vect(1)-sorx)/sqrt(rr)
         abc(2)=-(vect(2)-sory)/sqrt(rr)
         abc(3)=0.
         if(cc.gt.1.e-2)write(*,*) 'rfpip1 cc=',cc         
      endif

      COST1 = -(ABC(1)*VECT(4) + ABC(2)*VECT(5) + ABC(3)*VECT(6))

      iloop=0
 20   continue
      COST1 = -(ABC(1)*VECT(4) + ABC(2)*VECT(5) + ABC(3)*VECT(6))
      
      if(cost1.lt.0.)then
         if(iloop.eq.0.and.abc(3).eq.-1.)then
            rr=(vect(1)-sorx)**2+(vect(2)-sory)**2
            if(sqrt(rr).gt.(rpip1o-rpip1i)/2.+rpip1i)then
               abc(1)=(vect(1)-sorx)/sqrt(rr)
               abc(2)=(vect(2)-sory)/sqrt(rr)
               abc(3)=0.            
               iloop=iloop+1
               goto 20            
            else
               abc(1)=-(vect(1)-sorx)/sqrt(rr)
               abc(2)=-(vect(2)-sory)/sqrt(rr)
               abc(3)=0.            
               iloop=iloop+1
               goto 20            
            endif
         else if(iloop.eq.0.and.abc(3).eq.0.)then
            abc(1)=0.
            abc(2)=0.
            abc(3)=-1.
            iloop=iloop+1
            goto 20
         endif
         write(*,*)'error in refpip1 cost1=',cost1
         write(*,*)'x,y,z=',vect(1)-sorx,vect(2)-sory,vect(3)-sorz
         write(*,*)vect(4),vect(5),vect(6),abc(1),abc(2),abc(3)
      endif

      call refstn(cost1,abc)
      
      return
      end

