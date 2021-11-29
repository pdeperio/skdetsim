******************************************************
*     -----------------
      SUBROUTINE SGGVUS_ESCALE(ierr)
*     -----------------
*
*     (Purpose)
*        User routine for Generating Vector for 
*        Energy Scale
*
*   Added J. Hignight July 16, 2012
*     Only called when built with setenv ESCALE 1
*
******************************************************

      implicit none
#include "vcwork.h"
#include "geotnk.h"

      INTEGER IPART
      REAL POS(3),PMOM(3)
      real PP
      real dir(3),posmu(3)
      real mom
      integer ierr
      integer icall
      data icall/0/
      character filename*100

      real mup_mun_ratio,prob
      parameter(mup_mun_ratio=1.37)

      call getenv('TXTFILE',filename)

      if (icall.eq.0) then
         call geoset
cc         open(88,file='vcmake.dat',status='old')
         open(88,file=filename,status='old')
         print *, 'open:',filename
         icall=1
      endif
         
      read(88,*,end=20) mom,dir,posmu

c      print *,'dirmu',dir(1),dir(2),dir(3)

      pos(1)=posmu(1)
      pos(2)=posmu(2)
      pos(3)=posmu(3)
      PP = mom
      PMOM(1)=dir(1) * PP
      PMOM(2)=dir(2) * PP
      PMOM(3)=dir(3) * PP

c      IPART =13
      call grndm(prob,1)

      ipart = 13
      if (prob .le. (mup_mun_ratio/(mup_mun_ratio + 1.0))) ipart = -13

      write(*,901) ipart
      write(*,902) pos
      write(*,903) pmom

 901    format('SGGVUS: Particle ID:',I6)
 902      format('SGGVUS: Position   :',F12.5,F12.5,F12.5)
 903        format('SGGVUS: Momentum   :',F12.5,F12.5,F12.5)

      CALL VCMAKE(POS,PMOM,IPART)

      ierr = 0
      RETURN

 20      ierr = 1
      close(88)
      RETURN
      END
