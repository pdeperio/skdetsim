      subroutine refstn(cost1,abc)
#include "geant321/gctrak.inc"
#include "geant321/gcvolu.inc"
#include "geant321/gcking.inc"
#include "geant321/gctmed.inc"

#include "skdetsim_parm.h"

#include "linpara.h"

      PARAMETER(FIXSTEP = 5.e-3)
      real vecnew(3),abc(3),pos(3)
	real ang_inc,wave_len
	real pi_const
	parameter(pi_const=3.14159265)
	real ref
	integer matid
	real norm
      
*  Calculate the ratio of reflection on stainless

      NGPHOT = 0

*  Reflection !
***********************************************************
*  This part is new
*  mesured reflectivity is calculated
	ang_inc = acos(cost1) 
	ang_inc = 180.*ang_inc/pi_const
	wave_len = 1.23e-6/xphot(7,1)
	if(numed.eq.60)then
	  matid=1 ! Ti window 
	  call getreflect_lin(matid,ang_inc,wave_len,ref)
	else if(numed.eq.52.or.
     &        numed.eq.54)then
	  matid = 2 ! sample B
	  call getreflect_lin(matid,ang_inc,wave_len,ref)
	else if(numed.eq.53.or.numed.eq.55.or.numed.eq.50)then
	  matid = 3 ! sample C
	  call getreflect_lin(matid,ang_inc,wave_len,ref)
	end if

	prrefl = ref/100.
*  Sei hansha
	call grndm(rand_rand,1)
c	write(*,*)ang_inc,wave_len,prrefl,rand_rand
      IF(rand_rand.LT.PRREFL) THEN
c	write(*,*)'hansha',ang_inc,wave_len,prrefl,rand_rand
         nsref=nsref+1
         NGPHOT =1
         VECNEW(1) = VECT(4) + 2.*ABC(1)*COST1
         VECNEW(2) = VECT(5) + 2.*ABC(2)*COST1
         VECNEW(3) = VECT(6) + 2.*ABC(3)*COST1
         norm = 
     &       sqrt(vecnew(1)**2+vecnew(2)**2+vecnew(3)**2)
         vecnew(1) = vecnew(1)/norm 
         vecnew(2) = vecnew(2)/norm 
         vecnew(3) = vecnew(3)/norm 
      ENDIF

      IF(NGPHOT.GE.1) THEN

         STEP = FIXSTEP
         iloop=0
 10      continue
         do i=1,3
            pos(i)=vect(i)+step*vecnew(i)
         enddo

         call gmedia(pos,num)
         if(num.ne.2)then
            if(numed.eq.num.and.iloop.le.5)then
               step=step*2.
               iloop=iloop+1
               goto 10
            else if(numed.ne.num.and.iloop.le.5)then
               step=step/2.
               iloop = iloop+1
               goto 10
            endif
            write(*,*)'photon is not reflected to the water' 
            write(*,*)'numed=',num,cost1,abc(3),step
            write(*,*)'x,y,z=',vect(1)-sorx,vect(2)-sory,vect(3)-sorz
            goto 99
         endif

         DO I = 1,3
            XPHOT(I+3,1) = VECNEW(I)
            XPHOT(I,1) = pos(i)
         ENDDO
         XPHOT(7,1)  = VECT(7)
c         XPHOT(11,1) = TOFG + STEP/CLIGHT*refwtr
         XPHOT(11,1) = TOFG

**** porlarization is not yet

         if (OUTPUT_SCATTABLE.eq.1) then
            CALL GSKPHO(-ISTORY-20000000) ! track photon scattering and reflections
         else
            CALL GSKPHO(0)
         endif

         NGPHOT=0
      ENDIF

 99   continue
      RETURN
      END









