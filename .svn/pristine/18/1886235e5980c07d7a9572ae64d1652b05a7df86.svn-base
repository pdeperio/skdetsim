***********************************************************************
*     ----------------------------------
      SUBROUTINE sgbsdist(cost1,PRREFL,scatter,cost2)
*     ----------------------------------
*     
*     Black Sheet Reflection Distribution
*
*     ( purpose )
*       generate scattering of black sheet data
*
*     ( input ) 
*       cost1 : incident angle (from bs normal)
*       PRREFL: Fresnel fraction of reflected light
*
*     ( output ) 
*       scatter: type of reflection (lambert,specular,phong)
*       cost2  : reflected angle 
*
*     ( creation date and author )
*       P. Cravens  2005/03/21
*
*     ( comment )
*
***********************************************************************

      implicit none

      integer scatter

      real rdummy
      real PRREFL
      real lam,spec,phong,fact
      real cost1,cost2

      real pi
      parameter(pi=3.14159)

      real a,b

      fact = 10.                ! Phong shiny factor
      fact = 1./(fact+1.)

      lam = 0.015/PRREFL
      if (cost1.le.(2./3.)) then
         spec = -1.2*cost1+1.
      else
         spec = 0.2
      endif
      if (lam+spec.le.1.) then
         phong = 1.-lam-spec
      else
         phong = 0.0
         lam = 1.-spec
      endif

      a = lam
      b = lam+phong

      call grndm(rdummy,1)
      if (rdummy.gt.0.0.and.rdummy.le.a) then
         call grndm(rdummy,1)
         cost2 = sqrt(rdummy)
         scatter = 1            ! lambert
      else if (rdummy.gt.a.and.rdummy.le.b) then
         call grndm(rdummy,1)
         cost2 = rdummy**fact
         scatter = 2            ! phong
      else if (rdummy.gt.b.and.rdummy.le.1.0) then
         cost2 = cost1
         scatter = 3            ! specular
      endif


      RETURN
      END

