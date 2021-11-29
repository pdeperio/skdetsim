***********************************************************************
*     -----------------
      SUBROUTINE SGACRYLIC(ACRTRAN)
*     -----------------
*
*     (Purpose)
*       Cherenkov photons entering and exiting acrylic volume.
*     Reflection and refraction at water/acrylic interface.
*
*     (Input)
*       ACRTRAN(3)::    Position of the acrylic sphere's 
*                        reference system in MARS
*                        from GTRAN(*,NLEVEL) when photon entered acrylic
*                        or created in acrylic
*     (Output)
*       none
*
*     (Creation Date and Author)
*       2002.07.11 ; First version by Chris Sterner
*
************************************************************************
      IMPLICIT NONE
*includes
#include "geant321/gctrak.inc"
#include "geant321/gcking.inc"
#ifdef UPMU
#include "gumed_u.h"
#else
#include "gumed.h"
#endif
#include "geant321/gctmed.inc"

#include "skhead.h"
#include "skdetsim_parm.h"

      INTEGER normalconv, ienter, reflect
      REAL normal(3), costh,  prob, newvect(3), wlen
      REAL THETA1, SREF, PREF, THETA2, length, PRREFL
      REAL FIXSTEP, BACKSTEP
      PARAMETER (FIXSTEP=.01000)
      PARAMETER (BACKSTEP=.01000)
      REAL ACRTRAN(3)
      
      INTEGER*4 IREFSTRY

      integer nloop/0/
      real prvvect(3),rprv
      save nloop,prvvect

C for polarization
      REAL COSPOL, pdir(3), tanth,
     &     pdir_new(3), sdir_new(3), e_1s, e_1p,
     &     s1, s2, p1, p2, e_2s, e_2p
      REAL N1, N2, acrnsg, REFSG
C *******

      IREFSTRY = 0

C     to avoid infinit loop between acrylic and water
      if(nloop.eq.0) then
         prvvect(1) = vect(1)
         prvvect(2) = vect(2)
         prvvect(3) = vect(3)
         nloop = 1
      else
         rprv = sqrt((prvvect(1)-vect(1))**2 + (prvvect(2)-vect(2))**2
     &        + (prvvect(3)-vect(3))**2)
         if(rprv.le.1.e-3) then
            nloop=nloop+1
         else
            nloop=0
         endif
      endif
      if(nloop.ge.10000) then
         print*,'more than 10000 steps abandoned'
         nloop = 0
         return
      endif

C     determine if entering (IENTER=1) or exiting (IENTER=0) acrylic
C     assumed in new medium, INWVOL.EQ.1 
      if(numed.eq.IMDACR) then
         IENTER=1
      else if(numed.eq.imdwti.OR.numed.eq.IMDWTS) then
         IENTER=0
      else
         WRITE(6,*) 'ERROR in sgacrylic 2', numed
      endif

*    Find normal to acrylic cover surface
      call SGACRNOR(VECT,NORMAL,ACRTRAN)

*     normalize vect and normal
      length = sqrt(vect(4)**2+vect(5)**2+vect(6)**2)
      if(length .ne. 1.0 .or. length .ne. 0.0) then
         vect(4)= vect(4)/length
         vect(5)= vect(5)/length
         vect(6)= vect(6)/length
      endif
      length = sqrt(normal(1)**2+normal(2)**2+normal(3)**2)
      if(length .ne. 1.0 .or. length .ne. 0.0) then
         normal(1)= normal(1)/length
         normal(2)= normal(2)/length
         normal(3)= normal(3)/length
      endif      

*     cos between normal and incident  light
      costh = vect(4)*normal(1)+vect(5)*normal(2)+vect(6)*normal(3)
      normalconv = 1
      if(costh.lt.0) then
         costh=abs(costh)
         normalconv = -1
      endif

      wlen=1237.8E-09 / GEKIN
      THETA1=ACOS(costh)
C     DETERMINE refl/refr prob. and angle of refraction
      call SGREFACR(IENTER,WLEN,THETA1,SREF,PREF,THETA2)

*     Reflection or refraction at interface

*     calculate the ratio of pwave (after SK-V)
      if(SK_GEOMETRY .ge. SK_V) then
         tanth = 1./costh**2 - 1.
         if(costh.eq.0..or.costh.eq.1..or.tanth.le.0)then
C     Refraction
C     determine direction of refracted photon in MRS, newvect(3)
            call sgacrefr(vect(4),normal,normalconv,costh,theta2,newvect)
            reflect=-1
            goto 1000
         else
            pdir(1)=(vect(4)-1.0*FLOAT(normalconv)*normal(1)/costh)/sqrt(1./costh**2-1.)
            pdir(2)=(vect(5)-1.0*FLOAT(normalconv)*normal(2)/costh)/sqrt(1./costh**2-1.)
            pdir(3)=(vect(6)-1.0*FLOAT(normalconv)*normal(3)/costh)/sqrt(1./costh**2-1.)
         endif
      endif

      if(THETA2.gt.900) then    !total internal reflection
         PRREFL=2.0             !greater than 1 to ensure reflection
      else
         if(SK_GEOMETRY.le.SK_IV) then
            PRREFL= PREF*cos(THETA1)**2 + SREF*sin(THETA1)**2 ! bug! (THETA1 is incident angle)

         else                   ! fix bug after SK-V
            IF (abs(vect(4)*polar(1)+vect(5)*polar(2)+vect(6)*polar(3))
     &           .gt.0.0001) then
               PRREFL= PREF/2. + SREF/2. ! if no polarization
            ELSE
               COSPOL = (polar(1)*pdir(1)+polar(2)*pdir(2)+polar(3)*pdir(3))**2
               CALL GRNDM(PROB,1)
               if(PROB.lt.COSPOL)then
                  PRREFL = PREF
               else
                  PRREFL = SREF
               endif
            ENDIF
         endif
      endif

      CALL GRNDM(PROB,1)
      IF(PROB.LE.PRREFL) THEN
C       Reflection
C        determine direction of reflected photon in MRS, newvect(3)
         call sgacrefl(vect(4),normal,normalconv,costh,newvect)
         reflect=1

      ELSE
C       Refraction
C        determine direction of refracted photon in MRS, newvect(3)
         call sgacrefr(vect(4),normal,normalconv,costh,theta2,newvect)
         reflect=-1
      ENDIF

 1000 continue

*store photon
      ngphot=1
C     Normalize new direction
      length = sqrt(newvect(1)**2+newvect(2)**2+newvect(3)**2)
      if(length .ne. 1.0 .or. length .ne. 0) then
         newvect(1)= newvect(1)/length
         newvect(2)= newvect(2)/length
         newvect(3)= newvect(3)/length
      endif
      if(reflect.eq.-1) then 
C     position is FIXSTEP amount advanced to ensure placed 
C       in correct, new, medium
         STEP = FIXSTEP
         xphot(1,1)=vect(1)+STEP*newvect(1)
         xphot(2,1)=vect(2)+STEP*newvect(2)
         xphot(3,1)=vect(3)+STEP*newvect(3)
         IREFSTRY=-ISTORY
      else if(reflect.eq.1) then
C     position is BACKSTEP amount reversed along original direction to
C       ensure placed in correct, old, medium
         STEP = BACKSTEP
         xphot(1,1)=vect(1)-STEP*vect(4)
         xphot(2,1)=vect(2)-STEP*vect(5)
         xphot(3,1)=vect(3)-STEP*vect(6)
         IREFSTRY=-ISTORY-20000000
      endif
C     new direction
      xphot(4,1)=newvect(1)
      xphot(5,1)=newvect(2)
      xphot(6,1)=newvect(3)
C     same wavelength
      xphot(7,1)=vect(7)
C
C     no polarization before SK-IV
C     polarization after process from SK-V
C
      if(SK_GEOMETRY .ge. SK_V) then
         sdir_new(1) = vect(5)*normal(3)-vect(6)*normal(2)
         sdir_new(2) = vect(6)*normal(1)-vect(4)*normal(3)
         sdir_new(3) = vect(4)*normal(2)-vect(5)*normal(1)
         length = sqrt(sdir_new(1)**2+sdir_new(2)**2+sdir_new(3)**2)
         if(length .ne. 1.0 .or. length .ne. 0) then
            sdir_new(1)= sdir_new(1)/length
            sdir_new(2)= sdir_new(2)/length
            sdir_new(3)= sdir_new(3)/length
         endif
         pdir_new(1) = newvect(2)*sdir_new(3)-newvect(3)*sdir_new(2)
         pdir_new(2) = newvect(3)*sdir_new(1)-newvect(1)*sdir_new(3)
         pdir_new(3) = newvect(1)*sdir_new(2)-newvect(2)*sdir_new(1)
         length = sqrt(pdir_new(1)**2+pdir_new(2)**2+pdir_new(3)**2)
         if(length .ne. 1.0 .or. length .ne. 0) then
            pdir_new(1)= pdir_new(1)/length
            pdir_new(2)= pdir_new(2)/length
            pdir_new(3)= pdir_new(3)/length
         endif

C     FROM WATER TO ACRYLIC
         IF ( IENTER .EQ. 1 ) THEN
            N1 = REFSG(WLEN)
            N2 = acrnsg(WLEN)
C     FROM ACRYLIC TO WATER
         ELSE IF ( IENTER .EQ. 0 ) THEN
            N1 = acrnsg(WLEN)
            N2 = REFSG(WLEN)
         ENDIF
         s1 = N1 * costh
         s2 = N2 * cos(theta2)
         p1 = N1 * cos(theta2)
         p2 = N2 * costh

         e_1s = polar(1)*sdir_new(1)+polar(2)*sdir_new(2)+polar(3)*sdir_new(3)
         e_1p = sqrt((polar(1)-e_1s*sdir_new(1))**2
     &        + (polar(2)-e_1s*sdir_new(2))**2
     &        + (polar(3)-e_1s*sdir_new(3))**2)

C     for refraction
         if(reflect .eq. -1) then
            if((s1+s2).ne.0. .and. (p1+p2).ne.0.) then
               e_2s = 2*s1*e_1s/(s1+s2)
               e_2p = 2*s1*e_1p/(p1+p2)
               length = sqrt(e_2s**2 + e_2p**2)
               e_2s = e_2s/length
               e_2p = e_2p/length
            endif

C     for reflection
         else if(reflect.eq.1) then
            if((s1+s2).ne.0. .and. (p1+p2).ne.0.) then
               e_2s = 2*s1*e_1s/(s1+s2) - e_1s
               e_2p = (2*s1*e_1p/(p1+p2))*(N2/N1) - e_1p
               length = sqrt(e_2s**2 + e_2p**2)
               e_2s = e_2s/length
               e_2p = e_2p/length
            endif

         endif

         xphot(8,1)  = e_2s * sdir_new(1) + e_2p * pdir_new(1)
         xphot(9,1)  = e_2s * sdir_new(2) + e_2p * pdir_new(2)
         xphot(10,1) = e_2s * sdir_new(3) + e_2p * pdir_new(3)
         length = sqrt(xphot(8,1)**2 + xphot(9,1)**2 + xphot(10,1)**2)
         xphot(8,1)  = xphot(8,1)/length
         xphot(9,1)  = xphot(9,1)/length
         xphot(10,1) = xphot(10,1)/length
      endif
cccccc

      xphot(11,1)=tofg

      if (OUTPUT_SCATTABLE.eq.1) then
         call gskpho(IREFSTRY) ! track photon scattering and reflections
      else
         CALL gskpho(0)
      endif

      ngphot=0

      RETURN
      END

