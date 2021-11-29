************************************************************************
*     -----------------
      SUBROUTINE SGBST
*     -----------------
*
*     (Purpose)
*
*     (Input)
*       none
*
*     (Output)
*       none
*
*
*     (Creation Date and Author)
*       1992.12.02 ; First version by Y.Koshio
*       2001.02.09 ; RNDM -> GRNDM
*       2002.12.01 ; C. Sterner skdetsim upgrade
*       2003.06.07 ; increase refleaction by factor 2.7 by K.Ishihara
************************************************************************
      implicit none
#include "skhead.h"
#include "skdetsim_parm.h"
#include "skcnst.h"
#include "geant321/gctrak.inc"
#include "geant321/gcvolu.inc"
#include "geant321/gcking.inc"
#include "./guref.h"

      real fixstep
      PARAMETER(FIXSTEP = .1)

      REAL VECNEW(3),ABC(3),pdir(3)
      LOGICAL ROTATE
      character*4 cname
      real cost1,norm,normv,normabc
      real ftun1,ftun2,ftun3,ftun4
      data ftun1,ftun2,ftun3,ftun4/1.,-1.78974,1.89044,-0.686483/
      real wstn, z,zstn,rdummy,rr,aaa,ps,prrefl,rfbssg,costh,sinth
      real cosph,sinph,theta,xcos,phi
      integer ipol,i

      integer scatter
      real cost2,cost3

      real ref_angle,sint,cost,sinp,cosp
      real norm_vecnew,plane(3),normplane,pp

      logical r2

C for polarization
      real pol, length

      ISTOP=1
c      write(6,*)'bla---ck in /gekin/istak',gekin,istak


* reflection with holizontal bar of wall
      IF(SK_GEOMETRY.eq.SK_I) then 
              !SK-I had horizontal metal stripes on the barrel
         ngphot=0
         wstn=4.8/2.
         do 10 z=-26.,25.,1.
            zstn=35.35+70.7*z
            if(vect(3).gt.zstn-wstn.and.vect(3).lt.zstn+wstn)then
               CALL GRNDM(RDUMMY,1)
               if(RDUMMY.lt.0.5)then
                  RR = SQRT(VECT(1)**2+VECT(2)**2)
                  ABC(1) = -VECT(1)/RR
                  ABC(2) = -VECT(2)/RR
                  ABC(3) = 0.
                  NGPHOT = 1
                  COST1 = -(ABC(1)*VECT(4) + ABC(2)*VECT(5) 
     &                 + ABC(3)*VECT(6))
                  VECNEW(1) = VECT(4) + 2.*ABC(1)*COST1
                  VECNEW(2) = VECT(5) + 2.*ABC(2)*COST1
                  VECNEW(3) = VECT(6) + 2.*ABC(3)*COST1
                  norm = 
     &                 sqrt(vecnew(1)**2+vecnew(2)**2+vecnew(3)**2)
                  vecnew(1) = vecnew(1)/norm
                  vecnew(2) = vecnew(2)/norm
                  vecnew(3) = vecnew(3)/norm
               endif
               goto 11
            endif
 10      continue
      ENDIF
         
*  Calculate the incident angle
      CALL UHTOC(NAMES(3),4,CNAME,4)
      IF(CNAME.EQ.'WMED') THEN
         RR = SQRT(VECT(1)**2+VECT(2)**2)
         ABC(1) = -VECT(1)/RR
         ABC(2) = -VECT(2)/RR
         ABC(3) = 0.
  
      ELSE IF(CNAME.EQ.'XMED') THEN
         ABC(1) = 0.
         ABC(2) = 0.
         ABC(3) = -SIGN(1.,VECT(3))

      ELSE
         WRITE(6,*) 'Warning in SGBST'
         RETURN
      ENDIF
** normalization      
      normv = sqrt(vect(4)**2+vect(5)**2+vect(6)**2)
      normabc = sqrt(abc(1)**2+abc(2)**2+abc(3)**2)
      COST1 = -(ABC(1)*VECT(4) + ABC(2)*VECT(5) + ABC(3)*VECT(6))

      if(normv.ne.0.and.normabc.ne.0.)then
        cost1 = cost1/(normv*normabc)
      else
        write(*,*)'ERROR in sgbst normv or normabc = 0'
        write(*,*)'vect',normv,vect(4),vect(5),vect(6),
     &            'abc',normabc,abc
	  stop
      end if
 
* calculate the ratio of pwave
	aaa = 1./cost1**2 - 1.
      if(cost1.eq.0..or.cost1.eq.1..or.aaa.le.0)then
         ipol=0
         return
      else
         IF (abs(vect(4)*polar(1)+vect(5)*polar(2)+vect(6)*polar(3))
     &        .gt.0.0001) then

C            write(*,*) 'VECT=' ,vect(4),vect(5),vect(6)
C            write(*,*) 'PLOAR=',polar(1),polar(2),polar(3)
            
            ipol=0
         else
            pdir(1)=(vect(4)+abc(1)/cost1)/sqrt(1./cost1**2-1.)
            pdir(2)=(vect(5)+abc(2)/cost1)/sqrt(1./cost1**2-1.)
            pdir(3)=(vect(6)+abc(3)/cost1)/sqrt(1./cost1**2-1.)
            
C             write(*,*) 'VECT=',vect(4),vect(5),vect(6) 
C             write(*,*) 'ABC=' ,abc(1),abc(2),abc(3)
C             write(*,*) 'pdir=',pdir(1),pdir(2),pdir(3)
C             write(*,*) 'PLOAR=',polar(1),polar(2),polar(3)

c     
c            if(pdir(1)**2+pdir(2)**2+pdir(3)**2.lt.0.999)then
c               write(*,*)'error in sgbst pdir=',pdir
c            endif
c     
            ps=(polar(1)*pdir(1)+polar(2)*pdir(2)+polar(3)*pdir(3))**2
            
C             write(*,*) 'ps=',ps

c            if(ps.gt.1.01.or.ps.lt.0.)then
c               write(*,*)'error in sgbst ps=',ps
c            endif
C            if(rndm().lt.ps)then
            CALL GRNDM(RDUMMY,1)
            if(RDUMMY.lt.ps)then
               if(SK_GEOMETRY .le. SK_IV) then
                  ipol=2
               else
                  ipol=1
               endif
            else
               if(SK_GEOMETRY .le. SK_IV) then
                  ipol=1
               else
                  ipol=2
               endif
            endif
         endif
      endif

*  Calculate the ratio of reflection on B.S.

      PRREFL= 1.0 * RFBSSG(COST1,ipol)*BSTREF
      
      NGPHOT = 0
      
C         write(*,*) 'COST1=',cost1,'  ipol=',ipol

*  Reflection !
*  Sei hansha on wavy plane

      IF (REFLECTION_TUNE .EQ. AP_TUNE) THEN
         CALL GRNDM(RDUMMY,1)
         IF(RDUMMY.LT.PRREFL) THEN
            NGPHOT = 1
            VECNEW(1) = VECT(4) + 2.*ABC(1)*COST1
            VECNEW(2) = VECT(5) + 2.*ABC(2)*COST1
            VECNEW(3) = VECT(6) + 2.*ABC(3)*COST1
            CALL GFANG(VECNEW,COSTH,SINTH,COSPH,SINPH,ROTATE)

 7          CALL GRNDM(RDUMMY,1)
            THETA = ACOS(1.-RDUMMY**4)
C            THETA = (pi/2.)*RDUMMY  ! final tunig1 -- random reflection 

            xcos=cos(pi/2.-theta)
c            if(rndm().gt.
c    &       ftun1+ftun2*xcos+ftun3*xcos**2+ftun4*xcos**3)goto 7
C            PHI = (PI*2.)*rndm()
            CALL GRNDM(RDUMMY,1)
            PHI = (PI*2.)*RDUMMY
            VECNEW(1)=SIN(THETA)*COS(PHI)
            VECNEW(2)=SIN(THETA)*SIN(PHI)
            VECNEW(3)=COS(THETA)
            IF(ROTATE) CALL GDROT(VECNEW,COSTH,SINTH,COSPH,SINPH)
         ENDIF

      ELSE                      ! LE, UNI_TUNE

         CALL GRNDM(RDUMMY,1)
         IF(RDUMMY.LT.PRREFL) THEN
            NGPHOT = 1

            IF(SK_GEOMETRY.eq.SK_I) then 
               VECNEW(1) = VECT(4) + 2.*ABC(1)*COST1
               VECNEW(2) = VECT(5) + 2.*ABC(2)*COST1
               VECNEW(3) = VECT(6) + 2.*ABC(3)*COST1
               CALL GFANG(VECNEW,COSTH,SINTH,COSPH,SINPH,ROTATE)

               CALL GRNDM(RDUMMY,1)
               THETA = ACOS(1.-RDUMMY**4)
               CALL GRNDM(RDUMMY,1)
               PHI = (PI*2.)*RDUMMY
               VECNEW(1)=SIN(THETA)*COS(PHI)
               VECNEW(2)=SIN(THETA)*SIN(PHI)
               VECNEW(3)=COS(THETA)
               IF(ROTATE) CALL GDROT(VECNEW,COSTH,SINTH,COSPH,SINPH)

            ELSE

               if (COST1.lt.0.) then
                  NGPHOT=0
                  GOTO 888
CC                  NGPHOT=0
               endif

 777           call sgbsdist(COST1,PRREFL,scatter,cost2)
               sint = sqrt(1.-cost2**2)
               cost = cost2
               CALL GRNDM(RDUMMY,1)
               phi = RDUMMY*pi*2.
               sinp = sin(phi)
               cosp = cos(phi)


               if (scatter.eq.1) then ! Lambert
                  VECNEW(1) = sint*cosp
                  VECNEW(2) = sint*sinp
                  VECNEW(3) = cost
                  CALL GFANG(ABC,COSTH,SINTH,COSPH,SINPH,ROTATE)
                  IF (ROTATE) CALL GDROT(VECNEW,COSTH,SINTH,COSPH,SINPH)
               endif
               if (scatter.eq.2) then ! Specular 
                  VECNEW(1) = VECT(4) + 2.*ABC(1)*COST1
                  VECNEW(2) = VECT(5) + 2.*ABC(2)*COST1
                  VECNEW(3) = VECT(6) + 2.*ABC(3)*COST1
                  CALL GFANG(VECNEW,COSTH,SINTH,COSPH,SINPH,ROTATE)
                  VECNEW(1) = sint*cosp
                  VECNEW(2) = sint*sinp
                  VECNEW(3) = cost
                  IF (ROTATE) CALL GDROT(VECNEW,COSTH,SINTH,COSPH,SINPH)
                  IF ( (ABC(1)*VECNEW(1)+ABC(2)*VECNEW(2)+
     &                 ABC(3)*VECNEW(3)) .le. 0.0) then 
                     goto 777
                  endif
               endif
               if (scatter.eq.3) then ! Phong
                  VECNEW(1) = VECT(4) + 2.*ABC(1)*COST1
                  VECNEW(2) = VECT(5) + 2.*ABC(2)*COST1
                  VECNEW(3) = VECT(6) + 2.*ABC(3)*COST1
               endif
 888           CONTINUE
            ENDIF
         ENDIF
      ENDIF


 11   continue     

      IF(NGPHOT.GE.1) THEN

         STEP = FIXSTEP
         DO I = 1,3
            XPHOT(I+3,1) = VECNEW(I)
            XPHOT(I,1) = VECT(I) + STEP*VECNEW(I)
         ENDDO
         XPHOT(7,1)  = VECT(7)
c         XPHOT(11,1) = TOFG + STEP/CLIGHT*refwtr
         XPHOT(11,1) = TOFG

**** porlarization is not until SK-IV
C
C polarization (refer from random reflection formular in Geant4)
C
         if(SK_GEOMETRY .ge. SK_V) then
            pol = polar(1)*abc(1) + polar(2)*abc(2) + polar(3)*abc(3)
            xphot(8,1)  = -polar(1) + 2.*pol*abc(1)
            xphot(9,1)  = -polar(2) + 2.*pol*abc(2)
            xphot(10,1) = -polar(3) + 2.*pol*abc(3)
            length = sqrt(xphot(8,1)**2 + xphot(9,1)**2 + xphot(10,1)**2)
            xphot(8,1)  = xphot(8,1)/length
            xphot(9,1)  = xphot(9,1)/length
            xphot(10,1) = xphot(10,1)/length
         endif
cccccc

         if (OUTPUT_SCATTABLE.eq.1) then
            CALL GSKPHO(-ISTORY-20000000) ! track photon scattering and reflections
         else
            CALL GSKPHO(0)
         endif

         NGPHOT=0

      ENDIF
      
c      write(6,*)'bla---ck out /gekin/istak',gekin,istak


      RETURN
      END
