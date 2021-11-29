************************************************************************
*     -------------------------------------------------------
      SUBROUTINE SGREFPMT(NOPMT,OLDVEC,QE,VECNEW,COST1,IFLAG)
*     -------------------------------------------------------
*
*     (Purpose)
*       Calculate catch or reflection
*
*     (Input)
*       NOPMT  : Cable number
*       OLDVEC : Position(1-3) and Direction(4-6) of hitting PMT
*       QE     : SK-I/II, quantum efficiency (only for original ledetsim)
*                SK-III, TRACK_ALL (absorbed photon) ->
*                        p.e. conversion ratio in photo-cathode
*                        TRACK_QE_OVER_ABS_PHOTON  qetable/max(qetable)
*       VECNEW : given reflection direction from skdonuts
*                                   (for original ledetsim and UNI_MODE)
*       COST1  : Incident angle to PMT (only for original apdetsim)
*
*     (Output)
*       IFLAG  : SK-I/II/III ATMPD REFLECTION TUNE
*                 Flag (0: catch, 1: Sei hansha, 2: Ran hansha)
*                SK-I/II LE REFLECTION TUNE
*                 Flag (0: disappear, 1: Sei hansha, 2: Ran hansha)
*                >=SK-III LE_,UNI_TUNE
*                 Flag (0: strange, 1: Sei hansha, 2: Ran hansha,
*                       3: p.e. production, 4: transmission,
*                       5: absorbed w/o p.e.)
*
*     (Creation Date and Author)
*       1996. 2.28 ; First version by Y.Koshio
*       1998. 3. 3 ; Modified by K.Okumura
*       2003. 2.20 ; Merged version by Y.Koshio
*       2008. 4.23 ; Algorithm of reducing the number of photons is
*                    improved by Byeongsu Yang
*
************************************************************************

      implicit none

#include "skdetsim_parm.h"
#include "skhead.h"
#include "guref.h"
#include "geant321/gctrak.inc"

      integer nopmt, iflag
      real qe, OLDVEC(6), VECNEW(3),COST1

      integer ipol
      real PRREFL, PRABSL,PRABSL0, PROB, RANREF
*      real PRREFL, PRABSL, PROB, RANREF
      real blamda

      real norm,normoldvec(3),normvecnew(3),product,sdir(3),asdir,ps
      real pdir(3), tan2t1, abc(3)

      If (REFLECTION_TUNE.ne.AP_TUNE) then

*  Calculate the incident angle
         norm = sqrt(oldvec(4)**2+oldvec(5)**2+oldvec(6)**2)
         if(norm.gt.0)then
            normoldvec(1) = oldvec(4)/norm
            normoldvec(2) = oldvec(5)/norm
            normoldvec(3) = oldvec(6)/norm
         else
            write(*,*)'ERROR in sgrefpmt oldvec = 0'
            stop
         end if
         norm = sqrt(vecnew(1)**2+vecnew(2)**2+vecnew(3)**2)
         if(norm.gt.0)then
            normvecnew(1) = vecnew(1)/norm
            normvecnew(2) = vecnew(2)/norm
            normvecnew(3) = vecnew(3)/norm
         else
            write(*,*)'ERROR in sgrefpmt vecnew = 0'
            stop
         end if

         product = normoldvec(1)*normvecnew(1)
     &        +normoldvec(2)*normvecnew(2)
     &        +normoldvec(3)*normvecnew(3)
 
         if(abs(product).gt.1)then
            write(*,*)'ERROR in sgrefpmt |product|>1',product
            if(product.ge.1)then
               product = 1. - 0.0001
            else if(product.le.-1.)then
               product = -1. + 0.0001
            end if
         end if
         cost1=cos(acos((-1)*product)/2.)

* calculate the ratio of pwave
         if(SK_GEOMETRY .le. SK_IV) then
            if(cost1.eq.0..or.cost1.eq.1.)then
               ipol=0
            else   
               IF (abs(vect(4)*polar(1)+vect(5)*polar(2)+vect(6)*polar(3))
     &              .gt.0.0001) then
                  ipol=0
               else
                  sdir(1)=oldvec(4)+vecnew(1)
                  sdir(2)=oldvec(5)+vecnew(2)
                  sdir(3)=oldvec(6)+vecnew(3)
                  asdir=sqrt(sdir(1)**2+sdir(2)**2+sdir(3)**2)
                  sdir(1)=sdir(1)/asdir
                  sdir(2)=sdir(2)/asdir
                  sdir(3)=sdir(3)/asdir
                  ps=(polar(1)*sdir(1)+polar(2)*sdir(2)+polar(3)*sdir(3))**2
                  CALL GRNDM(PROB,1)
                  if(PROB.lt.ps)then
                     ipol=2
                  else
                     ipol=1
                  endif
               endif
            endif

         else                   ! bug(?) fix after SK-V
            tan2t1 = 1./cost1**2 - 1.
            if(cost1.eq.0..or.cost1.eq.1..or.tan2t1.le.0)then
               ipol=0
            else
               IF (abs(vect(4)*polar(1)+vect(5)*polar(2)+vect(6)*polar(3))
     &              .gt.0.0001) then
                  ipol=0
               else
                  abc(1) = -vect(4)+vecnew(1)
                  abc(2) = -vect(5)+vecnew(2)
                  abc(3) = -vect(6)+vecnew(3)
                  norm = sqrt(abc(1)**2+abc(2)**2+abc(3)**2)
                  abc(1) = abc(1)/norm
                  abc(2) = abc(2)/norm
                  abc(3) = abc(3)/norm

                  pdir(1)=(vect(4)+abc(1)/cost1)/sqrt(1./cost1**2-1.)
                  pdir(2)=(vect(5)+abc(2)/cost1)/sqrt(1./cost1**2-1.)
                  pdir(3)=(vect(6)+abc(3)/cost1)/sqrt(1./cost1**2-1.)
                  norm = sqrt(pdir(1)**2+pdir(2)**2+pdir(3)**2)
                  pdir(1) = pdir(1)/norm
                  pdir(2) = pdir(2)/norm
                  pdir(3) = pdir(3)/norm
            
                  ps=(polar(1)*pdir(1)+polar(2)*pdir(2)+polar(3)*pdir(3))**2
                  CALL GRNDM(PROB,1)
                  if(PROB.lt.ps)then
                     ipol=1
                  else
                     ipol=2
                  endif
               endif
            endif
         endif

      Else                      ! AP_TUNE
         ipol = 0
      Endif

*  Check if the photon is reflected or not
*  rfpmsg was changed to subroutines to get prrefl, prabsl, prabsl0
****      PRREFL= RFPMSG(COST1,ipol)

C-- Algorithm of reducing the number of photons is improved
C   by Byeongsu Yang 2008/04/23
C      call RFPMSG(COST1,ipol,PRREFL,PRABSL,PRABSL0)
      call RFPMSG(COST1,ipol,PRREFL,PRABSL)
      IF (SK_GEOMETRY.ge.SK_III .and. REFLECTION_TUNE.ne.AP_TUNE) THEN
        CALL GRNDM(PROB,1)

        IFLAG = 0
*  Check if the photon is absorbed without p.e. production
CCC        IF(PROB.LT.(PRABSL-QE*PRABSL/PRABSL0)) THEN
        IF(PROB.LT.PRABSL*(1.0-QE)) THEN
           IFLAG = 5
*  Check if the p.e. is produced or not
        ELSE IF(PROB.LT.PRABSL) THEN
           IFLAG = 3
*  Check if the photon do specular reflection
        ELSE IF(PROB.LT.(PRABSL+PRREFL*PMSREF)) THEN
           IFLAG = 1
*  Check if the photon is transmitted
        ELSE IF(PROB.GE.(PRABSL+PRREFL)) THEN
           IFLAG = 4
*  Check if the photon do diffused reflection
CCC        ELSE IF(PROB.LT.(PRABSL+PRREFL*(PMSREF+PMRREF))) THEN
        ELSE
           IFLAG = 2
C-- Remove unnecessary part. Byeongsu Yang 2008/04/23
C*  when PRABSL+PRREFL <=PROB< PRABSL+PRREFL, but this must be strange case,
C        ELSE
C           IFLAG = 0
        ENDIF

******* above checking order is a little bit faster
**  Check if the photon do specular reflection
*        IF(PROB.lt.(PRREFL*PMSREF)) THEN
*           IFLAG = 1
**  Check if the photon do diffused reflection
*        ELSE IF(PROB.lt.(PRREFL*(PMSREF+PMRREF))) THEN
*           IFLAG = 2
** should be strange case, assume a photon is absorbed without p.e. production
*        ELSE IF(PROB.lt.PRREFL) THEN
*           IFLAG = 0
**  Check if the p.e. is produced or not
*        ELSE IF(PROB.lt.(PRREFL+QE*PRABSL/PRABSL0)) THEN
*           IFLAG = 3
**  Check if the photon is absorbed without p.e. production
*        ELSE IF(PROB.lt.(PRREFL+PRABSL)) THEN
*           IFLAG = 5
**  Check if the photon is transmitted
*        ELSE
*           IFLAG = 4
*        ENDIF

*        WRITE(*,*)'Ref:',COST1,ACOS(COST1)*180./3.141593,ipol
*        WRITE(*,*)'PROB:',PRREFL,PRABSL,PRABSL0,QE*PRABSL/PRABSL0
        RETURN
      ENDIF	! SK_GEOMETRY condition

      If (REFLECTION_TUNE.eq.LE_TUNE) prrefl=prrefl/(1.-qe)
      CALL GRNDM(PROB,1)

      IF(PROB.GE.PRREFL) THEN
         IFLAG = 0
         RETURN
      ELSE
         CALL GRNDM(RANREF,1)
         IF(RANREF.LE.PMSREF) THEN
            IFLAG = 1
C-- Remove unnecessary part. Byeongsu Yang 2008/04/23
CC         ELSE IF(RANREF.LE.PMSREF+PMRREF) THEN
         ELSE
            IFLAG = 2
CC         ELSE
CC            IFLAG = 0
CC            RETURN
         ENDIF
      ENDIF

      RETURN
      END
