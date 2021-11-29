
*CMZ :  3.14/16 05/09/90  11.40.55  by  Rene Brun
*-- Author :
      SUBROUTINE GTCKOV
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Tracking routine for Cherenkov type tracks               *
C.    *                                                                *
C.    *    ==>Called by : GTRACK                                       *
C.    *       Author    R.Brun, F.Bruyant                              *
C.    *                                                                *
C.    *       Modified for SuperKamiokande : Geantino for              *
C.    *                                          Cherenkov light       *
C.    *        11/27/92 by Y.Koshio/Y.Suzuki                           *
C.    *        8/30/02 modified for acrylic in SK-II by C.Sterner      *
C.    *                                                                *
C.    *        2006.05.10 Bug correction for OD segmentation           *
C.    *                              Fanny Dufour                      *
C.    *                                                                *
C.    *        2009.1.06 impliment the asym  by H.Sekiya / Y.Koshio    *
C     *        2009.11.18 Make top bottom asymmetrical MC              * 
C     *                                        (beta version)          *
C     *                    Fix bug-OD's water parameter become         *
C     *                                      uniformly (T.Yokozawa)    *
C     *        2020.02.04 Bug correction for TBA calculation           *
C     *                    Y.Nakajima, Y.Koshio                        *
C.    *                                                                *
C     *        2020.05 Add photon survival probability penalty in      *
C     *                top & bottom OD endcap to account for steel     *
C     *                structure (T. Wester)                           *
C.    *                                                                *
C.    ******************************************************************
      IMPLICIT NONE
#include "geant321/gctmed.inc"
#include "geant321/gctrak.inc"
#include "geant321/gconst.inc"
#include "geant321/gckine.inc"
#include "gcbank.h"
#include "epsmac.h"
#include "rayfac.h" 
#include "guwtr.h"
#include "gumed.h"
#include "skhead.h"
#include "skdetsim_parm.h"

#include "dsparm.h"
#include "dsparma.h"

C local defines for implicit none
      integer iproc,i
      real blamda,rdummy,zintuw,stepuw,wtrsg,suwtr,vectmp,effnsg
      real sluwtr,zintac,stepac,acrsg,suacr,acrnsg
      real lesym, leasym, leabs ,rnext1,rnext2, ynext1, ynext2
      REAL WSYMSCSG,WASYSCSG,WABSSG
      integer icase, alphacount
      real trans0, alpha, zbound, tracklsk, trans,
     &     rend, rbound, transend, transbound, yend, ybound
      real factora, factorb, factorc,b24ac, lewatert
      real probss
      
      SAVE alphacount

      integer bound_flag

C.    ------------------------------------------------------------------


*-------------------WATER MEDIA TRANSPORT
*
**  Check if Media is Water
*
      IF(NUMED.EQ.1.OR.NUMED.EQ.2.OR.NUMED.EQ.8.or.NUMED.eq.9) THEN
*
* *** Compute current step size
*     
         STEP   = BIG
         IPROC  = 0
         BLAMDA = 1237.8E-09 / GEKIN
*
* *** Create interaction probabilities
*
         CALL GRNDM(RDUMMY,1)
         ZINTUW = LOG(RDUMMY)
* 
* **  Step limitation due to interaction with water
*
         STEPUW = WTRSG(BLAMDA)
         trans0 = 100./STEPUW
         lesym = WSYMSCSG(BLAMDA)/100.
         leasym = WASYSCSG(BLAMDA)/100.
         leabs = WABSSG(BLAMDA)/100.
c         lewatert = 1./leabs
c         write(21,*) vect(3),lewatert

c **  Set asymmetric factor
         alpha = TBADS
         zbound = -9999.

c        New flag to check if the photon reached the Z=-11m or detector boundary
         bound_flag = 0

         if(SK_GEOMETRY.eq.SK_III .or. SK_GEOMETRY.eq.SK_IV) zbound = -1100.
         if(alphacount.eq.0)then
            if(SK_GEOMETRY.ge.SK_III) then
               print *,'*********************************************'
               print *,'top-bottom asymmetry parameter by Yokozawa'
               if(SK_GEOMETRY.eq.SK_III .or. SK_GEOMETRY.eq.SK_IV)
     &              print *,'alpha,zbound = ',alpha,zbound
               if(SK_GEOMETRY.eq.SK_V)
     &              print *,'alpha = ',alpha,' No z boundary in SK-V'
               print *,'*********************************************'
               print *,'yokozawa debug: leabs'
            endif
            alphacount = 1
         end if

         IF(STEPUW.LE.EPSMAC)THEN
            SUWTR = BIG
         ELSE
            if(alpha.eq.0.)then
               SUWTR = -trans0*zintuw
               goto 1009
            end if
c
c  For SK-III and SK-IV (zboundary is defined)
c
c            if(SK_GEOMETRY.ge.SK_III .and. WATER_TUNE.eq.UNI_TUNE) then ! t-b asym
c            if( SK_GEOMETRY.ge.SK_III .and. ( (WATER_TUNE.eq.UNI_TUNE)
            if( (SK_GEOMETRY.eq.SK_III .or. SK_GEOMETRY.eq.SK_IV)
     &           .and. ( (WATER_TUNE.eq.UNI_TUNE)
     &           .or. (WATER_TUNE.eq.10) .or. (WATER_TUNE.EQ.11) )
     &           .and. ((NUMED .ne. 9) .and. (NUMED .ne. 1)) ) then ! only ID
c*** shift mean abs value to data by Yokozawa Jun.24th,2011
               if((1.0-69.63*alpha).ne.0 .and. IWATERP.eq.200904) then ! keep bug for ATMPD
                  leabs = WABSSG(BLAMDA)/(1.0-69.63*alpha)/100.
               else if((1.0+69.63*alpha).ne.0 .and. IWATERP.ne.200904) then !bug fix for LOWE
c     Fixed the sign in from of 69.63! Y. Nakajima, Mar. 17, 2020
                  leabs = WABSSG(BLAMDA)/(1.0+69.63*alpha)/100.
               else
                  leabs = WABSSG(BLAMDA)/100.
               end if
               rend = TRACKLSK(vect(1),vect(4))

               if(vect(3).ge.zbound) then
                  if(vect(6).ge.0.) then
                     icase = 1
                  else
                     rbound = (zbound-vect(3))/vect(6)
                     if(rend.le.rbound) then
                        icase = 1
                     else
                        icase = 2
                     endif
                  endif
               else
                  if(vect(6).ge.0.) then
                     rbound = (zbound-vect(3))/vect(6)
                     if(rend.le.rbound) then
                        icase = 3
                     else
                        icase = 4
                     endif
                  else
                     icase = 3
                  endif
               endif
               
c  top bottom asymmetry is only for absorption
               if(IWATERP.eq.200904) then
                  factora=leabs*alpha*vect(6)
               else
c Bug fix on Feb. 4, 2020! for LOWE
                  factora=0.5*leabs*alpha*vect(6)
               endif
               factorb=leabs+lesym+leasym+leabs*alpha*vect(3)
               factorc=zintuw
               b24ac = factorb*factorb-4*factora*factorc
               
               if(icase.eq.1) then
                  trans = 1/(leabs*(1+alpha*vect(3))+lesym+leasym)
                  transend = 1/(leabs*(1+alpha*(vect(3)+rend*vect(6)))+lesym+leasym)
                  if(trans.gt.0. .and. transend.gt.0.) then
                     if(IWATERP.eq.200904) then
                        yend = exp(-rend*(leabs*(1+alpha*(vect(3)+rend*vect(6)))+lesym+leasym))
                     else
                        yend = exp(-rend*(leabs*(1+alpha*(vect(3)+0.5*rend*vect(6)))+lesym+leasym))
                     endif
                  else
                     if(IWATERP.eq.200904) then
                        transend = trans0
                        yend = 9999.
                     else
                        write(6,*) 'Strange transparecy value in icase 1 !!! ',trans, transend
                        write(6,*) 'Call Nakajima or Koshio'
                        stop
                     endif
                  endif
                  if(RDUMMY.ge.yend) then
                     if(b24ac.lt.0)then
                        print *,'error!!!b**2-4ac<0'
                        SUWTR = -transend*zintuw
                        goto 1009
                     end if
                     if(factora.eq.0)then
                        SUWTR = (-factorc/factorb)
                     else
                        SUWTR = (-factorb+sqrt(b24ac))/(2*factora)
                     end if
                  else
                     if(IWATERP.eq.200904) then
                        SUWTR = -transend*zintuw
                     else
                        bound_flag = 1
                        SUWTR = rend + PREC
                     endif
                  endif

               else if(icase.eq.2) then

                  if(IWATERP.eq.200904) then ! keep bug for ATMPD
                     if(trans.gt.0. .and. transend.gt.0.) then
                        ybound = exp(-rbound*(leabs*(1+alpha*(vect(3)+rbound*vect(6)))+lesym+leasym))
                        yend = exp(-rend*(leabs*(1+alpha*(vect(3)+rend*vect(6)))+lesym+leasym))
                     else
                        transend = trans0
                        yend = 9999.
                        ybound = 9999.
                     endif
                     transbound = 1/(leabs*(1+alpha*zbound)+lesym+leasym)
                  else ! bug fix for LOWE
                     trans = 1/(leabs*(1+alpha*vect(3))+lesym+leasym)
                     transbound = 1/(leabs*(1+alpha*zbound)+lesym+leasym)
                     if(trans.gt.0. .and. transbound.gt.0.) then
                        ybound = exp(-rbound*(leabs*(1+alpha*(vect(3)+0.5*rbound*vect(6)))+lesym+leasym))
                     else
                        write(6,*) 'Strange transparecy value in icase 2 !!! ',trans, transbound
                        write(6,*) 'Call Nakajima or Koshio'
                        stop
                     endif
                  endif

                  if(RDUMMY.ge.ybound) then
                     if(b24ac.lt.0)then
                        print *,'error!!!b**2-4ac<0'
                        if(IWATERP.eq.200904) then ! keep bug for ATMPD
                           SUWTR = -transend*zintuw
                        else ! bug fix for LOWE
c                     May be better to use trandbound here?
                           SUWTR = -transbound*zintuw
                        endif
                        goto 1009
                     end if
                     if(factora.eq.0)then
                        SUWTR = (-factorc/factorb)
                     else
                        SUWTR = (-factorb+sqrt(b24ac))/(2*factora)
                     end if
                  else
                     if(IWATERP.eq.200904) then ! keep bug for ATMPD
                        SUWTR = -transbound*zintuw
                     else ! bug fix for LOWE
                        bound_flag = 1
                        SUWTR = rbound + PREC
                     endif
                  endif

               else if(icase.eq.3) then
                  transbound = 1/(leabs*(1+alpha*zbound)+lesym+leasym)
                  SUWTR = -transbound*zintuw

               else if(icase.eq.4) then
                  if(IWATERP.eq.200904) then ! keep bug for ATMPD
                     trans = 1/(leabs*(1+alpha*vect(3))+lesym+leasym)
                     transend = 1/(leabs*(1+alpha*(vect(3)+rend*vect(6)))+lesym+leasym)
                  endif
                  transbound = 1/(leabs*(1+alpha*zbound)+lesym+leasym)

                  if(IWATERP.eq.200904) then ! keep bug for ATMPD
                     if(trans.gt.0. .and. transbound.gt.0. .and. transend.gt.0.) then
                        ybound = exp(-rbound*(leabs*(1+alpha*(vect(3)+rbound*vect(6)))+lesym+leasym))
                        yend = exp(-rend*(leabs*(1+alpha*(vect(3)+rend*vect(6)))+lesym+leasym))
                     else
                        transend = trans0
                        yend = 9999.
                        ybound = 9999.
                     endif
                  else ! bug fix for LOWE
                     if(transbound.gt.0.) then
                        ybound = exp(-rbound/transbound)
                     else
                        write(6,*) 'Strange transparecy value in icase 4 !!! ', transbound
                        write(6,*) 'Call Nakajima or Koshio'
                        stop
                     endif
                  endif
                  if(RDUMMY.ge.ybound) then
                     SUWTR = -transbound*zintuw
                  else if(RDUMMY.ge.yend .and. IWATERP.eq.200904) then
                     if(b24ac.lt.0)then
                        print *,'error!!!b**2-4ac<0'
                        SUWTR = -transend*zintuw
                        goto 1009
                     end if
                     if(factora.eq.0)then
                        SUWTR = (-factorc/factorb)
                     else
                        SUWTR = (-factorb+sqrt(b24ac))/(2*factora)
                     end if
                  else
                     if(IWATERP.eq.200904) then ! keep bug for ATMPD
                        SUWTR = -transbound*zintuw
                     else ! bug fix for LOWE
                        bound_flag = 1
                        SUWTR = rbound + PREC
                     endif
                  endif
               endif
               if(SUWTR.lt.0.) SUWTR = -trans0*zintuw


*** first trial : top bottom asymmetry is for all the water parameters
c
c               if(icase.eq.1) then
c                  trans = trans0+alpha*vect(3)
c                  transend = trans0 + alpha*rend*vect(6)
c                  if(trans.gt.0. .and. transend.gt.0.) then
c                     yend = exp(- rend/(trans0+alpha*(vect(3)+rend*vect(6))))
c                  else
c                     transend = trans0
c                     yend = 9999.
c                  endif
c                  if(RDUMMY.ge.yend) then
c                     SUWTR = -(trans0*zintuw+alpha*vect(3)*zintuw)/
c     &                    (1+alpha*vect(6)*zintuw)
c                  else
c                     SUWTR = -transend*zintuw
c                  endif
c               else if(icase.eq.2) then
c                  ybound = exp(- rbound/(trans0+alpha*(vect(3)+rbound*vect(6))))
c                  yend = exp(- rend/(trans0+alpha*(vect(3)+rend*vect(6))))
c                  transbound = trans0 + alpha*zbound
c                  if(RDUMMY.ge.ybound) then
c                     SUWTR = -(trans0*zintuw+alpha*vect(3)*zintuw)/
c     &                    (1+alpha*vect(6)*zintuw)
c                  else
c                     SUWTR = -transbound*zintuw
c                  endif
c               else if(icase.eq.3) then
c                  transbound = trans0 + alpha*zbound
c                  SUWTR = -transbound*zintuw
c               else if(icase.eq.4) then
c                  trans = trans0+alpha*vect(3)
c                  transbound = trans0 + alpha*zbound
c                  transend = trans0 + alpha*rend*vect(6)
c                  if(trans.gt.0. .and. transbound.gt.0. .and. transend.gt.0.) then
c                     ybound =
c     &                    exp(- rbound/(trans0+alpha*(vect(3)+rbound*vect(6))))
c                     yend = exp(- rend/(trans0+alpha*(vect(3)+rend*vect(6))))
c                  else
c                     transend = trans0
c                     ybound = 9999.
c                     yend = 9998.
c                  endif
c               if(ybound.lt.yend) then
c                  print*,'hen2 in gtckov',ybound,yend,rbound,rend
c                  stop
c               endif
c                  if(RDUMMY.ge.ybound) then
c                     SUWTR = -transbound*zintuw
c                  else if(RDUMMY.ge.yend) then
c                     SUWTR = -(trans0*zintuw+alpha*vect(3)*zintuw)/
c     &                    (1+alpha*vect(6)*zintuw)
c                  else
c                     SUWTR = -transend*zintuw
c                  endif
c               endif
c               if(SUWTR.lt.0.) SUWTR = -trans0*zintuw
*** 

c  For SK-V (no zboundary)
            else if ( SK_GEOMETRY.eq.SK_V .and. WATER_TUNE.eq.UNI_TUNE
     &              .and. (NUMED.ne.9 .and. NUMED.ne.1) ) then ! only ID

               rend = TRACKLSK(vect(1),vect(4))
               
c  top bottom asymmetry is only for absorption
               factora=0.5*leabs*alpha*vect(6)
               factorb=leabs+lesym+leasym+leabs*alpha*vect(3)
               factorc=zintuw
               b24ac = factorb*factorb-4*factora*factorc

               trans = 1/(leabs*(1+alpha*vect(3))+lesym+leasym)
               transend = 1/(leabs*(1+alpha*(vect(3)+rend*vect(6)))+lesym+leasym)
               if(trans.gt.0. .and. transend.gt.0.) then
                  yend = exp(-rend*(leabs*(1+alpha*(vect(3)+0.5*rend*vect(6)))+lesym+leasym))
               else
                  print *,'error!!! trans or transend less than zero in gtckov'
                  transend = trans0
                  yend = 9999.
               endif
               if(RDUMMY.ge.yend) then
                  if(b24ac.lt.0)then
                     print *,'error!!!b**2-4ac<0 in gtckov'
                     SUWTR = -transend*zintuw
                     goto 1009
                  end if
                  if(factora.eq.0)then
                     SUWTR = (-factorc/factorb)
                  else
                     SUWTR = (-factorb+sqrt(b24ac))/(2*factora)
                  end if
               else
                  bound_flag = 1
                  SUWTR = rend
               endif
               if(SUWTR.lt.0.) SUWTR = -trans0*zintuw

            else
               SUWTR  = -ZINTUW * 100./STEPUW
            endif

         ENDIF
 1009    continue
         IF(SUWTR.LT.STEP) THEN
            STEP = SUWTR
            if (bound_flag .eq. 0) then
               IPROC= 1
            endif
         ENDIF

* **  Step limitaion due to geometry?
*
         CALL GTNEXT
         IF (STEP.GT.SNEXT) THEN
            STEP   = SNEXT + PREC
            INWVOL= 2
            IPROC = 0
            NMEC =1
            LMEC(1)=1
         ENDIF

c        T. Wester add travel distance penalty in top od region to 
c        mimic effect of steel structure
        
         if (SK_GEOMETRY.ge.SK_V) then
             if (numed .eq. imdwoe) then
                 probss = 1.0
                 if (vect(3) .gt. 0) then
                     probss = 1.0 - PROBODSST * step / 100.
                 else
                     probss = 1.0 - PROBODSSB * step / 100.
                 endif
                 if (probss .lt. 0.05) then
                     probss = 0.05
                 endif

                 call grndm(rdummy)
                 if (rdummy .gt. probss) then
                     iproc = 0
                     istop = 1
                 endif
             endif
         endif


**  Update local pointers if medium or particle code has changed

         IF(IUPD.EQ.0) THEN
            IUPD = 1
         ENDIF
      


c      IF(NUMED.EQ.7.OR.numed.eq.2.or.numed.eq.8) THEN
c         WRITE(*,*) 'ACRsteps', numed,STEP,SNEXT,SUWTR,STEPUW, 
c     &        blamda,nstep
c     ENDIF

*
* *** Linear transport
*
         DO 20 I = 1,3


            vectmp=vect(i)+step*vect(i+3)
            if(vectmp.eq.vect(i))then
               if(vect(i+3).ne.0.) then
                  vectmp=vect(i)+abs(vect(i))*sign(1.,vect(i+3))*
     +                 epsmac
               endif
            endif
            vect(i)=vectmp
 20      CONTINUE



         SLENG = SLENG + STEP
*
* *** Update time of flight
*     
         TOFG = TOFG +STEP/CLIGHT*EFFNSG(BLAMDA)
*
* *** Apply the selected process if any
*
         IF (IPROC.EQ.0) GO TO 999
         NMEC = 1
         LMEC(NMEC) = IPROC
*     
* **  interact with water ? 
*     
         IF (IPROC.EQ.1) THEN
*     
*     *** Calculate the situation of Cherencov light after interaction
*     (absorption or scattering ?)
*     

            CALL SGABSC(BLAMDA,SUWTR) 
            if (OUTPUT_SCATTABLE.eq.1) ISTORY=ISTORY+1
            

C---- Update Constant value
            SLUWTR=SLENG
            STEPUW=BIG
         ENDIF
*
*------------------------END WATER MEDIA TRANSPORT
*

*
*------------------------ACRYLIC MEDIA TRANSPORT
*
**      Check if Medium is Acrylic Cover     
      ELSE IF(NUMED.EQ.7) THEN
*
* *** Compute current step size
*      
         STEP   = BIG
         IPROC  = 0
         BLAMDA = 1237.8E-09 / GEKIN
*
* *** Create interaction probabilities
*
         CALL GRNDM(RDUMMY,1)
         ZINTAC = -LOG(RDUMMY)
* 
* **  Step limitation due to absorption by acrylic
*
         STEPAC = ACRSG(BLAMDA)
         IF(STEPAC.LE.EPSMAC)THEN
            SUACR = BIG
         ELSE
            SUACR  = ZINTAC * 100./STEPAC
         ENDIF
         IF(SUACR.LT.STEP) THEN
            STEP = SUACR
            IPROC= 1
         ENDIF
*
* **  Step limitaion due to geometry?
*
         CALL GTNEXT
         IF (STEP.GT.SNEXT) THEN
            STEP   = SNEXT + PREC
            INWVOL= 2
            IPROC = 0
            NMEC =1
            LMEC(1)=1
         ENDIF

**  Update local pointers if medium or particle code has changed

         IF(IUPD.EQ.0) THEN
            IUPD = 1
         ENDIF
      
c      IF(NUMED.EQ.7.OR.numed.eq.2.or.numed.eq.8) THEN
c         WRITE(*,*) 'ACRsteps', numed,STEP,SNEXT,SUWTR,STEPUW, 
c     &        blamda,nstep
c     ENDIF

*
* *** Linear transport
*
         DO 120 I = 1,3
            vectmp=vect(i)+step*vect(i+3)
            if(vectmp.eq.vect(i))then
               if(vect(i+3).ne.0.) then
                  vectmp=vect(i)+abs(vect(i))*sign(1.,vect(i+3))*
     +                 epsmac
               endif
            endif
            vect(i)=vectmp
 120      CONTINUE
         
         SLENG = SLENG + STEP
*     
* *** Update time of flight
*
         TOFG = TOFG +STEP/CLIGHT*ACRNSG(BLAMDA)
*
* *** Apply the selected process if any
*
         IF (IPROC.EQ.0) GOTO 999
         NMEC = 1
         LMEC(NMEC) = IPROC
*
* **  interact with Acrylic? 
*
         IF (IPROC.EQ.1) THEN
*
* *** Absorbed Cherencov light is stopped
*        (absorption only, no scattering in acrylic)
            istop=1
C            WRITE(*,*) 'ABSORPTION in ACRYLIC now'
          

C---- Update Constant value
            STEPAC=BIG
         ENDIF
*
*------------------------END ACRYLIC MEDIA TRANSPORT
*
         
*
*------------------------PHOTON IN WRONG MEDIUM
*
      ELSE

         IF(INWVOL.EQ.0) THEN
            PRINT *,"ERROR IN GTCKOV"
            PRINT *,NUMED,INWVOL
            STOP
         ELSE IF(INWVOL.EQ.1)THEN
            PRINT *,'generated in PMT'
         ENDIF

      ENDIF

 999  CONTINUE      
      RETURN
      END


