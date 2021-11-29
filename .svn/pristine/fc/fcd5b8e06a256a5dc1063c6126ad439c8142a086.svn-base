*CMZ :  3.14/16 17/10/90  15.14.39  by  Rene Brun
*-- Author :
      SUBROUTINE SGABSC(BLAMDA,SUWTR)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *  Calculate the situation of Cherencov light after interaction  *
C.    *  (absorption or scattering ?)                                  *
C.    *                                                                *
C.    *    ==>Called by : GTCKOV                                       *
C.    *       Author    Y.Koshio                                       *
C.    *                                                                *
C     *        2020.02.04 Bug correction for TBA calculation           *
C     *                    Y.Nakajima, Y.Koshio                        *
C.    ******************************************************************
C.
C     Feb.- 5-96  bug on scattered direction was fixed. (M.Shiozawa)
C     Feb.-13-96  MIE scattering is introduced. (M.Shiozawa)
C     July.-1-01  Total Mie Scattering is introduced. (T.Shibata)
C                     Rayleigh = Real Rayleigh + Symmmetric Mie 
C                     Mie = Asymmetric Mie
C                   Total Mie = Symmmetric + Asymmetric scattering   
C                     ( Total Mie = WTMISG(BLAMDA) )
C     Apr.- 4-06  Scattering type name changed for SK-III (B. Yang)
C                  Symmetric scattering (Rayleigh type of scattered angle distribution)
C                    WSYMSCSG
C                  Asymmetric scattering (Asymmetric Mie type of scattered angle distribution)
C                    WASYSCSG
C                 All Absorption coefficient are merged at WABSSG.F
C                    WABSSG
C     Apr.- 8-07  Water scattering information of SK-III LOWE printed out to standard out modified    (B. Yang)
C                 Meaningless Safe guard commented out (B. Yang)
C
C     2009.11.18  Make top bottom asymmetrical MC (beta version)
C                 Fix bug-OD's water parameter become uniform (T.Yokozawa) 
c     2011.10.27  bug fixed in backward rayleigh scattering by jp/mine
c
c     2020.01.15  Make for SK-V
c
#include "skdetsim_parm.h"
#include "skhead.h"
#include "geant321/gctmed.inc"
#include "geant321/gconst.inc"
#include "geant321/gckine.inc"
#include "geant321/gctrak.inc"
#include "gcbank.h"
#include "epsmac.h"
#include "rayfac.h" 
#include "guwtr.h"

#include "dsparm.h"

      LOGICAL ROTATE
      Integer icall
      REAL WSYMSCSG,WASYSCSG,rend
      Data icall/0/
C---- counting error for all events
      SAVE error

      real abs_corr
      real abs_corr_c

C---- Calculate the Probability of scattering

      alpha = TBADS
c      if(SK_GEOMETRY.ge.SK_III .and. WATER_TUNE.eq.UNI_TUNE) then ! t-b asym
c      if((SK_GEOMETRY.ge.SK_III .and. (WATER_TUNE.eq.UNI_TUNE
      if(((SK_GEOMETRY.eq.SK_III .or. SK_GEOMETRY.eq.SK_IV) .and.
     &     (WATER_TUNE.eq.UNI_TUNE .or. WATER_TUNE.eq.10.or.WATER_TUNE.eq.11))
     &     .and. ((NUMED .ne. 9) .and. (NUMED .ne. 1)) ) then ! only ID

         if(IWATERP.eq.200904) then ! keep bug for ATMPD
            rend = TRACKLSK(vect(1),vect(4))
            if(vect(3)+SUWTR*vect(6).ge.-1200)then
               if(SUWTR.lt.rend)then
                  COEFF = WTRSG(BLAMDA)+WABSSG(BLAMDA)*(alpha*(vect(3)+SUWTR*vect(6)))
               else
                  COEFF = WTRSG(BLAMDA)+WABSSG(BLAMDA)*(alpha*(vect(3)+rend*vect(6)))
               end if
            else
               if(SUWTR.ge.rend)then
                  COEFF = WTRSG(BLAMDA)+WABSSG(BLAMDA)*(alpha*(vect(3)+rend*vect(6)))
               else
                  COEFF = WTRSG(BLAMDA)+WABSSG(BLAMDA)*(alpha*(-1200.))
               end if
            end if
         else ! bug fix for LOWE

c     Add correction to the absorption factor to be consistent with the treatment in gtckov.F
c     In gtckov.F abs -> abs * 1/(1+c)
c     Here,  abs -> abs - abs * c/(1+c)         
            abs_corr_c = 69.63*alpha
            if (abs_corr_c .ne. 0) then
               abs_corr = WABSSG(BLAMDA) * abs_corr_c /(1.0+abs_corr_c)
            else
               abs_corr = 0.
            end if
            if(vect(3).ge.-1100)then
               COEFF = WTRSG(BLAMDA)+WABSSG(BLAMDA)*alpha*vect(3) - abs_corr
            else
               COEFF = WTRSG(BLAMDA)+WABSSG(BLAMDA)*alpha*(-1100.) - abs_corr
            end if
         endif

      else if(SK_GEOMETRY.eq.SK_V .and. WATER_TUNE.eq.UNI_TUNE ! for SK_V (no zboundary)
     &        .and. (NUMED.ne.9 .and. NUMED.ne.1) ) then ! only ID
         COEFF = WTRSG(BLAMDA)+WABSSG(BLAMDA)*alpha*vect(3)

      else
         COEFF = WTRSG(BLAMDA) 
      end if 
         
      IF ((SK_GEOMETRY.EQ.SK_III) .or. (SK_GEOMETRY.EQ.SK_IV) .or. (SK_GEOMETRY.EQ.SK_V) ) THEN
C----------------------------------------------------------------------
C----------------------------------------------------------------------
CC  SK-III scattering probability calculation : Apr.- 4-07 B. Yang  CC
C----------------------------------------------------------------------
C----------------------------------------------------------------------
        IF (WATER_TUNE .NE. AP_TUNE) THEN
          COEFRS = WSYMSCSG(BLAMDA)
          COEFMI = WASYSCSG(BLAMDA)
          PROBSC = (COEFRS+COEFMI)/COEFF
          PROBRS = COEFRS/COEFF
        ELSE                     ! original ap_tune
          COEFRS = WSYMSCSG(BLAMDA)
          COEFMI = WASYSCSG(BLAMDA)
          COEFSC = COEFRS + COEFMI
          PROBSC = COEFSC/COEFF
          PROBRS = COEFRS/COEFSC
C        ELSE
C          WRITE(0,*) 'ERROR : WATER_TUNE is not usual,', WATER_TUNE
C          STOP
        ENDIF
        IF (icall.eq.0) then
          ATTEL = 1.0/COEFF
          print *
          print *, '******************** water parameter **********'
     /            //'**********'
          IF (WATER_TUNE .NE. AP_TUNE) THEN
	     if(WATER_TUNE.eq.11) then
               print *, 'Gd effect as absorption'
             endif
             if(WATER_TUNE.eq.10) then
               print *,'Gd effect as scattering'
             endif 	
             if (SK_GEOMETRY.EQ.SK_III) then
                print *, 'SHORT abs. =',absfac
                print *, 'Sym. Sctt. =',rayfac
                print *, 'Asy. Sctt. =',miefac
             else if (SK_GEOMETRY.GE.SK_IV) then
                print *, '<Absorption>'
                print *, '   absfac  =',absfac
                print *, '   abssfac =',abssfac
                print *, '   abspow  =',abspow
                print *, '   absfit  =',absfit
                print *, '<Symmetric Scattering>'
                print *, '   rayfac  =',rayfac
                print *, '   raysfac =',raysfac
                print *, '<Asymmetric Scattering>'
                print *, '   miefac  =',miefac
                print *, '   miesfac =',miesfac
                print *, '   asyfit  =',asyfit
             endif

          ELSE                   ! original ap_tune
            print *, 'ABSORPTION AABL = ', AABL
            print *, 'SYMM. SCTT ARAS = ', ARAS
            print *, 'ASYM. SCTT AMIS = ', AMIS
          ENDIF
          print *, '***********************************************'
     /          //'**********'
          print *
          print *, '*********************** Probability **************'
          print *, 'LAMBDA     = ',BLAMDA
          print *, 'ABSORPTION  wabssg  = ',WABSSG(BLAMDA)
          print *, 'SYMM. SCTT wsymscsg = ',WSYMSCSG(BLAMDA)
          print *, 'ASYM. SCTT wasyscsg = ',WASYSCSG(BLAMDA)
          print *, 'ATTENUATION LENGHTH = ', ATTEL, 'm'
          print *, '**************************************************'
          print *
          icall = 1
        Endif


C---- Scattering or Absorption ? : ABSORPTION -> GOTO 100
        CALL GRNDM(RNSC,1)
        IF(RNSC.GE.PROBSC) GOTO 100

        IF (WATER_TUNE .EQ. AP_TUNE) THEN
          CALL GRNDM(RNRS,1)
          ROTATE = .FALSE.
        ENDIF
C---------------------------------------------------------------------
CC     SK-III scattering probability calculation over               CC
C---------------------------------------------------------------------

      ELSE IF (SK_GEOMETRY.EQ.SK_I .OR. SK_GEOMETRY.EQ.SK_II) THEN
C---------------------------------------------------------------------
C---------------------------------------------------------------------
CC     SK-I & II scattering probability calculation                 CC
C---------------------------------------------------------------------
C---------------------------------------------------------------------
        IF (WATER_TUNE .EQ. LE_TUNE) THEN
          COEFRS = WRAYSG(BLAMDA)
          COEFMI = WMIESG(BLAMDA)
        ELSE
          COEFSC = ASCL*WLNGSG(BLAMDA) + (ARAS+APSS)*WSRTSG(BLAMDA)
     &          + (1.0)*WTMISG(BLAMDA)
        END IF
        IF (icall.eq.0) THEN
           print *
           print *, '******************** water parameter **********'
     /            //'**********'
           IF (WATER_TUNE .EQ. LE_TUNE) THEN
              print *, 'sgabsc.F: Water_Tune = Le_Tune'
              print *, 'SHORT: abs. =',absfac
              print *, '     : scat.(RAY)=',rayfac
              print *, '     : scat.(MIE)=',miefac
           ELSE
              print *, 'LONG : abs., scat.= AABL, ASCL: ',AABL,ASCL
              print *, 'SHORT: abs., ray., mie, scat.= '
              print *, 'APAS,ARAS,AMIS,APSS',
     ,                  APAS,ARAS,AMIS,APSS
           ENDIF
           print *, '***********************************************'
     /          //'**********'
           print *
           IF (WATER_TUNE .EQ. AP_TUNE) THEN
              ATTEL = 1.0/COEFF
              print *, '*********************** Probability ***********'
              print *, 'sgabsc.F: Water_Tune = AP_Tune'
              print *, 'ABSORPTION =',WLNGSG(BLAMDA)
              print *, 'RAYLEIGH = ',WSRTSG(BLAMDA) 
              print *, 'TOTAL Mie =',WTMISG(BLAMDA)
              print *, 'ASYMMETRIC MIE =',WMIESG(BLAMDA)
              print *, 'ATTENUATION LENGHTH =',ATTEL
              print *, '***********************************************'
              print *
           END IF
           icall = 1
        ENDIF
        IF (WATER_TUNE .EQ. LE_TUNE) THEN
          PROBSC = (COEFRS+COEFMI)/COEFF
          PROBRS = COEFRS/COEFF
        ELSE
          PROBSC = COEFSC/COEFF
        END IF

C---- Scattering or Absorption ?

        CALL GRNDM(RNSC,1)
        IF(RNSC.GE.PROBSC) GOTO 100

C---- Rayreigh Scattering:
C      not all photons are tracked:
C      PROBRS is not Rayleigh . PROBRS is symmmetric scattering,
C--    Rayreigh or NOT
        IF (WATER_TUNE .EQ. AP_TUNE) THEN
           COEFRS = ARAS*WSRTSG(BLAMDA)
           COEFMI = AMIS*WMIESG(BLAMDA)
           COEFTM = 1.0*WTMISG(BLAMDA)
           COEFSM = COEFTM - COEFMI 
           PROBRS = (COEFRS+ COEFSM)/COEFSC  
           PROBMIE = (COEFRS+COEFTM)/COEFSC
           CALL GRNDM(RNRS,1)
           ROTATE = .FALSE.
C--after ASCL & APSS was rejected, PROBMIE is meaningless. That is 1
           IF (RNRS .GT. PROBMIE) THEN
              CALL GRNDM(RDUMMY,1)
              COSTH = (RDUMMY-0.5)*2.0
              TH = ACOS(COSTH)
              CALL GRNDM(RDUMMY,1)
              PHI = RDUMMY*TWOPI
              COSTH=COS(TH)
              SINTH=SIN(TH)
              COSPH=COS(PHI)
              SINPH=SIN(PHI)
              VECT(4) = SINTH*COSPH
              VECT(5) = SINTH*SINPH
              VECT(6) = COSTH
              polar(1) = costh*cosph
              polar(2) = costh*sinph
              polar(3) = -sinth
              RETURN
           END IF
        END IF
C----------------------------------------------------------------------
CC     SK-I & II  scattering probability calculation over           CC
C----------------------------------------------------------------------

      ELSE
        WRITE(0,*)
     &        'ERROR : SK_GEOMETRY is not supported yet,',SK_GEOMETRY 
        STOP
      ENDIF

C----------------------------------------------------------------------
C----------------------------------------------------------------------
CC     Handlie Scattering type and Scatterd angle                   CC
C----------------------------------------------------------------------
C----------------------------------------------------------------------
      totev = totev + 1
      IF (abs(vect(4)*polar(1)+vect(5)*polar(2)+vect(6)*polar(3))
     &     .gt.0.0001) then
         error = error+1
         If (ifix(error/10000)*10000.eq.error) then
            print *, 'Polar. error=',error,'/',totev
         Endif
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     averaged rayleigh scattering.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         IF (((WATER_TUNE .NE. AP_TUNE) .AND.
     .        (RNSC .LT. PROBRS)) .OR.
     .       ((WATER_TUNE .EQ. AP_TUNE) .AND.
     .        (RNRS .LE. PROBRS))) THEN  
            RAYMAX = 2.
            CALL GFANG(VECT(4),COSTH,SINTH,COSPH,SINPH,ROTATE)
 11         CONTINUE
            CALL GRNDM(COSTHRAY,1)
ccc jp/mine on 27oct,2011
            COSTHRAY = (COSTHRAY-.5)*2.
            RAY = 1.+ COSTHRAY**2
            CALL GRNDM(RDUMMY,1)
            IF(RDUMMY*RAYMAX.GT.RAY) GOTO 11
            TH  = ACOS(COSTHRAY)
            CALL GRNDM(RDUMMY,1)
            PHI = TWOPI*RDUMMY
            VECT(4) = SIN(TH)*COS(PHI)
            VECT(5) = SIN(TH)*SIN(PHI)
            VECT(6) = COS(TH)
            IF(ROTATE)CALL GDROT(VECT(4),COSTH,SINTH,COSPH,SINPH)
         ELSE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     averaged MIE scattering(no care for polarization)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
 13         CALL GRNDM(RDUMMY,1)
            COSTH = (RDUMMY-0.5)*2.0
            CALL GRNDM(RDUMMY,1)
            if (RDUMMY.gt.SGMIES(costh)) Goto 13

            CALL GFANG(VECT(4),rotate_cost,
     ,                 rotate_sint,rotate_cosp,rotate_sinp,ROTATE)
            TH = ACOS(COSTH)
            CALL GRNDM(RDUMMY,1)
            PHI   = TWOPI*RDUMMY
            COSTH=COS(TH)
            SINTH=SIN(TH)
            COSPH=COS(PHI)
            SINPH=SIN(PHI)
            vect(4) = sinth*cosph
            vect(5) = sinth*sinph
            vect(6) = costh
            IF (ROTATE) CALL GDROT(VECT(4),rotate_cost,
     ,          rotate_sint,rotate_cosp,rotate_sinp)
         Endif
      ELSE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     polarized rayleigh scattering.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C rotate
         vectxy = sqrt(vect(4)**2.0+vect(5)**2.0)
         If (vectxy.gt.0.) then
            rotate_sinp_p = vect(5)/vectxy
            rotate_cosp_p = vect(4)/vectxy
            Call GDROT(vect(4),1.0,0.,rotate_cosp_p,-rotate_sinp_p)
            Call GDROT(polar(1),1.0,0.,rotate_cosp_p,-rotate_sinp_p)
         Endif
         vectxz = sqrt(vect(4)**2.0+vect(5)**2.0+vect(6)**2.0)
         If (vectxz.gt.0.) then
            rotate_sint = sqrt(vect(4)**2.0+vect(5)**2.0)/vectxz
            rotate_cost = vect(6)/vectxz
            Call GDROT(vect(4),rotate_cost,-rotate_sint,1.0,0.)
            Call GDROT(polar(1),rotate_cost,-rotate_sint,1.0,0.)
         Endif
         polarxy = sqrt(polar(1)**2.0+polar(2)**2.0)
         If (polarxy.gt.0.) then
            rotate_sinp = polar(2)/polarxy
            rotate_cosp = polar(1)/polarxy
            Call GDROT(polar(1),1.0,0.,rotate_cosp,-rotate_sinp)
            Call GDROT(vect(4),1.0,0.,rotate_cosp,-rotate_sinp)
         Endif
         If (polar(1).lt.0.9999.or.vect(6).lt.0.9999) then
            print *, 'sgabsc: cannot rotate photon'
     &              //' vector/polarization.'
            print *, '   gamma(xyz)=',vect(4),vect(5),vect(6)
            print *, '   polar(xyz)=',polar(1),polar(2),polar(3)
            print *, '                           naiseki='
     &             ,vect(4)*polar(1)+vect(5)*polar(2)+vect(6)*polar(3)
         Endif
            
C     scatter
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     polarized rayleigh scattering.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         IF (((WATER_TUNE .NE. AP_TUNE) .AND.
     .        (RNSC .LT. PROBRS)) .OR.
     .       ((WATER_TUNE .EQ. AP_TUNE) .AND.
     .        (RNRS .LE. PROBRS))) THEN
            RAYMAX = 1.
 10         CONTINUE
            CALL GRNDM(RDUMMY,1)
            COSTH = (RDUMMY-0.5)*2.0
            TH = ACOS(COSTH)
            CALL GRNDM(RDUMMY,1)
            PHI   = TWOPI*RDUMMY
            RAY= 1.- sin(TH)**2.0*cos(PHI)**2.0
            CALL GRNDM(RDUMMY,1)
            IF(RDUMMY*RAYMAX.GT.RAY) GOTO 10
         ELSE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     polarized MIE scattering.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
 14         CALL GRNDM(RDUMMY,1)
            COSTH = (RDUMMY-0.5)*2.0
            CALL GRNDM(RDUMMY,1)
            if (RDUMMY.gt.SGMIES(costh)) Goto 14
            TH = ACOS(COSTH)
            CALL GRNDM(RDUMMY,1)
            PHI   = TWOPI*RDUMMY
         Endif
         COSTH=COS(TH)
         SINTH=SIN(TH)
         COSPH=COS(PHI)
         SINPH=SIN(PHI)
         vect(4) = sinth*cosph
         vect(5) = sinth*sinph
         vect(6) = costh
         polar(1) = sqrt(1.-sinth**2.0*cosph**2.0)
         polar(2) = -(sinth**2.0*sinph*cosph)/polar(1)
         polar(3) = -(sinth*costh*cosph)/polar(1)
            
C     rotate bak
         If (polarxy.gt.0.) then
            Call GDROT(polar(1),1.0,0.,rotate_cosp,rotate_sinp)
            Call GDROT(vect(4),1.0,0.,rotate_cosp,rotate_sinp)
         Endif
         If (vectxz.gt.0.) then
            Call GDROT(vect(4),rotate_cost,rotate_sint,1.0,0.)
            Call GDROT(polar(1),rotate_cost,rotate_sint,1.0,0.)
         Endif
         If (vectxy.gt.0.) then
            Call GDROT(vect(4),1.0,0.,rotate_cosp_p,rotate_sinp_p)
            Call GDROT(polar(1),1.0,0.,rotate_cosp_p,rotate_sinp_p)
         Endif
            
         IF (abs(vect(4)*polar(1)+vect(5)*polar(2)+vect(6)*polar(3))
     &        .gt.0.0001) then
            print *, 'sgabsc: cannot rotate back photon'
     &              //'vector/polar.'
            print *, '      gamma(xyz)=',vect(4),vect(5),vect(6)
            print *, '      polar(xyz)=',polar(1),polar(2),polar(3)
            print *, '                              naiseki='
     &              ,vect(4)*polar(1)+vect(5)*polar(2)+vect(6)*polar(3)
         Endif
      Endif

      GOTO 200

C---- Absorption

 100  ISTOP=1

 200  RETURN
      END
