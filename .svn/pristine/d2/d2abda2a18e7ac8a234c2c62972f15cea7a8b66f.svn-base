***********************************************************************
*     -----------------
      SUBROUTINE tr2fac(tr)
*     -----------------
* 
*     Called from dscard.F and gukine.F
*
*     (Purpose)
*     (SK-I) translate attenuation length to Rayleigh coefficient
*     (SK-II&III) translate to Rayleigh coefficient
*
*     (Creation Date and Author)
*     2003/01/16; Merge to skdetsim from ledetsim
*     2003/05/19; J.H.Yoo : Tune SK-II by SK-I data
*     2003.06.07; rfact from exp(23.57) to 152E8 by K.Ishihara
*     2003.06.28; C. Sterner : Replace SK-I tune removed by J.H.Yoo
*     2006.04.18 ; K Bays : Updated for SK-III
*     2007.04.08 ; B. Yang : SK-III scattering tuned by water scattering laser data
*     2007.04.13 ; B. Yang : SK-III absorption tuned by water scattering laser data
*     2007.09.28 ; B. Yang : SK-III water parameter retuned due to MC bug
*     2008.03.14 ; B. Yang : SK-III water parameter retuned with Feb, 08 data
*     2009.03.14 ; B. Yang : SK-III water parameter retuned with Feb, 08 data
*     2009.05.16 ; B. Yang : SK-III water parameter as a function of time
*     2009.09.04 ; T. Mclachlan : SK-IV monthly water parameters installed.
*     2010.09.07 ; S.H. Cho : SK-IV monthly water parameters installed with new laser
*                       Monthly variation set by variables in this file.
*     2011.09.23 ; Y.Koshio : SK-IV parameters for LINAC 2010 by T.Yokozawa
***********************************************************************

      IMPLICIT NONE

#include "skdetsim_parm.h"
#include "skhead.h"
#include "rayfac.h"
#include "dsparm.h"

c   absfac about water transparency was introduced by B. Yang 2009/05/16
C   related variables were also defined
      integer NMAX, I
      PARAMETER(NMAX = 23)
      real*4 tr, absfac0, wt2afac, sctt_part,absfac1,absfac2
      PARAMETER(absfac0 = 0.2904585E+08, wt2afac = 0.1208289E+13, sctt_part = 0.5898086E+08)

c----

      IF (WATER_TUNE .EQ. AP_TUNE) THEN
         write(0,*) 'TR2FAC is not used for AP_TUNE'
         stop
      ENDIF

CCS-----------------------------
CCS         SK-I from ledetsim
CCS
CCS-----------------------------
      IF(SK_GEOMETRY.eq.SK_I) THEN
         call trabs(tr,absfac,rayfac,miefac)
         print *,'abs,ray,mie',absfac,rayfac,miefac
         return
      
CCS-----------------------------
CCS         SK-II Temporary tuning
CCS
CCS-----------------------------

c========================================
c      absfac=ABS  ! (Absorption parameter)
c      rayfac=RAY  ! (Rayleigh parameter)
c      miefac=MIE  ! (Mie correction)
c========================================

      ELSE IF(SK_GEOMETRY.eq.SK_II) THEN

c   The scattering parameters tuned by Parker 2003/11/01

         absfac = 0.0
         rayfac = 1.70E8
         miefac = 0.05E8

CCS------------------------END SK-II---------------------------------------
      
c====== SK-III case=======

      ELSE IF(SK_GEOMETRY.eq.SK_III) THEN
c   The scattering parameters except for short wavlength absorption parameter tuned by B. Yang 2007/04/08
c   But the definition of scattering attenuation length is changed. Refer wsymscsg.F and wasyscsg.F

c   The short wavlength absorption parameter tuned by B. Yang 2007/04/13
c   But the long wavelength absorption parameter is defined at wabssg.F

C      absfac = 0.6586E8
C      rayfac = 1.137E8
C      miefac = 0.15E8

c   Water parameter retuned due to bug about black sheet by B. Yang 2007/09/28

C      absfac = 0.7014E8
C      rayfac = 0.8669E8
C      miefac = 0.1369E8

c   Water parameter retuned after reflection tuning on PMTs by B. Yang 2008/03/14

C         absfac = 0.2904585E+08	! 2008/06/XX
         absfac = absfac0
         if(WATER_TUNE.eq.UNI_TUNE) then
            rayfac = 0.8582E+08*1.05
         else
            rayfac = 0.8582E+08
         endif
         miefac = 0.05476010E+08

         if( tr .le. 0 ) then	! set them default value
            write(*,*) 'Warning: absfac =', absfac, 'at water transparency = ', tr,'cm'
            absfac = absfac0
            call set_wabscoeff(tr) 
            return
         endif

CCCC   absfac and LELONGFAC about water transparency was introduced by B. Yang 2008/06/XX
CCCC      absfac = wt2afac*(1.0/tr-1.0/tr0) + absfac0
CCCC      absfac = wt2afac/tr - sctt_part
         absfac = -.3383244E+08+exp(16.82205+14827.33/tr)

         if( tr .ge. 12500. ) then
            absfac1 = 0.3106475E+07+exp(14.72945+30522.71/tr)
            if (absfac1 .ge. absfac ) absfac = (absfac+2*absfac1)/3.
         elseif( tr .le. 10000. ) then
            absfac2 = -.6063145E+08 + 0.1166306E+13/tr
            if (absfac2 .le. absfac ) absfac = (absfac+2*absfac2)/3.
         endif

         if( absfac .lt. 0.0) then
            write(*,*) 'Warning: absfac =', absfac, 'at water transparency = ', tr,'cm'
            write(*,*) '         forced to be 0.0'
            absfac = 0.0
         endif

         call set_wabscoeff(tr)

c===== end SK-III =======

c====== SK-IV case=======

      ELSE IF(SK_GEOMETRY.eq.SK_IV) THEN

         absfac  = 0.2961E+8
         abssfac = 1.164
         abspow  = 10.94
         miefac  = 0.4623E+7
         miesfac = 392.4

         call set_monthly_waterp_sk4(tr) ! set the time dependent parameters

c===== end SK-IV =======


c====== SK-V case (given by J.Yang / calibration group based on Aug.-Sep.2019 data) =======
c====== SK-V case (given by Ji-Woong Seo /calibration group in Nov.26 based on Aug.-Sep.2019 data) =======
       

      ELSE IF(SK_GEOMETRY.eq.SK_V) THEN

c         absfac  = 0.2961E+8  ! fixed from last tuning in SK-IV.
c         abssfac = 2.703172   ! 0.0754185 / 0.0279
c         abspow  = 10.94      ! fixed from last tuning in Sk-IV
c         miefac  = 0.276092E+7  
c         miesfac = 392.4      ! fixed

c         absfac  = 0.196324E+8  ! P1. Fitted by Seo.
c         abssfac = 0.984        ! P2./0.0279 = 0.0274536 / 0.0279 . 0.0279 is from A(1) in wabssg.F
c         abspow  = 9.79045      ! P3. Fitted by Seo.
c         miefac  = 0.0299065E+7 ! P7
c         miesfac = 877.109      ! P8

c        Updated by Seo san Jan.14 2021
         absfac  = 0.199448E+8  ! P1. Fitted by Seo.
         abssfac = 1.21539      ! P2./0.0279 = 0.0339093 / 0.0279 . 0.0279 is from A(1) in wabssg.F
         abspow  = 12.3540      ! P3. Fitted by Seo.
         miefac  = 0.0299065E+7 ! P7
         miesfac = 877.109      ! P8         

         call set_monthly_waterp_sk4(tr) ! set the time dependent parameters

c===== end SK-IV =======
         
      ENDIF
      return
      end


      SUBROUTINE set_wabscoeff(tr)

      IMPLICIT NONE

#include "rayfac.h"

      integer NMAX, I
      PARAMETER(NMAX = 23)
      real*4 absfac0, LELONGFAC0,tr
      real*4 LELFL(NMAX),BLSPL(NMAX), BLSPLFAC

      PARAMETER(absfac0 = 0.2904585E+08)
      PARAMETER(LELONGFAC0 = 0.8013058)
      PARAMETER(BLSPLFAC = 450.5400710)	! = 500. * (POPEFRYSG(450.003)/A(1))**(1/LELONGPOW)

      DATA LELFL/0.533415, 0.549915, 0.556915, 0.588415, 0.601415, 0.637415, 0.666415,
     & 0.673416, 0.674416, 0.684916, 0.702416, 0.704416, 0.705416, 0.712916, 0.733416,
     & 0.746916, 0.8013058,0.825416, 0.850416, 0.881916, 0.940416, 1.003916, 1.01279 /

      DATA BLSPL/507.484,  504.683,  500.04,   497.446,  495.023,  492.424,  489.997,
     & 488.092,  484.877,  482.6655, 490.01,   479.118,  474.825,  472.526,  469.997,
     & 467.5255, 463.918,  462.459,  460.084,  457.504,  454.9675, 452.3955, 450.003 /

      if( tr .le. 0 ) then	! set them default value
        absfac = absfac0
        LELONGFAC = LELONGFAC0
        BLAMDASP = BLSPL(17)
        return
      endif

      LELONGFAC = LELONGFAC0 * absfac / absfac0

CCCC   An cross point(BLAMDASP) between POPE/FRY's absorption coefficients and TUNED LONGER wavelength fitting curve
CCCC   about LELONGFACwas introduced by B. Yang 2008/06/XX

      if ( LELONGFAC .GT. LELFL(NMAX) ) then
        BLAMDASP = BLSPLFAC * (1.0/LELONGFAC)**(1.0/LELONGPOW)
CC around 2.108794, 420nm
        if ( BLAMDASP .lt. 420. ) write(*,*) 'Warning(tr2fac) : BLAMDASP is shorter than 420nm',BLAMDASP,'nm, at', tr,'cm'
      elseif ( LELONGFAC .GE. LELFL(1) ) then
        DO I=1,NMAX-1
          if(LELONGFAC.le.LELFL(I+1)) then
            BLAMDASP = BLSPL(I) + (BLSPL(I+1)-BLSPL(I))/(LELFL(I+1)-LELFL(I))*(LELONGFAC-LELFL(I))
            GOTO 20
          endif
        ENDDO
      else
        BLAMDASP = BLSPL(1)
      endif
   20 CONTINUE

      return
      end

      SUBROUTINE set_monthly_waterp_sk4(tr)

*     Monthly Water Parameter Tuning. 
*     Thomas McLachlan. Last Update: 20/04/2010
*     SungHoon Cho. Last Update: 07/09/2010
*     Yokozawa Last Update untl Mar.2011
*     Yokozawa add IWATERP=999999 : day by day for SK-IV watert->absfit 
*     Koshio add IWATERP=888888 : water dep absfit

      IMPLICIT NONE

      real tr,absfac1,absfac2,yabsfac,absfitfac


#include "rayfac.h"
#include "dsparm.h"

      if (IWATERP .eq. 200810) then ! 10/2008 data

         absfit  = 0.9053885
         rayfac  = 0.9229056E+08
         raysfac = 95622.83
         asyfit  = 0.1879478E-03


      elseif (IWATERP .eq. 200811) then ! 11/2008 data

         absfit  = 1.009911
         rayfac  = 0.9233614E+08
         raysfac = 96264.66
         asyfit  = 0.1295007E-03


      elseif (IWATERP .eq. 200812) then ! 12/2008 data

         absfit  = 1.064494
         rayfac  = 0.8506918E+08
         raysfac = 115243.5
         asyfit  = 0.1287528E-03


      elseif (IWATERP .eq. 200901) then ! 1/2009 data

         absfit  = 0.7194061
         rayfac  = 0.9025358E+08
         raysfac = 100854.4
         asyfit  = 0.1469163E-03


      elseif (IWATERP .eq. 200902) then ! 2/2009 data

         absfit  = 0.7888419
         rayfac  = 0.8915427E+08
         raysfac = 102780.9
         asyfit  = 0.1470308E-03


      elseif (IWATERP .eq. 200903) then ! 3/2009 data

         absfit  = 0.6696396
         rayfac  = 0.9344493E+08
         raysfac = 90823.89
         asyfit  = 0.1376437E-03


      elseif (IWATERP .eq. 200904) then ! 4/2009 data

         absfit  = 0.6429981
         rayfac  = 0.8615809E+08
         raysfac = 111136.5
         asyfit  = 0.9859612E-04


      elseif (IWATERP .eq. 200905) then ! 5/2009 data

         absfit  = 0.6651415
         rayfac  = 0.8607781E+08
         raysfac = 108982.0
         asyfit  = 0.1005731E-03


      elseif (IWATERP .eq. 200906) then ! 6/2009 data

         absfit  = 0.7091202
         rayfac  = 0.8465350E+08
         raysfac = 110403.9
         asyfit  = 0.1090088E-03


      elseif (IWATERP .eq. 200907) then ! 7/2009 data

         absfit  = 0.6982910
         rayfac  = 0.8911509E+08
         raysfac = 101877.2
         asyfit  = 0.1035123E-03


      elseif (IWATERP .eq. 200908) then ! 8/2009 data

         absfit  = 0.5973664
         rayfac  = 0.8987661E+08
         raysfac = 100035.3
         asyfit  = 0.1086015E-03

      elseif (IWATERP .eq. 20090728) then ! Linac 
         absfit  = 0.6649       ! Yokozawa test
         rayfac  = 0.9109565E+08
         raysfac = 96678.71
         asyfit  = 0.1057341E-03
         write(*,*) '***************************************'
         write(*,*) 'linac elapse day 9'
         write(*,*) '***************************************'

      elseif (IWATERP .eq. 20090729) then ! Linac 
         absfit  = 0.657        ! Yokozawa test
         rayfac  = 0.9109565E+08
         raysfac = 96678.71
         asyfit  = 0.1057341E-03
         write(*,*) '***************************************'
         write(*,*) 'linac elapse day 10'
         write(*,*) '***************************************'

      elseif (IWATERP .eq. 20090730) then ! Linac 
         absfit  = 0.6491       ! Yokozawa test
         rayfac  = 0.9109565E+08
         raysfac = 96678.71
         asyfit  = 0.1057341E-03
         write(*,*) '***************************************'
         write(*,*) 'linac elapse day 11'
         write(*,*) '***************************************'

      elseif (IWATERP .eq. 20090731) then ! Linac 
         absfit  = 0.6412       ! Yokozawa test
         rayfac  = 0.9109565E+08
         raysfac = 96678.71
         asyfit  = 0.1057341E-03
         write(*,*) '***************************************'
         write(*,*) 'linac elapse day 12'
         write(*,*) '***************************************'

      elseif (IWATERP .eq. 20090801) then ! Linac 
         absfit  = 0.6333       ! Yokozawa test
         rayfac  = 0.9109565E+08
         raysfac = 96678.71
         asyfit  = 0.1057341E-03
         write(*,*) '***************************************'
         write(*,*) 'linac elapse day 13'
         write(*,*) '***************************************'

      elseif (IWATERP .eq. 20090802) then ! Linac 
         absfit  = 0.6254       ! Yokozawa test
         rayfac  = 0.9109565E+08
         raysfac = 96678.71
         asyfit  = 0.1057341E-03
         write(*,*) '***************************************'
         write(*,*) 'linac elapse day 14'
         write(*,*) '***************************************'

      elseif (IWATERP .eq. 20090803) then ! Linac 
         absfit  = 0.6175       ! Yokozawa test
         rayfac  = 0.9109565E+08
         raysfac = 96678.71
         asyfit  = 0.1057341E-03
         write(*,*) '***************************************'
         write(*,*) 'linac elapse day 15'
         write(*,*) '***************************************'

      elseif (IWATERP .eq. 20090804) then ! Linac 
         absfit  = 0.6096       ! Yokozawa test
         rayfac  = 0.9109565E+08
         raysfac = 96678.71
         asyfit  = 0.1057341E-03
         write(*,*) '***************************************'
         write(*,*) 'linac elapse day 16'
         write(*,*) '***************************************'

      elseif (IWATERP .eq. 20090805) then ! Linac 
         absfit  = 0.6017       ! Yokozawa test
         rayfac  = 0.9109565E+08
         raysfac = 96678.71
         asyfit  = 0.1057341E-03
         write(*,*) '***************************************'
         write(*,*) 'linac elapse day 17'
         write(*,*) '***************************************'

      elseif (IWATERP .eq. 20090806) then ! Linac 
         absfit  = 0.5938       ! Yokozawa test
         rayfac  = 0.9109565E+08
         raysfac = 96678.71
         asyfit  = 0.1057341E-03
         write(*,*) '***************************************'
         write(*,*) 'linac elapse day 18'
         write(*,*) '*************************************** '

      elseif (IWATERP .eq. 20090807) then ! Linac 
         absfit  = 0.5859       ! Yokozawa test
         rayfac  = 0.9109565E+08
         raysfac = 96678.71
         asyfit  = 0.1057341E-03
         write(*,*) '***************************************'
         write(*,*) 'linac elapse day 19'
         write(*,*) '***************************************'

      elseif (IWATERP .eq. 20090808) then ! Linac 
         absfit  = 0.578  
         rayfac  = 0.9109565E+08
         raysfac = 96678.71
         asyfit  = 0.1057341E-03
         write(*,*) '***************************************'
         write(*,*) 'linac elapse day 20'
         write(*,*) '***************************************'




      elseif (IWATERP .eq. 200909) then ! 9/2009 data

         absfit  = 0.6396149
         rayfac  = 0.8990191E+08
         raysfac = 99664.21
         asyfit  = 0.1046891E-03


      elseif (IWATERP .eq. 200910) then ! 10/2009 data

         absfit  = 0.6834394
         rayfac  = 0.8745879E+08
         raysfac = 106454.2
         asyfit  = 0.1095626E-03


      elseif (IWATERP .eq. 200911) then ! 11/2009 data

         absfit  = 0.5977888
         rayfac  = 0.9133789E+08
         raysfac = 96468.29
         asyfit  = 0.1137638E-03


      elseif (IWATERP .eq. 200912) then ! 12/2009 data

         absfit  = 0.5748957
         rayfac  = 0.9297710E+08
         raysfac = 91879.08
         asyfit  = 0.1144302E-03

C     1/2010, 2/2010 Water parameters seem to be wrong. Changed by S.H. Cho 10/09/07
C         elseif (IWATERP .eq. 201001) then  ! 1/2010 data
C
C          absfit  = 0.9347053
C          rayfac  = 0.4076335E+08
C          raysfac = 344329.0
C          asyfit  = 0.1268748E-03
C
C
C         elseif (IWATERP .eq. 201002) then  ! 2/2010 data
C
C          absfit  = 1.018283
C          rayfac  = 0.2273096E+08
C          raysfac = 714595.5
C          asyfit  = 0.1166678E-03

      elseif (IWATERP .eq. 201001) then ! 1/2010 data

         absfit  = 0.6358912
         rayfac  = 0.9332435E+08
         raysfac = 90964.41
         asyfit  = 0.1099797E-03

      elseif (IWATERP .eq. 201002) then ! 2/2010 data

         absfit  = 0.6487944
         rayfac  = 0.8239963E+08
         raysfac = 123917.6
         asyfit  = 0.1045012E-03


      elseif (IWATERP .eq. 201003) then ! 3/2010 data

         absfit  = 0.6148528
         rayfac  = 0.8915349E+08
         raysfac = 102369.9
         asyfit  = 0.1018177E-03


      elseif (IWATERP .eq. 201004) then ! 4/2010 data

         absfit  = 0.6548554
         rayfac  = 0.9369887E+08
         raysfac = 93688.78
         asyfit  = 0.1056951E-03


      elseif (IWATERP .eq. 201005) then ! 5/2010 data

         absfit  = 0.6710849
         rayfac  = 0.8714606E+08
         raysfac = 106780.0
         asyfit  = 0.1108840E-03


      elseif (IWATERP .eq. 201006) then ! 6/2010 data

         absfit  = 0.6703039
         rayfac  = 0.9370E+08
         raysfac = 88800.0
         asyfit  = 0.1008570E-03


      elseif (IWATERP .eq. 201007) then ! 7/2010 data

         absfit  = 0.6963342
         rayfac  = 0.8396788E+08
         raysfac = 115395.1
         asyfit  = 0.1135528E-03


      elseif (IWATERP .eq. 201008) then ! 8/2010 data

         absfit  = 0.6953708
         rayfac  = 0.7788437E+08
         raysfac = 137648.2
         asyfit  = 0.1192991E-03


      elseif (IWATERP .eq. 201009) then ! 9/2010 data

         absfit  = 0.7787179
         rayfac  = 0.8227479E+08
         raysfac = 121765.4
         asyfit  = 0.1221935E-03


      elseif (IWATERP .eq. 201010) then ! 10/2010 data

         absfit  = 0.7258545
         rayfac  = 0.8471648E+08
         raysfac = 108818.1
         asyfit  = 0.1250875E-03


      elseif (IWATERP .eq. 201011) then ! 11/2010 data

         absfit  = 0.7017593
         rayfac  = 0.8413286E+08
         raysfac = 110839.6
         asyfit  = 0.1270272E-03


      elseif (IWATERP .eq. 201012) then ! 12/2010 data

         absfit  = 0.7090331
         rayfac  = 0.8650171E+08
         raysfac = 104610.1
         asyfit  = 0.1295903E-03


      elseif (IWATERP .eq. 201101) then ! 1/2011 data

         absfit  = 0.7203000
         rayfac  = 0.8435E+08
         raysfac = 110400.
         asyfit  = 0.1233325E-03


      elseif (IWATERP .eq. 201102) then ! 2/2011 data

         absfit  = 0.7024496
         rayfac  = 0.8571E+08
         raysfac = 105500.
         asyfit  = 0.1357234E-03


      elseif (IWATERP .eq. 201103) then ! 3/2011 data

         absfit  = 0.7210646
         rayfac  = 0.8395+08
         raysfac = 111100.
         asyfit  = 0.1307218E-03

      elseif (IWATERP .eq. 888888) then ! for water transparency dep
         call trabs_sk4(tr,absfitfac)
         absfit  = 0.7142 * absfitfac
         rayfac  = 0.8702E+08   ! average over after 2009 April
         raysfac = 0.1059E+06
         asyfit  = 0.1205E-03   ! Linac period (it should be variable, but may not affect)

*     developed by yokozawa : for day by day absorption parameter.
      elseif (IWATERP .eq. 999999) then ! for day by day param.

         abssfac = 1.164
         abspow  = 10.94
         miefac  = 0.4623E+7
         miesfac = 392.4

         absfit  = 0.6240657
         rayfac  = 0.8513785E+08
         raysfac = 113809.3
         asyfit  = 0.1002079E-03

         if(tr.lt.30000) then
            absfac = -730.1*tr +0.2929E+8

            yabsfac = -.3383244E+08+exp(16.82205+14827.33/tr)

            if( tr .ge. 12500. ) then
               absfac1 = 0.3106475E+07+exp(14.72945+30522.71/tr)
               if (absfac1 .ge. yabsfac ) yabsfac = (yabsfac+2*absfac1)/3.
            elseif( tr .le. 10000. ) then
               absfac2 = -.6063145E+08 + 0.1166306E+13/tr
               if (absfac2 .le. yabsfac ) yabsfac = (yabsfac+2*absfac2)/3.
            endif

         else
            print *,'Worning:water transparency is too large'
            print *,'set absafac to defolt value'
            absfac = 0.2961E+8
            yabsfac = 0.2961E+8
         end if
         print *,'yokozawa check:absfac=',absfac
         print *,'yokozawa for debug: yamada absfac=',yabsfac

C for LINAC 2010 period by T.Yokozawa (laser data during 20100726 - 20100815)

      elseif (IWATERP .eq. 20100726) then ! Linac Period

         absfit  = 0.7142       ! Linac period
         rayfac  = 0.8702E+08   ! average over after 2009 April
         raysfac = 0.1059E+06
         asyfit  = 0.1205E-03   ! Linac period

      elseif (IWATERP .eq. 201908) then ! Aug. 2019 data

c     Yang-san's result on Aug.6 2020
c         absfit  = 0.260922
c         rayfac  = 1.1964E+08
c         raysfac = 54848.2
c         asyfit  = 0.989802E-04

cc     Seo's result on Nov.26 2020
c         absfit  = 0.53109       ! P0
c         rayfac  = 1.12716E+08   ! P4
c         raysfac = 65169.3       ! P5
c         asyfit  = 0.500047E-04  ! P6

c     Seo's result on Jan.14 2021
         absfit  = 0.635126      ! P0
         rayfac  = 1.12716E+08   ! P4
         raysfac = 65169.3       ! P5
         asyfit  = 0.500047E-04  ! P6
         
      else                      ! Apr.2009 is used if IWATERP is not set.

         write(*,*) '***************************************'
         write(*,*) 'TR2FAC: WATERP is not set in card file.'
         write(*,*) 'Water Parameter tuning is set to'
         write(*,*) 'Apr.2009 results (IWATERP=200904)'
         write(*,*) '***************************************'

         absfit  = 0.6240657
         rayfac  = 0.8513785E+08
         raysfac = 113809.3
         asyfit  = 0.1002079E-03

      endif

**    End of Monthly Water Parameter Tuning. 

      return
      end

