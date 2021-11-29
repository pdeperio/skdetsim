************************************************************************
*     -------------------------------------
      SUBROUTINE DSINIP
*     -------------------------------------
*
*     (Purpose)
*       Initialization of parameters for detector Simulator
*
*       Called by : SGINIT
*
*     (Creation Date and Author)
*       1992.05.18 ; First version by K.S.Hirata
*       2003.01.17 ; Merge by Y.Koshio
*       July 2012  : Changes following Rayleigh scattering bug by T. Akiri
*       2013.06.18 ; Gd Effect by P.Weatherly
*
************************************************************************
      implicit none

#include "guwtr.h"
#include "guref.h"
#include "skparm.h"
#include "dsparm.h"
#include "skhead.h"
#include "skdetsim_parm.h"
c MF
#include "reflections.h"

      integer i

      Real  gain_factor
      Parameter (gain_factor=1.02)

C     WATER_TUNE           Water Parameter tuning
C     (INTEGER*4)          = LE_TUNE (= 1)  LE group tuning         
C                          = AP_TUNE (= 2)  ATMPD group tuning 
C     REFLECTION_TUNE      Reflection tuning
C     (INTEGER*4)          = LE_TUNE (= 1)  LE group tuning     
C                          = AP_TUNE (= 2)  ATMPD group tuning 

*     Set PMT qetable version
      cqetab = '3_0'            ! default
      if (SK_GEOMETRY .eq. SK_V) then
         cqetab = '5_2'         ! default
      endif
      
      if( (REFLECTION_TUNE .eq. LE_TUNE) .or.
     &        (REFLECTION_TUNE .eq. UNI_TUNE)) THEN	

         if(IWATERP.eq.201208) then ! modified by Tarek
*
* values updated on october 3, 2007 by M Fechner
* after the black sheet bug fix
* the old values (with bug) are listed for reference 
*
* Vaues updated by Tarek Akiri, July 2012
* Retuning of PMT reflectivity because of bug in Rayleigh
* scattering for non polarized light (affects ref. light)
* Used same data as Max for tuning
* Called by rfpmsg.F and rfbssg.F
*
            if( (SK_GEOMETRY .eq. SK_III) .or. (SK_GEOMETRY .eq. SK_IV) .or. (SK_GEOMETRY .eq. SK_V) )then
               ntuning =4           
               at(1) = 337.
               at(2) = 365.
               at(3) = 400.
               at(4) = 420.
* Black sheet reflectivity
* Old values with black-sheet bug : 0.085,0.090,0.080,0.087
* Old values with Rayleigh bug : 0.071,0.065,0.061,0.069
               tuned_bs(1) = 0.079 
               tuned_bs(2) = 0.079 
               tuned_bs(3) = 0.077 
               tuned_bs(4) = 0.077 
* PMT reflectivity
* Old values with black-shet bug : 0.27,0.37,0.48,0.55
* Old values with Rayleigh bug : 0.29,0.42,0.49,0.53
               tuned_pm(1) = 0.37
               tuned_pm(2) = 0.43
               tuned_pm(3) = 0.53
               tuned_pm(4) = 0.59
               print *,"************************************************************************************"
               print *,"Reflections: setting auto-calib tuning values for PM and BS reflectivity in dsinip.F"
               print *,"************************************************************************************"
            endif
            print *,"# tuning values :" , ntuning
            print *,at
            print *,tuned_bs
            print *,tuned_pm
            print *,"***********************************"         
         endif

*     Set PMT correction efficiency

         if(SK_GEOMETRY.eq.SK_I) then
            COREPMT = 0.913
            REDUCFCT = 1.0
         else if(SK_GEOMETRY.eq.SK_II) then
cc            COREPMT = 0.7623      !  should be tuned by LINAC
cc            COREPMT = 0.8713        !  Tuned by 12,8 MeV nhit 3/11
cc            COREPMT = 0.7925        !  Tuned by 12,8 MeV nhit 3/11
            COREPMT = 0.8431921 !  add phong,lambert scatter 05/04/14
            REDUCFCT = 1.2
         else if(SK_GEOMETRY.eq.SK_III) then
C    2006.4.18 M. Smy set to SK-II for now
C            COREPMT = 0.8431921
C    2007.5.25 B. Yang Tuned by LINAC 8MeV taken May 2007
C            COREPMT = 0.9285
C    2007.10.24 B. Yang  Because of bug in skdetsim, reuned by LINAC 8MeV taken May 2007
C            COREPMT = 0.979707
C    2008.06.23 B. Yang  PMT reflection, QE, etc. was retuned. so retuned by LINAC 8MeV taken Feb 2008
C            COREPMT = 0.88481
C    2008.11.01 B. Yang   LINAC 8MeV taken Feb 2008
C            COREPMT = 0.8826642
C    2009.05.19 Y.Koshio  LINAC 8MeV taken Feb 2008
C            COREPMT = 0.885843
C    2009.09.08 T.Yokozawa  LINAC 8MeV taken Feb 2008
C            COREPMT = 0.887600449
C    2016.03.02 J.Kameda  Additional reduction factor estimated by CRmu (SK4) after Bug Fixes.
            COREPMT = 0.887600449*0.9949

C-- Algorithm of reducing the number of photons is improved and can be applied also LE mode. Byeongsu Yang 2008/04/23
            maxfct = 1.3
            REDUCFCT = maxfct * COREPMT
         else if(SK_GEOMETRY.eq.SK_IV) then
C            COREPMT = 0.8826642 * 1.011   ! Tuned by J.Kameda using through-mu. 2010/04/26
c            COREPMT = 0.8826642 * 1.011 * 0.9925  ! Tuned by J.Kameda using through-mu. 2013/12/17
            COREPMT = 0.8826642 * 1.011 * 0.9925 * 0.9949 ! Tuned by J.Kameda after bug fix of BST collar (sgacredge.F).

            maxfct = 1.3
            REDUCFCT = maxfct * COREPMT

         else if(SK_GEOMETRY.eq.SK_V) then 

c            COREPMT = 0.8826642 * 1.011 * 0.9925 * 0.9949 ! Tuned by J.Kameda after bug fix of BST collar (sgacredge.F).
c            COREPMT = 0.8826642 * 1.011 * 0.9925 * 0.9949 * 0.9966 ! final value is given by thr-going muon by J.K (Mar.2021).
            COREPMT = 0.8826642 * 1.011 * 0.9925 * 0.9949 * 0.9973 ! final value is given by thr-going muon by J.K (Jun.2021).

            maxfct = 2.6  ! change to cover HK PMTs
            REDUCFCT = maxfct * COREPMT

         endif

*     Set PMT reflection parameter

         PMSREF = 1.
C         PMRREF = 0.

*     Set Black Sheet reflection parameter
         if(SK_GEOMETRY.le.2) then
            BSTREF = 1.         ! Parker 05/04/15
         else
            BSTREF = 0.5        ! Kameda 09/05/01
         endif
cc         BSTREF = 1.            ! Parker 05/04/15
cc         if(SK_GEOMETRY.eq.1) then
cc            BSTREF = 1.
cc         else if(SK_GEOMETRY.eq.2) then
cc            BSTREF = 2.7
cc         endif

      else if(REFLECTION_TUNE.EQ.AP_TUNE) THEN

C--------------------

*     Set PMT correction efficiency

         COREPMT = 0.88*gain_factor
#ifdef TUNED_1PE
C
C     eff. is modified for tuning MC of LINAC and decay-e
C     May-6-97  by M.Shiozawa
C
         COREPMT = COREPMT/0.95

C
C     Oct-17-98 Y.Hatakeyama
C
         COREPMT = COREPMT*0.975
C
#endif TUNED_1PE

C-- Algorithm of reducing the number of photons is improved and can be applied also LE mode. Byeongsu Yang 2008/04/23

CCC moved from sgpmt.F and ggckov.F 2008.04.24 Byeongsu Yang
cc  0.75 = 1.0 - 0.25
         if(SK_GEOMETRY.eq.SK_I) then
            REDUCFCT = 1.0
         else if(SK_GEOMETRY.eq.SK_II) then
C increasing # of photons by 20% for the SK-II MC tuning
C                          Shoei NAKAYAMA (15-APR-2003)
            REDUCFCT = 1.2
            COREPMT = COREPMT*0.929
            maxfct = 1.3
         else if(SK_GEOMETRY.eq.SK_III) then

CCC     An amount of photon increase is changed to 1.25 (by G.Mitsuka, Oct 2007)

            REDUCFCT = 1.25
* tuned by cosmic ray muon
CCC            COREPMT = COREPMT*0.929
CCC     SK-III charge scale is tuned by G.Mitsuka(Apr11,2007)
CCC            COREPMT = COREPMT*0.978
CCC     SK-III charge scale is re-tuned after BLST bug-fix by G.Mitsuka(Oct,2007)

            COREPMT = COREPMT*0.989
            maxfct = 1.3

         else if (SK_GEOMETRY.eq.SK_IV) then ! just copy from SK-III
            REDUCFCT = 1.25
            COREPMT = COREPMT*0.989
            maxfct = 1.3
         else if (SK_GEOMETRY.eq.SK_V) then ! just copy from SK-IV for temporal
            REDUCFCT = 1.25
            COREPMT = COREPMT*0.989
            maxfct = 1.3
         else
            write(*,*) '******************************************************'
            write(*,*) 'DSINIP: Not Supported SK_GEOMETRY version:',SK_GEOMETRY
            write(*,*) '******************************************************'
            stop
         endif

*     Set PMT reflection parameter
C    
C    2001-7-24 T.Shibata ; BSTREF changed 
C
         PMSREF = 1.0
C         PMRREF = 0.

*     Set Black Sheet reflection parameter

         BSTREF = 0.5

* Maximilien Fechner, march 2007
* ATMPD reflection tuning -- see reflections.h and doc/ directory
* set wavelength dependant PMT & BS reflection parameters here
* laser injector wavelengths

* values updated on october 3, 2007 by M Fechner
* after the black sheet bug fix
* the old values (with bug) are listed for reference 
         if( (SK_GEOMETRY .eq. SK_III) .or. (SK_GEOMETRY .eq. SK_IV) .or. (SK_GEOMETRY .eq. SK_V) )then
            ntuning =4           
            at(1) = 337.
            at(2) = 365.
            at(3) = 400.
            at(4) = 420.
*     black sheet reflectivity / tuning value
* Old values with black-sheet bug : 0.085,0.090,0.080,0.087
            tuned_bs(1) = 0.071
            tuned_bs(2) = 0.065
            tuned_bs(3) = 0.061
            tuned_bs(4) = 0.069
*     PMT reflectivity
* Old values with black-sheet bug : 0.27,0.37,0.48,0.55
            tuned_pm(1) = 0.29
            tuned_pm(2) = 0.42
            tuned_pm(3) = 0.49
            tuned_pm(4) = 0.53
            print *,"***********************************"
            print *,"Reflections: setting SK-III ATMPD tuning values in dsinip.F"
         else       ! SK-1 or 2
c use previous tuning parameters
            ntuning = 5
            at(1) = 365.
            at(2) = 400.
            at(3) = 420.
            at(4) = 460.
            at(5) = 500.
            tuned_bs(1) = 0.045
            tuned_bs(2) = 0.045
            tuned_bs(3) = 0.055
            tuned_bs(4) = 0.060
            tuned_bs(5) = 0.055
            tuned_pm(1) = 0.467
            tuned_pm(2) = 0.570
            tuned_pm(3) = 0.650
            tuned_pm(4) = 1.133
            tuned_pm(5) = 1.067
            print *,"***********************************"
            print *,"Reflections: setting SK-I / SK-II ATMPD tuning values in dsinip.F"
         endif
            print *,"# tuning values :" , ntuning
            print *,at
            print *,tuned_bs
            print *,tuned_pm
            print *,"***********************************"

      endif

      if ( (WATER_TUNE .eq. LE_TUNE) .or. (WATER_TUNE .eq. UNI_TUNE)
     &	.or. (WATER_TUNE.eq.10) .or. (WATER_TUNE.eq.11) ! related to Gd
     &     ) THEN

*     Water transparency parameter
#ifndef LINAC
         WATRDS = 8300.
c
c move to dscard.F by Y.Koshio 2021/11/06
c         call skrunday          ! fill the common of skday.h
c         if(SK_GEOMETRY.ge.4) then ! fill the common of skwt.h
c            call skwt_gain_corr
c         else
c            call skwt
c         endif
c
#endif
CCC         write(6,*) 'using water transparency of ',watrds

      else if(WATER_TUNE.EQ.AP_TUNE) THEN

*     Set water transparency parameter

         ASCL = 0.
         APAS = 0.0
         APSS = 0.0

C    
C    July-3-01 T.Shibata ; AABL,ARAS,AMIS changed 
C      
         AABL = 1.0
         IF(SK_GEOMETRY.eq.SK_I) then
            ARAS = 0.72
            AMIS = 0.8
         ELSE IF(SK_GEOMETRY.eq.SK_II) then
C
C    May-11-03 HKSeo (S. Moriyama committed) ; ARAS, AMIS changed by -10%
C     
            ARAS = 0.65
            AMIS = 0.72
         ELSE IF(SK_GEOMETRY.eq.SK_III) then
C
C    2006.4.18 M. Smy set to SK-II for now
C
C    Apr-04-07 B. Yang ; The expression of absorption and scattering changed.
C              refer to wtrsg.F, wsymscsg.F, wasyscsg.F, wabssg.F
C
C            ARAS = 1.0
C            AMIS = 0.1683
C
C    Sep-28-07 B. Yang ; SK-III water parameter and coefficients expression retuned
C               because of MC bug. So symmetric scattering and absorption coefficient
C               function changed again.  Please refer to wsymscsg.F and wabssg.F
C
            ARAS = 1.025
            AMIS = 0.1362

         ELSE IF(SK_GEOMETRY.eq.SK_IV) then  ! same as SK_III
            ARAS = 1.025
            AMIS = 0.1362
         ELSE IF(SK_GEOMETRY.eq.SK_V) then   ! just copy from SK_IV
            ARAS = 1.025
            AMIS = 0.1362
         else
            write(*,*) '******************************************************'
            write(*,*) 'DSINIP: Not Supported SK_GEOMETRY version:',SK_GEOMETRY
            write(*,*) '******************************************************'
            stop
         ENDIF

      endif

*     Electornics parameter
*     If you set them in card, the number in card are taken priority.

C     DS_TUNE              Tuning in DS routines
C     (INTEGER*4)          = LE_TUNE (= 1)  LE group tuning         
C                          = AP_TUNE (= 2)  ATMPD group tuning 
C meaning less
C      if(DS_TUNE.EQ.LE_TUNE) THEN
C         SGGATE = 600.
C      else if(DS_TUNE.EQ.AP_TUNE) THEN
C         SGGATE = 300.
C      endif

      TCNVSK = 2.5
      QCNVSK = 5.5

      TRESDS = 3.
      QRESDS = 0.85

      GATEDS = 1.3
      GLTRDS = 800.

      IF(SK_GEOMETRY .ge. SK_IV) then ! QBEE
        SIGWDS = 400.
      elseif (sk_geometry .le. SK_III) then
        SIGWDS = 300.
      else
        write(*,*) 'Not supported SK_GEOMETRY. Stop.', SK_GEOMETRY
        stop
      ENDIF

      SGGATE = SIGWDS

      DEADDS = 900.

      TWINDS = 200.

      IF(SK_GEOMETRY .LE. SK_III) THEN 
         IF(TRACK_PHOTON.EQ.TRACK_ALL_PHOTON) THEN
CCC dsnois.F
            OFFSET0 = 540.
         else if ( (TRACK_PHOTON.EQ.TRACK_QE_PHOTON) .OR.
     &           (TRACK_PHOTON.EQ.TRACK_QE_OVER_ABS_SAVE_MEMORY)) THEN
CCC dshigh.F
            OFFSET0 = 960.
         ELSEIF (TRACK_PHOTON.EQ.TRACK_QE_OVER_ABS_PHOTON) THEN
CCC dstrig_multi.F
            OFFSET0 = 960.
            OFFSET0 = 840.
            TRGDS = 260.
            DTRGOFF = 260.
         ENDIF
      Elseif ((SK_GEOMETRY.GE.SK_IV) .AND. 
     &        ((TRACK_PHOTON.EQ.TRACK_QE_PHOTON) .OR.
     &        (TRACK_PHOTON.EQ.TRACK_QE_OVER_ABS_SAVE_MEMORY))
     &        ) THEN
         OFFSET0 = 960.
      ENDIF

      write(*,*) 'DSINIP::::COREPMT is set to :', COREPMT
      RETURN
      END
