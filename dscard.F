***********************************************************************
*     -------------------------------------
      SUBROUTINE DSCARD
*     -------------------------------------
*
*     (Purpose)
*       Read DS card and setup common blocks.
*
*       Called by : SGINIT
*
*     (Input)
*
*     (Creation Date and Author)
*       1995.02.01 ; First version by Koshio
*
************************************************************************

      implicit none

#include "skhead.h"
#include "skparm.h"
#include "dsparm.h"
#include "dsparma.h"
#include "dsbnkf.h"
#include "guflag.h"
#include "guref.h"

#ifdef FADC
#include "fadcwrbnk.h"
#endif

      CHARACTER*40 TMPFIL, COMMAND
      INTEGER I,IRET, j
      INTEGER LRUNW, leladay, IQETAB
      REAL WATRCD, COREPMTC

      integer LRUNTBA
      real TBACD

** for gain simulation
      integer lenchr
      CHARACTER*6 tmpchar

** for TBA with gain 2020/03/26 M.Shinoki
      real slope_ms, intercept_ms

** for run by run trigger threshold
#include "skruninf.h"
      integer nrunsk_for_trigthr

#include "skday.h"

#ifdef LINAC
#include "lininfo2.h"
#endif

#include "skdetsim_parm.h"

C --- initialize parameters
      IFTRIG = 0
      IFSWTRG = 0
      IT0_OFFSET = 0
      IMWATR = 0
      WATRCD = 8300
      LRUNW = -9999
      IWATERP = -9999
      GD_SCAT_FRAC = 1.
      IMTBA = 0
      TBACD = 0.
      LRUNTBA = -9999

      COREPMTC = COREPMT
      IQETAB = 0

C     PMT-by-PMT dark-rate toggle
      DOPMTDARKRT = 0
      PMTNOISEFILE='./DarkRates_sk4.dat' ! need this for initialize
      
      DARKDS = 3.
      DTHRSK = 0.2
      TRIGDS = 30.

      IMCPHI = .TRUE.
      IMCRI  = .TRUE.
      IMCRAI = .TRUE.
      IMCTQI = .TRUE.
      ITQI   = .TRUE.
      ITQAI  = .TRUE.
      ITQREALI  = .TRUE.
      ITQAREALI  = .TRUE.
      ITQZI  = .TRUE.

C      FULWIDA = 32.*1e3
      FULWIDA = 7.*1e3
C      DSANTT0 = 16.*1e3
      DSANTT0 = 4*1e3
      DARKDSA(1) = 5.0
      DARKDSA(2) = 5.0
      DARKDSA(3) = 5.0
      IAGEQAL = 1
C-Y.H      DSTYVKR = 1.
      DSTYVKR = 0.9
c-J.G.

      DSSEGRTOP = 1.
      DSSEGRBOT = 1.

      DSTYVKRTOP(1) = 0.68
      DSTYVKRBOT(1) = 0.92
      DSTYVKRWAL(1) = 0.9
      ! tyvek reflectivity on the outer walls
      DSTYVKRTOP(2) = -1.0
      DSTYVKRBOT(2) = -1.0
      DSTYVKRWAL(2) = -1.0
      DSTYVKM = 0.

c-F.D., K.S.
      PROBTT = 0.
      PROBTB = 0.
      PROBODSST = 0.
      PROBODSSB = 0.

C-Y.H      APMTCOR = 1.
      APMTCOR = 0.22
      APMTEFR = 1.
      APMTEFT = 1.
      APMTEFB = 1.
      APMTEFW = 1.
C      AWLSCOR = 1.
      AWLSCOR = 0.007
      ADSCTHR = 0.25
      ADRANGE = 50.
      ADSAPO(1) = 1.
      ADSAPO(2) = 0.002
      ADSAPO(3) = 0.0003
      ADSAPN(1) = 1.
      ADSAPN(2) = 0.0007
      ADSAPN(3) = 0.0003
      ADSTHRF = 0.01

     
! --- Plateau + Gaussian response paraters for old OD Tubes
!     20080213 rvw 
      GAIN_PLATEAU = 0.6
      GAINFACT     = 0.9
 
C-- by Baran Bodur 201910
      ODDARKRUN = 0
C      ODINDSPE = 0
      ODINDSPE = .false.
      ODPEDSMEAR = 0
      ODIMBSMEAR = 0

! --- Default setting is to count photons inside of sgpmta
      OD_PHOTON = 0

! --- Default not simulating the after pulse       
      ADD_AFTER_PULSE = 0

C --- read card

      CALL GETUNI(TMPFIL)
      write (*,*) 'tmpfil = ',TMPFIL
      write (*,*) 'fcard  = ',fcard
      OPEN(UNIT=90, FILE=FCARD, STATUS='OLD')
      OPEN(UNIT=91, FILE=TMPFIL, STATUS='UNKNOWN')
      CALL EXCARD('DS', 90, 91, IRET)
      CLOSE(UNIT=91)
      CLOSE(UNIT=90)

      OPEN(UNIT=94, FILE=TMPFIL, STATUS='OLD')
#ifdef absoft
C     Ugly work-around of problem passing LUN from Absoft to f2c or g77.
      call system('/bin/ln -sf '//TMPFIL//' fort.94')
#endif
      CALL FFINIT(0)
      CALL FFSET('LINP', 94)
      CALL FFSET('SIZE', 9)
C for water
      CALL FFKEY('MWAT',   IMWATR, 1, 'INTEGER')
      CALL FFKEY('WATRA',  WATRCD, 1, 'REAL')
      CALL FFKEY('WATRUN', LRUNW,  1, 'INTEGER')

C for SK-IV tuning (20090806 Shiozawa)
      CALL FFKEY('WATERP',  IWATERP, 1, 'INTEGER')

C for top-bottom asymmetry
      CALL FFKEY('MTBA',   IMTBA,   1, 'INTEGER')
      CALL FFKEY('TBA' ,   TBACD,   1, 'REAL')
      CALL FFKEY('TBARUN', LRUNTBA, 1, 'INTEGER')

C for PMT QE
      CALL FFKEY('NORMQE', COREPMTC, 1, 'REAL')
      CALL FFKEY('QETAB', IQETAB, 1, 'INTEGER')

C for PMT gain (Added on Apr. 24, by Y. Nakajima)
      nrunsk4_for_gain = 0
      CALL FFKEY('GAINRUN', nrunsk4_for_gain, 1, 'INTEGER')
      
C for trigger threshold (Added on Nov. 3, 2021, by Y. Koshio)
c     Modified by Y. Nakano Nov 3, 2021.
      nrunsk_for_trigthr = 0
      CALL FFKEY('TRIGRUN', nrunsk_for_trigthr, 1, 'INTEGER')
      if(SK_GEOMETRY.ge.SK_IV) then
         if(nrunsk_for_trigthr.ne.0) then
            call softtrg_get_cond(softtrg_detector, softtrg_thr, 
     $           softtrg_t0_offset, softtrg_pre_t0, softtrg_post_t0)
            call get_run_softtrg(nrunsk_for_trigthr,
     $           softtrg_detector, softtrg_thr,
     $           softtrg_t0_offset, softtrg_pre_t0, softtrg_post_t0)
            call softtrg_set_cond(softtrg_detector, softtrg_thr, 
     $           softtrg_t0_offset, softtrg_pre_t0, softtrg_post_t0)
         endif
      endif

C for simulation
      CALL FFKEY('ANTI',   IFANTI, 1, 'INTEGER')
      CALL FFKEY('TRIGSIM',IFTRIG, 1, 'INTEGER')
      CALL FFKEY('NONSWTRG',IFSWTRG, 1, 'INTEGER')
      CALL FFKEY('IT0OFF',IT0_OFFSET, 1, 'INTEGER')

C for Gadolinium Scattering Fraction -Pi
      CALL FFKEY('GDSCATFRAC',GD_SCAT_FRAC, 1, 'REAL')

C for ID
      CALL FFKEY('DOPMTDARKRT',   DOPMTDARKRT, 1, 'INTEGER')
      CALL FFKEY('PMTNOISEFILE',PMTNOISEFILE,256,'MIXED')
      CALL FFKEY('PMTNOISERUN',PMTNOISERUN,1,'INTEGER')
      CALL FFKEY('DARK',   DARKDS, 1, 'REAL')
      CALL FFKEY('DISTHR', DTHRSK, 1, 'REAL')
      CALL FFKEY('TRGTHR', TRIGDS, 1, 'REAL')
      CALL FFKEY('GEGATE', SGGATE, 1, 'REAL')
      CALL FFKEY('ATMVET', SIGWDS, 2, 'REAL')
      CALL FFKEY('GLTTIM', GLTRDS, 1, 'REAL')
      CALL FFKEY('OFFSET', OFFSET0, 1, 'REAL')
      CALL FFKEY('TRGDS',  TRGDS, 1, 'REAL')
      CALL FFKEY('TRGOFF', DTRGOFF, 1, 'REAL')
c      CALL FFKEY('TIMGATE',GATEDS, 1, 'REAL')
c      CALL FFKEY('TWINDW', TWINDS, 1, 'REAL')


CCCC for OD
      CALL FFKEY('FULWIDA', FULWIDA , 1 , 'REAL')
      CALL FFKEY('DSANTT0', DSANTT0 , 1 , 'REAL')
      CALL FFKEY('DARKDSA', DARKDSA , 3 , 'REAL')
      CALL FFKEY('IAGEQAL', IAGEQAL , 1 , 'INTEGER')
CCC(obsolote)      CALL FFKEY('DSTYVKR', DSTYVKR , 1 , 'REAL')
      CALL FFKEY('DSTYVKM', DSTYVKM , 1 , 'REAL')
      CALL FFKEY('DSTYVKRTOP', DSTYVKRTOP , 2 , 'REAL')
      CALL FFKEY('DSTYVKRBOT', DSTYVKRBOT , 2 , 'REAL')
      CALL FFKEY('DSTYVKRWAL', DSTYVKRWAL , 2 , 'REAL')

c - Factor for reflection off segmentation Tyvek
      CALL FFKEY('DSSEGRTOP', DSSEGRTOP , 1 , 'REAL')
      CALL FFKEY('DSSEGRBOT', DSSEGRBOT , 1 , 'REAL')
c - add for transmission      
      CALL FFKEY('PROBTT', PROBTT , 1 , 'REAL')
      CALL FFKEY('PROBTB', PROBTB , 1 , 'REAL')
      CALL FFKEY('PROBODSST', PROBODSST , 1 , 'REAL')
      CALL FFKEY('PROBODSSB', PROBODSSB , 1 , 'REAL')
c - ANTI detector PMT's factors      
      CALL FFKEY('APMTCOR', APMTCOR , 1 , 'REAL')
      CALL FFKEY('APMTEFR', APMTEFR , 1 , 'REAL')
      CALL FFKEY('APMTEFT', APMTEFT , 1 , 'REAL')
      CALL FFKEY('APMTEFB', APMTEFB , 1 , 'REAL')
      CALL FFKEY('APMTEFW', APMTEFW , 1 , 'REAL')
c - wave length shifter
      CALL FFKEY('AWLSCOR', AWLSCOR , 1 , 'REAL')
c - discriminator
      CALL FFKEY('ADSCTHR', ADSCTHR , 1 , 'REAL')
c - Q dynamic range
CCC(obsolete)      CALL FFKEY('ADRANGE', ADRANGE , 1 , 'REAL')
      CALL FFKEY('ADSAPO', ADSAPO , 3 , 'REAL')
      CALL FFKEY('ADSAPN', ADSAPN , 3 , 'REAL')
c - PMT threshold
      CALL FFKEY('ADSTHRF', ADSTHRF , 1 ,'REAL')
c - gain
      CALL FFKEY('GAINPLATEAU', GAIN_PLATEAU , 1 , 'REAL')
      CALL FFKEY('GAINFACT'   , GAINFACT     , 1 , 'REAL')
c - More for OD by Baran B. in 201910 for sk5
      call FFKEY('ODDARKRUN',ODDARKRUN,1,'INTEGER')
      call FFKEY('ODINDSPE',ODINDSPE,1,'LOGICAL')
      call FFKEY('ODPEDSMEAR',ODPEDSMEAR,1,'REAL')
      call FFKEY('ODIMBSMEAR',ODIMBSMEAR,1,'REAL')
      

! photon counting method for OD 
      CALL FFKEY('ODPHOTON'   , OD_PHOTON    , 1 , 'INTEGER')

c - set timing window of saving hit information after SK-IV
      CALL FFKEY('SHORTWIN',IFSHORTWIN, 1, 'INTEGER')

c - flag to simulate after-pulse
      CALL FFKEY('AFTERPULSE',ADD_AFTER_PULSE, 1, 'INTEGER')

c - writing bank
C      CALL FFKEY('WRTBANK',IMCPHI, 6, 'LOGI')
C      CALL FFKEY('WRTBANK',IMCPHI, 7, 'LOGI')
      CALL FFKEY('WRTBANK',IMCPHI, 9, 'LOGI')
#ifdef FADC
      CALL FFKEY('WRTFADCBNK',IFADCI, 1, 'LOGI')
#endif
      CALL FFGO
      CLOSE(UNIT=94)
#ifdef absoft
C     Ugly work-around of problem passing LUN from Absoft to f2c or g77.
      call system('/bin/rm fort.94')
#endif
      COMMAND = '/bin/rm '//TMPFIL
      CALL SYSTEM(COMMAND)

      write(*,*) 'Water parameters are set for one period (SK-IV tuning)'
      write(*,*) 'IWATERP = ',IWATERP

C --- Conversion for LE_and_UNI_WATER mode

#ifndef LINAC
c move from dsinip.F by Y.Koshio 2021/11/06
      call skrunday             ! fill the common of skday.h
      if(SK_GEOMETRY.ge.4) then ! fill the common of skwt.h
         if(SK_GEOMETRY.eq.4) then
            if(IMTBA.eq.2) then 
               call skwt
            else
               call skwt_gain_corr
            endif
         else if(SK_GEOMETRY.eq.5) then
            call skwt
         else
            write(6,*) 'Not supported yet SK_GEOMETRY : ',SK_GEOMETRY
         endif
      endif
#endif

C      IF (WATER_TUNE .EQ. LE_TUNE) THEN
      IF (WATER_TUNE .NE. AP_TUNE) THEN
         IF(IMWATR.EQ.1) THEN
           IF(WATRCD .GT. 0) THEN
             WATRDS = WATRCD
           ELSE
             WRITE(0,*)'DSCARD : Error!! Bad DS-WATRA ', WATRCD
             STOP 1
           ENDIF
         ELSEIF(IMWATR.EQ.2) THEN
           IF(LRUNW .GT. 0) THEN
             leladay = relapse(LRUNW) ! get elapse day from 96/1/1
             call lfwater(leladay,WATRDS) ! calculate water tranceparancy
           ELSE
             WRITE(0,*)'DSCARD : Error!! Bad DS-WATRUN ', LRUNW
             STOP 1
           ENDIF
         ELSEIF(IMWATR.NE.0.AND.IMWATR.NE.-1) THEN
           WRITE(0,*)'DSCARD : Error!! Bad DS-MWAT ', IMWATR
           STOP 1
         ENDIF
         call tr2fac(WATRDS)
      ENDIF

*** for top-bottom asym.
*     add day by day tbads for SK-IV by Yokozawa
      IF(IMTBA.EQ.1) THEN
         TBADS = TBACD * 0.00001
      ELSE IF(IMTBA.EQ.2) THEN
         IF(LRUNTBA .ge. 30000.and.LRUNTBA.lt.60000.and.SK_GEOMETRY.EQ.SK_III) THEN
            leladay = relapse(LRUNTBA) ! get elapse day from 96/1/1
            call lftba(leladay,TBADS) ! calculate top bottom asym for SK=III
            TBADS = -2.04035 * (TBADS*100) * 0.00001
         ELSE IF(LRUNTBA.ge.60000.and.LRUNTBA.le.77958.and.SK_GEOMETRY.EQ.SK_IV) then
            leladay = relapse(LRUNTBA) ! get elapse day from 96/1/1
            call lftba(leladay,TBADS) ! calculate top bottom asym for SK=IV
c     modify atode
            TBADS = (-0.163*(TBADS*100)*(TBADS*100) - 3.676*(TBADS*100) )*0.00001

c            TBADS = -2.04035 * (TBADS*100) * 0.00001
         else if(LRUNTBA.gt.77958.and.SK_GEOMETRY.EQ.SK_IV) then
            print *,'time varience for SK-IV is supported until 2018/5/30'
            print *,'temporary, input tba param of Apr.2009'
            TBADS = 0.0001003
            STOP 1              ! is this required ? by Y. Nakano
cccccccFor SK-V (2021/11/03 M.Shinoki)
         else if(LRUNTBA.gt.80000.and.LRUNTBA.le.85985.and.SK_GEOMETRY.GE.SK_V) then
            leladay = relapse(LRUNTBA) ! get elapse day from 96/1/1
            call lftba(leladay,TBADS) ! calculate top bottom asym for SK-V/SK-VI
            WRITE(6,*)'Check: Elapsed day = ', leladay
            WRITE(6,*)'Check: Official TBA value = ', TBADS
c            TBADS = (-5.584 * (TBADS*100) - 11.31) * 0.00001 ! Fixed Y-int watert
            TBADS = (-5.263 * (TBADS*100) - 10.92) * 0.00001 ! Free Y-int watert
         else if(LRUNTBA.gt.85985)then
             print *,'time varience for SK-V/VI is supported until R#85985'
             print *,'Please upload TBA file!'
             print*, 'Ask expert, Y. Koshio, Y. Nakano, L. Luis'
             STOP 1
         ELSE
            WRITE(0,*)'DSCARD : Error!! Bad DS-TBARUN ', LRUNTBA
            WRITE(0,*)'Support after SK-III'
            STOP 1
         ENDIF
      ELSE IF(IMTBA.EQ.3) THEN  ! gain correction only for SK-IV
         IF(LRUNTBA.ge.60000.and.LRUNTBA.le.77958.and.SK_GEOMETRY.EQ.SK_IV) then
            leladay = relapse(LRUNTBA) ! get elapse day from 96/1/1
            call lftba_gain_corr(leladay,TBADS) ! calculate top bottom asym for SK=IV
            WRITE(6,*)'Check: Elapsed day = ', leladay ! 2020/03/17 M.Shinoki
            WRITE(6,*)'Check: Official TBA value = ', TBADS ! 2020/03/17 M.Shinoki
c     modify atode
c            TBADS = (-0.163*(TBADS*100)*(TBADS*100) - 3.676*(TBADS*100) )*0.00001
c            TBADS = -2.04035 * (TBADS*100) * 0.00001
c            slope_ms = (-1.1175e-10)*WATRDS*WATRDS + (2.35166e-6)*WATRDS - 0.0164681 ! 2020/03/26 M.Shinoki
c            slope_ms = (-9.88564e-11)*WATRDS*WATRDS + (2.14864e-6)*WATRDS - 0.0163051 ! 2020/05/09 M.Shinoki
            slope_ms = (1.04159e-10)*WATRDS*WATRDS - (3.30078e-6)*WATRDS + 0.0203749 ! 2020/05/17 M.Shinoki
c            intercept_ms = 0.953854*slope_ms + 0.00195722 ! 2020/03/26 M.Shinoki
c            intercept_ms = -0.332257*slope_ms - 0.0110225 ! 2020/05/09 M.Shinoki
            intercept_ms = 2.42604*slope_ms + 0.003059 ! 2020/05/17 M.Shinoki
            TBADS =(slope_ms*(TBADS*100)+intercept_ms)*0.01 ! 2020/03/26 M.Shinoki
         else if(LRUNTBA.gt.77958.and.SK_GEOMETRY.EQ.SK_IV) then
            print *,'time varience for SK-IV is supported until R#77958'
            print *,'Please upload TBA file!'
            print *,'temporary, input tba param of Apr.2009'
            TBADS = 0.0001003
            STOP 1
         ELSE
            WRITE(0,*)'DSCARD : Error!! Bad DS-TBARUN ', LRUNTBA
            WRITE(0,*)'Support (gain corr tba) after SK-IV'
            STOP 1
         ENDIF
      ELSEIF(IMTBA.NE.0.AND.IMTBA.NE.-1) THEN
         WRITE(0,*)'DSCARD : Error!! Bad DS-MTBA ', IMTBA
         STOP 1
      ENDIF
      WRITE(6,*)''
      WRITE(6,*)'*******************************************************'
      WRITE(6,*)' The top-bottom asym. param. is set to', TBADS,'.'
      WRITE(6,*)'*******************************************************'
      WRITE(6,*)''

      WRITE(6,*)''
      WRITE(6,*)'******************************************************************************************'

c      if(IQETAB.eq.2) cqetab='4_1'
      if(IQETAB.eq.2) cqetab='4_3' ! changed in 2020 June
      if(IQETAB.eq.50) cqetab='5_0'

      if ( COREPMT .NE. COREPMTC ) then
         WRITE(6,*)' The default absolute normalizatino factor for PMT QE is ', COREPMT,'.'
         COREPMT = COREPMTC
         if(SK_GEOMETRY.ge.SK_III) then
           if(REFLECTION_TUNE.NE.AP_TUNE) THEN       ! LE_TUNE or UNI_TUNE
               REDUCFCT = maxfct * COREPMT
           else if(TRACK_PHOTON.ne.TRACK_QE_PHOTON) then
               REDUCFCT = maxfct * COREPMT
           endif
         endif
      endif
      WRITE(6,*)' The absolute normalizatino factor for PMT QE is set to ', COREPMT,'.'
      WRITE(6,*)'******************************************************************************************'
      WRITE(6,*)''


      IF(IFANTI.EQ.0) THEN
        WRITE(6,*)'OD OFF.'
      ELSEIF(IFANTI.EQ.1) THEN
        WRITE(6,*)'OD ON.'
      ELSE
        WRITE(0,*)'dscard: Error, DS-ANTI ',IFANTI,'.'
      ENDIF

      IF(IFTRIG.EQ.0) THEN
        WRITE(6,*)'Trigger simualation OFF.'
      ELSEIF(IFTRIG.EQ.1) THEN
        WRITE(6,*)'Trigger simualation ON.'
      ELSE
        WRITE(0,*)'dscard: Error, DS-IFTRIG ',IFTRIG,'.'
      ENDIF
      
      IF(IFSWTRG.EQ.0) THEN
        WRITE(6,*)'SW Trigger simualation ON'
      ELSEIF(IFSWTRG.EQ.1) THEN
        WRITE(6,*)'SW Trigger simualation OFF. IT0 offset is ',IT0_OFFSET
      ELSE
        WRITE(0,*)'dscard: Error, DS-NONSWTRG ',IFSWTRG,'.'
      ENDIF
      
      WRITE(6,*)
      WRITE(6,*) '*** Trigger simulation parameter '
      WRITE(6,*)

      WRITE(6,600) SGGATE
 600  FORMAT(' Gate width for high energy simulation was set to '
     &     ,F10.3,' (nsec).')

#ifdef LINAC
***  Dark rate should be set each run..
      if(DARKDS.lt.0) then
         if(lin_dark.eq.-1) then
            write(6,*) 'Dark rate cannot be set in LINAC run...'
            write(6,*) 'Check Run number and'
            write(6,*) '/home/sklowe/linac/const/linac_sk(SK_GEOMETRY)_runsum.dat'
            stop 1
         else

            if(SK_GEOMETRY.eq.SK_I) then
              darkfac = 1.0
            elseif(SK_GEOMETRY.eq.SK_II) then
              darkfac = 1.201
            elseif(SK_GEOMETRY.eq.SK_III) then
              darkfac = 1.286
            elseif(SK_GEOMETRY.eq.SK_IV) then
              darkfac = 1.269
            elseif(SK_GEOMETRY.eq.SK_V) then
c  temporal value 
              darkfac = 1.269 
            else
               write(*,*) 'SK_GEOMETRY not supported number ', sk_geometry
               stop
            endif

            darkds = lin_dark * darkfac
            write(6,*) 'Apply darkfac = ',darkfac
         endif
      endif
#endif

      WRITE(6,610) DARKDS
 610  FORMAT(' Uniform dark noise rate was set to ',F10.3,' (kHz).')
      DARKDS = DARKDS*1.e3

C     PMT-by-PMT dark noise rate for SK-IV
      if (DOPMTDARKRT.ne.0) then
         if(SK_GEOMETRY.ge.SK_IV) then
            if(DOPMTDARKRT.eq.1) then
               write(6,*) 'Using individual PMT-by-PMT dark rate for SK-IV from file:'
               write(6,*) PMTNOISEFILE

            else if(DOPMTDARKRT.eq.2) then
               write(tmpchar,'(I6.6)')  PMTNOISERUN
               PMTNOISEFILE='/home/sklowe/const/darkr_for_skdetsim/output/darkr.'//tmpchar//'.txt'
               print *,'================================='
               print *,PMTNOISERUN
               print *,'================================='
            else
               write(6,*) 'Not support darkrate option ',DOPMTDARKRT
               stop
            endif

            write(6,*) 'Using individual PMT-by-PMT dark rate for SK-IV from file:'
            write(6,*) PMTNOISEFILE(1:lenchr(PMTNOISEFILE))
            open (unit=83,file=PMTNOISEFILE,status="old")
            do i = 1, maxpm
               read (83,*) j, darkrt(i)
               darkrt(i) = darkrt(i)*1000
            end do
            close (83)

         else
            write(6,*) 'Individual PMT-by-PMT dark rate not prepared for non-SK-IV.'
            write(6,*) 'Using uniform dark rate.'
            DOPMTDARKRT=0
            do i = 1, maxpm
               darkrt(i) = DARKDS
            end do
         endif
      else
         do i = 1, maxpm
            darkrt(i) = DARKDS
         end do
      endif

      WRITE(6,620) DTHRSK
 620  FORMAT(' Discri. threshold was set to ',F10.3,' (p.e.).')

      WRITE(6,630) INT(TRIGDS)
 630  FORMAT(' Trigger threshold was set to ',I5,' (hit).')

      WRITE(6,640) SIGWDS,DEADDS
 640  FORMAT(' ATM veto was set to ',F10.3,' ->',F10.3,' (nsec)')

      WRITE(6,650) GATEDS
 650  FORMAT(' Dynamic range (unchangeable) was set to '
     &     ,F10.3,' (micro sec).')
      GATEDS = GATEDS*1.e3

      WRITE(6,660) GLTRDS
 660  FORMAT(' the time of global trigger was set to ',F10.3,' (nsec).')

      WRITE(6,670) TWINDS
 670  FORMAT(' Hitsum pulse window (unchangeable) was set to '
     &     ,F10.3,' (nsec).')

      WRITE(6,*) ' the default place of geant T0 was set to ',OFFSET0,' (nsec).'

      IF(TRACK_PHOTON.EQ.TRACK_QE_OVER_ABS_PHOTON) THEN
        WRITE(6,*) 'TRACK_QE_OVER_ABS_PHOTON and SK-I to III modes need additional parameters'
        WRITE(6,*) ' Geant trigger timing in no trigger simulation case was set to',TRGDS ,' (nsec).'
        WRITE(6,*) ' difference of trigger timing and t0 was set to ',DTRGOFF,' (nsec).'
      ENDIF

CCC for OD
      write(6,680) IAGEQAL
 680  FORMAT(' Quality of ANTI detector Simulation(1:Quick/2:Exact) ',
     $       I7)

      write(6,681) FULWIDA
 681  FORMAT(' Anti detector TDC FULL width  ',F10.3,' (nsec).')

      write(6,682) DSANTT0
 682  FORMAT(' Anti detector TDC T0 offset   ',F10.3,' (nsec).')

      write(6,683) (DARKDSA(i)*1.e3, i=1,3)
 683  FORMAT(' Anti detector Dark Rate   ',3F10.3,' (KHz).')      
      DARKDSA(1) = DARKDSA(1)*1.e3
      DARKDSA(2) = DARKDSA(2)*1.e3
      DARKDSA(3) = DARKDSA(3)*1.e3

CCC(obsolete)      write(6,684) DSTYVKR
CCC(obsolete) 684  FORMAT(' Relative ratio of Tyvek reflection  ',F10.3)

      if( dstyvkrtop(2) .lt. 0. ) then
          dstyvkrtop(2) = dstyvkrtop(1)
      endif
      if( dstyvkrbot(2) .lt. 0. ) then 
          dstyvkrbot(2) = dstyvkrbot(1)
      endif
      if( dstyvkrwal(2) .lt. 0. ) then 
          dstyvkrwal(2) = dstyvkrwal(1)
      endif
 
      write(6,284) DSTYVKRTOP
 284  FORMAT(' Relative ratio of Tyvek reflection TOP  ',2F10.3)

      write(6,285) DSTYVKRBOT
 285  FORMAT(' Relative ratio of Tyvek reflection BOT ',2F10.3)

      write(6,286) DSTYVKRWAL 
 286  FORMAT(' Relative ratio of Tyvek reflection WAL ',2F10.3)


      write(6,281) DSSEGRTOP
 281  FORMAT(' Relative ratio of segmentation Tyvek reflection TOP  ',F10.3)

      write(6,282) DSSEGRBOT
 282  FORMAT(' Relative ratio of segmentation Tyvek reflection BOT ',F10.3)



c      if(DSTYVKRTOP(1).ne.0.80.or. DSTYVKRBOT(1).ne.0.80 
c     $.or.DSTYVKRWAL(1).ne.0.80)then
c         write(6,*)' **************************************************'
c         write(6,*)'   '
c         write(6,*)' WARNING ! WARNING ! WARNING ! WARNING ! WARNING ! '
c         write(6,*)'   >>>>> You are using unofficial TYVEK parameters '
c         write(6,*)'   '
c         write(6,*)' **************************************************'
c      endif

      write(6,287) PROBTT
 287  FORMAT(' TRANSMISSION PROBABILITY, TOP ',F10.3)
      write(6,288) PROBTB
 288  FORMAT(' TRANSMISSION PROBABILITY, BOTTOM ',F10.3)

      write(6,685) APMTCOR
 685  FORMAT(' Global Anti detector PMT collection efficiency ',F10.3)
      

      write(6,687) AWLSCOR 
 687  FORMAT(' Global Anti detector WLS collection efficiency ',F10.3)


      write(6,694) APMTEFR
 694  FORMAT(' Efficiency factor for old tubes',F10.3)

      write(6,695) APMTEFT
 695  FORMAT(' Efficiency factor for top tubes',F10.3)

      write(6,696) APMTEFB
 696  FORMAT(' Efficiency factor for bottom tubes',F10.3)

      write(6,697) APMTEFW
 697  FORMAT(' Efficiency factor for wall tubes',F10.3)


      write(6,688) ADSCTHR 
 688  FORMAT(' Global Anti detector discriminator threshold ',F10.3,
     $       '(p.e.)')

      write(6,690) ADSTHRF
 690  FORMAT('Anti detector discriminator threshold factor ',F10.3)

CCC      write(6,689) ADRANGE 
CCC 689  FORMAT(' Anti detector dynamic range, new tubes ',F10.3,
CCC     $       '(p.e.)')

      write(6,692) ADSAPN 
 692  FORMAT(' Anti detector saturation parameters, new tubes ',3F13.7
     $      )

      write(6,693) ADSAPO
 693  FORMAT(' Anti detector saturation parameters, old tubes ',3F13.7
     $       )

      write(6,698) GAIN_PLATEAU
 698  FORMAT(' Plateau portion of response for old od-tubes ',F10.3)

      write(6,699) GAINFACT
 699  FORMAT(' Gaussian portion of response for old od-tubes ',F10.3)

      write(6,700) OD_PHOTON
 700  FORMAT(' OD_PHOTON counting  ',I5)

C -- bank control

      WRITE(6,*)
      WRITE(6,*) '*** Bank control '
      WRITE(6,*)

      DO 10 I=1,2
         NBNKDS(I)=1
 10   CONTINUE

      IF (IMCPHI) THEN
         NBNKDS(1) = 1
         WRITE(6,701)
 701     FORMAT(' Writing MCPHOTON bank : ON')
      ELSE
         NBNKDS(1) = 0
         WRITE(6,702)
 702     FORMAT(' Writing MCPHOTON bank : OFF')
      END IF
            
      IF (IMCRI) THEN
         NBNKDS(2)=1
         WRITE(6,703)
 703     FORMAT(' Writing MCRAWTQ bank : ON')
      ELSE
         NBNKDS(2)=0
         WRITE(6,704)
 704     FORMAT(' Writing MCRAWTQ bank : OFF')
      END IF

      IF (IMCRAI) THEN
         NBNKDS(3)=1
         WRITE(6,705)
 705     FORMAT(' Writing MCRAWTQA bank : ON')
      ELSE
         NBNKDS(3)=0
         WRITE(6,706)
 706     FORMAT(' Writing MCRAWTQA bank : OFF')
      END IF

      IF (IMCTQI) THEN
         NBNKDS(4)=1
         WRITE(6,707)
 707     FORMAT(' Writing MCTQ bank : ON')
      ELSE
         NBNKDS(4)=0
         WRITE(6,708)
 708     FORMAT(' Writing MCTQ bank : OFF')
      END IF

      IF (ITQI) THEN
         NBNKDS(5)=1
         WRITE(6,709)
 709     FORMAT(' Writing TQ bank : ON')
      ELSE
         NBNKDS(5)=0
         WRITE(6,710)
 710     FORMAT(' Writing TQ bank : OFF')
      END IF

      IF (ITQAI) THEN
         NBNKDS(6)=1
         WRITE(6,711)
 711     FORMAT(' Writing TQA bank : ON')
      ELSE
         NBNKDS(6)=0
         WRITE(6,712)
 712     FORMAT(' Writing TQA bank : OFF')
      END IF

      IF (ITQREALI) THEN
         NBNKDS(7)=1
         WRITE(6,715)
 715     FORMAT(' Writing TQREAL bank : ON')
      ELSE
         NBNKDS(7)=0
         WRITE(6,716)
 716     FORMAT(' Writing TQREAL bank : OFF')
      END IF

      IF (ITQAREALI) THEN
         NBNKDS(8)=1
         WRITE(6,717)
 717     FORMAT(' Writing TQAREAL bank : ON')
      ELSE
         NBNKDS(8)=0
         WRITE(6,718)
 718     FORMAT(' Writing TQAREAL bank : OFF')
      END IF

      IF (ITQZI) THEN
         NBNKDS(9)=1
         WRITE(6,719)
 719     FORMAT(' Writing TQZ bank : ON')
      ELSE
         NBNKDS(9)=0
         WRITE(6,720)
 720     FORMAT(' Writing TQZ bank : OFF')
      END IF



#ifdef FADC
      IF (IFADCI) THEN
         WRITE(6,713)
 713     FORMAT(' Writing FADC bank : ON')
      ELSE
         WRITE(6,714)
 714     FORMAT(' Writing FADC bank : OFF')
      END IF
#endif

      RETURN
      END
