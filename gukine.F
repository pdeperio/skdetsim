************************************************************************
*     -----------------
      SUBROUTINE GUKINE
*     -----------------
*
*     (Purpose)
*        Define the kinematics of all tracks in the current event.
*
*    Modified:  4-Aug-98 to handle muon decays immediately
*               (D.Casper)
*    Add FADC (only if you defined): 30-Apr-2004 (M.Ishitsuka)
*    Add ESCALE (only if defined at build): 16-July-2012 (J. Hignight)
************************************************************************

      implicit none

#include "geant321/gcflag.inc"
#include "geant321/gckine.inc"
#include "geant321/gctrak.inc"
#include "skmc.h"
#include "gupht.h"
#include "dsbnkf.h"
#include "apmultmp.h"      
#include "apcard.h"
#include "aptrkseg.h"
#include "muvtx.h"

#include "skday.h"
#include "dsparm.h"

#include "skhead.h"
#include "skdetsim_parm.h"
#include "skparm.h"
#include "skruninf.h"  

#include "guflag.h"

#ifdef FADC
#include "fadcwrbnk.h"
#endif

      integer maxvec
      PARAMETER(MAXVEC=100)
      integer IVTX(MAXVEC)
      real VERTEX(3),PLAB(3)
      CHARACTER*32 CNAME
      integer IBUF(20), ibuf_head(6)
      real RBUF(20), rbuf_head(2)
      EQUIVALENCE (IBUF(1),RBUF(1))
      EQUIVALENCE (IBUF_head(1),RBUF_head(1))

      REAL vmod, xmlep

      integer lrun,elapsed,lrunsav
      real watert
      data lrunsav/-1/
      save lrunsav

      integer nsube, isube,
     &     ier, nw, nvtxvc, i, nv, iring, nt
      integer iseed_rdm1,iseed_rdm2,iseed_rdm3,iseed_rdm4
     &     ,iseed_rdm5,lux


** for TBA with gain 2020/03/26 M.Shinoki (2021/11/04 Y. Nakano)
      real slope_ms, intercept_ms

c** for gain simulation
      integer lenchr, j
      CHARACTER*6 tmpchar
      logical*4 ex

#ifdef FADC
      if (IFADCI) then
** FADC initialization per event
            Call fadcwrinit
      endif
#endif

C-- output random number
      call grndmq(iseed_rdm1,iseed_rdm2,0,'G')
      CALL RLUXAT(LUX,ISEED_RDM3,ISEED_RDM4,ISEED_RDM5)
      if(IFMESS.eq.0) write(6,*) 'EVENT SEEDS',iseed_rdm1,iseed_rdm2,iseed_rdm3
     &     ,iseed_rdm4,iseed_rdm5

C---- Decide if there is a muon decay to process.  If there is an
C---- event in memory right now, it means there is.

      call kzenum(nsube,isube)
      if (nsube .gt. 0) then    !  Yes, there is a delayed trigger
         call dssubvc(ier)
         if (ier.ne.0) then
            write(6,*) 'GUKINE: failed to create sub-event; stopping'
            call kzeclr         !Current event is invalid
            ieorun = 1
            ieotri = 1
            goto 99
         endif
         call iniseg
      else
         
C --- Generate Vector

         IF(IFLGV.EQ.0) THEN    !ZBS kinematics file

            DO WHILE (ISKIPEV.GE.1)
               CALL KZECLR
               CALL KZREAD(LIN,IER)
               IF(IER.NE.0) THEN
                  IEORUN = 1
                  IEOTRI = 1
                  WRITE(6,*) 'The data in your file is finished'
                  WRITE(6,*) 'num of event = ', IEVENT-1
                  GOTO 99
               ENDIF
               call skverchk
               ISKIPEV = ISKIPEV - 1
            ENDDO
            call iniseg
            CALL KZECLR
            CALL KZREAD(LIN,IER)
            IF(IER.NE.0) THEN
               IEORUN = 1
               IEOTRI = 1
               WRITE(6,*) 'The data in your file is finished'
               WRITE(6,*) 'num of event = ', IEVENT-1
               GOTO 99
            ENDIF
            call skverchk

*** check MCHEAD bank and add the random seed into the 1st segment
            call kzbloc('MCHEAD',ier)
            if(ier.ne.0) then
               CALL KZGET1('MCHEAD',0,NW,IBUF)
               call kzbdel('MCHEAD')
               if (ibuf(1) .ne. (SK_GEOMETRY+1000)) then
                write(*,*) 'GUKINE: SK_GEOMETRY in MCHEAD in input file was:',ibuf(1)
                write(*,*) 'GUKINE: It is overwrite to:',SK_GEOMETRY+1000
                ibuf(1) = SK_GEOMETRY+1000
               endif
            else
               NW = 13
               IBUF(1) = SK_GEOMETRY+1000
               IBUF(2) = 1
               IBUF(3) = 1
               IBUF(4) = 1
               IBUF(5) = 1
               IBUF(6) = 1
            endif
            CALL KZBCR1('MCHEAD',NW,'-I',2,6,'6I',IER)
            IF (IER.NE.0) then
               write(6,*) 'MCHEAD bank creation failed'
               stop
            endif
            CALL KZREP1('MCHEAD',0,NW,IBUF)
            IBUF(1) = iseed_rdm1 ! Geant seed
            IBUF(2) = iseed_rdm2
            IBUF(3) = LUX       ! RLUX seed
            IBUF(4) = iseed_rdm3
            IBUF(5) = iseed_rdm4
            IBUF(6) = iseed_rdm5
            CALL KZREP1('MCHEAD',1,6,IBUF)
***

ccyt            IF(IMWATR.EQ.-1 .or. IMTBA.EQ.-1) THEN
ccyn            IF(IMWATR.EQ.-1 .or. IMTBA.EQ.-1 .or. IMTBA.EQ.3) THEN
            IF(IMWATR.EQ.-1 .or. IMTBA.EQ.-1 .or. IMTBA.EQ.3 .or. IMTBA.EQ.2) THEN
               CNAME = 'LOWMC'
               call kzbloc(CNAME,IER)
               IF(IER.NE.0) THEN
                  CALL KZGET1(CNAME,1,NW,IBUF)
                  lrun = IBUF(1)
                  IF (WATER_TUNE .NE. AP_TUNE) THEN ! LE, UNI_TUNE
                     if(lrun.ne.lrunsav) then
                        lrunsav = lrun
                        nrunsk4_for_gain = lrun  ! added by y.t. (2018/04/15)
                        elapsed = relapse(lrun) ! get elapse day from 96/1/1
                        call lfwater(elapsed,watert) ! cal. water tranceparancy
                        WATRDS = watert
                        call tr2fac(WATRDS)
                        print*,'Run # ',lrun,' ; water tra.=',watert

c**** trigger threshold change 34 hit -> 31 hit by y.nakano/y.takeuchi (2018/03/06)
c                        call softtrg_get_cond(softtrg_detector, softtrg_thr, 
c     $                       softtrg_t0_offset, softtrg_pre_t0, softtrg_post_t0)
c                        if (lrun.ge.73583 .and. softtrg_thr(2).ne.31) then
c                           print *, 'gukine: SLE softtrg_thr(2) is changed from ',
c     $                               softtrg_thr(2), ' to ', 31, ' at Run = ', lrun
c                           softtrg_thr(2) = 31
c                           call softtrg_set_cond(softtrg_detector, softtrg_thr, 
c     $                          softtrg_t0_offset, softtrg_pre_t0, softtrg_post_t0)
c                        endif

                        if(SK_GEOMETRY.ge.SK_IV) then
                           if(nrunsk4_for_gain.ne.0) then
                             call softtrg_get_cond(softtrg_detector, softtrg_thr,
     $                            softtrg_t0_offset, softtrg_pre_t0, softtrg_post_t0)
                             call get_run_softtrg(nrunsk4_for_gain,
     $                            softtrg_detector, softtrg_thr,
     $                            softtrg_t0_offset, softtrg_pre_t0, softtrg_post_t0)
                             call softtrg_set_cond(softtrg_detector, softtrg_thr,
     $                            softtrg_t0_offset, softtrg_pre_t0, softtrg_post_t0)
                           endif
                        endif

                        if(IMTBA.eq.-1) then ! top-bottom asym.

ccyt                           IF(lrun .ge. 30000) THEN
ccyt                              call lftba(elapsed,TBADS)
ccyt                              TBADS = (-0.163*(TBADS*100)*(TBADS*100) - 3.676*(TBADS*100) )*0.00001
ccyt                           ELSE
ccyt                              WRITE(0,*)'Error!!'
ccyt                              WRITE(0,*)'TBA is supported after SK-III'
ccyt                              STOP 1
ccyt                           ENDIF

ccyt TBA copied from dscard.F (2018/03/02)
                           IF(LRUN .ge. 30000.and.LRUN.lt.60000) THEN
                              call lftba(elapsed,TBADS) ! calculate top bottom asym for SK=III
                              TBADS = -2.04035 * (TBADS*100) * 0.00001
                           ELSE IF(LRUN.ge.60000.and.LRUN.le. 77958)then
                              call lftba(elapsed,TBADS) ! calculate top bottom asym for SK=IV
                              TBADS = (-0.163*(TBADS*100)*(TBADS*100) - 3.676*(TBADS*100) )*0.00001
                           else if(LRUN.gt.77958)then
                              print *,'time varience for SK-IV is supported until 2018/5/30'
                              print *,'temporary, input tba param of Apr.2009'
                              TBADS = 0.0001003
                           ELSE
                              WRITE(0,*)'gukine : Error!! Bad LRUN ', LRUN
                              WRITE(0,*)'Support after SK-III'
                              STOP 1
                           ENDIF
                           WRITE(6,*)' The top-bottom asym. param. is set to', TBADS,'.'
                         ELSE IF(IMTBA.EQ.2) THEN
                           IF(LRUN.ge. 30000.and.LRUN.lt.60000) THEN
                              call lftba(elapsed,TBADS) ! calculate top bottom asym for SK=III
                              TBADS = -2.04035 * (TBADS*100) * 0.00001
                           ELSE IF(LRUN.ge.60000.and.LRUN.le. 77958) then
                              call lftba(elapsed,TBADS) ! calculate top bottom asym for SK=IV
cccccccmodify atode
cccccccModified by Y. Nakano and M. Harada 20211104 Mid-night                              
                              TBADS = (-0.163*(TBADS*100)*(TBADS*100) - 3.676*(TBADS*100) )*0.00001
c                              TBADS = -2.04035 * (TBADS*100) * 0.00001
                           else if(LRUN.gt.77958.and.SK_GEOMETRY.EQ.SK_IV)then
                              print *,'time varience for SK-IV is 
     &                                    supported until 2018/5/30'
                              print *,'temporary, input tba param
     &                                   of Apr.2009'
                              TBADS = 0.0001003
                              STOP 1
                           else if (LRUN.gt.80000.and.LRUN.le.85985) then
                              call lftba(elapsed,TBADS) ! calculate top bottom asym for SK-V/SK-IV
                              WRITE(6,*)'Check: Elapsed day = ', elapsed
                              WRITE(6,*)'Check: Official TBA value = ', TBADS
c                              TBADS = (-5.584 * (TBADS*100) - 11.31) * 0.00001 ! Fixed Y-int watert
                              TBADS = (-5.263 * (TBADS*100) - 10.92) * 0.00001 ! Free Y-int watert
                           else if(LRUN.gt.85985)then
                              print *,'time varience for SK-V/IV is supported until R#85985'
                              print *,'Please upload TBA file!'
                              STOP 1
                           ELSE
                              WRITE(0,*)'DSCARD : Error!! Bad DS-TBARUN',
     $                                   LRUN
                              WRITE(0,*)'Support after SK-III'
                              STOP 1
                           ENDIF
                        ELSE IF(IMTBA.EQ.3.and.SK_GEOMETRY.EQ.SK_IV) THEN ! gain correction only for SK-IV
                           IF(LRUN.ge.60000.and.LRUN.le.77958)then

                              call lftba_gain_corr(elapsed,TBADS) ! calculate top bottom asym for SK=IV
c     removed by Y. Nakano and M. Harada 20211104
c                              TBADS = (-0.163*(TBADS*100)*(TBADS*100) - 3.676*(TBADS*100) )*0.00001
                              WRITE(6,*)'Check: Elapsed day = ', elapsed ! 2020/03/17 M.Shinoki
                              WRITE(6,*)'Check: Official TBA value = ', TBADS ! 2020/03/17 M.Shinoki
c                              TBADS = (-0.163*(TBADS*100)*(TBADS*100) - 3.676*(TBADS*100) )*0.00001
c                              TBADS = -2.04035 * (TBADS*100) * 0.00001
c                              slope_ms = (-1.1175e-10)*WATRDS*WATRDS + (2.35166e-6)*WATRDS - 0.0164681 ! 2020/03/26 M.Shinoki
c                              slope_ms = (-9.88564e-11)*WATRDS*WATRDS + (2.14864e-6)*WATRDS - 0.0163051 ! 2020/05/09 M.Shinoki
                              slope_ms = (1.04159e-10)*WATRDS*WATRDS - (3.30078e-6)*WATRDS + 0.0203749 ! 2020/05/17 M.Shinoki
c                              intercept_ms = 0.953854*slope_ms + 0.00195722 ! 2020/03/26 M.Shinoki
c                              intercept_ms = -0.332257*slope_ms - 0.0110225 ! 2020/05/09 M.Shinoki
                              intercept_ms = 2.42604*slope_ms + 0.003059 ! 2020/05/17 M.Shinoki
                              TBADS =(slope_ms*(TBADS*100)+intercept_ms)*0.01 ! 2020/03/26 M.Shinoki 
                           else if(LRUN.gt.77958.and.SK_GEOMETRY.EQ.SK_IV) then
                              print *,'time varience for SK-IV is supported until R#77958'
                              print *,'Please upload TBA file!'
                              STOP 1
                           ELSE
                              WRITE(0,*)'gukine : Error!! Bad LRUN ', LRUN
                              WRITE(0,*)'Support (gain corr tba) after SK-IV'
                              STOP 1
                           ENDIF
                           WRITE(6,*)' The gain corr top-bottom asym. param. is set to', TBADS,'.'
                        ENDIF

ccyt PMT-by-PMT dark noise rate for SK-IV   copied from dscard.F (2018/03/02)
                        if (DOPMTDARKRT.ne.0) then
                           if(SK_GEOMETRY.ge.SK_IV) then
                              if(DOPMTDARKRT.eq.1) then
c                                do nothing
                              else if(DOPMTDARKRT.eq.2) then
                                 write(tmpchar,'(I6.6)')  lrun
                                 PMTNOISEFILE='/home/sklowe/const/darkr_for_skdetsim/output/darkr.'//tmpchar//'.txt'
                              else
                                 write(6,*) 'gukine: Not support darkrate option ',DOPMTDARKRT
                                 stop
                              endif

c*** open dark rate file, if exist
                              inquire(file=PMTNOISEFILE, exist=ex)
                              if (ex) then
                                 write(6,*) 'Using individual PMT-by-PMT dark rate for SK-IV from file:'
                                 write(6,*) PMTNOISEFILE(1:lenchr(PMTNOISEFILE))
                              else
                                 do i = 1, 100 
                                    write(tmpchar,'(I6.6)')  lrun-i
                                    PMTNOISEFILE='/home/sklowe/const/darkr_for_skdetsim/output/darkr.'//tmpchar//'.txt'
                                    inquire(file=PMTNOISEFILE, exist=ex)
                                    if (ex) then
                                       write(6,*) 'Alternative PMT-by-PMT dark rate for SK-IV from file:'
                                       write(6,*) PMTNOISEFILE(1:lenchr(PMTNOISEFILE))
                                       goto 7
                                    endif
                                 enddo
                                 write(6,*) 'gukine: Error cannot find PMT-by-PMT dark rate for SK-IV'
                                 stop
                              endif
 7                            continue
                              open (unit=83,file=PMTNOISEFILE,status="old")
                              do i = 1, maxpm
                                 read (83,*) j, darkrt(i)
                                 darkrt(i) = darkrt(i)*1000
                              end do
                              close (83)
                           else
                              write(6,*) 'Individual PMT-by-PMT dark rate not prepared for non-SK-IV.'
                              stop
                           endif
                        endif

                     endif
                  ENDIF
               ELSE
                  print*,'missing LOWMC bank...'
                  stop
               ENDIF
            ENDIF

         ELSE IF(IFLGV.EQ.1) THEN !Get kinematics from steering cards

            CALL SGGVCD
            CALL VCMKVC
            CALL VCMKVX
            call iniseg

         ELSE IF(IFLGV.EQ.2) THEN !Get kinematics from text file ala NUANCE

            do while (iskipev.ge.0)
               CALL KZECLR
               CALL SGGTXT(IER)
               IF (IER.NE.0) THEN
                  IEORUN = 1
                  IEOTRI = 1
                  WRITE(6,*) 'End of kinematic input file after ',IEVENT, ' events'
                  GOTO 99
               ENDIF
               ISKIPEV = ISKIPEV-1
            enddo
            call iniseg
            ISKIPEV = 0
         ELSE IF(IFLGV.EQ.-2) THEN

 11         CONTINUE
            CALL KZECLR
            If (ISKIPEV.gt.1) then
               ISKIPEV=ISKIPEV-1
               CALL KZREAD(LIN,IER)
               IF(IER.NE.0) THEN
                  IEORUN = 1
                  IEOTRI = 1
                  WRITE(6,*) 'The data in your file is finished'
                  WRITE(6,*) 'num of event = ', IEVENT-1
                  GOTO 99
               ENDIF
               call skverchk
               Goto 11
            Endif
            CALL KZECLR
            call iniseg
            CALL APDSSBREAD(LIN,IER)
         
            IF(IER.NE.0) THEN
               IEORUN = 1
               IEOTRI = 1
               WRITE(6,*) 'The data in your file is finished'
               WRITE(6,*) 'num of event = ', IEVENT-1
               GOTO 99
            ENDIF

            CALL VCMKVC
            CALL VCMKVX
         ELSE                   !Call user-provided kinemtic generator
! Added by J. Hignight on July 16, 2012
! For use by energy-scale group
! Previously only had CALL SGGVUS
#ifdef ESCALE    
            CALL SGGVUS_ESCALE(ier)
            if(ier.ne.0) then
               ieorun = 1
               ieotri = 1
               goto 99
            endif
#else 
            CALL SGGVUS
#endif
! End of part added by J. Hignight
            CALL VCMKVC
            CALL VCMKVX
            call iniseg
         ENDIF
      endif                     !If processing secondary event

CCS      write(*,*) "SKDETSIM::gukine: IVENT=",IEVENT

C --- Read 'MCVERTEX' bank      

      CNAME = 'MCVERTEX'
      CALL KZGET1(CNAME,0,NW,IBUF)
      NVTXVC = IBUF(1)

      DO 100 I=1,NVTXVC

         CALL KZGET1(CNAME,I,NW,IBUF)

         VERTEX(1) = RBUF(1)
         VERTEX(2) = RBUF(2)
         VERTEX(3) = RBUF(3)
         TOFG      = RBUF(6)

         CALL GSVERT(VERTEX,0,0,0,0,NV)
         IF(NV.EQ.0) THEN
            PRINT*,'ERROR IN GSVERT'
            STOP
         ENDIF
         IVTX(I) = NV

 100  CONTINUE

C --- Read 'MCVECT' bank      

C --- Add 1998/03/17 Y.hayato
      CALL VCRDVC

      CNAME = 'MCVECT'
      CALL KZGET1(CNAME,0,NW,IBUF)
      NRING = IBUF(1)
      If (nring.gt.MAXRNG) then
         print *, 'Error in gukine: too much particles. reduced '
     &        ,nring,' to ',MAXRNG
         nring=maxrng
      Endif

      DO 200 IRING=1,NRING

         CALL KZGET1(CNAME,IRING,NW,IBUF)
         
C         if (iring.eq.1) then
C            ip_nu = ibuf(1)
C            call ucopy(rbuf(3),p_nu,3)
C            p_nu(4) = vmod(rbuf(3),3)
C            p_nu(5) = p_nu(4)
C         else if (iring.eq.3) then
C            ip_lep = ibuf(1)
C            call ucopy(rbuf(3),p_lep,3)
C            call mcmass(ibuf(1),xmlep)
C            p_lep(4) = sqrt(vmod(rbuf(3),3)**2+xmlep**2)
C            p_lep(5) = vmod(rbuf(3),3)
C            p_had(1) = p_nu(1) - p_lep(1)
C            p_had(2) = p_nu(2) - p_lep(2)
C            p_had(3) = p_nu(3) - p_lep(3)
C            p_had(4) = p_nu(4) + 939. - p_lep(4)
C            p_had(5) = vmod(p_had,3)
C         endif


C     Check the neutrino energy. This is a special treatment 
C     of high energy atmospheric neutrino events.

      if (FLAG_TRACK_SW .eq. 1) then
       if (TRACK_PHOTON_ORG .eq. 3) then
        if (
     & ((abs(ibuf(1)) .eq. 12) .or. (abs(ibuf(1)) .eq. 14) .or.
     &  (abs(ibuf(1)) .eq. 16) ) .and. 
     &  (rbuf(2) .gt. 1000000.0) ) then
           write(*,*) '******************************************************'
           write(*,*) 'GUKINE:'
           write(*,*) 'Neutrino Energy is higher than 1TeV.'
           write(*,*) 'Tracking method is changed to 4 from',TRACK_PHOTON_ORG
           write(*,*) 'ID and Enu:',ibuf(1),rbuf(2)
           write(*,*) '******************************************************'

          TRACK_PHOTON = 4

            IBUF_head(1) = TRACK_PHOTON
            IBUF_head(2) = TRACK_PHOTON_ORG
            CALL KZREP1('MCHEAD',2,6,ibuf_head)

        endif
       endif
      endif


C     This next line checks the "tracking" flag for this particle.  If
C     it is set to zero it should not be put on the stack and passed to
C     geant.  In this line, IBUF(8) corresponds to ICRNVC.

         IF(IBUF(8).EQ.0) GOTO 200

         CALL PD2GT(IBUF(1),IPRT(IRING))
         IF(IPRT(IRING).EQ.0) THEN
            WRITE(6,*) 'Your Particle is NOT here!'
     &           ,'  ring,pid=',iring,IBUF(1)
            WRITE(6,*) 'Please check the particle number'
            STOP
         ENDIF

         PLAB(1) = RBUF(3)*(1.E-3)
         PLAB(2) = RBUF(4)*(1.E-3)
         PLAB(3) = RBUF(5)*(1.E-3)

         NV = IVTX(IBUF(7))

         CALL GSKINE(PLAB,IPRT(IRING),NV,0,0,NT)
         IF(NT.EQ.0) THEN
            PRINT*,'ERROR IN GSKINE'
            STOP

         ENDIF

C Check for initial muon, so upmu information generated later on is
C not confused by secondary muons
C not perfect, but works most of the time
         if ((IPRT(IRING).eq.5).or.(IPRT(IRING).eq.6)) then
c            print *,'We have a muon - nt ',nt
            first_id_mu = nt
         endif

         itrcon(nt)=iring
         
 200  CONTINUE   

 99   RETURN
      END
