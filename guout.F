*     CMZ :  3.12/27 06/09/88  14.33.45  by  Rene Brun
*--   Author :
      SUBROUTINE GUOUT
C     .
C     .    ******************************************************************
C     .    *                                                                *
C     .    *       User routine called at the end of each event             *
C     .    *                                                                *
C     .    *    ==>Called by : GTRIG                                        *
C     .    *                                                                *
C     .    *        First Version : 1993.02.25    by Y.Koshio               *
C     .    *                                                                *
C     .    *   Modified: 4 August 1998 to process subevents concurrently    *
C     .    *             (D.Casper)                                         *
C     .    ******************************************************************
C     .
C     .
C     .    ------------------------------------------------------------------
C     .
      
      
      implicit none
      
#include "skhead.h"
#include "geant321/gcflag.inc"
#include "dsbnkf.h"
#include "aptrkseg.h"
#include "skdetsim_parm.h"
#include "skparm.h"
#include "skdetsim_test.h"
#include "guflag.h"

#ifdef FADC
#include "fadcwrbnk.h"
#endif
      
#include "apscndry.h"

      logical more
      integer icount, icall, MAXNPHASE
      save icount, icall
      DATA icount/0/
      DATA icall/0/
      parameter(MAXNPHASE=5)
      integer ntot, icur, nmue
      character*256 mc_subnum
      integer eventsave
      save eventsave
      integer lbank
      
C     Added by T. Wester
      integer maxnbuf
      parameter (maxnbuf=8)
      integer ibuf(maxnbuf)
      integer ierr
      
      if(OUTPUT_TYPE.eq.2) then ! For ATMPD group
         
         If (icall.eq.0) then
C     read subrun number from environment variable 
C     MC_SUBNUM
C     added by M.Shiozawa  99/10/04
            call getenv('MC_SUBNUM',mc_subnum)
            read(mc_subnum,'(i3)') nsubsk
            icall=1
         Endif
         call skfakehd
         call skmkruninf
CC to save the storage in laser MC
#ifndef LASER
         call sgcardbnk
         Call skmkmcver
#endif
         
c        Update Header bank for SK-IV. Trigger informations are added.

         if (SK_GEOMETRY .ge. SK_IV) then
            call sktrginfo
         endif
c

      if (FLAG_TRACK_SW .eq. 1) then
         if (TRACK_PHOTON .ne. TRACK_PHOTON_ORG) then

          write(*,*) '******************************************************'
          write(*,*) 'GUOUT:'
          write(*,*) 'TRACKING method is changed back to',TRACK_PHOTON_ORG  
          write(*,*) '******************************************************'
          track_photon = track_photon_org
         endif
      endif

C     Write out muon vertex bank and reset common
#ifndef LASER
         call mkmubnk
#endif
         call clrmuver
         
c     support manual specification of badch file
c     added by t.wester. 2019
c    ---------------------------------------
         write(*,*) 'BADCHRUN set to',BADCHRUN
         if (BADCHRUN.ne.0) then
c     write APBADCHINF bank which is used by skread
             ibuf(1) = BADCHRUN 
             ibuf(2) = 1
             ibuf(3) = 999
             ibuf(4) = 0
             call kzbcr1('APBADCHINF',maxnbuf,'-I',0,0,' ',ierr)
             if(ierr.ne.0)stop 'SKMCBADCH-I-cannot create bank'
             call kzrep1('APBADCHINF',0,maxnbuf,ibuf)
         endif
c     -------------------------------------


         if (SK_GEOMETRY .le. SK_III) then ! Consider sub-events

C     Check for particles outside the current event gate
         
            call dssubchk(more, nmue)
         
            if (more) then
               nevent = nevent + 1
#ifdef FADC
               flag_fadc_sube = 1
#endif
            endif
      
*     if (flgseg .ne. seg_noseg) then
            if (OUTPUT_HBOOK.eq.OUT_SEG)
     &           then

               call packseg
               call packprt
               ndk = nmue
               call kzenum(ntot,icur)
               subevt = icur
               call hfnt(seg_idh)
            endif
      
            if (more) then
               call kzenum(ntot,icur)
               call kzemrk      ! Mark the current (sub-)event for merging/output
            else            
               call kzenum(ntot,icur)
               call kzemrk      ! Mark the current (sub-)event for merging/output
               call kzefst      ! Return to the primary event
               call kzemrg      ! Merge all sub-events with it
               if (validate) call do_validate ! do this just before writing out, 
                                ! to cover all attempts to break things
               if ((TRACK_PHOTON .ne. NOT_TRACK_PHOTON).and.(OUTPUT_ZBS)) then
                  CALL KZWRIT(LOUT)
               endif
               eventsave=eventsave+1
               if(IFMESS.eq.0) WRITE(6,*) 'SKDETSIM::guout: Processed RUN: ',nrunsk,' SUB:',nsubsk
     &              ,' EVENT:',eventsave
               CALL KZECLR
               call skapmcclr   ! reset count of apmcstop bank segments
               call dssubcl     !Clear running trigger time variable
               icount = icount + 1 !Don't change subrun between subevents!
            endif

         elseif (SK_GEOMETRY .ge. SK_IV) then

            if (validate) call do_validate ! do this just before writing out, 
                                ! to cover all attempts to break things
            if ((TRACK_PHOTON .ne. NOT_TRACK_PHOTON).and.(OUTPUT_ZBS)) then
               if(IFSHORTWIN.eq.1) call delete_outside_hits()
               CALL KZWRIT(LOUT)
            endif

            eventsave=eventsave+1
            if(IFMESS.eq.0) WRITE(6,*) 'SKDETSIM::guout: Processed RUN: ',nrunsk,' SUB:',nsubsk
     &           ,' EVENT:',eventsave
            CALL KZECLR
            call skapmcclr      ! reset count of apmcstop bank segments
            icount = icount + 1 !Don't change subrun between subevents!
         endif

C     
      else                      ! For LOWE group

c        Update Header bank for SK-IV. Trigger informations are added.

         if (SK_GEOMETRY .ge. SK_IV) then
            call skmkruninf
            call sktrginfo
         endif

         if (validate) call do_validate ! do this just before writing out, 
                                ! to cover all attempts to break things
         if ((TRACK_PHOTON .ne. NOT_TRACK_PHOTON).and.(OUTPUT_ZBS)) then
            if (SK_FILE_FORMAT.eq.0) then ! ZBS format
               if(IFSHORTWIN.eq.1) call delete_outside_hits()
               CALL KZWRIT(LOUT)
            else                ! skroot

c*** load TQ info. into standard variables
               goto 20
 10            continue
               print *, 'skread error!:guout: RUN: ',nrunsk,' SUB:',
     $              nsubsk,' EVENT:',eventsave  
               stop
 20            continue
               SK_FILE_FORMAT=0 ! read from zbs
               call skread(-LOUT,*10,*10,*10,*10)
               SK_FILE_FORMAT=1 ! output to skroot in only this program

c*** generate dummy header info. 
               nrunsk = 999999
               nsubsk = 999999
               nevsk  = eventsave + 1  ! to make 1st event = 1 (2018/03/04 y.takeuchi)
               nevhwsk = nevsk*100+2   ! (tentative) "ntrg+2" in sktrginfo.F   y.t
               gatewsk = 76800         !  9600+ 67200 (same as sktrginfo.F  y.t)

c*** apply 1.3 micro sec gate
               if(IFSHORTWIN.eq.1) call delete_outside_hits()
c*** store MC info
               call skroot_set_mcinfo(LOUT)
c*** store secondary info
               call skroot_set_secondary(LOUT,
     $              nscndprt, itrkscnd, istakscnd, iprtscnd,
     $              iprntprt, lmecscnd, iprnttrk, iorgprt,
     $              iflgscnd, iprntidx, nchilds, ichildidx,
     $              iflgscnd_temp, iprntid_temp, ichildid_temp,
     $              vtxscnd, pscnd, tscnd, pprnt, pprntinit,
     $              vtxprnt, npvcscnd, nchildsvc, ichildidxvc)
c*** store header & tq info
               call skroot_set_tree(LOUT)
c*** write an event
               call skroot_fill_tree(LOUT)
            endif
         endif
         eventsave=eventsave+1
         if(IFMESS.eq.0) WRITE(6,*) 'SKDETSIM::guout: Processed RUN: ',nrunsk,' SUB:',nsubsk
     &        ,' EVENT:',eventsave         
         CALL KZECLR
      endif
      

C
C+ Standard Diagnostic Test Histograms
C
      if (OUTPUT_HBOOK .eq. OUT_TEST) then
         TEST_Events=TEST_Events+1
         call HFNT(22)
         call HF1(40,FLOAT(TEST_Gen_phot),1.)
         call HF1(41,FLOAT(TEST_Coll_phot),1.)
         call HF1(42,FLOAT(TEST_Colla_phot),1.)
         call HF1(43,TEST_NHit,1.)
         TEST_gen_phot=0
         TEST_Coll_phot=0
         TEST_Colla_phot=0
         TEST_NHit=0.
      endif


      END
