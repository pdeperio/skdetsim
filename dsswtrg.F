*
      subroutine dsswtrg
*
*     software trigger and related procedures
*
*     2008/07 First version by Jun Kameda
*     2013/02 T. Wongjirad
*             Modified to store trigger information in ZBS.
*             Now stores event true t0 in TDC units.
*             Also stores elapsed time from t0 to trigger.
*             The variable idtgsk is set in similar manner to data,
*             in which overlapping triggers all set the bits of idtgsk.
*             Stores all triggers found by software trigger as well.
*             Has been checked to not interfere with simulation.
*             These mods. only turn on if DSSWTRIGSTOREINFO compiler flag
*             is included. Is NOT on by default.
*

      implicit none

#include "skparm.h"
#include "skhead.h"
#include "skpdst.h"
#include "geopmt.h"
#include "sktq_raw.h"
#include "softtrg_tblF.h"
#include "skmc.h"
#include "dsdata.h"
#include "dsparm.h"
#include "geant321/gcflag.inc"
#include "skruninf.h"

#include "guflag.h"

      integer ncand, in_gate_only
      integer softtrg_inittrgtbl
      external softtrg_inittrgtbl
      integer iii, flag_t0
      integer bankerr
      integer i,j
      integer dstrgbuf_master(5)
      integer dstrgibuf_seg1(4) ! Store all trigger info
      real dstrgrbuf_seg1(4)
      integer dsevtrig_ibuf(4) ! Store primary trigger info
      real dsevtrig_rbuf(4)
      equivalence (dstrgibuf_seg1(1), dstrgrbuf_seg1(1))
      equivalence (dsevtrig_ibuf(1), dsevtrig_rbuf(1))
      integer dsswtrgid, trigbit, dsswbankerr
      integer numdsswtrgs
      integer it0sk_geantt0, it0sk_temp
      real geantt0
      integer prim_trg
      real prim_pret0

#ifdef DSSWTRIGSTOREINFO
c     Clear the number of triggers
      numdsswtrgs = 0
      dsswtrgid = -1
c     First clear the software trigger buffers
      do 22 I=1,4
         dstrgibuf_seg1(i) = -1 ! used to fill each trigger
         dsevtrig_ibuf(i) = -1 ! used to store primary trigger
 22   continue
      dsevtrig_rbuf(3) = 0.0
      dstrgrbuf_seg1(3) = 0.0
c     define and clear bank
      dsswbankerr = 0
      CALL KZBCR1('DSTRGOFF', 5, '5I', 10, 4, '2I1F1I', dsswbankerr)
      do 23 I=1,10
         CALL KZREP1('DSTRGOFF', I, 4, dstrgibuf_seg1)
 23   continue

c     Store the geant0, in clock units
      it0sk_geantt0 = (ntrg+2)*2**15 + int(trg_offset/tdc_unit_sk4) ! (data block+2)*(16-bits) + geantt0 pos (in TDC)
c     geantt0 in ns (0 is set by arb. clock position)
      geantt0  = real(it0sk_geantt0*tdc_unit_sk4)
c     Define the gate width (taken from sofkl/src/softtrg/make_swtrg_list.c)
#define SAME_GATE_WIDTH 768

#endif

      if(IFSWTRG.eq.0) then
c     USE THE SOFTWARE TRIGGER
         flag_t0 = 0
         dsswtrgid = -1
c     GET A LIST OF SOFTWARE TRIGGERS
         ncand = softtrg_inittrgtbl(nrunsk, 0, in_gate_only, 1280)

c     DETERMINE THE CANDIDATE
         prim_trg = -1          ! the current primary trigger (-1=none found)
         prim_pret0 = 0.0       ! t0 of trigger
         if (IFMESS.eq.0) then
            write(*,*) '===== DSSWTRG: RUNNING SOFTWARE TRIGGER ======'
            write(*,*) '  SKRUNINF Trigger mask',softtrg_mask
            write(*,*) '  SW Trigger. # of trigger Cand.:',ncand
            write(*,*) '  # of hit PMT in full range:',nqisk_raw
         endif
         
c     -------------------------------------------------------------------------------------
c     original code: loop. look for LE/HE/OD trigger and store it0sk. ignores SLE triggers.
c     
         if (IFMESS.eq.0) then
            write(*,*) '  ------- Result of Software Trigger ----------'
         endif
         do 20 I=1,ncand
            
            if ( (swtrgtype(i) .eq. 0 ) .or. (swtrgtype(i) .eq. 1 )
     &           .or.(swtrgtype(i) .eq. 3)) then ! LE or HE or OD
               
               if ( flag_t0 .eq. 0 ) then
c     IF THIS THE FIRST TRIGGER (flag_t0==0), THEN STORE AND SET FLAG
                  it0sk = swtrgt0ctr(i)
                  flag_t0 = 1
                  prim_trg = i  ! new code
                  idtgsk = ISHFT(1, swtrgtype(i)) ! new code
               else             ! end of if first trigger (flag_t0==0)
c     flag already set
                  if (swtrgt0ctr(i).le.it0sk) then
c     trigger time is ealier, change the it0sk value
                     it0sk = swtrgt0ctr(i)
                     prim_trg = i ! new code
                     idtgsk = ISHFT(1, swtrgtype(i)) ! new code
                  endif         ! if trigger time before current trigger
                  
               endif            ! first trigger or not
               
            endif               ! if LE/HE/OD

c     output info for each trigger
c            if (IFMESS.eq.0) then
c               write(*,'(A22,I12)')
c     $              '   SWTRG Cand     #   ',i
c               write(*,'(A22,Z12.12)') 
c     $              '   SWTRG TYPE      :0x',ishft(1,swtrgtype(i))
c               write(*,'(A22,I12)') 
c     $              '   SWTRG INDEX     :  ',swtrgidx(i)
c               write(*,'(A22,I12)') 
c     $              '   SWTRG T0 HW-CTR :  ',swtrgt0hwctr(i)
c               write(*,'(A22,I12)') 
c     $              '   SWTRG T0 LOC    :  ',swtrgt0loc(i)
c               write(*,'(A22,I12)') 
c     $              '   SWTRG T0 TDC CTR:  ',swtrgt0ctr(i)
c               write(*,'(A22,I12)') 
c     $              '   SWTRG GS HW-CTR :  ',swtrggshwctr(i)
c               write(*,'(A22,I12)')
c     $              '   SWTRG GS TDC CTR:  ',swtrggsctr(i)
c               write(*,*)
c     $              '   SWTRGTYPE:  ',swtrgtype(i),it0sk,swtrgt0ctr(i),
c     $              it0sk-swtrgt0ctr(i)
c            endif
c     End of trigger search (original code)
 20      continue               ! loop over candidate triggers

c-------------------------------------------------------------------------------------------         
c***  add SLE/LE/HE/OD trigger information around T0+/-200 TDC count (2018/03/30 y.takeuchi)
         do i = 1, ncand
            if ( (swtrgtype(i) .eq. 0 ) .or. (swtrgtype(i) .eq. 1 )
     &           .or.(swtrgtype(i) .eq. 2) 
     &           .or.(swtrgtype(i) .eq. 3)) then 

               if ( flag_t0 .eq. 0 ) then
c     IF THIS THE FIRST TRIGGER (flag_t0==0), THEN STORE AND SET FLAG (may add 1st SLE trigger)
                  it0sk = swtrgt0ctr(i)
                  flag_t0 = 1
                  prim_trg = i  ! new code
                  idtgsk = ISHFT(1, swtrgtype(i)) ! new code
               endif             ! end of if first trigger (flag_t0==0)

               if (abs(swtrgt0ctr(i) - it0sk).lt.384) then   ! 1.92count/ns -> 200ns = 384 count
                  idtgsk = ior(idtgsk, ISHFT(1, swtrgtype(i)))
               endif

            endif
         enddo
c-------------------------------------------------------------------------------------------         


         if (flag_t0 .eq. 0) then ! no LE/HE/OD trigger, t0 is set to GEANT t0 (in units of hardware clock position)
            it0sk = (ntrg+2)*2**15 + int(trg_offset/tdc_unit_sk4)
            idtgsk = 0
         endif

c         if (IFMESS.eq.0) then
c            write(*,*) 'SET IT0SK: ',it0sk
c            write(*,*) 'SET IDTGSK: ',idtgsk
c         endif
c     ===============================================================================

         
#ifdef DSSWTRIGSTOREINFO
c     ===============================================================================
c     NEW CODE: IMMITATE BEHAVIOR OF TRIGGER FOR ONLINE DATA, STORE TRIGGER INFO
c     Basically want the trigger id to match the way the data works.
c     If overlapping triggers, then bits are combined using OR

c     [Following Removed for backwards compatibility]
c     ************************************************************************
c     Search for primary trigger.
c     According to Hayato-san's Collaboration talk,
c     the primary trigger will be the one with the earliest pre_t0 time.
c     see skofl/src/softtrg/softtrg_cond.h for details of algorithm.
c     code here mimics that code.
c     However, going to use original code's version of selecting trigger.
c
c         do 25 i=1,ncand
c            if (prim_trg==-1) then
c     first trigger, just fill the variables
c               prim_trg = i
c               prim_pret0 = swtrggsctr(i)
c            else
c     subsequent triggers, replace if pret0 is earlier
c               if (swtrggsctr(i).lt.prim_pret0) then
c                  prim_trg = i
c                  prim_pret0 = swtrggsctr(i)
c               endif
c            endif
c 25      continue
c     ************************************************************************************

c     STORGARE OF TRIGGER INFO
c     Store primary trigger info

         if (prim_trg.ne.-1) then
            it0sk_temp = swtrgt0ctr(prim_trg) ! trigger time, in TDC counts
            trigbit = ISHFT(1, swtrgtype(prim_trg)) ! trigger bit
            dsswtrgid = trigbit
         else
            it0sk_temp = it0sk  ! should be set to geantt0
            trigbit = 0
            dsswtrgid = -1
         endif
         dsevtrig_ibuf(1) = trigbit ! Trigger Type
         dsevtrig_ibuf(2) = it0sk_temp ! Trigger time
         dsevtrig_rbuf(3) = real((it0sk_temp-(ntrg+2)*2**15)
     &        *tdc_unit_sk4)-trg_offset ! time of event trigger relative to geant t0
         dsevtrig_ibuf(4) = prim_trg

c     Look for overlap triggers (mimics output of data)
         if (prim_trg.ne.-1) then
            do 26 i=1,ncand
c     loop over triggers again to get overlap triggers
               if ( i.ne.dsevtrig_ibuf(4).and.
     &              (ABS(it0sk-swtrgt0ctr(i))< SAME_GATE_WIDTH)) then
                  trigbit = ISHFT(1, swtrgtype(i)) ! trigger bit
                  dsevtrig_ibuf(1) = IOR( dsevtrig_ibuf(1), trigbit )
               endif
 26         continue
         endif

c     STORE EACH TRIGGER, UP TO MAX (9 total)
         do 27 i=1,ncand
            numdsswtrgs = numdsswtrgs + 1
            if (numdsswtrgs < 10 ) then
               trigbit = ISHFT(1, swtrgtype(i))
               dstrgibuf_seg1(1) = trigbit
               dstrgibuf_seg1(2) = swtrgt0ctr(i)
               dstrgrbuf_seg1(3) = real( ((ntrg+2)*2**15-swtrgt0ctr(i))
     &              *tdc_unit_sk4)+trg_offset ! ns from event trigger to geant t0 (should be negative)
               dstrgibuf_seg1(4) = numdsswtrgs
               CALL KZREP1('DSTRGOFF', 1+numdsswtrgs, 4, dstrgibuf_seg1)
            endif
 27      continue               ! loop over candidate triggers
#endif                          ! end of #ifdef DSSWTRIGSTOREINFO
            
c     END OF STATEMENTS IF USE SOFTWARE TRIGGER FLAG IS ON (IFSWTRG.eq.0)
      else if(IFSWTRG.eq.1) then
c     SET THE SOFTWARE TRIGGER to position of geant t0 + user offset
         it0sk = (ntrg+2)*2**15+int(trg_offset/tdc_unit_sk4)+IT0_OFFSET
      else
         write(6,*) 'Not support; DS-NONSWTRG ',IFSWTRG
         stop
      endif
      it0xsk=it0sk

#ifdef DSSWTRIGSTOREINFO
c     STORE TRIGGER INFO TO BANK. SUMMARIZE TRIGGER SEARCH RESULTS
      dsswtrgid = dsevtrig_ibuf(1)
      if (IFMESS.eq.0) then
         write(*,*) '  -- Triger Summary --'
         write(*,*) '   Trigger time is set to : ',it0sk,
     &        ' (= ',real(it0sk*tdc_unit_sk4),' )'
         write(*,*) '   IT0_OFFSET = ',IT0_OFFSET
         write(*,*) '   dsswtrgid: ',dsswtrgid
         write(*,*) '   primary trigger type: ',dsevtrig_ibuf(1)
         write(*,*) '   primary trigger it0sk: ',dsevtrig_ibuf(2)
         write(*,*) '   primary trigger pre_t0: ',prim_pret0
         write(*,*) '   primary trigger dgeantt0: ',dsevtrig_rbuf(3)
         write(*,*) '   primary trigger number: ',dsevtrig_ibuf(4)
         write(*,*) '   geant t0 in hardware time: ',it0sk_geantt0
         write(*,*) '   geant t0 relative to software trigger:',
     &        geantt0,' ns'
      endif

c     FILL MASTER SEGMENT
      dstrgbuf_master(1) = 4    ! length of master segment
      dstrgbuf_master(2) = 4 ! length of other segments
      dstrgbuf_master(3) = IT0_OFFSET
      dstrgbuf_master(4) = it0sk_geantt0
      dstrgbuf_master(5) = numdsswtrgs
      CALL KZREP1('DSTRGOFF', 0, 5, dstrgbuf_master)

c FILL PRIMARY TRIGGER SEGMENT: SEGMENT 1
      if (IFMESS.eq.0) then
         write(*,*) '  STORING TRIGGER ID:',dsevtrig_ibuf(1)
         write(*,*) '  STORING IT0SK:',dsevtrig_ibuf(2)
         write(*,*) '  STORING dGEANT0:',dsevtrig_rbuf(3)
         write(*,*) '  STORING TRIGNUM:',dsevtrig_ibuf(4)
         write(*,*) '  SET IDTGSK to ',dsevtrig_ibuf(1)
      endif
      idtgsk = dsevtrig_ibuf(1)
      CALL KZREP1('DSTRGOFF', 1, 4, dsevtrig_ibuf)
      
      if (IFMESS.eq.0) then
      write(*,*) 'end of trigger summary'
      write(*,*) '-----------------------------------'
      endif
#endif
      return
      end
