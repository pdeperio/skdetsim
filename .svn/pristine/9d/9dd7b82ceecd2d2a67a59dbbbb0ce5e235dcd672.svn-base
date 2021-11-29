***********************************************************************
*     -------------------------------------
      SUBROUTINE dscnt(nphoton, ph_t_list, ph_pmt_list, ngate,
     & iclist_s,tlist_s,iqlist_s,qlist_s)
*     -------------------------------------
*
*     (Purpose)
*
*        Count # of photons in each QTC gate and hit timing.
*        Multi hit of the PMT is taken into account.
*
*     (Input) 
*        nphoton     : number of photons
*        ph_t_list   : list of the timing of the photons
*        ph_pmt_list : list of the hit PMT ID
*
*     (Output)
*        ngate_s           : number of hits (gates)  
*        iclist_s          : list of the hit PMT number     (time order)
*        iqlist_s          : number of photons in each gate (time order)
*
*     (Creation Date and Author)
*
*       2008.07.08 ; First version by Jun Kameda
*
************************************************************************

      implicit none

#include "skparm.h"
#include "skhead.h"
#include "skdetsim_parm.h"
#include "dsparm.h"
#include "skmc.h"
#include "dsdata.h"

#include "guflag.h"

      integer photon_id, nphoton
      integer ph_pmt_list(MXPHDS)
      integer index(MXPHDS)
      real    ph_t_list(MXPHDS)
      integer i,j,k,ii
      real pe1,eff,dseff_qbee
      external dseff_qbee

      real    tlist(MXGATE) ,QLIST(MXGATE),qlist_s(MXGATE)
      integer iclist(MXGATE),iqlist(MXGATE)
      

      !!!! sorted variables
      real    tlist_s(MXGATE)
      integer iclist_s(MXGATE),iqlist_s(MXGATE)
      integer sort_index(MXGATE)

      real    now_tstart(MAXPM)
      integer gate_flag(MAXPM),now_nphoton(MAXPM),now_gate_id(MAXPM)
      integer hit_flag(MAXPM),list_hit_pmt(MAXPM)
      integer ngate,npmt,nhit_pmt

      real now_charge(MAXPM)

      integer IBUF(3)
      real RBUF(3)
      EQUIVALENCE (IBUF(1),RBUF(1))



      integer iring,iph,ndata,ic,iq,ibef,ierr,ifthr,icab,
     &     nwmast,nsegm,lsegm,isseg,isegm
      real t,dt,absdt,pe,thr,tsig,trnd,qpc


*** For gain increasing (2017 June, by LOWE group)

      real gain_param
      real gain_nom(5), gain_nom5(5)
c     average 2008, Oct. for SK-IV
      data gain_nom/0.992057, 1.00184, 0.989083, 0.98216, 0.987678/
c     average 2019, Oct. for after SK-V
      data gain_nom5/0.99161, 0.98805, 0.99163, 0.99177, 0.99919/

      real gain5(5)
      save gain5

      integer readrun
      character*200 rskip(15)
      real pmtgain(MAXPM)
      save pmtgain
      integer gain_run
      data gain_run/0/
      save gain_run

      real prodyear
      integer pmtid
      character*200 pskip(4)

      real gainfac

      character*200 file_name

ccyt added to suppress too many logs (y.takeuchi)
      integer nprint
      save nprint
      data nprint/0/
***

cc  Sorting

      call sortzv(ph_t_list,index,nphoton,1,0,0)

cc  Initialize variables

      nhit_pmt = 0
      ngate = 0

      call vzero(now_gate_id,MAXPM) ! current gate ID for each PMT
      call vzero(now_nphoton,MAXPM) ! # of photons in current gate for each PMT
      call vzero(now_charge,MAXPM)  ! # of p.e's in current gate for each PMT
      call vzero(hit_flag, MAXPM)   ! 
      call vzero(now_tstart,MAXPM)  ! starting time of the current gate

      call vzero(gate_flag,MAXPM)   ! flag for gate

      call vzero(iclist,MXGATE)     ! PMTs ID of each hit
      call vzero(iqlist,MXGATE)     ! # of photons in each QTC gate
      call vzero(qlist,MXGATE)      ! # of p.e.s in each QTC gate
      call vzero(tlist,MXGATE)      ! Timing of each hit


*** For gain increasing (2017 June, by LOWE group)

ccccc artificial gain adjustment factor to test behaviour of gain vs nhit
      gainfac = 1.0
c      write (6,*) 'WARNING!!! gains are artificially multiplied by factor', gainfac
c      write (6,*) 'This was introduced to study gain vs nhit'
c      write (6,*) 'Make sure to remove this for any other studies!'

      if (nrunsk4_for_gain .ne. gain_run) then
         write (6,*) "Current run", nrunsk4_for_gain
         write (6,*) "Saved gain run", gain_run
         write (6,*) 'saving the gain run'
         gain_run = nrunsk4_for_gain
cc read pmt gain file
         call findconsts(file_name, 'wt_table', 'water.ave10.2', 2)
         open(101, file=file_name,status='old',form='formatted')
         readrun=0;
         write(6,*) "run",nrunsk4_for_gain

         do while(readrun.lt.nrunsk4_for_gain)
            read(101,*) readrun, rskip, gain5
         enddo

         write (6,*) "gain", gain5

         close(101)
cc read PMT production year and assign pmt gain
         if (SK_GEOMETRY.eq.SK_IV) then
            call findconsts(file_name,'pmt_prod_year','pmt_prod_year_sk4.dat',2)
         else
            call findconsts(file_name,'pmt_prod_year','pmt_prod_year.dat',2)
         endif
         
         open(102, file=file_name)
         do i=1,11146
            read(102,*) pmtid, prodyear, pskip
c         write(6,*) pmtid, prodyear,pskip
c         do j=1,4
c            read(102,*) dskip
c         enddo
            if (SK_GEOMETRY.eq.SK_IV) then

               if(prodyear.lt.1000)then
c               pmtgain(pmtid) = 0. !Bug!!!
                  pmtgain(pmtid) = 1.0
               else if(prodyear.lt.1996)then
                  pmtgain(pmtid) = gain5(1)/gain_nom(1)
               else if(prodyear.lt.1998)then
                  pmtgain(pmtid) = gain5(2)/gain_nom(2)
               else if(prodyear.lt.2004)then
                  pmtgain(pmtid) = gain5(3)/gain_nom(3)
               else if(prodyear.lt.2005)then
                  pmtgain(pmtid) = gain5(4)/gain_nom(4)
               else
                  pmtgain(pmtid) = gain5(5)/gain_nom(5)
               endif
            else if (SK_GEOMETRY.ge.SK_V) then
               if(prodyear.lt.1000)then
                  pmtgain(pmtid) = 1.0
               else if(prodyear.lt.1996)then
                  pmtgain(pmtid) = gain5(1)/gain_nom5(1)
               else if(prodyear.lt.1998)then
                  pmtgain(pmtid) = gain5(2)/gain_nom5(2)
               else if(prodyear.lt.2004)then
                  pmtgain(pmtid) = gain5(3)/gain_nom5(3)
               else if(prodyear.lt.2005)then
                  pmtgain(pmtid) = gain5(4)/gain_nom5(4)
               else
                  pmtgain(pmtid) = gain5(5)/gain_nom5(5)
               endif
            else
               write(*,*) 'dscnt: sk_geometry is not correct ',SK_GEOMETRY, ' stop.'
               stop
            endif

            pmtgain(pmtid) = gainfac * pmtgain(pmtid)

            gain_corr(pmtid) = pmtgain(pmtid)

c     write (6,*) 'pmtid/year/gain5', pmtid, prodyear, pmtgain(pmtid)
         enddo
         close(102)
      else
         if (nprint.lt.10) then
            write (6,*) 'skipping to read gain values'
            write (6,*) "gain", gain5
            nprint=nprint+1
         else if (nprint.eq.10)then
            write (6,*) 'dscnt: suppress gain message'
            nprint=nprint+1
         endif
      endif
***

cc  scan photons
!            write(79,*)'                  ', 0, 0,0, 0

       do 100 i = 1, nphoton

         photon_id = index(i)
                 t = ph_t_list(photon_id)
              npmt = ph_pmt_list(photon_id)

cc  fill hit PMT list

         if (hit_flag(npmt) .eq. 0) then
            nhit_pmt = nhit_pmt + 1
            hit_flag(npmt) = 1
            list_hit_pmt(nhit_pmt) = npmt
         endif

         !write(6,*) ' dscnt: ', nphoton

cc  Process photons 

         if (gate_flag(npmt) .eq. 0) then  ! ATM gate is not opened now.

             gate_flag(npmt)   =  1           ! open new gate
             now_tstart(npmt)  =  t           
             now_nphoton(npmt) =  1

             if ((sk_geometry .eq. sk_i).or.(sk_geometry .eq. sk_ii).or.(sk_geometry .eq. sk_iii)) then
                pe1 = 1.0
             else if (sk_geometry .eq. sk_iv) then
                call skrn1pe_sk3(npmt,pe1)
             else if ((sk_geometry .eq. sk_v).or.(sk_geometry .eq. sk_vi)) then
                call skrn1pe_sk5(npmt,pe1)                
             else
                write(*,*) 'dscnt: sk_geometry is not correct. stop.'
                stop
             endif
           
             if(gain_run .ne. 0) then ! gain increasing
                gain_param = pmtgain(npmt)
                pe1 = pe1 * gain_param
             endif
             now_charge(npmt) = pe1

             if (ngate .eq. MXGATE) then
              write(*,*) 'DSRES_SK4:Number of QBEE gates exceeds MXGATE'
              stop
             endif

             ngate =  ngate + 1               
             iclist(ngate)     =  npmt        
             tlist(ngate)      =  t           
             now_gate_id(npmt) =  ngate       

         else if (gate_flag(npmt) .eq. 1) then ! ATM gate is now opened.

            DT = t - now_tstart(npmt)

          if (DT .le. SIGWDS) then  ! in gate
           
             now_nphoton(npmt) = now_nphoton(npmt) + 1

             if ((sk_geometry .eq. sk_i).or.(sk_geometry .eq. sk_ii).or.(sk_geometry .eq. sk_iii)) then
                pe1 = 1.0
             else if (sk_geometry .eq. sk_iv) then
                call skrn1pe_sk3(npmt,pe1)
             else if ((sk_geometry .eq. sk_v) .or. (sk_geometry .eq. sk_vi)) then
                call skrn1pe_sk5(npmt,pe1)                
             else
                write(*,*) 'dscnt: sk_geometry is not correct. stop.'
                stop
             endif
             
             if(gain_run .ne. 0) then ! gain increasing
                gain_param = pmtgain(npmt)
                pe1 = pe1 * gain_param
             endif
             eff = dseff_qbee(dt)
             now_charge(npmt) = now_charge(npmt)+pe1*eff

          else if ( (DT .gt. SIGWDS)
     &      .and. (DT .le. DEADDS)) then  ! QBEE VETO

           ! write(6,*)' --- ignored  ', t, dt, npmt, sigwds, deadds

          else if (DT .gt. DEADDS) then  ! New gate

             gate_flag(npmt)   = 0   ! close gate
             iqlist(now_gate_id(npmt))
     &                         = now_nphoton(npmt)

             qlist(now_gate_id(npmt))
     &                         = now_charge(npmt)

             gate_flag(npmt)   = 1   ! open new gate
             now_tstart(npmt)  = t   ! gate opening timing
             now_nphoton(npmt) = 1   ! number of photons in this gate

             if ((sk_geometry .eq. sk_i).or.(sk_geometry .eq. sk_ii).or.(sk_geometry .eq. sk_iii)) then
                pe1 = 1.0
             else if (sk_geometry .eq. sk_iv) then
                call skrn1pe_sk3(npmt,pe1)
             else if ((sk_geometry .eq. sk_v) .or. (sk_geometry .eq. sk_vi)) then
                call skrn1pe_sk5(npmt,pe1)                
             else
                write(*,*) 'dscnt: sk_geometry is not correct. stop.'
                stop
             endif
             
             if(gain_run .ne. 0) then ! gain increasing
                gain_param = pmtgain(npmt)
                pe1 = pe1 * gain_param
             endif
             now_charge(npmt) = pe1

             ngate = ngate + 1       ! add new hit to list

             if (ngate .eq. MXGATE) then
              write(*,*) 'DSRES_SK4:Number of QBEE gates exceeds MXGATE'
              stop
             endif

             iclist(ngate)     = npmt
             tlist(ngate)      = t
             now_gate_id(npmt) = ngate

          endif
         endif

 100     continue

CC     close QBEE gates which are opened until the last photon.

         do 200 j = 1,nhit_pmt

            npmt = list_hit_pmt(j)

            if (gate_flag(npmt) .eq. 1) then    ! still opened gate
               gate_flag(npmt) = 0
               iqlist(now_gate_id(npmt)) = now_nphoton(npmt)
                qlist(now_gate_id(npmt)) = now_charge(npmt)
               now_nphoton(npmt) = 0
            endif

 200     continue

       if(IFMESS.eq.0) write(*,*) 'dscnt:  ngate', ngate, nhit_pmt
        
       call sortzv(tlist,sort_index,ngate,1,0,0)
       do 201 j = 1, ngate
          tlist_s(j)   =  tlist( sort_index(j)  )
          iqlist_s(j)  = iqlist( sort_index(j)  )
          qlist_s(j)   = qlist( sort_index(j)  )
          iclist_s(j)  = iclist( sort_index(j)  )
  201  continue


C ---
      RETURN
      END
