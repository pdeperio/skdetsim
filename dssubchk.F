      subroutine dssubchk(more, nmue)
*******************************************************
*  Search the current event for decay electrons and generate
*  sub-events for them.  
*  
*     Output:  more    .true. if sub-events remain to be processed,
*                      .false. otherwise.
*
*  D.Casper 4-08-98
*******************************************************
      implicit none
#include "skparm.h"
C
      logical more
      integer nmue
      integer ier
C
      real gate_width
      parameter (gate_width = 950.)
      integer convvect_goodlen
      parameter (convvect_goodlen = 11)
      integer max_particles
      parameter (max_particles = 100)
      real big_number
      parameter (big_number = 1.e30)
      integer electron_code, photon_code
      parameter (electron_code = 11)
      parameter (photon_code   = 22)
      real p_cut
      parameter (p_cut = 5.)  !Minimum secondary particle energy to consider
C
      integer ierror
      integer event_number, run_number
      integer current_subevt, total_subevt
      real tgate, tnext
      real convvect_data(convvect_goodlen)
      integer convvect_ihdr(convvect_goodlen)
      integer convvect_idata(convvect_goodlen)
      equivalence (convvect_data, convvect_idata)
      integer convvect_len, convvect_hdrlen
      integer np, ip, n_2nd
      real vtxadd(3,max_particles), padd(3,max_particles),
     &     tadd(max_particles)
      integer iadd(max_particles)
      real wallsk, vmod
C
      data tgate /gate_width/
      
C
      more = .false.
      nmue = 0
C Find the current sub-event and the total number.  These should agree.

      call kzenum(total_subevt,current_subevt)
      if (current_subevt.ne.total_subevt) then
         write(6,*)'DSSUBCHK:evts disagree:',current_subevt,total_subevt
         tgate = gate_width
         return
      endif

C Look for the first particle not yet treated 

      call kzefst
      tnext = big_number
      call kzget1('CONVVECT',0,convvect_hdrlen,convvect_ihdr)
      if (convvect_hdrlen.gt.0) then
         n_2nd = convvect_ihdr(1)
         do ip = 1, n_2nd
            call kzget1('CONVVECT',ip,convvect_len,convvect_idata)
            if (convvect_len.eq.convvect_goodlen) then
C  Do we care about this particle?
C               if ((convvect_data(9).lt.tnext) .and.
C     &                (convvect_data(9).gt.tgate) .and.
                if (
     &                ((abs(convvect_idata(2)).eq.electron_code .or.
     &                 abs(convvect_idata(2)).eq.photon_code)) .and.
     &                (vmod(convvect_data(6),3).ge.p_cut) .and.
     &                (wallsk(convvect_data(3)).gt.0.)) then
C                     print *,'dssubchk: found something'
                     if (abs(convvect_idata(2)).eq.electron_code .and. 
     &                   convvect_data(9).gt.300.) then
C                        print *,vmod(convvect_data(6),3),convvect_data(9)
                        nmue = nmue + 1
                     endif
                     if (convvect_data(9).lt.tnext .and.
     &                 convvect_data(9).gt.tgate) then
                        tnext = convvect_data(9)
                     endif
               endif  !particle is candidate for subevent?
            else  !CONVVECT not valid
               write(6,*) 'DSSUBCHK: invalid CONVVECT bank'
            endif  !CONVVECT is valid
         enddo  !Particles in subevt
      endif !Subevt has CONVVECT

C  Restore current event

      call kzecur(current_subevt,ierror)
      if (ierror.ne.0) then
         write (6,*)'DSSUBCHK: error restoring current event'
         tgate = gate_width
         return
      endif

C  Did we find any particles to process?

      if (tnext.lt.big_number) then !Yes, 
         more = .true.
         tgate = tnext
      else                          !No, reset tgate
         tgate = gate_width
      endif
      return

      entry dssubvc(IER)
****************************************************
*  Create a new event with all particles within gate_width
*  of tgate, and update tgate.
*
*  IER returns 0 if successful.
*  The program must stop if IER.NE.0
****************************************************
C  get CONVVECT bank for secondary particle data
      ier = -1
      call kzefst
      call kzget1('CONVVECT',0,convvect_hdrlen,convvect_ihdr)
      if (convvect_hdrlen.le.0) then
         write(6,*) 'DSSUBVC: no CONVVECT bank'
         return
      endif
      n_2nd = convvect_ihdr(1)

      np = 0
      do ip = 1, n_2nd
         call kzget1('CONVVECT',ip,convvect_len, convvect_idata)
         if (convvect_len .ne. convvect_goodlen) then
            write(6,*) 'DSSUBVC: invalid CONVVECT segment length'
            return
         endif
         if (convvect_data(9).ge.tgate .and.
     &       convvect_data(9).le.tgate+gate_width) then
            if (np.lt.max_particles) then
              np = np+1
              iadd(np) = convvect_idata(2)
              call ucopy(convvect_data(3),vtxadd(1,np),3)
              call ucopy(convvect_data(6),padd(1,np),3)
              tadd(np) = convvect_data(9) * (1.e-9)
            else if (np.eq.max_particles) then
               write(6,*) 'DSSUBVC: warning - too many secondaries'
            endif
         endif
      enddo

C Update tgate

      tgate = tgate + gate_width

C Create a new event

      call kzelst
      call kzenew
      call kzcrmh
      call kzenum(total_subevt,current_subevt)
      call apvcmpmk(np,vtxadd,padd,tadd,iadd)
      call vcmkvx
      call vcmkvc
C
      ier = 0
      return

      entry dssubcl
**************************************************
*  Provided for safety to explicitly reset tgate
**************************************************
      tgate = gate_width

      return
      end
