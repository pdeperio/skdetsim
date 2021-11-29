      subroutine packprt
      implicit none
*
*  Pack the list of MC particles for segment data
*
#include "aptrkseg.h"
      integer ibuf(20), np, iseg, ndata, i, j, lbank, parent_row, nrow_mcvect
      real rbuf(20), mass, vmod
      equivalence (ibuf,rbuf)
      logical nu_evt, tau_evt, tau_dk
      integer nue_code, num_code, nut_code, tau_code, proton_code, neutron_code,
     &        oxygen_code
      parameter (nue_code = 12)
      parameter (num_code = 14)
      parameter (nut_code = 16)
      parameter (tau_code = 15)
      parameter (proton_code = 2212)
      parameter (neutron_code = 2112)
      parameter (oxygen_code = 8016)
C
C+ entry point
C
      npart = 0
      iseg = 0
      call kzget1('MCVECT',iseg,ndata,ibuf)
      np = ibuf(1) ! Number of particles in MCVECT bank
      nrow_mcvect = np
      nu_evt = .false.
      tau_evt = .false.
      tau_dk = .false.
      do i = 1, np
         call kzget1('MCVECT',i,ndata,ibuf)
         if (ndata.gt.0 .and. npart.lt.mxpart) then
            npart = npart + 1
            ip_part(npart) = ibuf(1)
            if (i.eq.1 .and. (abs(ibuf(1)).eq.nue_code.or.
     &                        abs(ibuf(1)).eq.num_code.or.
     &                        abs(ibuf(1)).eq.nut_code).and.
     &          ibuf(8).eq.0) then   ! First particle is an informational neutrino
               nu_evt = .true.
            endif
            if (nu_evt .and. (i.eq.2) .and. (((ibuf(1).ne.proton_code) .and.
     &                                       (ibuf(1).ne.neutron_code).and.
     &                                       (ibuf(1).ne.oxygen_code)).or.
     &                                       ibuf(8).eq.1)) then
               nu_evt = .false. ! Second particle is not an informational hadronic target
            endif 
            call ucopy(rbuf(3),p_part(1,npart),3)
            p_part(4,npart) = rbuf(2)
            p_part(5,npart) = vmod(rbuf(3),3)
            if (rbuf(3).ne.0 .or. rbuf(4).ne.0) then
               ph_part(npart) = atan2(rbuf(4),rbuf(3))
            else
               ph_part(npart) = 0.
            endif
            if (p_part(5,npart).gt.0) then
               th_part(npart) = acos(rbuf(5)/p_part(5,npart))
            else
               th_part(npart) = 0.
            endif
            c_part(npart) = 0 ! For non-nu event, primary particles are class 0
            if (nu_evt) then
               if (i.le.2) then
                  c_part(npart) = 0
               else if (i.eq.3) then
                  c_part(npart) = 1
                  if (abs(ip_part(npart)).eq.tau_code) tau_evt = .true.
               else
                  c_part(npart) = 2 !Assume hadronic 
                  if (tau_evt) then
                     if (abs(ip_part(npart)).eq.nut_code) then
                        if (tau_dk) then ! Oops, marked some particles incorrectly
                           do j = 1, npart-1
                              if (c_part(j).eq.3) c_part(j) = 2
                           enddo
                        endif
                        tau_dk = .true.
                     endif
                     if (tau_dk) then
                        c_part(npart) = 3
                     endif
                  endif ! tau lepton produced by neutrino
               endif ! 3 particle or subsequent
            endif ! neutrino event
         endif !track exists
      enddo

      call kzbloc('SCNDPRT',lbank)
      if (lbank.le.0) return
      call kzget1('SCNDPRT',0,ndata,ibuf)
      np = ibuf(1)
      do i = 1, np
         call kzget1('SCNDPRT',i,ndata,ibuf)
         if (ndata.gt.0 .and. npart.lt.mxpart) then
            npart = npart + 1
            ip_part(npart) = ibuf(2)
            call ucopy(rbuf(6), p_part(1,npart),3)
            p_part(5,npart) = vmod(rbuf(6),3)
            call mcmass(ip_part(npart),mass)
            p_part(4,npart) = sqrt(p_part(5,npart)**2 + mass**2)
            if (rbuf(6).ne.0 .or. rbuf(7).ne.0) then
               ph_part(npart) = atan2(rbuf(7),rbuf(6))
            else
               ph_part(npart) = 0.
            endif
            if (p_part(5,npart).gt.0) then
               th_part(npart) = acos(rbuf(8)/p_part(5,npart))
            else
               th_part(npart) = 0.
            endif
            c_part(npart) = 1 ! for non-neutrino events, secondaries are class 1
            if (nu_evt) then
               parent_row = ibuf(12)
               if (parent_row.gt.nrow_mcvect) then
                  write(*,*) 'seg_packprt: parent row is outside MCVECT range'
                  stop
               endif
               if (ibuf(13).ne.ip_part(parent_row)) then
                  write(*,*) 'seg_packprt: parent code is ',ip_part(parent_row),
     &                       ' but expected ',ibuf(13)
                  stop
               endif
               c_part(npart) = c_part(parent_row) + 2
            endif ! neutrino event
         endif
      enddo

C
C+  Last but not least, store an informational track representing the total momentum and net hadronic vector
C
      if (nu_evt) then
         npart = npart + 1
         ip_part(npart) = 0
         call vadd(p_part(1,1),p_part(1,2),p_part(1,npart),4)
         p_part(5,npart) = vmod(p_part(1,npart),3)
         if (p_part(1,npart).ne.0 .or. p_part(2,npart).ne.0) then
            ph_part(npart) = atan2(p_part(2,npart),p_part(1,npart))
         else
            ph_part(npart) = 0.
         endif
         if (p_part(5,npart).gt.0) then
            th_part(npart) = acos(p_part(3,npart)/p_part(5,npart))
         else
            th_part(npart) = 0.
         endif
         c_part(npart) = 6  ! Total momentum

         npart = npart + 1
         ip_part(npart) = 0
         call vsub(p_part(1,npart-1),p_part(1,3),p_part(1,npart),4)
         p_part(5,npart) = vmod(p_part(1,npart),3)
         if (p_part(1,npart).ne.0 .or. p_part(2,npart).ne.0) then
            ph_part(npart) = atan2(p_part(2,npart),p_part(1,npart))
         else
            ph_part(npart) = 0.
         endif
         if (p_part(5,npart).gt.0) then
            th_part(npart) = acos(p_part(3,npart)/p_part(5,npart))
         else
            th_part(npart) = 0.
         endif
         c_part(npart) = 7 ! Hadronic momentum
      endif

      return
      end
