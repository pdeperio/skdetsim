        subroutine sgcnt_aphoton(nopmt)
!
!       Count # of photo-electrons in a gate.
!
!       input   nopmt : PMT ID
!       output  dsmulti.h   ( OD variables )
!
!       2008/09/30 Jun Kameda for SK-IV memory-saved version
!       2008.11.18 adapted by rvw
!

      implicit none

#include "skdetsim_parm.h"
#include "geant321/gctrak.inc"
#include "geant321/gckine.inc"
#include "geant321/gcking.inc"
#include "skparm.h"
#include "skmc.h"
#include "dsdata.h"
#include "dsparm.h"
#include "dsmulti.h"

      integer nopmt
      real tofgsk,rdummy,del
      integer flag_gate(2),gate_id,igate, current_gate_id

               ! now in [ns] 
               !tofgsk = tofg * (1.e+09)
               tofgsk = tofg
               !write(6,*)'sgcnt_aphoton: begin- ', nopmt 

                if (num_gate_anti(nopmt) .eq. 0)  then ! first gate

                  num_gate_all_anti =  num_gate_all_anti + 1    ! open gate
                  num_gate_anti(nopmt) =  num_gate_anti(nopmt) +1
                  nphoto_gate_anti(num_gate_all_anti) = nphoto_gate_anti(num_gate_all_anti) + 1
                  timing_gate_anti(num_gate_all_anti) = tofgsk
                  pmt_id_anti(num_gate_all_anti) = nopmt
                  width_photons_anti(num_gate_all_anti) = 0.0
                  list_gate_anti(nopmt,num_gate_anti(nopmt)) = num_gate_all_anti
                  current_gate_id = num_gate_all_anti

                else if (num_gate_anti(nopmt) .ge. 1) then ! already exist gate

                    flag_gate(1) = 0
                    flag_gate(2) = 0

                  do igate = 1, num_gate_anti(nopmt)

                     gate_id = list_gate_anti(nopmt,igate)

CC               ADC veto will be considered in merging part

                    if (
     &                (tofgsk .gt. timing_gate_anti(gate_id)) .and.
     &                (tofgsk .le. (timing_gate_anti(gate_id)+SIGWDS) )
     &                 ) then 
                      if (flag_gate(1) .eq. 0) then 
                         flag_gate(1) = gate_id
                      else
                      write(*,*) 'SGPMT:Someting wrong',nopmt,gate_id,tofgsk
     &                           ,timing_gate_anti(gate_id),num_gate_anti(nopmt)
                      endif
                     endif


                     if ( (tofgsk .le. timing_gate_anti(gate_id)) .and.
     &                (timing_gate_anti(gate_id) .le. (tofgsk+SIGWDS) )
     &                  ) then 
                       if (flag_gate(2) .eq. 0) then 
                          flag_gate(2) = gate_id
                       else
                        write(*,*) 'SGPMT_SK4:Someting wrong',nopmt,gate_id
                       endif
                      endif

                    enddo

                  if (flag_gate(1) .ne. 0) then  

                     gate_id = flag_gate(1)
                     current_gate_id = gate_id
                     nphoto_gate_anti(gate_id) = nphoto_gate_anti(gate_id) + 1
                     if ( (timing_gate_anti(gate_id)-tofgsk) .ge. width_photons_anti(gate_id))
     &                  width_photons_anti(gate_id)  = abs(tofgsk - timing_gate_anti(gate_id))


                  else if (flag_gate(2) .ne. 0) then                       

                     gate_id = flag_gate(2)
                     current_gate_id = gate_id
                     nphoto_gate_anti(gate_id) = nphoto_gate_anti(gate_id) + 1
                     width_photons_anti(gate_id)  = width_photons_anti(gate_id) + 
     &                         abs(tofgsk - timing_gate_anti(gate_id))
                     timing_gate_anti(gate_id) = tofgsk
 
                  else  ! open new gate

                     num_gate_all_anti =  num_gate_all_anti + 1  ! open gate
                     pmt_id_anti(num_gate_all_anti) = nopmt
                     num_gate_anti(nopmt) =  num_gate_anti(nopmt) +1
                     current_gate_id = num_gate_all_anti 
                     nphoto_gate_anti(num_gate_all_anti) = nphoto_gate_anti(num_gate_all_anti) + 1
                     timing_gate_anti(num_gate_all_anti) = tofgsk
                     width_photons_anti(num_gate_all_anti) = 0.0
                     list_gate_anti(nopmt,num_gate_anti(nopmt)) = num_gate_all_anti
 
                   endif

                 endif


          !write(6,*) 'sgcnt_aphoton: end- ', num_gate_all_anti, timing_gate_anti(gate_id)
      return
      end
