	subroutine sgcnt_photon(NOPMT)
c
c	Count # of photo-electrons in a gate.
c
c       input   NOPMT : PMT ID
c       output  dsmulti.h
c
c       2008/09/30 Jun Kameda for SK-IV memory-saved version
c

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

      integer NOPMT
      real tofgsk,rdummy,del
      integer flag_gate(2),gate_id,igate, current_gate_id

               TOFGSK = TOFG * (1.E+09)

cc  Comment out. 
cc               call rngausk(100.,10.,del)               
cc               CALL GRNDM(RDUMMY,1)
cc               if(RDUMMY.lt.0.02) tofgsk=tofgsk+del

                if (num_gate(NOPMT) .eq. 0)  then ! first gate

                  num_gate_all =  num_gate_all + 1    ! open gate
                  num_gate(NOPMT) =  num_gate(NOPMT) +1
                  nphoto_gate(num_gate_all) = nphoto_gate(num_gate_all) + 1
                  timing_gate(num_gate_all) = TOFGSK
                  pmt_id(num_gate_all) = NOPMT
                  width_photons(num_gate_all) = 0.0
                  list_gate(NOPMT,num_gate(NOPMT)) = num_gate_all
                  current_gate_id = num_gate_all

                else if (num_gate(NOPMT) .ge. 1) then ! already exist gate

                    flag_gate(1) = 0
                    flag_gate(2) = 0

                  do igate = 1, num_gate(NOPMT)

                     gate_id = list_gate(NOPMT,igate)

CC               ADC veto will be considered in merging part

                    if (
     &                (TOFGSK .gt. timing_gate(gate_id)) .and.
     &                (TOFGSK .le. (timing_gate(gate_id)+SIGWDS) )
     &                 ) then 
                      if (flag_gate(1) .eq. 0) then 
                         flag_gate(1) = gate_id
                      else
                      write(*,*) 'SGPMT:Someting wrong',NOPMT,gate_id,TOFGSK
     &  ,timing_gate(gate_id),num_gate(NOPMT)
                      endif
                     endif


                     if ( (TOFGSK .le. timing_gate(gate_id)) .and.
     &                (timing_gate(gate_id) .le. (TOFGSK+SIGWDS))
     &                  ) then 
                       if (flag_gate(2) .eq. 0) then 
                          flag_gate(2) = gate_id
                       else
                        write(*,*) 'SGPMT_SK4:Someting wrong',NOPMT,gate_id
                       endif
                      endif

                    enddo

                  if (flag_gate(1) .ne. 0) then  

                   gate_id = flag_gate(1)
                   current_gate_id = gate_id
                   nphoto_gate(gate_id) = nphoto_gate(gate_id) + 1
                   if ( (timing_gate(gate_id)-TOFGSK) .ge. width_photons(gate_id))
     &             width_photons(gate_id)  = abs(TOFGSK - timing_gate(gate_id))


                  else if (flag_gate(2) .ne. 0) then                       

                   gate_id = flag_gate(2)
                   current_gate_id = gate_id
                   nphoto_gate(gate_id) = nphoto_gate(gate_id) + 1
                   width_photons(gate_id)  = width_photons(gate_id) + 
     &                         abs(TOFGSK - timing_gate(gate_id))
                   timing_gate(gate_id) = TOFGSK
 
                   else  ! open new gate

                num_gate_all =  num_gate_all + 1  ! open gate
                pmt_id(num_gate_all) = NOPMT
                num_gate(NOPMT) =  num_gate(NOPMT) +1
                current_gate_id = num_gate_all 
                nphoto_gate(num_gate_all) = nphoto_gate(num_gate_all) + 1
                timing_gate(num_gate_all) = TOFGSK
                width_photons(num_gate_all) = 0.0
                list_gate(NOPMT,num_gate(NOPMT)) = num_gate_all
 
                   endif

                 endif


      return
      end
