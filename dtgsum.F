cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine dtgsum
c
c	(input) none
c
c	(output) none
c
c	(purpos) fill the run information to dtginfo.h
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	implicit none

#include "dtginfo.h"

	character*100 fname
	integer run,mode,trg
	real x,y,z
	integer badrun
	integer i

cc	fname = '/home/jcravens/lowe/dtg/fitter/RefTab/'//
cc     &        'dtg_sk2_runsum.dat'
cc	fname = '/home/sklowe/dtg/const/'//
cc     &          'dtg_sk3_runsum.dat'
cc     &        'dtg_sk2_runsum.dat.noshift'
cc      fname = '/home/sekiya/skrep/skdetsim_dtg/'//
cc     &        'dtg_sk3_runsum.dat'

cc      fname = '/home/sekiya/skrep/'//
      fname = '/home/sito/skrep/'//
     &        'dtg_sk4_runsum.dat'

	open(10,file=fname,status = 'old',form='formatted')

	do i=1,dtg_num
	  dtg_mode(i) = -1
	  dtg_trg(i) = -1
	  dtg_x(i) = -20000.
	  dtg_y(i) = -20000.
	  dtg_z(i) = -20000.
	  dtg_badrun(i) = -1
	end do
 
	i = 0
 2	read(10,*,end=200)run,mode,trg,x,y,z,badrun

cc	print*,'dtgsum',x,y,z

	i = i+ 1
	dtg_mode(run) = mode
	dtg_trg(run) = trg
	dtg_x(run) = x
	dtg_y(run) = y
	dtg_z(run) = z
	dtg_badrun(run) = badrun
	goto 2
 200	continue
	dtg_runnum = i
	close(10)

	return
	end

