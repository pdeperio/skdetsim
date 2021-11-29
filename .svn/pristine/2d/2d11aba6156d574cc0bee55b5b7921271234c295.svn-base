ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	dtginfo.h
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	integer dtg_num
c	parameter(dtg_num=20000)
c        parameter(dtg_num=40000)
	parameter(dtg_num=80000)
	integer dtg_mode(dtg_num),dtg_trg(dtg_num)
	real dtg_x(dtg_num),dtg_y(dtg_num),dtg_z(dtg_num)
	integer dtg_badrun(dtg_num)
	integer dtg_runnum
        real dtg_dark(dtg_num)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c	dtg_trg    0    normal run (2m pull-up)
c	           1    trigger run (no pull-ups)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	common/dtginfo/dtg_mode,dtg_x,dtg_y,dtg_z,
     &      dtg_badrun,dtg_trg,dtg_dark
	common/dtginfox/dtg_runnum
