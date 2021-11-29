	program txtzbs
*
*  Convert a USMC/Nuance text file of vectors to a ZBS file
*
*  D.Casper 5 Dec 1999
*
	implicit none
#include "dsbnkf.h"
	integer luno, ierr, nzbs, kz, iend
	integer ierr1, ierr2
	integer nout
        parameter (luno = 20)

        parameter (nzbs = 1500000)
        common /kzbs/kz(nzbs)


        call kzlogl(-2)
        call kzlimit(nzbs)

        call skopenf(10,1,'F',ierr1)
        call skopenf(20,1,'Z',ierr2)
	lin = 10
	iend = 0
	nout = 0
	do while (iend.eq.0)
		call sggtxt(iend)
		if (iend.eq.0) then
C			print *,'writing'
			call kzwrit(20)
			nout = nout + 1
		endif
	enddo

	call skclosef(10)
	call skclosef(20)
	print *,nout,' events written'

	end

