c	Igloo 3:6:3 - Equal Latitude   V-1.0  	Last modified: 29-6-98
c
c 	This contains subroutines for the 3:6:3 Equal Latitude
c	pixelization.  There are subroutines to go between three 
c	different coordinates: theta-phi, azimuthal pixel number 
c	and hierarchical pixel number.   To create and transform 
c	temperature maps, the temperatures are store in the 
c	azimuthal ordering.  For coarsening and nearest neighbors(?), 
c 	the hierarchical ordering should be used.  

c	This program uses six central routines:   
c		ang2pix_azim(iazim,ncoars,theta,phi)
c			: Theta-Phi => Azimuthal pixel number 
c		pix2ang_azim(iazim,ncoars,theta,phi)
c			: Azimuthal pixel number => Theta-Phi 
c		ang2pix_hier(ihier,ncoars,theta,phi)
c			: Theta-Phi => Hierarchical pixel number 
c		pix2ang_hier(ihier,ncoars,theta,phi)
c			: Hierarchical pixel number => Theta-Phi 
c		azim2hier(iazim,ihier,ncoars)
c			: Azimuthal pixel number => Hierarchical 
c		hier2azim(ihier,iazim,ncoars)
c			: Hierarchical pixel number => Azimuthal

c	There are also two subroutines to initialize the tree arrays, 
c	which should be globally defined and common blocked as they are 
c	in the driver program:
c		init_tree 
c		init_inv_tree 
c	These need to be intialized prior to running the hierarchical 
c	routines, but need to be run only once. 

c	NCOARS is the variable which defines the resolution of the 
c	pixelization.  It can range from 0 (12 pixels) to 12 which is 
c	12 x 4^12 = 200,000,000 pixels.  While this can easily be increased, 
c	this should be sufficient for our purposes.   

c	The azimuthal base pixel numbering is arranged as follows: 
c	________________________________________________   
c	|		|		|		|  North Pole
c	|	0	|	1	|	2	|
c	________________________________________________
c	|	|	|	|	|	|	|
c	|   3	|   4	|   5	|   6	|   7	|   8	|  Equator
c	________________________________________________
c	|		|		|		|
c	|	9	|	10	|	11	|
c	________________________________________________   South Pole		
c	
c	I've since thought of a slightly better way to do it, but it's 
c	probably not worth coding it up. 
c 
c	I offer no guarantees for this software.  
c	I've tested it somewhat.  It seems reasonably fast compared 
c	to other algorithms I've seen.  Obviously, it still needs to 
c	be fully tested, and quite likely the routines can be improved.
c	Let me know if you have any suggestions or find any bugs. 
c	   
c 			Rob Crittenden  29-6-98
c			R.G.Crittenden@damtp.cam.ac.uk


      subroutine pix2ang_azim(ipix, ncoars, theta, phi) 
      implicit none
      real*4 theta, phi, twopi, thirdpi, dtheta 
      parameter (twopi = 6.2831853) 
      parameter (thirdpi = 1.0471976) 
      integer*4 ipix, irow, icol, izone, ipixbar
      integer*4 ncoars, ntot, nside, nsidesq, nbelow, nrow, nquart  

      nside = ncoars 
      ntot = 12*nside*nside 
      dtheta = thirdpi/float(nside) 
      nquart = ntot/4
      
      if (ipix .ge. nquart .and. ipix .lt. 3*nquart) then 	!Equator 
        nbelow = nquart
	irow = nside 
	nrow = 6*nside

      elseif (ipix .lt. nquart .and. ipix .ge. 0) then    	!North cap 
	nbelow = ntot/16 
        irow = nside/2
        do while (nbelow .gt. ipix) 
	  nbelow = nbelow/4 
          irow = irow/2
        enddo 
        nrow = 9*irow 
        if (nrow .eq. 0) nrow = 3 
      
      elseif (ipix .ge. 3*nquart .and. ipix .lt. ntot) then	!South cap
	nbelow = ntot/16 
        irow = nside/2
        ipixbar = ntot - ipix - 1
        do while (nbelow .gt. ipixbar) 
	  nbelow = nbelow/4 
          irow = irow/2
        enddo 
        nrow = 9*irow 
        nbelow = ntot - nbelow*4
        irow = 3*nside - irow*2 
        if (nrow .eq. 0) then  
	  nrow = 3 
          irow = 3*nside - 1 
          nbelow = ntot - 3
        endif 

      else 
        Pause 'Pixel out of range - Azimuthal'  
      endif 

      izone = ipix - nbelow 
      irow = irow + izone/nrow 
      icol = mod(izone,nrow) 
      theta = (float(irow) + 0.5)*dtheta
      phi = (float(icol) + 0.5)*twopi/float(nrow) 
	
      return 
      end



      subroutine ang2pix_azim(ipix, ncoars, theta, phi) 
      implicit none
      real*4 theta, phi, twopi, thirdpi, dtheta 
      parameter (twopi = 6.2831853) 
      parameter (thirdpi = 1.0471976) 
      integer*4 ipix, ncoars, ntot, nside, nsidesq 
      integer*4 irow, icol, nbelow, nrow, ibelow, irowbar 

      nside = ncoars !6
      ntot = 12*nside*nside !432
      dtheta = thirdpi/float(nside) !0.17453293333
      irow = int(theta/dtheta)      !12

      if (irow .ge. nside .and. irow .lt. 2*nside) then 	!Equator 
        nbelow = ntot/4
	ibelow = nside 
	nrow = 6*nside

      elseif (irow .ge. 0 .and. irow .lt. nside) then  		!North cap 
        ibelow = nside/2
        do while (ibelow .gt. irow) 
          ibelow = ibelow/2
        enddo 
        nrow = 9*ibelow
	nbelow = ibelow*ibelow*3 
        if (ibelow .eq. 0) nrow = 3 
      
      elseif (irow .ge. 2*nside .and. irow .lt. 3*nside) then	!South cap
        ibelow = nside/2 !3
        irowbar = 3*nside - irow -1 !5
        do while (ibelow .gt. irowbar) 
          ibelow = ibelow/2
        enddo 
        nrow = 9*ibelow !27
        nbelow = ntot - ibelow*ibelow*12 !432 - 3*3*12=216
        ibelow = 3*nside - ibelow*2 !12
        if (nrow .eq. 0) then  
	  nrow = 3 
          ibelow = ibelow - 1 
          nbelow = ntot - 3
        endif 

      else 
        Pause 'Angle out of range - Azimuthal'  
      endif

      icol = mod(int(float(nrow)*(phi/twopi)),nrow)
      ipix = nbelow + (irow - ibelow)*nrow + icol 

      return 
      end 


      subroutine init_tree
      integer*4 it, ip, i 
      integer*4 treet(0:255), treep(0:255) 
      common /tree/ treet, treep
      real twopi, thirdpi
      common/piconst/ twopi, thirdpi

      twopi = 8.*atan(1.)
      thirdpi = twopi/6.

      do i = 0, 255
        ip =  iand(i,85)		
        it = iand(i,170)
        treep(i) = mod(ip,2)+2*mod(ip/4,2)+4*mod(ip/16,2)+8*mod(ip/64,2)
        treet(i) = mod(it,4)/2+mod(it/4,4)+2*mod(it/16,4)+4*mod(it/64,4)
      enddo 
      write(*,*) 'Pixel tree initialized'
      return 
      end 


      subroutine init_inv_tree
      integer*4 it, ip, i, mult, sum 
      integer*4 invtree(0:15,0:15)
      common /inv_tree/ invtree 
      do it = 0,15 
        do ip = 0,15 
          mult = 1
          sum = 0
          do i = 0,4 
            sum = sum + (2*ibits(it,i,1) + ibits(ip,i,1))*mult 
            mult = mult*4
	  enddo 
          invtree(it,ip) = sum
        enddo 
      enddo 
      return 
      end


      subroutine pix2ang_hier(ipix, ncoars, theta, phi)
      implicit none
      real*4 theta, phi, twopi, thirdpi, dtheta 
      parameter (twopi = 6.2831853) 
      parameter (thirdpi = 1.0471976) 
      integer*4 ipix, ncoars, ntot, nside 
      integer*4 irow, icol, ndiv, nsidesq, nrow
      integer*4 word1, word2, word3, iremain, iquad, ibase
      integer*4 treet(0:255), treep(0:255) 
      common /tree/ treet, treep

      nside = ncoars !6
      nsidesq = nside*nside  !36
      dtheta = thirdpi/float(nside) !0.174532933...

      iremain = mod(ipix,nsidesq) !0
      ibase = ipix/nsidesq !12
	
      if (ibase .ge. 3 .and. ibase .lt. 9) then 	!Equator 
        word1 = mod(iremain, 256) 
        word2 = mod(iremain/256,256) 
        word3 = iremain/65536
        irow = 256*treet(word3) + 16*treet(word2) + treet(word1) 
        icol = 256*treep(word3) + 16*treep(word2) + treep(word1) 
      
        irow = irow + nside 
        icol = icol + nside*(ibase - 3)
        nrow = 6*nside
      
      else 
        if (iremain .eq. 0) then 			!Polar pixel
          nrow = 3 
          irow = 0 
          icol = mod(ibase,3) !0
        else  						!Others
	  ndiv = nside/2
          iquad = iremain/ndiv/ndiv  
          do while (iquad .eq. 0) 
            ndiv= ndiv/2 
            iquad = iremain/ndiv/ndiv  
          enddo 
          iremain = mod(iremain, ndiv*ndiv)
          word1 = mod(iremain, 256) 
          word2 = mod(iremain/256,256) 
          word3 = iremain/65536
          irow = 256*treet(word3) + 16*treet(word2) + treet(word1) 
          icol = 256*treep(word3) + 16*treep(word2) + treep(word1) 
          irow = ndiv + irow
          icol = icol + mod(ibase,3)*ndiv*3 + (iquad-1)*ndiv
          nrow = 9*ndiv
        endif  

	if (ibase .ge. 0 .and. ibase .lt. 3) then	!North cap 
        elseif (ibase .ge. 9 .and. ibase .lt. 12) then 	!South cap 
          irow = 3*nside - irow - 1
	else
          Pause 'Pixel out of range - Heirarchical'
	endif 

      endif 

      theta = (float(irow) + 0.5)*dtheta
      phi = (float(icol) + 0.5)*twopi/float(nrow) 
      
      return
      end 


      subroutine ang2pix_hier(ipix, ncoars, theta, phi) 
      implicit none
      real*4 theta, phi, twopi, thirdpi, dtheta 
C      parameter (twopi = 6.2831853) 
C      parameter (thirdpi = 1.0471976) 
      integer*4 ipix, ncoars, ntot, nside, ibase, iquad, irowbar 
      integer*4 irow, icol, nbelow, ndiv, ibelow, nsidesq 
      integer*4 invtree(0:15,0:15)
      common /inv_tree/ invtree 
      common/piconst/ twopi, thirdpi

      nside = ncoars 
      ntot = 12*nside*nside 
      dtheta = thirdpi/float(nside) 
      
      irow = theta/dtheta
      if (irow.ge.3*nside) then
         print *,'ang2pix_hier: theta is too large: ',theta
         print *,'ang2pix_hier: truncating theta'
      endif
      irow = max(min(irow,3*nside-1),0)  !Avoid round off errors

      if (irow .ge. nside .and. irow .lt. 2*nside) then 	!Equator 
        icol = mod(int(float(nside)*phi/thirdpi),6*nside)
        ibase = 3 + icol/nside 
	icol = mod(icol,nside) 
        irow = irow - nside 
        ipix = 65536*invtree(irow/256,icol/256) 
        ipix = ipix + 256*invtree(mod(irow/16,16),mod(icol/16,16)) 
	ipix = ipix + invtree(mod(irow,16),mod(icol,16))

      elseif (irow .lt. nside .and. irow .ge. 0) then    	!North cap 
        ibelow = nside/2
        do while (ibelow .gt. irow) 
          ibelow = ibelow/2
        enddo 
        if (ibelow .eq. 0)  then 
          ipix = 0
          ibase = mod(int(3.*phi/twopi),3)
        else 
          irow = irow - ibelow
          icol = mod(int(9.*float(ibelow)*(phi/twopi)),9*ibelow)
          ibase = icol/ibelow/3 
          iquad = mod(icol/ibelow,3) + 1
          icol = mod(icol,ibelow)
          ipix = 65536*invtree(irow/256,icol/256) 
          ipix = ipix + 256*invtree(mod(irow/16,16),mod(icol/16,16)) 
      	  ipix = ipix + invtree(mod(irow,16),mod(icol,16))
          ipix = ipix + ibelow*ibelow*iquad
        endif 
      
      elseif (irow .ge. 2*nside .and. irow .lt. 3*nside) then	!South cap
        ibelow = nside/2 !3
        irowbar = 3*nside - irow -1 !5
        do while (ibelow .gt. irowbar) 
          ibelow = ibelow/2
        enddo 
        if (ibelow .eq. 0)  then 
          ipix = 0
          ibase = mod(int(3.*phi/twopi),3) + 9
        else 
          irow = irowbar - ibelow !2
          icol = mod(int(9.*float(ibelow)*(phi/twopi)),9*ibelow) !25
          ibase = icol/ibelow/3 + 9 !11
          iquad = mod(icol/ibelow,3) + 1 !3
          icol = mod(icol,ibelow) !1
          ipix = 65536*invtree(irow/256,icol/256) 
          ipix = ipix + 256*invtree(mod(irow/16,16),mod(icol/16,16)) 
      	  ipix = ipix + invtree(mod(irow,16),mod(icol,16))
          ipix = ipix + ibelow*ibelow*iquad
        endif 

      else 
        print *,theta,phi
        Pause 'Angle out of range - Hierarchical'  
      endif

      ipix = ipix + nside*nside*ibase  

      return 
      end

      subroutine azim2hier(iazim, ihier, ncoars) 
      implicit none
      integer*4 iazim, ihier, nbelow, ibelow, iazimbar
      integer*4 ipix, ncoars, ntot, nside 
      integer*4 irow, icol, ndiv, nsidesq, nrow, nquart
      integer*4 word1, word2, word3, iremain, iquad, ibase
      integer*4 invtree(0:15,0:15)
      common /inv_tree/ invtree 

      nside = ncoars 
      ntot = 12*nside*nside 
      nquart = 3*nside*nside
      
      if (iazim .ge. nquart .and. iazim .lt. 3*nquart) then 	!Equator 
        iremain = iazim - nquart
	irow = iremain/6/nside  
	icol = mod(iremain,nside)
        ibase = 3 + mod(iremain/nside,6) 
        ipix = 65536*invtree(irow/256,icol/256) 
        ipix = ipix + 256*invtree(mod(irow/16,16),mod(icol/16,16)) 
	ipix = ipix + invtree(mod(irow,16),mod(icol,16))

      elseif (iazim .lt. nquart .and. iazim .ge. 0) then   	!North cap 
        if (iazim .lt. 3)  then 
          ipix = 0
          ibase = iazim
        else 
          nbelow = ntot/16 
          ibelow = nside/2
          do while (nbelow .gt. iazim) 
	    nbelow = nbelow/4 
            ibelow = ibelow/2
          enddo 
          ndiv = 3*ibelow 
          iremain = iazim - nbelow
          irow =  iremain/3/ndiv
	  icol = mod(iremain,ibelow)
          ipix = 65536*invtree(irow/256,icol/256) 
          ipix = ipix + 256*invtree(mod(irow/16,16),mod(icol/16,16)) 
      	  ipix = ipix + invtree(mod(irow,16),mod(icol,16))
          icol = mod(iremain,3*ndiv)
          ibase = icol/ndiv
          iquad = mod(icol/ibelow,3) + 1
          ipix = ipix + ibelow*ibelow*iquad
	endif
      
      elseif (iazim .ge. 3*nquart .and. iazim .lt. ntot) then	!South cap
        iazimbar = ntot - 1 -iazim
        if (iazimbar .lt. 3)  then 
          ipix = 0
          ibase = 9+  mod(iazim ,3)
        else 
          nbelow = ntot/16 
          ibelow = nside/2
          do while (nbelow .gt. iazimbar) 
	    nbelow = nbelow/4 
            ibelow = ibelow/2
          enddo 
          ndiv = 3*ibelow 
          iremain =  iazim + 4*nbelow - ntot  
          irow =  (iazimbar - nbelow)/3/ndiv
	  icol = mod(iremain,ibelow)
          ipix = 65536*invtree(irow/256,icol/256) 
          ipix = ipix + 256*invtree(mod(irow/16,16),mod(icol/16,16)) 
      	  ipix = ipix + invtree(mod(irow,16),mod(icol,16))
          icol = mod(iremain,3*ndiv)
          ibase = icol/ndiv + 9
          iquad = mod(icol/ibelow,3) + 1
          ipix = ipix + ibelow*ibelow*iquad
	endif
      
      else 
        Pause 'Pixel out of range - Azimuthal'  
      endif 

      ihier = ipix + nside*nside*ibase  
      return 
      end 


      subroutine hier2azim(ihier, iazim, ncoars)
      implicit none
      integer*4 ihier, iazim, nbelow, ibelow
      integer*4 ipix, ncoars, ntot, nside 
      integer*4 irow, icol, ndiv, nsidesq, nrow
      integer*4 word1, word2, word3, iremain, iquad, ibase
      integer*4 treet(0:255), treep(0:255) 
      common /tree/ treet, treep

      nside = ncoars 
      nsidesq = nside*nside 
      ntot = 12*nsidesq

      iremain = mod(ihier,nsidesq)
      ibase = ihier/nsidesq 
	
      if (ibase .ge. 3 .and. ibase .lt. 9) then 	!Equator 
        word1 = mod(iremain, 256) 
        word2 = mod(iremain/256,256) 
        word3 = iremain/65536
        irow = 256*treet(word3) + 16*treet(word2) + treet(word1) 
        icol = 256*treep(word3) + 16*treep(word2) + treep(word1) 
      
        icol = icol + nside*(ibase - 3)
	nbelow = 3*nsidesq
	ndiv = 6*nside
      
      else 
        if (iremain .eq. 0) then 			!Polar pixel
          nbelow = 0 
	  ndiv = 3 
          irow = 0 
          icol = mod(ibase,3)
        else  						!Others
	  ndiv = nside/2
          iquad = iremain/ndiv/ndiv  
          do while (iquad .eq. 0) 
	    ndiv= ndiv/2 
            iquad = iremain/ndiv/ndiv  
          enddo 
          iremain = mod(iremain, ndiv*ndiv)
          word1 = mod(iremain, 256) 
          word2 = mod(iremain/256,256) 
          word3 = iremain/65536
          irow = 256*treet(word3) + 16*treet(word2) + treet(word1) 
          icol = 256*treep(word3) + 16*treep(word2) + treep(word1) 
          icol = icol + mod(ibase,3)*ndiv*3 + (iquad-1)*ndiv
	  nbelow = ndiv*ndiv*3 
	  ndiv = ndiv*9
        endif  

	if (ibase .ge. 0 .and. ibase .lt. 3) then	!North cap 
        elseif (ibase .ge. 9 .and. ibase .lt. 12) then 	!South cap 
	  nbelow = ntot - 4*nbelow
          irow = ndiv/9 - irow -1
	else
          Pause 'Pixel out of range - Heirarchical'
	endif 

      endif 

      iazim = nbelow + irow*ndiv + icol 
      
      return
      end 
