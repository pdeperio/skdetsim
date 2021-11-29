ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine array_locate_lin(array,array_num,x,ix,od_flag)
c
c	(input)	array   :array (real)
c	            array_num:dimension of array
c	            x       :value
c	            
c	(output)	ix      :index which satisfy the relation as follow,
c	                     array(ix) < x < array(ix+1)
c	            od_flag :-1  if x < array(1) ;under flow
c	                    :0   if array(1) < x <array(array_num)
c	                    :1   if array(array_num) < x ;over flow
c	(purpos)    get the index 'ix' 
c	            
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	implicit none

	integer array_num
	real array(array_num), x
	integer ix,od_flag
	integer i,imax,imin,imidle


** check whether array was sorted or not
	do i=1,array_num-1
	  if(array(i).gt.array(i+1))then
	    write(*,*)'array is not sorted!!!!!'
	    stop
	  end if
	end do

	od_flag = 0

	if(x.ge.array(array_num))then
	  od_flag = 1
	  ix = array_num
	  return
	end if

	if(x.lt.array(1))then
	  od_flag = -1
	  ix = 1
	  return
	end if

	imax = array_num+1
	imin = 0

	do while((imax-imin).gt.1)
	  imidle = int(float(imax+imin)/2.)
	  if(x.ge.array(imidle))then
	    imin = imidle
	  else 
	    imax = imidle
	  end if
	end do

	ix = imin

	return
	end
	  
	
