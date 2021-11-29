ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine getreflect_lin(matid,angle,lamda,ref)
c
c	(input)	matid:	material id
c					      =1:Ti window
c					      =2:stainless beam pipe 
c					      =3:stainless beam pipe arround Ti window
c			angle: 	an angle of incidence ( degree )
c			lamda:	wave length of light ( nano meter)
c
c	(output)	ref:		reflectivitiy
c
c	(purpose)  get percent reflectance of material 
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	implicit none
#include "linreflect.h"

	integer matid
	real angle,lamda,ref,wlen
	integer iangle,ilamda,od_angle,od_lamda
	real t,u,z1,z2,z3,z4
	real denomi,interp2

	wlen = lamda*1.e-3

	if(matid.eq.1)then ! if Ti window
	  call array_locate_lin(ti_angle,ti_ang_num,angle,iangle,od_angle)
	  call array_locate_lin(ti_wavelen,ti_wlen_num,wlen,ilamda,od_lamda)
	  if(od_angle.eq.-1.or.od_angle.eq.1)then 
	    if(od_lamda.eq.-1.or.od_lamda.eq.1)then
	      ref = ti_ref(ilamda,iangle)
	      return
	    else if(od_lamda.eq.0)then
	      u = 0.
	      denomi =ti_wavelen(ilamda+1)-ti_wavelen(ilamda) 
	      if(denomi.ne.0)then
	        t = (wlen-ti_wavelen(ilamda))/denomi
	      else
	        t =0.
	      endif
	      z1 = ti_ref(ilamda,iangle)
	      z2 = ti_ref(ilamda+1,iangle)
	      z3 = ti_ref(ilamda+1,iangle)
	      z4 = ti_ref(ilamda,iangle)
	      ref = interp2(z1,z2,z3,z4,t,u)
	      return
	    end if
	  else if(od_angle.eq.0)then
	     if(od_lamda.eq.-1.or.od_lamda.eq.1)then
	      t = 0.
	      denomi = ti_angle(iangle+1)-ti_angle(iangle)
	      if(denomi.ne.0)then
	        u = (angle-ti_angle(iangle))/denomi
	      else 
	        u = 0.
	      end if
	      z1 = ti_ref(ilamda,iangle)
	      z2 = ti_ref(ilamda,iangle)
	      z3 = ti_ref(ilamda,iangle+1)
	      z4 = ti_ref(ilamda,iangle+1)
	      ref = interp2(z1,z2,z3,z4,t,u)
	      return
	    else
	      denomi = ti_wavelen(ilamda+1)-ti_wavelen(ilamda)
	      if(denomi.ne.0)then
	        t = (wlen-ti_wavelen(ilamda))/denomi
	      else
	        t =0.
	      endif
	      denomi = ti_angle(iangle+1)-ti_angle(iangle)
	      if(denomi.ne.0)then
	        u = (angle-ti_angle(iangle))/denomi
	      else 
	        u = 0.
	      end if
	      z1 = ti_ref(ilamda,iangle)
	      z2 = ti_ref(ilamda+1,iangle)
	      z3 = ti_ref(ilamda+1,iangle+1)
	      z4 = ti_ref(ilamda,iangle+1)
	      ref = interp2(z1,z2,z3,z4,t,u)
	      return
	    end if
	  end if
	else if(matid.eq.2)then ! if  stainless beam pipe sample B
	  call array_locate_lin(stain_angle_b,stain_ang_num_b
     &      ,angle,iangle,od_angle)
	  call array_locate_lin(stain_wavelen_b,stain_wlen_num_b
     &      ,wlen,ilamda,od_lamda)
	  if(od_angle.eq.-1.or.od_angle.eq.1)then 
	    if(od_lamda.eq.-1.or.od_lamda.eq.1)then
	      ref = stain_ref_b(ilamda,iangle)
	      return
	    else if(od_lamda.eq.0)then
	      u = 0.
	      denomi =stain_wavelen_b(ilamda+1)
     &             -stain_wavelen_b(ilamda) 
	      if(denomi.ne.0)then
	        t = (wlen-stain_wavelen_b(ilamda))/denomi
	      else
	        t =0.
	      endif
	      z1 = stain_ref_b(ilamda,iangle)
	      z2 = stain_ref_b(ilamda+1,iangle)
	      z3 = stain_ref_b(ilamda+1,iangle)
	      z4 = stain_ref_b(ilamda,iangle)
	      ref = interp2(z1,z2,z3,z4,t,u)
	      return
	    end if
	  else if(od_angle.eq.0)then
	     if(od_lamda.eq.-1.or.od_lamda.eq.1)then
	      t = 0.
	      denomi = stain_angle_b(iangle+1)
     &              -stain_angle_b(iangle)
	      if(denomi.ne.0)then
	        u = (angle-stain_angle_b(iangle))/denomi
	      else 
	        u = 0.
	      end if
	      z1 = stain_ref_b(ilamda,iangle)
	      z2 = stain_ref_b(ilamda,iangle)
	      z3 = stain_ref_b(ilamda,iangle+1)
	      z4 = stain_ref_b(ilamda,iangle+1)
	      ref = interp2(z1,z2,z3,z4,t,u)
	      return
	    else
	      denomi = stain_wavelen_b(ilamda+1)
     &              -stain_wavelen_b(ilamda)
	      if(denomi.ne.0)then
	        t = (wlen-stain_wavelen_b(ilamda))/denomi
	      else
	        t =0.
	      endif
	      denomi = stain_angle_b(iangle+1)
     &              -stain_angle_b(iangle)
	      if(denomi.ne.0)then
	        u = (angle-stain_angle_b(iangle))/denomi
	      else 
	        u = 0.
	      end if
	      z1 = stain_ref_b(ilamda,iangle)
	      z2 = stain_ref_b(ilamda+1,iangle)
	      z3 = stain_ref_b(ilamda+1,iangle+1)
	      z4 = stain_ref_b(ilamda,iangle+1)
	      ref = interp2(z1,z2,z3,z4,t,u)
	      return
	    end if
	  end if
	else if(matid.eq.3)then ! if stainless arround Ti window
	  call array_locate_lin(stain_angle_c,stain_ang_num_c
     &      ,angle,iangle,od_angle)
	  call array_locate_lin(stain_wavelen_c,stain_wlen_num_c
     &      ,wlen,ilamda,od_lamda)
	  if(od_angle.eq.-1.or.od_angle.eq.1)then 
	    if(od_lamda.eq.-1.or.od_lamda.eq.1)then
	      ref = stain_ref_c(ilamda,iangle)
	      return
	    else if(od_lamda.eq.0)then
	      u = 0.
	      denomi =stain_wavelen_c(ilamda+1)
     &             -stain_wavelen_c(ilamda) 
	      if(denomi.ne.0)then
	        t = (wlen-stain_wavelen_c(ilamda))/denomi
	      else
	        t =0.
	      endif
	      z1 = stain_ref_c(ilamda,iangle)
	      z2 = stain_ref_c(ilamda+1,iangle)
	      z3 = stain_ref_c(ilamda+1,iangle)
	      z4 = stain_ref_c(ilamda,iangle)
	      ref = interp2(z1,z2,z3,z4,t,u)
	      return
	    end if
	  else if(od_angle.eq.0)then
	     if(od_lamda.eq.-1.or.od_lamda.eq.1)then
	      t = 0.
	      denomi = stain_angle_c(iangle+1)
     &              -stain_angle_c(iangle)
	      if(denomi.ne.0)then
	        u = (angle-stain_angle_c(iangle))/denomi
	      else 
	        u = 0.
	      end if
	      z1 = stain_ref_c(ilamda,iangle)
	      z2 = stain_ref_c(ilamda,iangle)
	      z3 = stain_ref_c(ilamda,iangle+1)
	      z4 = stain_ref_c(ilamda,iangle+1)
	      ref = interp2(z1,z2,z3,z4,t,u)
	      return
	    else
	      denomi = stain_wavelen_c(ilamda+1)
     &              -stain_wavelen_c(ilamda)
	      if(denomi.ne.0)then
	        t = (wlen-stain_wavelen_c(ilamda))/denomi
	      else
	        t =0.
	      endif
	      denomi = stain_angle_c(iangle+1)
     &              -stain_angle_c(iangle)
	      if(denomi.ne.0)then
	        u = (angle-stain_angle_c(iangle))/denomi
	      else 
	        u = 0.
	      end if
	      z1 = stain_ref_c(ilamda,iangle)
	      z2 = stain_ref_c(ilamda+1,iangle)
	      z3 = stain_ref_c(ilamda+1,iangle+1)
	      z4 = stain_ref_c(ilamda,iangle+1)
	      ref = interp2(z1,z2,z3,z4,t,u)
	    end if
	  end if
	else
	  write(*,*)'material not defined'
	end if
	return
	end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	real function interp2(z1,z2,z3,z4,t,u)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	implicit none

	real z1,z2,z3,z4,t,u

	interp2 = (1.-t)*(1.-u)*z1 + t*(1.-u)*z2 +
     &          t*u*z3 + (1.-t)*u*z4

	return 
	end

	

	

	

