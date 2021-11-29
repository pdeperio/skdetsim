cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine muminus(z,pp,h)
c
c	(input)  z   Z value of atom that mu muinus was captured
c	         pp  momentum(MeV/c)
c
c	(output) h  probability of energy spectrum
c
c	energy spectrum for the decay electron from mu minus captured
c	by atom which has Z proton
c	from Phys Lett 22(1974)119 P.Hanggi et.al
c
c	(creation date and author)
c	1997.6.24     H.Ishino
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#include "weakparm.h"
	real pp,e_e,h
	integer z
	real w0
	real x,jougen,kagen
	real espe
#ifdef DEBUG
	real rplus,rminus,delta,c_parm,a_parm,e_parm,b_parm
	real sp_a,sp_b,sp_c,sp_d,sp_e,sp_f,sp_g
#endif

	e_e = pp
	x = w0(z) - e_e
	jougen = espe(x,z,e_e)
	x = 0.
	kagen = espe(x,z,e_e)

	h = (256./pi_const)*(float(z)*fine_struct)**5*e_e*
     &    sqrt(e_e**2-e_mass**2)*
     &    (jougen - kagen)
#ifdef DEBUG
	write(*,*)'espe',espe(x,z,e_e)
	write(*,*)'rplus',rplus(z,e_e,x)
	write(*,*)'rminus',rminus(z,e_e,x)
	write(*,*)'delta',delta(z)
	write(*,*)'c_parm',c_parm(z,e_e)
	write(*,*)'a_parm',a_parm(z,e_e)
	write(*,*)'e_parm',e_parm(z,e_e)
	write(*,*)'b_parm',b_parm(e_e)
	write(*,*)'sp_a',sp_a(x,z,e_e)
	write(*,*)'sp_b',sp_b(x,z,e_e)
	write(*,*)'sp_c',sp_c(x,z,e_e)
	write(*,*)'sp_d',sp_d(x,z,e_e)
	write(*,*)'sp_e',sp_e(x,z,e_e)
	write(*,*)'sp_f',sp_f(x,z,e_e)
	write(*,*)'sp_g',sp_g(x,z,e_e)
#endif

	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	real function espe(x,z,e_e)
#include "weakparm.h"
	real x,e_e
	integer z
      real sp_a,sp_b,sp_c,sp_d,sp_e,sp_f,sp_g


	espe = sp_a(x,z,e_e)*sp_b(x,z,e_e)*sp_c(x,z,e_e)
     &     + sp_d(x,z,e_e) + sp_e(x,z,e_e) + sp_f(x,z,e_e)
     &     + sp_g(x,z,e_e)

	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	real function sp_a(x,z,e_e)
#include "weakparm.h"
	real x,e_e
	integer z
	real b_parm,delta

	sp_a = atan((2.*x - b_parm(e_e))/sqrt(delta(z))) + 
     &       atan((2.*x + b_parm(e_e))/sqrt(delta(z)))

	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	real function sp_b(x,z,e_e)
#include "weakparm.h"
	real x,e_e
	integer z
	real delta

	sp_b = 1./(delta(z)*sqrt(delta(z)))

	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	real function sp_c(x,z,e_e)
#include "weakparm.h"
	real x,e_e
	integer z
	real c_parm,a_parm,e_parm,b_parm,delta
	
	sp_c = 2.*(c_parm(z,e_e)-a_parm(z,e_e)+
     &           (2.*e_parm(z,e_e)*a_parm(z,e_e))/b_parm(e_e))
     &         /(delta(z))
     &       -e_parm(z,e_e)/b_parm(e_e)

	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	real function sp_d(x,z,e_e)
#include "weakparm.h"
	real x,e_e
	integer z
	real e_parm,b_parm,a_parm,rplus,rminus,delta

	sp_d = (e_parm(z,e_e)/b_parm(e_e))*
     &      (a_parm(z,e_e)*x - x**3)/
     &      (rplus(z,e_e,x)*rminus(z,e_e,x)*delta(z))

	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	real function sp_e(x,z,e_e)
#include "weakparm.h"
	real x,e_e
	integer z
	real p,q,r,s
	real c_parm,a_parm,e_parm,b_parm,delta,rplus,rminus

	p = (c_parm(z,e_e)-a_parm(z,e_e))/2. + 
     &     e_parm(z,e_e)*a_parm(z,e_e)/b_parm(e_e)
	q = delta(z)**2
	r = (2.*x + b_parm(e_e))/rplus(z,e_e,x)
 	s = (2.*x - b_parm(e_e))/rminus(z,e_e,x)

	sp_e = (p/q)*(r+s)

	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	real function sp_f(x,z,e_e)
#include "weakparm.h"
	real x,e_e
	integer z
	real p,q,r,s
	real c_parm,b_parm,e_parm,a_parm,delta,rplus,rminus

	p = c_parm(z,e_e)/b_parm(e_e) + 
     &    (e_parm(z,e_e)*a_parm(z,e_e))/b_parm(e_e)**2
	q = 6.*delta(z)
	r = (2.*a_parm(z,e_e) + b_parm(e_e)*x)/rplus(z,e_e,x)**2
c	s = (-2.*a_parm(z,e_e) + b_parm(e_e)*x)/rplus(z,e_e,x)**2
	s = (-2.*a_parm(z,e_e) + b_parm(e_e)*x)/rminus(z,e_e,x)**2

	sp_f = (p/q)*(r+s)

	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	real function sp_g(x,z,e_e)
#include "weakparm.h"
	real x,e_e
	integer z
	real p,q,r,s,t,u
	real b_parm,e_parm,a_parm,delta,rplus,rminus

	p = (e_parm(z,e_e)/b_parm(e_e)**2)-(1./b_parm(e_e))
	q = x**2+(2.*a_parm(z,e_e)**2+a_parm(z,e_e)*b_parm(e_e)*x)/
     &    delta(z)
	r = x**2+(2.*a_parm(z,e_e)**2-a_parm(z,e_e)*b_parm(e_e)*x)/
     &    delta(z)
	s = rplus(z,e_e,x)**2
	t = rminus(z,e_e,x)**2
	u = q/s - r/t

	sp_g = (p/6.)*u

	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	real function rplus(z,e_e,x)
#include "weakparm.h"
	real e_e,x,a_parm,b_parm
	integer z

	rplus = a_parm(z,e_e) + b_parm(e_e)*x + x*x

	return 
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	real function rminus(z,e_e,x)
#include "weakparm.h"
	real e_e,x,a_parm,b_parm
	integer z

	rminus = a_parm(z,e_e) - b_parm(e_e)*x + x*x

	return 
	end
	
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	real function delta(z)
#include "weakparm.h"
	integer z

	delta = 4.*(float(z)*mu_mass*fine_struct)**2

	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	real function c_parm(z,e_e)
#include "weakparm.h"
	real e_e,w0
	integer z

	c_parm = 3.*(w0(z)-e_e)**2

	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	real function  a_parm(z,e_e)
#include "weakparm.h"
	real e_e
	integer z

	a_parm =4.*(float(z)*mu_mass*fine_struct)**2 
     &        + e_e**2 - e_mass**2

	return 
	end 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	real function  e_parm(z,e_e)
#include "weakparm.h"
	real e_e,w0,b_parm
	integer z

	e_parm = b_parm(e_e)*(w0(z)-e_e)/e_e

	return 
	end 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	real function  b_parm(e_e)
#include "weakparm.h"
	real e_e

	b_parm = 2*sqrt(e_e**2-e_mass**2)

	return 
	end 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	real function w0(z)
#include "weakparm.h"
	integer z

	w0 = mu_mass - (mu_mass*(float(z)*fine_struct)**2)/2.

	return
	end

