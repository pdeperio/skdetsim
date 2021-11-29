ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c	linreflect.h
c
c	This header file define the reflectivity of the Ti window and
c	stainless beam pipe as the function of wave length of cherenkov
c	light and angle of incidence of that.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	integer ti_ang_num,ti_wlen_num
	parameter(ti_ang_num=8)
	parameter(ti_wlen_num=15)
	
	integer stain_ang_num_b,stain_wlen_num_b
	parameter(stain_ang_num_b=7)
	parameter(stain_wlen_num_b=15)

	integer stain_ang_num_c,stain_wlen_num_c
	parameter(stain_ang_num_c=8)
	parameter(stain_wlen_num_c=15)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	ti_angle angle for Ti window
c	stain_angle_b stainless angle for sample B(beam pipe component)
c	stain_angle_c stainless angle for sample C(arround the Ti window)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	real ti_angle(ti_ang_num)
     &    ,stain_angle_b(stain_ang_num_b)
     &    ,stain_angle_c(stain_ang_num_c)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	ti_wavelen wavelength for Ti window
c	stain_wavelen_b stainless wavelength for sample B(beam pipe component)
c	stain_wavelen_c stainless wavelength for sample C(arround the Ti window)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	real ti_wavelen(ti_wlen_num)
     &    ,stain_wavelen_b(stain_wlen_num_b)
     &    ,stain_wavelen_c(stain_wlen_num_c)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	ti_ref reflectivity for Ti window
c	stain_ref_b stainless reflectivity for sample B(beam pipe component)
c	stain_ref_c stainless reflectivity for sample C(arround the Ti window)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	real ti_ref(ti_wlen_num,ti_ang_num)
     &    ,stain_ref_b(stain_wlen_num_b,stain_ang_num_b)
     &    ,stain_ref_c(stain_wlen_num_c,stain_ang_num_c)

      data ti_angle /20.,30.,40.,50.,60.,70.,75.,80./
      data stain_angle_b /20.,30.,40.,50.,60.,70.,75./
      data stain_angle_c /20.,30.,40.,50.,60.,70.,75.,80./

      data ti_wavelen /0.30,0.36,0.37,0.38,0.40
     &                ,0.42,0.44,0.45,0.46,0.47
     &                ,0.48,0.50,0.53,0.54,0.60/
      data stain_wavelen_b /0.30,0.36,0.37,0.38,0.40
     &                ,0.42,0.44,0.45,0.46,0.47
     &                ,0.48,0.50,0.53,0.54,0.60/
      data stain_wavelen_c /0.30,0.36,0.37,0.38,0.40
     &                ,0.42,0.44,0.45,0.46,0.47
     &                ,0.48,0.50,0.53,0.54,0.60/

      data ti_ref /38.10,38.20,38.50,39.60,39.00,39.90,41.70
     &      ,42.60,42.80,45.60,47.80,48.40,48.60,49.80,49.90
     &      ,38.00,38.00,39.20,39.80,39.40,40.20,41.50,42.90
     &      ,43.70,43.70,44.70,47.70,47.67,47.66,47.60
     &      ,38.20,38.30,39.30,40.20,39.80,41.80,43.90,44.30
     &      ,44.40,44.50,44.60,45.10,45.73,45.76,47.50
     &      ,39.00,39.10,40.40,41.60,41.90,42.80,43.10,43.10  
     &      ,43.10,43.10,44.50,46.40,46.30,46.42,46.80
     &      ,41.60,42.75,43.90,45.20,45.40,44.80,45.90,46.90
     &      ,46.92,46.94,46.96,47.00,46.97,46.96,46.90
     &      ,40.00,40.10,44.60,44.80,45.40,44.40,45.13,45.06
     &      ,45.00,45.15,45.30,45.60,46.20,46.40,47.30
     &      ,45.00,45.60,46.10,46.60,47.10,47.20,47.66,47.90
     &      ,47.96,48.03,48.10,48.40,49.00,49.20,50.20
     &      ,51.70,51.40,51.46,51.53,51.00,51.20,52.10,52.10
     &      ,52.10,52.00,52.16,52.50,53.80,53.80,53.80/

      data stain_ref_b/
     &       51.00,51.10,54.20,57.20,57.00,55.40,58.12,58.25
     &      ,58.37,58.50,58.66,59.00,59.48,59.64,60.10
     &      ,51.10,51.20,51.90,54.70,57.00,57.60,58.00,57.90
     &      ,58.00,58.10,58.03,57.90,58.56,58.78,59.90
     &      ,51.00,51.17,51.20,52.30,53.40,58.00,58.40,58.48
     &      ,58.56,58.65,58.73,58.90,59.14,59.22,59.70
     &      ,51.10,51.18,51.20,51.80,53.90,53.90,56.30,56.40
     &      ,57.15,57.90,58.00,58.20,59.52,59.65,60.30
     &      ,52.60,52.68,52.70,53.60,53.50,55.50,56.70,57.00
     &      ,57.30,57.30,57.50,58.00,59.20,59.30,59.60
     &      ,55.40,55.48,55.50,54.40,57.10,57.00,58.70,58.66
     &      ,58.63,58.60,58.66,58.80,59.04,59.12,59.60
     &      ,57.70,57.78,57.80,58.40,59.40,59.80,59.55,59.90
     &      ,60.13,60.36,60.60,60.50,60.50,60.70,60.50/

      data stain_ref_c/
     &       52.60,54.40,58.30,58.15,57.70,56.60,58.50,58.50
     &      ,58.70,59.50,59.63,59.90,60.20,60.65,61.70
     &      ,52.10,55.60,58.90,59.20,57.00,56.80,58.25,58.40
     &      ,59.00,59.10,59.33,59.80,59.80,59.70,61.30
     &      ,52.00,53.20,58.40,58.30,56.70,58.60,58.40,58.50
     &      ,58.80,59.10,59.20,59.40,60.03,60.06,61.10
     &      ,49.90,50.20,53.40,54.50,54.30,57.90,59.35,59.80
     &      ,59.95,60.10,59.80,60.40,60.64,60.72,61.00
     &      ,50.70,51.10,51.60,52.10,53.10,54.20,57.10,57.40
     &      ,57.50,57.80,57.80,57.80,58.40,58.50,59.10
     &      ,51.60,52.20,52.90,52.96,53.10,54.36,55.30,55.60
     &      ,55.70,55.80,56.03,56.50,57.12,57.24,57.90
     &      ,52.00,52.10,53.30,53.30,53.30,54.60,55.46,55.73
     &      ,56.00,56.80,56.80,56.80,57.70,57.70,57.70
     &      ,52.60,52.70,53.40,53.66,54.20,56.50,58.20,58.26
     &      ,58.33,58.40,58.50,58.70,58.85,58.90,59.20/

