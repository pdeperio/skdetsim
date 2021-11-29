*     DOME PARAMETERS
*     Note, must include "geotnk.h" before this header

      INTEGER npar_dome_cyl
      REAL dome_cyl_z
      REAL dome_cyl_halfz
      REAL dome_cyl_centerz
      REAL dome_cyl_center(3), dome_cyl_par(3)

      INTEGER npar_dome_sphe
      REAL dome_sphe_radius
      REAL dome_sphe_zpos
      REAL dome_sphe_par(6), dome_sphe_center(3)
      
      PARAMETER( npar_dome_cyl = 3 )
      PARAMETER( dome_cyl_z = 3200 )
      PARAMETER( dome_cyl_halfz = dome_cyl_z/2.0 )
      PARAMETER( dome_cyl_centerz = ZPTKTK+dome_cyl_halfz )
      DATA dome_cyl_par/0.0,RTKTK,dome_cyl_halfz/ ! /rmin, rmax, half-height/
      DATA dome_cyl_center/0.0,0.0,dome_cyl_centerz/

      PARAMETER( npar_dome_sphe = 6 )
      PARAMETER( dome_sphe_radius = 2208.84375 ) ! from sqrt(1200.0**2+RTKTK**2)/(2400.0) )
      PARAMETER( dome_sphe_zpos = -1*dome_sphe_radius )
      DATA dome_sphe_par/0.0,dome_sphe_radius, 0., 90.0, 0.0, 360.0/ ! /rmin, rmax, the1, the2, phi1, phi2/
      DATA dome_sphe_center/0.0,0.0,dome_sphe_zpos/ ! origin set in the middle of the dome cylinder
