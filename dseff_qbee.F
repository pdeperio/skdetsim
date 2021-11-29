      real function dseff_qbee(tdiff)
*
*     input  :
*            tdiff : time from gate start time (nsec)
*     output :
*            dseff_qbee : efficiency
*            The efficiency is normalized to 1.0 at t=0.
*            Leak of the charge is alread corrected by measured 1p.e. distrubtion.
*
      implicit none
      
      integer j
      real tdiff,t_to_end

      real tbl_eff(22),tbl_t(22)
      data tbl_t/0. ,15.5038757 ,21.7054272 ,31.0077515 ,46.511631 ,
     & 62.0155029 ,68.2170563 ,74.418602 ,86.8217087 ,96.1240311 ,
     & 111.627914 ,117.829453 ,145.736435 ,164.34108 ,186.046524 ,
     & 217.054276 ,248.062012 ,282.170532 ,319.379852 ,350.387604 ,
     & 372.093048 ,400./

      data tbl_eff/0. ,0.00272665289 ,0.0059986366 ,0.0256305337 ,0.112883449 ,
     & 0.227402836 ,0.275937289 ,0.325017005 ,0.422085851 ,0.489161521 ,
     & 0.58868438 ,0.625766814 ,0.750920236 ,0.813087881 ,0.868439019 ,
     & 0.917518735 ,0.950238585 ,0.971233785 ,0.984321713 ,0.992501676 ,
     & 0.996319056 ,1. /


      t_to_end = 400.0 - tdiff
      dseff_qbee = 0.0

      if ( ( t_to_end .lt. 0) .or. ( t_to_end .gt. 400.0)) return
      if ( (t_to_end .lt. tbl_t(1)) .or.  (t_to_end .gt. tbl_t(22)) ) return

      if ( t_to_end .eq. tbl_t(22)) then
         dseff_qbee = 1.00000
*         dseff_qbee = dseff_qbee * 1.0 / 0.987469018
         return
      endif

      do j = 2,22
         if ( (t_to_end .ge. tbl_t(j-1) ) .and. (t_to_end .lt. tbl_t(j) ) ) then

           dseff_qbee = tbl_eff(j)*(t_to_end-tbl_t(j-1)) +
     &                  tbl_eff(j-1)*(tbl_t(j)-t_to_end)
           dseff_qbee = dseff_qbee / (tbl_t(j)-tbl_t(j-1))
*           dseff_qbee = dseff_qbee * 1.0 / 0.987469018

           return
         endif
      enddo


      return
      end



      
