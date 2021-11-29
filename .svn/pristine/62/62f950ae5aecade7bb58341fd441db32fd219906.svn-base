************************************************************************:

      subroutine dssatu_sk4(q_true, q_obs)

*
*     Consider saturation effect of PMT and QBEE in SK-IV.
*     Input 
*             Q_true : True Q
*     Output :
*             Q_obs  : Q after saturation effect
*
*     2009/02/27: Jun Kameda
*                 based on study by Y.Yokosawa & S.Moriyama
*          03/04: Modify original function to fix discontinuity
*                 at 110.0p.e.s
*

      implicit none
      real q_true, q_obs
      real par(5)
      data par/1.0009,-1.49716E-05,-2.83942E-07,4.6018E-10,-2.38401E-13/

      if ( q_true .le. 110.0) then  ! no effect

         q_obs = q_true

      elseif ( (q_true .gt. 110.0) .and. (q_true .lt. 1135.19324))  then  

         q_obs = par(1)*q_true+par(2)*q_true**2+par(3)*q_true**3
     &          +par(4)*q_true**4+par(5)*q_true**5

         q_obs = q_obs*1.003618  ! to adjust at 110.0

      elseif ( q_true .ge. 1135.19324) then

         q_obs = 1020.0

      endif

      return

      end
