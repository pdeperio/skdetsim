c Maximilien Fechner -- reflection tuning for ATMPD
c March 2007

      real at(5)
      real tuned_bs(5)
      real tuned_pm(5)      
      integer ntuning
      common/refatmpd/at,tuned_bs,tuned_pm,ntuning

c at : wavelengths of the laser injector (4 in SK-3 , 5 previously)
c tuned_bs : black sheet relfectivity, tuned values 
c tuned_pm : pmt reflectivity, tuned values
c ntuning = 4 or 5, important for linear interpolation purposes 

c PREVIOUS values (before march 2007)
c Black sheet -- from rfbssg
c      dimension tuned(5),at(5)
c      data tuned/0.045,0.045,0.055,0.060,0.055/
c      data at/365.,400.,420.,460.,500./

c PMT reflectivity -- from rfpmssg
c     data tuned/0.467,0.570,0.650,1.133,1.067/  
c     data at/365.,400.,420.,460.,500./

c The new values are set in DSINIP.F
c The main lsrsim routine cam override these values
