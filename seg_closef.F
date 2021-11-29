      subroutine seg_closef
      implicit none
#include "aptrkseg.h"

      integer icycle

      data icycle /0/

      call hrout(0,icycle,'NT')
      call hrend('segments')

      return
      end
