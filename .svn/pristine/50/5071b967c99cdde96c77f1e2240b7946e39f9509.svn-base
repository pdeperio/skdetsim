************************************************************************
*     -------------------------------------
      SUBROUTINE DSPRSV
*     -------------------------------------
*
*     (Purpose)
*        Save /DSPARM/ into 'MCPARMCONV' bank
*
*     (Input)
*        
*
*     (Output)
*       
*       
*       
*
*     (Creation Date and Author)
*       1992.05.21 ; First version by K.S.Hirata
*
************************************************************************

#include "skparm.h"
#include "dsparm.h"

      CHARACTER*32 CNAME
      DIMENSION RBUF(12)

C ---
       
      CNAME = 'MCPARMCONV'
      ISEG=0
      NDATA=12

      CALL KZBCR1(CNAME,NDATA,'-F',0,0,' ',IERR) 
       
      RBUF(1) = DARKDS
      RBUF(2) = TCNVSK
      RBUF(3) = QCNVSK
      RBUF(4) = TRESDS
      RBUF(5) = QRESDS
      RBUF(6) = DTHRSK
      RBUF(7) = TWINDS
      RBUF(8) = TRIGDS
      RBUF(9) = GATEDS
      RBUF(10) = BEFTDS
      RBUF(11) = DEADDS
      RBUF(12) = SIGWDS

      CALL KZREP1(CNAME,0,NDATA,RBUF)

C ---
      RETURN
      END
