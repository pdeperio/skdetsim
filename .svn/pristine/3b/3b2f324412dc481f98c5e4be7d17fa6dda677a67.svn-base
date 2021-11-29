************************************************************************
*     -------------------------------------
      SUBROUTINE DSBDEL
*     -------------------------------------
*
*     (Purpose)
*       Delete the banks specified in "dsbnkw".
*
*     (Input)
*
*
*     (Creation Date and Author)
*       1993.04.07 ; First version by K.S.Hirata
*
************************************************************************

#include "dsbnkf.h"
      
C --- Delete the specified banks

      CALL KZBLOC('MCPHOTON',LBIR)
      IF(NBNKDS(1).EQ.0.AND.LBIR.NE.0) CALL KZBDEL('MCPHOTON')
      CALL KZBLOC('MCPHOTAN',LBIR)
      IF(NBNKDS(1).EQ.0.AND.LBIR.NE.0) CALL KZBDEL('MCPHOTAN')
      CALL KZBLOC('MCRAWTQ',LBIR)
      IF(NBNKDS(2).EQ.0.AND.LBIR.NE.0) CALL KZBDEL('MCRAWTQ')
      CALL KZBLOC('MCRAWTQA',LBIR)
      IF(NBNKDS(3).EQ.0.AND.LBIR.NE.0) CALL KZBDEL('MCRAWTQA')
      CALL KZBLOC('MCTQ',LBIR)
      IF(NBNKDS(4).EQ.0.AND.LBIR.NE.0) CALL KZBDEL('MCTQ')
      CALL KZBLOC('TQ',LBIR)
      IF(NBNKDS(5).EQ.0.AND.LBIR.NE.0) CALL KZBDEL('TQ')
      CALL KZBLOC('TQA',LBIR)
      IF(NBNKDS(6).EQ.0.AND.LBIR.NE.0) CALL KZBDEL('TQA')
      CALL KZBLOC('TQREAL',LBIR)
      IF(NBNKDS(7).EQ.0.AND.LBIR.NE.0) CALL KZBDEL('TQREAL')
      CALL KZBLOC('TQAREAL',LBIR)
      IF(NBNKDS(8).EQ.0.AND.LBIR.NE.0) CALL KZBDEL('TQAREAL')
      CALL KZBLOC('TQZ',LBIR)
      IF(NBNKDS(9).EQ.0.AND.LBIR.NE.0) CALL KZBDEL('TQZ')
!     CALL KZBLOC('TQAZ',LBIR)
!     IF(NBNKDS(9).EQ.0.AND.LBIR.NE.0) CALL KZBDEL('TQAZ')
       

C ---

      RETURN
      END
