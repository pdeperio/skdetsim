*******************************************************************
*
      SUBROUTINE ROTATEVECT(PVECIN,PVECOUT,COSP,SINP,COST,SINT)
*     -----------------------------------------------------------



      IMPLICIT NONE

      REAL PVECIN(3),PVECOUT(3)
      REAL COSP,SINP,COST,SINT

      REAL TEMPV(3)

      TEMPV(1)=COSP*PVECIN(1)-SINP*PVECIN(2)
      TEMPV(2)=SINP*PVECIN(1)+COSP*PVECIN(2)
      TEMPV(3)=PVECIN(3)

      PVECOUT(1)=COST*TEMPV(1)+SINT*TEMPV(3)
      PVECOUT(2)=TEMPV(2)
      PVECOUT(3)=-SINT*TEMPV(1)+COST*TEMPV(3)

      RETURN
      END
