************************************************************************
*     -----------------
      SUBROUTINE SGGVCD
*     -----------------
*
*     (Purpose)
*        Generate Vector
*
************************************************************************

      implicit none

      integer imvtx,imdir,ipart,immom
      real rvtx(3),rwall,rdir(3),rmom(2),POS(3),PMOM(3),rdummy,pp
      INTEGER RNDCARD(5)
      COMMON/VECTCARD/IMVTX,RVTX,RWALL,IMDIR,RDIR,IPART,IMMOM,RMOM,RNDCARD

C --  GENERATE VECTOR
C

C --  SET ID,POSITION,MOMMENTUM
C
      IF(IMVTX.EQ.1) THEN
         POS(1)=RVTX(1)
         POS(2)=RVTX(2)
         POS(3)=RVTX(3)
      ELSE
         IF(IMVTX.EQ.2) CALL RNPOS(POS,RWALL)
      ENDIF

      IF(IMDIR.EQ.1) THEN
         PMOM(1)=RDIR(1)
         PMOM(2)=RDIR(2)
         PMOM(3)=RDIR(3)
      ELSE
         IF(IMDIR.EQ.2) CALL RNDIR(PMOM)
      ENDIF

      IF(IMMOM.EQ.1) THEN
         PP = RMOM(1)
      ELSE
C         IF(IMMOM.EQ.2) PP = RMOM(1)+RNDM()*(RMOM(2)-RMOM(1))
         CALL GRNDM(RDUMMY,1)
         IF(IMMOM.EQ.2) PP = RMOM(1)+RDUMMY*(RMOM(2)-RMOM(1))
      ENDIF
      PMOM(1)=PMOM(1) * PP
      PMOM(2)=PMOM(2) * PP
      PMOM(3)=PMOM(3) * PP
C
C --     SET PARAMETER TO COMMON 
C

      CALL VCMAKE(POS,PMOM,IPART)

      RETURN
      END


