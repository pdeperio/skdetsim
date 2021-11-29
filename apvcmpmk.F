**********************************************************************
*     ----------------------------------------
      SUBROUTINE APVCMPMK(NVCT,POS,PMOM,TVTX,IP)
*     ----------------------------------------
*     
*     ( purpose )
*       VECTOR GENERATION ( MULTI PARTICLE VERSION )
*
*     ( input ) 
*       NVCT          : NUMBER OF PARTICLE < MAXVC ( 101 )
*       POS(3,MAXVC)  : VERTEX POSITION 
*       PMOM(3,MAXVC) : MOMENTUM (MEV/C)
*       TVTX(MAXVC)   : TIME ( SEC )      
*       IP(MAXVC)     : PARICLE CODE
*
*     ( output ) 
*       COMMON VCWORK,VCVRTX
*
*     ( creation date and author )
*       1993.11.25 ; First version by K.Kaneyuki
*
**********************************************************************
#include "vcwork.h"
#include "vcvrtx.h"
      INTEGER IP(MAXVC)
      REAL POS(3,MAXVC),PMOM(3,MAXVC),TVTX(MAXVC)

C
C --  CLEAR COMMON
C
      CALL VCCLCM
C
C --  MAKE HEADER
C
      CALL MCMKHD
      CALL MCMKMH
      CALL MCMKWT
C
C --  SET PARAMETER TO VCWORK
C
      NVC=NVCT
C
C --  GENERAL VERTEX = POSVC(I,1)
C      
      CALL UCOPY(POS,POSVC,3)

      DO 10 I=1,NVCT
         IPVC(I)=IP(I)
         CALL MCMASS(IP(I),AMASVC(I))
         CALL UCOPY(PMOM(1,I),PVC(1,I),3)
         IORGVC(I)=0
         IFLGVC(I)=0
         ICRNVC(I)=1
         TIMVC(I)=TVTX(I)
         CALL UCOPY(POS(1,I),POSIVC(1,I),3)
ccc         IVTIVC(I)=1
         IVTIVC(I)=i
 10   CONTINUE
C     
C --  SET PARAMETER TO VCVRTX
C
      NVTXVC=NVCT
      DO 20 I=1,NVCT   
         CALL UCOPY(POS(1,I),PVTXVC(1,I),3)
         IFLVVC(I)=1
         IPARVC(I)=0
         TIMVVC(I)=TVTX(I)
 20   CONTINUE
      
      RETURN
      END







