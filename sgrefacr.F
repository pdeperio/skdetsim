************************************************************************
*     ----------------------------------------------------
      SUBROUTINE SGREFACR(IENTER,WLEN,THETA1,SREF,PREF,THETA2)
*     ----------------------------------------------------
*
*     (Purpose)
*       Calculate reflectibity on the surface of acrylic/water and the angle
*        of refraction
*
*     (Input)
*       IENTER     : flag (1: water to acryl, 0: acryl to water,
*                      otherwise: error)
*       WLEN   : wave length of the light (nm)
*       THETA1 : input angle (0 to pi/2) = reflection angle
*                  incident angle with respect to normal
*
*     (Output)
*       SREF   : reflection as S-polarization
*       PREF   : reflection as P-polarization
*       THETA2 : refraction angle (0 to pi/2)
*
*          If it is total reflection,
*             SREF = 1.
*             PREF = 1.
*             THETA2 = 999.
*
*     (Comment)
*       1) In the case of linear polarization, the total probability is
*            PREF*cos(THETA1)**2 + SREF*sin(THETA1)**2
*
*       2) We don't consider the change of polarization
*          after total reflection yet.
*
*       3) We don't consider the refraction as complex number yet
*
*
*     (Creation Date and Author)
*       2002.08.01 ; First version by R.Nambu (modified by Y.Koshio, C.Sterner)
*
************************************************************************

      IMPLICIT NONE

C REFRACTION INDEX
      real acrnsg
      real REFSG

C INPUT:
      REAL WLEN,THETA1
      INTEGER IENTER
C OUTPUT:
      REAL SREF,PREF,THETA2
C LOCAL VARIABLE:
      REAL N1,N2
      REAL C_ANGL
C BEGIN:

C FROM WATER TO ACRYLIC
      IF ( IENTER .EQ. 1 ) THEN
         N1 = REFSG(WLEN)
         N2 = acrnsg(WLEN)
         IF(N2.LT.N1) THEN
            C_ANGL = ASIN(N2/N1)
            IF ( THETA1 .GE. C_ANGL ) THEN
               SREF = 1.
               PREF = 1. 
               THETA2 = 999.
            ELSE
               CALL CAL_ACRYL(N1,N2,THETA1,SREF,PREF,THETA2)
            ENDIF
         ELSE
            CALL CAL_ACRYL(N1,N2,THETA1,SREF,PREF,THETA2)
         ENDIF
C FROM ACRYLIC TO WATER
      ELSE IF ( IENTER .EQ. 0 ) THEN
         N1 = acrnsg(WLEN)
         N2 = REFSG(WLEN)
         IF(N2.LT.N1) THEN
            C_ANGL = ASIN(N2/N1)
            IF ( THETA1 .GE. C_ANGL ) THEN
               SREF = 1.
               PREF = 1. 
               THETA2 = 999.
            ELSE

               CALL CAL_ACRYL(N1,N2,THETA1,SREF,PREF,THETA2)
            END IF
         ELSE
            CALL CAL_ACRYL(N1,N2,THETA1,SREF,PREF,THETA2)
         ENDIF
      ELSE
         WRITE(*,*) 'FLAG "IENTER" ERROR in sgrefacr.F'
         stop
      END IF
      RETURN
      END

C*************CALCULATION*************************
      SUBROUTINE CAL_ACRYL(N1,N2,THETA1,SREF,PREF,THETA2)
      IMPLICIT NONE
C INPUT 
      REAL THETA1
C OUTPUT
      REAL THETA2
      REAL SREF,PREF
C LOCAL VARIABLE:
      REAL N1,N2
      REAL COS1,SIN1,COS2
C FUNCTION
      REAL SREFLECT,PREFLECT
C BEGIN
      COS1 = COS(THETA1)
      SIN1 = SIN(THETA1)
      COS2 = SQRT(1-(N1/N2*SIN1)**2)
      THETA2 = ACOS(COS2)
      SREF = SREFLECT(N1,N2,COS1,COS2)
      PREF = PREFLECT(N1,N2,COS1,COS2)
      RETURN
      END

C**** REFLECTIVITY OF S_POLARIZATION *****
      FUNCTION SREFLECT(N1,N2,C1,C2)
C INPUT
      REAL N1,N2,C1,C2
C BEGIN
      SREFLECT = ((N1*C1-N2*C2)/(N1*C1+N2*C2))**2
      RETURN
      END

C**** REFLECTIVITY OF P_POLARIZATION *****
      FUNCTION PREFLECT(N1,N2,C1,C2)
C INPUT
      REAL N1,N2,C1,C2
C BEGIN
      IF((C1.EQ.0).OR.(C2.EQ.0)) THEN
         PREFLECT = 1
      ELSE
         PREFLECT = ((N1/C1-N2/C2)/(N1/C1+N2/C2))**2
      END IF
      RETURN
      END
