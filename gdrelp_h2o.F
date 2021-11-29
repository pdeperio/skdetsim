CDECK  ID>, GDRELP.
*CMZ :  3.21/02 29/03/94  15.41.21  by  S.Giani
*-- Author :
      SUBROUTINE GDRELP_H2O(A,Z,DENS,T,DEDX)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *  Calculates the mean 1/DENS*dE/dx of a  MUON  with kinetic     *
C.    *  energy T in an element of atomic number Z, atomic weight A    *
C.    *  and density DENS ( the density is just used for the           *
C.    *  calculation of the density effect in the case of high T).     *
C.    *                                                                *
C.    *   ------------------------------------------------------       *
C.    *   It calls GDRELP with the substitution EMMU -- > PMASS !      *
C.    *    ( and put back the right value of PMASS after)       !      *
C.    *   -----------------------------------------------------        *
C.    *                                                                *
C.    *    ==>Called by : GDRELA                                       *
C.    *       Author    L.Urban    *********                           *
C.    *                                                                *
C.    ******************************************************************
C.
      DOUBLE PRECISION PI,TWOPI,PIBY2,DEGRAD,RADDEG,CLIGHT,BIG,EMASS
      DOUBLE PRECISION EMMU,PMASS,AVO
*
      PARAMETER (PI=3.14159265358979324D0)
      PARAMETER (TWOPI=6.28318530717958648D0)
      PARAMETER (PIBY2=1.57079632679489662D0)
      PARAMETER (DEGRAD=0.0174532925199432958D0)
      PARAMETER (RADDEG=57.2957795130823209D0)
      PARAMETER (CLIGHT=29979245800.D0)
      PARAMETER (BIG=10000000000.D0)
      PARAMETER (EMASS=0.0005109990615D0)
      PARAMETER (EMMU=0.105658387D0)
      PARAMETER (PMASS=0.9382723128D0)
      PARAMETER (AVO=0.60221367D0)
*
*
      CALL GDRELX_H2O(A,Z,DENS,T,REAL(PMASS),DEDX)
*
      END
