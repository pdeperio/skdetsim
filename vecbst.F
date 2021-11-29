      SUBROUTINE VECBST(P,RM,EV,GAM)
C     S. ODAKA    11-OCT-'79
C     ******************************************************************
C     *    TO BOOST ANY MOMENTUM VECTOR 'P' OF THE PARTICLE WITH THE   *
C     * MASS 'RM' TO THE DIRECTION 'EV'(UNIT VECTOR) WITH THE GIVEN    *
C     * GAMMA FACTOR 'GAM' ,AND RESPOND THE RESULT INTO THE VECTOR 'P' *
C     * WHICH IS THE SAME ARREY AS THE GIVEN VECTOR TO BE BOOSTED      *
C     ******************************************************************
      DIMENSION P(3),EV(3)
C     
      PP=P(1)*EV(1)+P(2)*EV(2)+P(3)*EV(3)
      ES=P(1)*P(1)+P(2)*P(2)+P(3)*P(3)+RM*RM
      CE=(GAM-1.)*PP+SQRT((GAM*GAM-1.)*ES)
C     
      DO 10 I=1,3
         P(I)=CE*EV(I)+P(I)
   10 CONTINUE
C     
      RETURN
      END
