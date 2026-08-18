      PROGRAM DRV
      IMPLICIT NONE
      REAL SEL
      INTEGER I
      DO 10 I=1,6
        CALL TMRANN(SEL)
        WRITE(*,'(Z8.8)') SEL
   10 CONTINUE
      END
