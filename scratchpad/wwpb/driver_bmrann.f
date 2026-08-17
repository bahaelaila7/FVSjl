      PROGRAM DRV
C     Golden generator for WWPB bmrann.f (MINSTD a=16807 m=2^31-1,
C     SEL=BMS1/2^31, default seed 55329). Seed via BMRNSD(.TRUE.,55329)
C     to mirror the BLOCK DATA default (BMS0=55329, odd), then draw 8.
      IMPLICIT NONE
      REAL SEL
      INTEGER I
      CALL BMRNSD(.TRUE., 55329.0)
      DO 10 I=1,8
        CALL BMRANN(SEL)
        WRITE(*,'(Z8.8)') SEL
   10 CONTINUE
      END
