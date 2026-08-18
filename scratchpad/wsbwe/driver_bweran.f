      PROGRAM DRV
C     Golden generator for WSBWE bweran.f.
C     BWERAN LCG: S1=DMOD(16807*S0,2147483647); SEL=S1/2^31; S0=S1.
C     DATA default S0=55329 (odd). Draw 8, print Float32 hex.
C     Also test BWERSD(.TRUE.,55329) reseed + BWERGT/BWERPT save/restore.
      IMPLICIT NONE
      REAL SEL
      DOUBLE PRECISION DS
      INTEGER I
C     (a) direct stream from default seed
      DO 10 I=1,8
        CALL BWERAN(SEL)
        WRITE(*,'(Z8.8)') SEL
   10 CONTINUE
      WRITE(*,*) '--- reseed 55329 ---'
      CALL BWERSD(.TRUE., 55329.0)
      DO 20 I=1,8
        CALL BWERAN(SEL)
        WRITE(*,'(Z8.8)') SEL
   20 CONTINUE
      WRITE(*,*) '--- get/put seed save-restore ---'
      CALL BWERSD(.TRUE., 55329.0)
      CALL BWERAN(SEL)
      CALL BWERGT(DS)
      CALL BWERAN(SEL)
      WRITE(*,'(Z8.8)') SEL
      CALL BWERPT(DS)
      CALL BWERAN(SEL)
      WRITE(*,'(Z8.8)') SEL
      END
