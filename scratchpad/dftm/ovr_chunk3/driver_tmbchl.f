      PROGRAM DRVBCHL
C     Standalone driver over pristine dftm/tmbchl.f + tmrann.f.
C     Seed default 55329.  Draw TMBCHL(9,2) x3 then TMBCHL(11,3) x3,
C     dumping each result + the raw-draw counter as Float32 hex.
      IMPLICIT NONE
      REAL TMBCHL, X
      INTEGER I
      DO 10 I=1,3
        X = TMBCHL(9.0, 2.0)
        WRITE(0,'(A,I3,1X,Z8.8)') 'BCHL_DF ',I,TRANSFER(X,1)
   10 CONTINUE
      DO 20 I=1,3
        X = TMBCHL(11.0, 3.0)
        WRITE(0,'(A,I3,1X,Z8.8)') 'BCHL_GF ',I,TRANSFER(X,1)
   20 CONTINUE
      END
