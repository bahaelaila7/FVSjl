C  Standalone driver over the EXACT brann.f MINSTD LCG (seed 55329), printing
C  Float32 draws as IEEE-754 hex (bit-exact goldens for wpbr_rand!/wpbr_seed!).
C  Mirrors scratchpad/dftm/driver_tmrann.f. Compile: gfortran-16 -std=legacy.
      PROGRAM DRVBRANN
      IMPLICIT NONE
      DOUBLE PRECISION BRS0, BRS1, BRSS
      REAL SEL, SEED
      INTEGER I, IHEX
      EQUIVALENCE (SEL, IHEX)
      BRSS = 55329.D0
      BRS0 = 55329.D0
C  --- 6-draw stream from the default seed 55329 ---
      DO I = 1, 6
         BRS1 = DMOD(16807.D0*BRS0, 2147483647.D0)
         SEL  = BRS1/2147483648.D0
         BRS0 = BRS1
         WRITE(*,'(A,Z8.8)') 'STREAM ', IHEX
      END DO
C  --- BRNSED reseed: even 100 -> forced odd 101, then 1 draw ---
      SEED = 100.0
      IF (AMOD(SEED,2.0).EQ.0.) SEED = SEED + 1
      BRSS = SEED
      BRS0 = SEED
      BRS1 = DMOD(16807.D0*BRS0, 2147483647.D0)
      SEL  = BRS1/2147483648.D0
      BRS0 = BRS1
      WRITE(*,'(A,Z8.8)') 'SEED101 ', IHEX
C  --- reseed odd 55, then 2 draws ---
      SEED = 55.0
      IF (AMOD(SEED,2.0).EQ.0.) SEED = SEED + 1
      BRSS = SEED
      BRS0 = SEED
      BRS1 = DMOD(16807.D0*BRS0, 2147483647.D0)
      SEL  = BRS1/2147483648.D0
      BRS0 = BRS1
      WRITE(*,'(A,Z8.8)') 'SEED55A ', IHEX
      BRS1 = DMOD(16807.D0*BRS0, 2147483647.D0)
      SEL  = BRS1/2147483648.D0
      BRS0 = BRS1
      WRITE(*,'(A,Z8.8)') 'SEED55B ', IHEX
C  --- LSET=false: reset BRS0 to SS(=55), then 1 draw (== SEED55A) ---
      SEED = BRSS
      BRS0 = SEED
      BRS1 = DMOD(16807.D0*BRS0, 2147483647.D0)
      SEL  = BRS1/2147483648.D0
      BRS0 = BRS1
      WRITE(*,'(A,Z8.8)') 'RESET55 ', IHEX
      END
