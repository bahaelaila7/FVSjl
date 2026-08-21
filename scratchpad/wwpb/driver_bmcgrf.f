      PROGRAM DRVGRF
      INCLUDE 'PRGPRM.F77'
      INCLUDE 'PPEPRM.F77'
      INCLUDE 'BMPRM.F77'
      INCLUDE 'BMCOM.F77'
      INCLUDE 'BMPCOM.F77'
      INTEGER ICLS
      REAL OLDGRF(NSCL)
      REAL BASUM
C     set up one stand (ISTD=1): host BA increasing by class, nonhost small
      PBSPEC = 1
      LCDENS = .TRUE.
      LBMDEB = .FALSE.
      BASUM = 0.0
      DO 10 ICLS=1,NSCL
        BAH(1,ICLS)  = 2.0 + ICLS * 1.5
        BANH(1,ICLS) = 0.5
        GRF(1,ICLS)  = 0.0
        SDMR(1,ICLS) = 0.0
        SRR(1,ICLS)  = 0.0
        SSR(1,ICLS)  = 0.0
        OTHATT(1,ICLS)= 0.0
        TOPKLL(1,ICLS)= 0.0
        SCORCH(1,ICLS)= 0.0
        STRIKE(1,ICLS)= 0.0
        RVDSC(1,ICLS) = 0.90
        RVDFOL(1,ICLS)= 0.0
        BASUM = BASUM + BAH(1,ICLS) + BANH(1,ICLS)
   10 CONTINUE
      BAH(1,NSCL+1)  = BASUM - 0.5*NSCL
      BANH(1,NSCL+1) = 0.5*NSCL
      CALL BMCGRF(1, 2000, OLDGRF)
      WRITE(*,'(A,Z8.8)') 'GRFSTD ', GRFSTD(1)
      WRITE(*,'(A,Z8.8)') 'RVDNST ', RVDNST(1)
      DO 20 ICLS=1,NSCL
        WRITE(*,'(A,I2,1X,Z8.8)') 'GRF ', ICLS, GRF(1,ICLS)
   20 CONTINUE
      END
