      PROGRAM DRVINIT
C     Replicates bminit.f MSBA/UPBA/INC computation verbatim (lines 71-120).
      INCLUDE 'PRGPRM.F77'
      INCLUDE 'PPEPRM.F77'
      INCLUDE 'BMPRM.F77'
      INCLUDE 'BMCOM.F77'
      INTEGER ISIZ, LOW, I
      REAL    MID, X
      DATA UPSIZ /3,6,9,12,15,18,21,25,30,50/
      X = PIE / (24. * 24.)
      DO 100 ISIZ = 1,NSCL
         IF (ISIZ .EQ. 1) THEN
            LOW = 0
         ELSE
            LOW = UPSIZ(ISIZ-1)
         ENDIF
         MID = FLOAT(UPSIZ(ISIZ) + LOW) * 0.5
         MSBA(ISIZ) = MID * MID * X
         UPBA(ISIZ) = UPSIZ(ISIZ) * UPSIZ(ISIZ) * X
  100 CONTINUE
      RSLOPE = 1.0/10.0
      REPMAX = 4.0
      DBHMAX = 36.0
      REPLAC = 6.0
      B = 1 - (RSLOPE * REPLAC)
      LOW = 0
      DO 110 I = 1,NSCL
         IF (I .GT. 1) LOW = UPSIZ(I-1)
         DBHMID = (UPSIZ(I) + LOW) / 2.0
         IF (DBHMID .LT. DBHMAX) THEN
            INC(1,I) = (RSLOPE * DBHMID) + B
         ELSE
            INC(1,I) = REPMAX
         ENDIF
         INC(2,I) = INC(1,I)
         INC(3,I) = INC(1,1) * 0.1
  110 CONTINUE
      DO 200 ISIZ=1,NSCL
        WRITE(*,'(A,I2,1X,Z8.8,1X,Z8.8,1X,Z8.8)') 'C',ISIZ,
     >     MSBA(ISIZ), UPBA(ISIZ), INC(1,ISIZ)
  200 CONTINUE
      WRITE(*,'(A,Z8.8)') 'INC3_1 ', INC(3,1)
      END
