      PROGRAM DRVNUM
      INCLUDE 'PRGPRM.F77'
      INCLUDE 'PPEPRM.F77'
      INCLUDE 'BMPRM.F77'
      INCLUDE 'BMCOM.F77'
      INTEGER ISIZ, LOW, I
      REAL    MID, X, BASUM
      DATA UPSIZ /3,6,9,12,15,18,21,25,30,50/
C     MSBA (bminit)
      X = PIE / (24. * 24.)
      DO 100 ISIZ=1,NSCL
        IF (ISIZ.EQ.1) THEN
          LOW=0
        ELSE
          LOW=UPSIZ(ISIZ-1)
        ENDIF
        MID=FLOAT(UPSIZ(ISIZ)+LOW)*0.5
        MSBA(ISIZ)=MID*MID*X
  100 CONTINUE
      PBSPEC=1
      IPSON=.FALSE.
      LBMDEB=.FALSE.
      ISCMIN(1)=3
      ISCMIN(2)=3
      ISCMIN(3)=2
      USERA(1)=1.0
      USERA(2)=1.0
      USERA(3)=1.0
      GRFSTD(1)=0.80
      REPPHE(1)=0.0
      SSBATK(1)=0.0
      BASUM=0.0
      DO 110 ISIZ=1,NSCL
        BAH(1,ISIZ)=2.0+ISIZ*1.5
        BANH(1,ISIZ)=0.5
        TREE(1,ISIZ,1)=50.0
        PITCH(1,ISIZ)=0.0
        STRIKE(1,ISIZ)=0.0
        TOPKLL(1,ISIZ)=0.0
        OTHATT(1,ISIZ)=0.0
        SCORCH(1,ISIZ)=0.0
        BASUM=BASUM+BAH(1,ISIZ)+BANH(1,ISIZ)
  110 CONTINUE
      ATRPHE(1)=0.0
      BAH(1,NSCL+1)=BASUM-0.5*NSCL
      BANH(1,NSCL+1)=0.5*NSCL
      DO 120 I=1,2
        DO 121 ISIZ=1,3
          DWPHOS(1,I,ISIZ)=0.0
  121   CONTINUE
  120 CONTINUE
      CALL BMCNUM(1,2000)
      WRITE(*,'(A,Z8.8)') 'NUMER1 ', NUMER(1,1)
      WRITE(*,'(A,Z8.8)') 'TFOOD1 ', TFOOD(1,1)
      END
