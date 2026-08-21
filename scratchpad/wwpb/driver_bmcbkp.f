      PROGRAM DRVBKP
      INCLUDE 'PRGPRM.F77'
      INCLUDE 'PPEPRM.F77'
      INCLUDE 'BMPRM.F77'
      INCLUDE 'PPCNTL.F77'
      INCLUDE 'BMCOM.F77'
      INTEGER ISIZ, LOW, I
      REAL    MID, X, OLDGRF(NSCL), RSLP, B
      DATA UPSIZ /3,6,9,12,15,18,21,25,30,50/
C     BMINIT coeffs (verbatim)
      X = PIE / (24. * 24.)
      DO 100 ISIZ=1,NSCL
        IF (ISIZ.EQ.1) THEN
          LOW=0
        ELSE
          LOW=UPSIZ(ISIZ-1)
        ENDIF
        MID = FLOAT(UPSIZ(ISIZ)+LOW)*0.5
        MSBA(ISIZ) = MID*MID*X
  100 CONTINUE
      RSLP = 1.0/10.0
      B = 1 - (RSLP*6.0)
      LOW=0
      DO 110 I=1,NSCL
        IF (I.GT.1) LOW=UPSIZ(I-1)
        MID = (UPSIZ(I)+LOW)/2.0
        IF (MID.LT.36.0) THEN
          INC(1,I) = (RSLP*MID)+B
        ELSE
          INC(1,I) = 4.0
        ENDIF
        INC(2,I)=INC(1,I)
        INC(3,I)=INC(1,1)*0.1
  110 CONTINUE
C     scenario
      PBSPEC = 1
      NBGEN  = 1
      NIBGEN = 2
      IPSON  = .FALSE.
      IPSMIN = 2
      LFDBK  = .FALSE.
      LBMDEB = .FALSE.
      BMEND  = 1
      ICNT   = 0
      DO 120 ISIZ=1,NSCL
        PBKILL(1,ISIZ)=0.0
        ALLKLL(1,ISIZ)=0.0
        TOPKLL(1,ISIZ)=0.0
        STRIP(1,ISIZ) =0.0
        TREE(1,ISIZ,1)=0.0
  120 CONTINUE
      PBKILL(1,3)=10.0
      PBKILL(1,5)=5.0
      FINAL(1,1)=0.0
      FINAL(1,2)=0.0
      FINAL(1,3)=0.0
      BKP(1)=0.0
      BKPIPS(1)=0.0
      OLDBKP(1)=0.0
      DO 130 I=1,MXDWHC
        DO 131 ISIZ=1,MXDWSZ
          PSLASH(1,I,ISIZ)=0.0
  131   CONTINUE
  130 CONTINUE
      CALL BMCBKP(1, 2000, OLDGRF)
      WRITE(*,'(A,Z8.8)') 'BKP    ', BKP(1)
      WRITE(*,'(A,Z8.8)') 'BKPIPS ', BKPIPS(1)
      WRITE(*,'(A,Z8.8)') 'OLDBKP ', OLDBKP(1)
      END
C     ---- stubs for the absent PPE scheduler ----
      SUBROUTINE GPGET2 (IACT,IYR,MAXP,NPRMS,PRMS,MXS,SCNT,MYLST,LOK)
      LOGICAL LOK
      DIMENSION PRMS(*), MYLST(*)
      LOK = .FALSE.
      NPRMS = 0
      RETURN
      END
      SUBROUTINE GPADD (KODE,IDT,IACT,NPRMS,PRMS,J,MYLST)
      DIMENSION PRMS(*), MYLST(*)
      RETURN
      END
