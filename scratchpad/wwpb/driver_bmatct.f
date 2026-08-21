      PROGRAM DRVATCT
      INCLUDE 'PRGPRM.F77'
      INCLUDE 'PPEPRM.F77'
      INCLUDE 'BMPRM.F77'
      INCLUDE 'PPCNTL.F77'
      INCLUDE 'BMCOM.F77'
      INCLUDE 'BMPCOM.F77'
      INTEGER ISIZ
      REAL X, MID
      INTEGER LOW
      DATA UPSIZ /3,6,9,12,15,18,21,25,30,50/
C     MSBA (needed for ATTBYK path? no, but set anyway)
      X = PIE/(24.*24.)
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
      OUTOFF=.TRUE.
      ACDONE=.FALSE.
      BMSTND=1
      BMSDIX(1)=1
      STOCK(1)=.TRUE.
      SDD=0.0
      USERC(1)=100.0
      SELFA(1)=1.0
      URMAX(1)=15.0
      USERA(1)=1.0
      NUMER(1,1)=1000.0
      NUMER(2,1)=0.0
      BKP(1)=100.0
      BKPIPS(1)=0.0
      TFOOD(1,1)=50.0
      TFOOD(1,2)=0.0
      CALL BMATCT(2000)
      WRITE(*,'(A,Z8.8)') 'BKP    ', BKP(1)
      WRITE(*,'(A,Z8.8)') 'BKPIPS ', BKPIPS(1)
      END
C     ---- stubs for absent PPE spatial + scheduler ----
      SUBROUTINE SPLAAR (ISTD, AREA, IRC)
      AREA = 10.0
      IRC = 0
      RETURN
      END
      SUBROUTINE SPLALO (ISTD, XP, YP, IRC)
      XP = 0.0
      YP = 0.0
      IRC = 0
      RETURN
      END
      SUBROUTINE SPLADS (I1, I2, DIST, IRC)
      DIST = 0.0
      IRC = 0
      RETURN
      END
      SUBROUTINE GPGET2 (IA,IY,MP,NP,PR,MX,SC,ML,LOK)
      LOGICAL LOK
      DIMENSION PR(*), ML(*)
      LOK=.FALSE.
      RETURN
      END
