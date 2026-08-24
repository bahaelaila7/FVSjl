      PROGRAM DRVATCTM
C     Multi-stand BMATCT golden: exercises the interstand spatial BKP
C     redistribution (SCORE/PROP/SPLADS) + the Outside-World path.
C     Geometry (XLOC,YLOC,AREA) supplied to SPLAAR/SPLALO/SPLADS via /GEOM/.
      INCLUDE 'PRGPRM.F77'
      INCLUDE 'PPEPRM.F77'
      INCLUDE 'BMPRM.F77'
      INCLUDE 'PPCNTL.F77'
      INCLUDE 'BMCOM.F77'
      INCLUDE 'BMPCOM.F77'
      REAL GXL(4),GYL(4),GAR(4)
      COMMON /GEOM/ GXL,GYL,GAR
      INTEGER I

C     ---- landscape geometry (meters, meters, acres) ----
      GXL(1)=0.0    ; GYL(1)=0.0    ; GAR(1)=100.0
      GXL(2)=1000.0 ; GYL(2)=0.0    ; GAR(2)=200.0
      GXL(3)=0.0    ; GYL(3)=1500.0 ; GAR(3)=150.0
      GXL(4)=1200.0 ; GYL(4)=1200.0 ; GAR(4)=80.0

C     ---- MSBA (size-class BA) ----
      CALL SETMSBA

C     ---- common BM setup ----
      PBSPEC=1
      IPSON=.FALSE.
      LBMDEB=.FALSE.
      BMSTND=4
      DO 10 I=1,4
        BMSDIX(I)=I
        STOCK(I)=.TRUE.
   10 CONTINUE
      SDD=0.0
      USERA(1)=1.0
      SELFA(1)=1.0
      USERC(1)=100.0
      URMAX(1)=1.0
      USERA(3)=1.0
      SELFA(3)=1.0
      USERC(3)=100.0
      URMAX(3)=1.0
C     per-stand attractiveness numerator + BKP + food
      NUMER(1,1)=1000.0 ; BKP(1)=100.0 ; TFOOD(1,1)=50.0
      NUMER(1,2)= 300.0 ; BKP(2)= 40.0 ; TFOOD(2,1)=80.0
      NUMER(1,3)= 700.0 ; BKP(3)= 10.0 ; TFOOD(3,1)=30.0
      NUMER(1,4)= 500.0 ; BKP(4)= 60.0 ; TFOOD(4,1)=20.0
      DO 11 I=1,4
        NUMER(2,I)=0.0
        BKPIPS(I)=0.0
        TFOOD(I,2)=0.0
        GRFSTD(I)=1.0
   11 CONTINUE
      LBAD=.FALSE.
      STOCKO=1.0
      RVOD=1.0
      IBADBB=1
      BADREP(1)=1.0
      BADREP(3)=1.0

C     ======== SCENARIO A: Outside World OFF (pure interstand) ========
      OUTOFF=.TRUE.
      ACDONE=.FALSE.
      CALL BMATCT(2000)
      CALL DUMP('A')

C     ======== SCENARIO B: Outside World ON, floating (UFLOAT=-1) ======
C     restore driving inputs (BMATCT overwrites BKP)
      BKP(1)=100.0 ; BKP(2)=40.0 ; BKP(3)=10.0 ; BKP(4)=60.0
      OUTOFF=.FALSE.
      UFLOAT=-1.0
      ACDONE=.FALSE.
      CALL BMATCT(2000)
      CALL DUMP('B')

C     ======== SCENARIO C: Outside World ON, fixed constants ==========
      BKP(1)=100.0 ; BKP(2)=40.0 ; BKP(3)=10.0 ; BKP(4)=60.0
      OUTOFF=.FALSE.
      UFLOAT=0.0
      CBAO=50.0
      CRVOND=1.0
      CBAHO(1)=20.0
      CBASPO(1)=5.0
      CSPO(1)=2.0
      CBKPO(1)=10.0
      ACDONE=.FALSE.
      CALL BMATCT(2000)
      CALL DUMP('C')

C     ==== SCENARIO D: MULTI-YEAR OW-off cascade (ACDONE persists yrs 2..5) ====
C     Reset BKP, hold NUMER fixed, let BKP evolve through repeated BMATCT calls
C     across 5 years (exercises the ACDONE spatial-cache persistence that the
C     master-cycle bmdrv loop relies on). Dump BKP per stand per year.
      OUTOFF=.TRUE.
      ACDONE=.FALSE.
      BKP(1)=100.0 ; BKP(2)=40.0 ; BKP(3)=10.0 ; BKP(4)=60.0
      DO 300 IYR=1,5
        CALL BMATCT(2000+IYR)
        DO 310 I=1,4
          WRITE(*,'(A2,I1,1X,I2,1X,A,1X,Z8.8)') 'D',IYR,I,'BKP',BKP(I)
  310   CONTINUE
  300 CONTINUE

C     ==== SCENARIO E: MULTI-YEAR OW-floating cascade ====
      OUTOFF=.FALSE.
      UFLOAT=-1.0
      ACDONE=.FALSE.
      BKP(1)=100.0 ; BKP(2)=40.0 ; BKP(3)=10.0 ; BKP(4)=60.0
      DO 400 IYR=1,5
        CALL BMATCT(2100+IYR)
        DO 410 I=1,4
          WRITE(*,'(A2,I1,1X,I2,1X,A,1X,Z8.8)') 'E',IYR,I,'BKP',BKP(I)
  410   CONTINUE
  400 CONTINUE

      END

      SUBROUTINE SETMSBA
      INCLUDE 'PRGPRM.F77'
      INCLUDE 'PPEPRM.F77'
      INCLUDE 'BMPRM.F77'
      INCLUDE 'PPCNTL.F77'
      INCLUDE 'BMCOM.F77'
      INTEGER ISIZ,LOW
      REAL X,MID
      DATA UPSIZ /3,6,9,12,15,18,21,25,30,50/
      X=PIE/(24.*24.)
      DO 100 ISIZ=1,NSCL
        IF (ISIZ.EQ.1) THEN
          LOW=0
        ELSE
          LOW=UPSIZ(ISIZ-1)
        ENDIF
        MID=FLOAT(UPSIZ(ISIZ)+LOW)*0.5
        MSBA(ISIZ)=MID*MID*X
  100 CONTINUE
      RETURN
      END

      SUBROUTINE DUMP(TAG)
      INCLUDE 'PRGPRM.F77'
      INCLUDE 'PPEPRM.F77'
      INCLUDE 'BMPRM.F77'
      INCLUDE 'PPCNTL.F77'
      INCLUDE 'BMCOM.F77'
      INCLUDE 'BMPCOM.F77'
      CHARACTER*1 TAG
      INTEGER I
      DO 200 I=1,BMSTND
        WRITE(*,'(A1,1X,I2,1X,A,1X,Z8.8,1X,A,1X,Z8.8)')
     >    TAG,I,'BKP',BKP(I),'ATTC',ATTC(1,I)
        WRITE(*,'(A1,1X,I2,1X,A,1X,Z8.8,1X,A,1X,Z8.8)')
     >    TAG,I,'ATB',ATTBYK(1,I),'BOU',BKPOUT(1,I)
        WRITE(*,'(A1,1X,I2,1X,A,1X,Z8.8,1X,A,1X,Z8.8)')
     >    TAG,I,'BIN',BKPIN(1,I),'SLF',SELFBKP(1,I)
        WRITE(*,'(A1,1X,I2,1X,A,1X,Z8.8)')
     >    TAG,I,'BPS',BKPS(I)
  200 CONTINUE
      RETURN
      END

C     ---- geometry stubs backed by /GEOM/ ----
      SUBROUTINE SPLAAR (ISTD, A, IRC)
      REAL GXL(4),GYL(4),GAR(4)
      COMMON /GEOM/ GXL,GYL,GAR
      INTEGER ISTD,IRC
      REAL A
      A=GAR(ISTD)
      IRC=0
      RETURN
      END
      SUBROUTINE SPLALO (ISTD, XP, YP, IRC)
      REAL GXL(4),GYL(4),GAR(4)
      COMMON /GEOM/ GXL,GYL,GAR
      INTEGER ISTD,IRC
      REAL XP,YP
      XP=GXL(ISTD)
      YP=GYL(ISTD)
      IRC=0
      RETURN
      END
      SUBROUTINE SPLADS (I1, I2, DIST, IRC)
      REAL GXL(4),GYL(4),GAR(4)
      COMMON /GEOM/ GXL,GYL,GAR
      INTEGER I1,I2,IRC
      REAL DIST,DX,DY
      DX=GXL(I1)-GXL(I2)
      DY=GYL(I1)-GYL(I2)
      DIST=SQRT(DX*DX+DY*DY)
      IRC=0
      RETURN
      END
      SUBROUTINE GPGET2 (IA,IY,MP,NP,PR,MX,SC,ML,LOK)
      LOGICAL LOK
      DIMENSION PR(*), ML(*)
      LOK=.FALSE.
      RETURN
      END
