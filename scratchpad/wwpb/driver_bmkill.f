      PROGRAM DRVKILL
      INCLUDE 'PRGPRM.F77'
      INCLUDE 'PPEPRM.F77'
      INCLUDE 'BMPRM.F77'
      INCLUDE 'PPCNTL.F77'
      INCLUDE 'PLOT.F77'
      INCLUDE 'OUTCOM.F77'
      INCLUDE 'ARRAYS.F77'
      INCLUDE 'CONTRL.F77'
      INCLUDE 'BMCOM.F77'
      INTEGER ISIZ, ISPC
      DATA UPSIZ /3,6,9,12,15,18,21,25,30,50/
      LBMDEB=.FALSE.
      IBMYR1=1
      BMSTND=1
      BMSDIX(1)=1
      STOCK(1)=.TRUE.
      ISTND=1
      IBMMRT=1
C     host designation: sp7 (LP) is MPB host
      PBSPEC=1
      DO 5 ISPC=1,MAXSP
        HSPEC(1,ISPC)=0
        HSPEC(2,ISPC)=0
        HSPEC(3,ISPC)=0
        ISCT(ISPC,1)=0
        ISCT(ISPC,2)=0
    5 CONTINUE
      HSPEC(1,7)=1
C     BM state: OTPA + TPBK for class 3 host
      DO 6 ISIZ=1,NSCL
        OTPA(1,ISIZ,1)=0.0
        OTPA(1,ISIZ,2)=0.0
        TPBK(1,ISIZ,1,1)=0.0
        TPBK(1,ISIZ,1,2)=0.0
        TPBK(1,ISIZ,1,3)=0.0
    6 CONTINUE
      OTPA(1,3,1)=50.0
      TPBK(1,3,1,3)=60.0
C     FVS treelist: 2 LP trees, class 3 (dbh 7), prob 25 each, wk2=2
      ITRN=2
      ISCT(7,1)=1
      ISCT(7,2)=2
      IND1(1)=1
      IND1(2)=2
      ISP(1)=7
      ISP(2)=7
      DBH(1)=7.0
      DBH(2)=7.0
      PROB(1)=25.0
      PROB(2)=25.0
      WK2(1)=2.0
      WK2(2)=2.0
      CFV(1)=10.0
      CFV(2)=10.0
      CALL BMKILL
      WRITE(*,'(A,Z8.8)') 'WK2_1 ', WK2(1)
      WRITE(*,'(A,Z8.8)') 'WK2_2 ', WK2(2)
      END
      SUBROUTINE OPBISR(N,LST,IS,IBM)
      DIMENSION LST(*)
      IBM=1
      RETURN
      END
      SUBROUTINE SVMORT(I,A,J)
      DIMENSION A(*)
      RETURN
      END
