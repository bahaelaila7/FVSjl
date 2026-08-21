      PROGRAM DRVISTD
      INCLUDE 'PRGPRM.F77'
      INCLUDE 'PPEPRM.F77'
      INCLUDE 'BMPRM.F77'
      INCLUDE 'BMCOM.F77'
      INCLUDE 'BMRCOM.F77'
      REAL X, MID
      INTEGER ISIZ, LOW, I, J
      DATA UPSIZ /3,6,9,12,15,18,21,25,30,50/
      DATA WPSIZ /10,20,60/
      BMS0=55329.0D0
      BMS1=0.0D0
      BMSS=55329.0
C     MSBA + L2D (bminit subset)
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
      J=1
      DO 105 ISIZ=1,NSCL
        IF (UPSIZ(ISIZ).LE.3) THEN
          L2D(ISIZ)=0
        ELSE
          IDIF=(UPSIZ(ISIZ)-WPSIZ(J))+(UPSIZ(ISIZ-1)-WPSIZ(J))
          IF (UPSIZ(ISIZ).LE.WPSIZ(J)) THEN
            L2D(ISIZ)=J
          ELSEIF (IDIF.LE.0) THEN
            L2D(ISIZ)=J
          ELSE
            J=J+1
            L2D(ISIZ)=J
          ENDIF
        ENDIF
  105 CONTINUE
      LBMDEB=.FALSE.
      IQPTYP(1,1)=1
      DO 110 ISIZ=1,NSCL
        TREE(1,ISIZ,1)=0.0
        GRF(1,ISIZ)=1.0
        PBKILL(1,ISIZ)=0.0
        PITCH(1,ISIZ)=0.0
        STRIP(1,ISIZ)=0.0
        SPCLT(1,ISIZ,1)=0.0
        TVOL(1,ISIZ,1)=10.0
  110 CONTINUE
      TREE(1,3,1)=50.0
      TREE(1,7,1)=40.0
      TREE(1,5,1)=30.0
      SPRAY(1,1)=0.0
      BKP(1)=40.0
      FINAL(1,1)=0.0
      FINAL(1,2)=0.0
      FINAL(1,3)=0.0
      CALL BMISTD(1, 2000)
      DO 120 ISIZ=1,NSCL
        WRITE(*,'(A,I2,1X,Z8.8,1X,Z8.8,1X,Z8.8)') 'K',ISIZ,
     >    PBKILL(1,ISIZ), PITCH(1,ISIZ), STRIP(1,ISIZ)
  120 CONTINUE
      WRITE(*,'(A,Z8.8)') 'BKP  ', BKP(1)
      WRITE(*,'(A,Z8.8,1X,Z8.8,1X,Z8.8)') 'FIN ',
     >    FINAL(1,1),FINAL(1,2),FINAL(1,3)
      END
      SUBROUTINE SPLAAR(ISTD,AREA,IRC)
      AREA=1.0
      IRC=0
      RETURN
      END
