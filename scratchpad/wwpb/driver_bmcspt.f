      PROGRAM DRVSPT
      INCLUDE 'PRGPRM.F77'
      INCLUDE 'PPEPRM.F77'
      INCLUDE 'BMPRM.F77'
      INCLUDE 'BMCOM.F77'
      INTEGER ISIZ
      PBSPEC = 1
      IPSON  = .FALSE.
      LBMDEB = .FALSE.
      DO 5 ISIZ=1,NSCL
        TREE(1,ISIZ,1)=100.0
        PITCH(1,ISIZ)=0.0
        STRIKE(1,ISIZ)=0.0
        TOPKLL(1,ISIZ)=0.0
        OTHATT(1,ISIZ)=0.0
        SCORCH(1,ISIZ)=0.0
    5 CONTINUE
      ATRPHE(1)=0.0
C     case A: all zero (class 1)
      CALL BMCSPT(1,1,1,2000)
C     case B: one SP=0.3 (pitch), class 2
      PITCH(1,2)=0.30
      CALL BMCSPT(1,2,1,2000)
C     case C: two SP (pitch 0.3, strike 0.4), class 3
      PITCH(1,3)=0.30
      STRIKE(1,3)=0.40
      CALL BMCSPT(1,3,1,2000)
C     case D: three SP (pitch .2, strike .3, scorch .4), class 4
      PITCH(1,4)=0.20
      STRIKE(1,4)=0.30
      SCORCH(1,4)=0.40
      CALL BMCSPT(1,4,1,2000)
C     case E: one SP = 1.0 (saturated), class 5
      SCORCH(1,5)=1.0
      CALL BMCSPT(1,5,1,2000)
      WRITE(*,'(A,Z8.8)') 'A ', SPCLT(1,1,1)
      WRITE(*,'(A,Z8.8)') 'B ', SPCLT(1,2,1)
      WRITE(*,'(A,Z8.8)') 'C ', SPCLT(1,3,1)
      WRITE(*,'(A,Z8.8)') 'D ', SPCLT(1,4,1)
      WRITE(*,'(A,Z8.8)') 'E ', SPCLT(1,5,1)
      END
