      PROGRAM GENPHOTO
      IMPLICIT NONE
      INCLUDE 'PRGPRM.F77'
      INCLUDE 'FMPARM.F77'
      REAL FOTOVAL(MXFLCL), FOTOVALS(9)
      INTEGER IREF, ICODE, I
      DO IREF = 1, 32
        DO ICODE = 1, 200
          DO I = 1, MXFLCL
            FOTOVAL(I) = -99.0
          ENDDO
          DO I = 1, 9
            FOTOVALS(I) = -99.0
          ENDDO
          CALL FMPHOTOVAL(IREF, ICODE, FOTOVAL, FOTOVALS)
          IF (FOTOVAL(1) .LT. 0.0 .AND. FOTOVAL(1) .GT. -98.0) GOTO 100
          IF (FOTOVAL(1) .LT. -98.0) GOTO 100
          WRITE (6,'(I3,",",I3,20(",",ES16.9))') IREF, ICODE,
     &      (FOTOVAL(I),I=1,MXFLCL), (FOTOVALS(I),I=1,9)
        ENDDO
  100   CONTINUE
      ENDDO
      END
