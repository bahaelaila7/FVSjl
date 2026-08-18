      PROGRAM DRIVER_DFB
      IMPLICIT NONE
C
C     Standalone golden-oracle driver for the DETERMINISTIC Douglas-fir
C     Beetle routines: DFBIND, DFBDBH, DFBER, DFBPRB.
C     Populates ARRAYS/CONTRL/PLOT/DFBCOM commons with a synthetic
C     all-Douglas-fir stand, calls the real (pristine) DFB code, and
C     dumps Float32 bit patterns (Z8 hex) so a Julia port can be checked
C     bit-exact.  Reads pristine Fortran only; writes only its dump.
C
      INCLUDE 'PRGPRM.F77'
      INCLUDE 'ARRAYS.F77'
      INCLUDE 'CONTRL.F77'
      INCLUDE 'PLOT.F77'
      INCLUDE 'DFBCOM.F77'

      INTEGER N, I, SZ, IB
      LOGICAL LMIN
      REAL PROTBK, TESTD(14)

C.... synthetic DF stand: DBH (in) and TPA (PROB) per tree record
      REAL DBHV(12), TPAV(12)
      DATA DBHV /  3.0,  5.0,  8.0,  8.7,  9.0, 10.0,
     &            12.4, 15.0, 18.0, 20.0, 24.0, 30.0 /
      DATA TPAV / 40.0, 30.0, 25.0, 20.0, 18.0, 15.0,
     &            12.0, 10.0,  8.0,  6.0,  4.0,  2.0 /

      N = 12

C.... ARRAYS: all records are Douglas-fir (IDFSPC set in block data)
      BA = 0.0
      DO 10 I = 1, N
         DBH(I)  = DBHV(I)
         PROB(I) = TPAV(I)
         CFV(I)  = DBHV(I) * DBHV(I) * 0.30
         IND1(I) = I
         BA = BA + 0.005454154 * DBH(I) * DBH(I) * PROB(I)
   10 CONTINUE

C.... CONTRL: tree count and species-class index table
      ITRN = N
      DO 20 I = 1, MAXSP
         ISCT(I,1) = 0
         ISCT(I,2) = 0
   20 CONTINUE
      ISCT(IDFSPC,1) = 1
      ISCT(IDFSPC,2) = N

C.... DFBCOM run-time flags (block data provides WINSUC/IFVSSP/PERDD/
C.... ROWDOM/IDFSPC/S0/SS; set the rest to the DFBINT defaults we need)
      DEBUIN = .FALSE.
      LBAMOD = .FALSE.
      JODFB  = 6

C.... ---- call the real deterministic routines ----
      CALL DFBER (LMIN)
      CALL DFBPRB (PROTBK)
      CALL DFBDBH

C.... ---- dump ----
      OPEN (51, FILE='dfb_golden.txt', STATUS='UNKNOWN')
      WRITE (51,'(A,I0)') 'N ', N
      WRITE (51,'(A,I0)') 'IDFSPC ', IDFSPC
      WRITE (51,'(A,I0)') 'LMIN ', MERGE(1,0,LMIN)
      CALL PR(51,'BA9',    BA9)
      CALL PR(51,'BADF9',  BADF9)
      CALL PR(51,'A45DBH', A45DBH)
      CALL PR(51,'PBADF4', PBADF4)
      CALL PR(51,'BA',     BA)
      CALL PR(51,'PROTBK', PROTBK)
      DO 30 I = 1, 20
         IB = TRANSFER(START(I),1)
         WRITE (51,'(A,I0,A,Z8.8)') 'START ', I, ' ', IB
   30 CONTINUE

C.... DFBIND size-class mapping over a spread of diameters
      TESTD(1)=0.4
      TESTD(2)=1.0
      TESTD(3)=1.9
      TESTD(4)=2.0
      TESTD(5)=8.7
      TESTD(6)=9.0
      TESTD(7)=10.0
      TESTD(8)=12.4
      TESTD(9)=19.0
      TESTD(10)=20.0
      TESTD(11)=24.0
      TESTD(12)=30.0
      TESTD(13)=40.0
      TESTD(14)=45.0
      DO 40 I = 1, 14
         CALL DFBIND (TESTD(I), SZ)
         IB = TRANSFER(TESTD(I),1)
         WRITE (51,'(A,Z8.8,A,I0)') 'DFBIND ', IB, ' ', SZ
   40 CONTINUE

      CLOSE (51)
      WRITE (6,*) 'DONE LMIN=', LMIN, ' PROTBK=', PROTBK
      END

      SUBROUTINE PR (U, NAME, V)
      INTEGER U
      CHARACTER*(*) NAME
      REAL V
      WRITE (U,'(A,1X,Z8.8,1X,1PE16.8)') NAME, TRANSFER(V,1), V
      RETURN
      END
