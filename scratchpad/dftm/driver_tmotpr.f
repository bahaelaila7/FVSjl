      PROGRAM DRV
      IMPLICIT NONE
      REAL ELEV,SLOPE,ASPECT,TOPO,RELDEN,RELDSP3,RELDSP4,TPROB
      REAL STDCLO,HOST,AVCRDI,PROTBK,TMASHD,PGFBA,BA
C  --- method 1 (Heller) inputs ---
      ELEV=35.0
      SLOPE=0.30
      ASPECT=1.2
      TOPO=2.0
      RELDEN=120.0
      RELDSP3=40.0
      RELDSP4=25.0
      TPROB=200.0
      STDCLO=100.0
      IF (RELDEN .LT. 150.0) STDCLO=(RELDEN/150.0)*100.0
      HOST=((RELDSP3+RELDSP4)/RELDEN)*100.0
      AVCRDI=2.0*SQRT((435.6*(RELDEN/TPROB))/3.141593)
      PROTBK = 1.0 / (1.0 + EXP(-(-0.852926 -
     >         (0.00023762 * ELEV  * 100.0) -
     >         (0.00052207 * SLOPE * 100.0) +
     >         (0.00348125 * SLOPE * 100.0 * COS(ASPECT)) +
     >         (0.00733782 * SLOPE * 100.0 * SIN(ASPECT)) -
     >         (0.37562400 * TOPO) +
     >         (0.02049930 * STDCLO) +
     >         (0.01866930 * HOST) +
     >         (0.01685520 * AVCRDI))))
      WRITE(*,'(A,Z8.8)') 'M1 ', PROTBK
C  --- method 2 (Mika-Moore w/ ash) ---
      TOPO=3.0
      TMASHD=15.93
      PGFBA=0.4
      BA=180.0
      PROTBK = 1.0 / (1.0 + EXP(-(-12.7746 -
     >         1.4764 * (TOPO - 1.0) -
     >         0.0992 * TMASHD +
     >         4.8334 * PGFBA +
     >         2.5529 * ALOG(BA))))
      WRITE(*,'(A,Z8.8)') 'M2 ', PROTBK
C  --- method 3 (Mika-Moore no ash) ---
      PROTBK = 1.0 / (1.0 + EXP(-(-12.3652 -
     >         1.4050 * (TOPO - 1.0) +
     >         3.8230 * PGFBA +
     >         2.2426 * ALOG(BA))))
      WRITE(*,'(A,Z8.8)') 'M3 ', PROTBK
      END
