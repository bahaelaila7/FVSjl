      PROGRAM GOLDEN_BMMAIN
C     Replicates bminit.f MSBA + bmout.f MAINOUT aggregation (DO-20 + IPS slash)
C     on a controlled WWPB stand state, to golden the jl wwpb_main_report port.
C     PIE=3.14159 (REAL param, matching jl WWPB_PI24). NSCL=10.
      IMPLICIT NONE
      INTEGER NSCL, MXDWHC, MXDWHZ, ISIZ, J, K, LOW
      PARAMETER (NSCL=10, MXDWHC=2, MXDWHZ=3)
      REAL PIE, X, MID
      PARAMETER (PIE=3.14159)
      INTEGER UPSIZ(NSCL)
      REAL MSBA(NSCL)
      REAL TREE(NSCL,2), BAH(NSCL), TVOL(NSCL,2)
      REAL PBKILL(NSCL), ALLKLL(NSCL), SPCLT(NSCL)
      REAL DWPHOS(MXDWHC,MXDWHZ)
      REAL BASTD, GRFSTD, OLDBKP, BKPA
      REAL TPAKLL, PROPHKLD, BAK_SC, TPA_SC, TVOL_SC, HVOL_SC, VOLK_SC
      REAL BAK_YR, TPA_YR, TPAK_YR, VOL_YR, VOLH_YR, VOLK_YR
      REAL BA_SP, SPCL_TPA, IPS_SLSH
      DATA UPSIZ /3,6,9,12,15,18,21,25,30,50/
C     ---- controlled state ----
      DO ISIZ=1,NSCL
        TREE(ISIZ,1)=0.0
        TREE(ISIZ,2)=0.0
        BAH(ISIZ)=0.0
        TVOL(ISIZ,1)=0.0
        TVOL(ISIZ,2)=0.0
        PBKILL(ISIZ)=0.0
        ALLKLL(ISIZ)=0.0
        SPCLT(ISIZ)=0.0
      ENDDO
      DO J=1,MXDWHC
        DO K=1,MXDWHZ
          DWPHOS(J,K)=0.0
        ENDDO
      ENDDO
C     class 3 (7-9in): 100 host + 20 nonhost, 50 host BA, vols 8/6, kills 5+2, 10% special
      TREE(3,1)=100.0
      TREE(3,2)=20.0
      BAH(3)=50.0
      TVOL(3,1)=8.0
      TVOL(3,2)=6.0
      PBKILL(3)=5.0
      ALLKLL(3)=2.0
      SPCLT(3)=0.10
C     class 5 (13-15in): 40 host, 30 BA, vol 15, kill 3, 5% special
      TREE(5,1)=40.0
      BAH(5)=30.0
      TVOL(5,1)=15.0
      PBKILL(5)=3.0
      SPCLT(5)=0.05
C     ips slash pools
      DWPHOS(1,1)=1.5
      DWPHOS(2,3)=0.75
      BASTD=95.0
      GRFSTD=1.7
      OLDBKP=12.3
      BKPA=4.2
C     ---- bminit MSBA ----
      X= PIE / (24. * 24.)
      DO 100 ISIZ = 1,NSCL
          IF (ISIZ .EQ. 1) THEN
              LOW= 0
          ELSE
              LOW= UPSIZ(ISIZ-1)
          ENDIF
          MID = FLOAT(UPSIZ(ISIZ) + LOW) * 0.5
          MSBA(ISIZ) = MID * MID * X
  100 CONTINUE
C     ---- bmout MAINOUT DO-20 ----
      BAK_YR=0.0
      TPA_YR=0.0
      TPAK_YR=0.0
      VOL_YR=0.0
      VOLH_YR=0.0
      VOLK_YR=0.0
      BA_SP=0.0
      SPCL_TPA=0.0
      IPS_SLSH=0.0
      DO 20 ISIZ=1,NSCL
        TPAKLL = PBKILL(ISIZ) + ALLKLL(ISIZ)
        IF (TREE(ISIZ,1) .GT. 1E-6) THEN
           PROPHKLD = TPAKLL / TREE(ISIZ,1)
        ELSE
           PROPHKLD = 0.0
        ENDIF
        BAK_SC = BAH(ISIZ) * PROPHKLD
        TPA_SC = TREE(ISIZ,1) + TREE(ISIZ,2)
        TVOL_SC = (TVOL(ISIZ,1)*TREE(ISIZ,1) + TVOL(ISIZ,2)*TREE(ISIZ,2))
        HVOL_SC = TVOL(ISIZ,1)*TREE(ISIZ,1)
        VOLK_SC = TPAKLL * TVOL(ISIZ,1)
        BAK_YR  = BAK_YR  + BAK_SC
        TPA_YR  = TPA_YR  + TPA_SC
        TPAK_YR = TPAK_YR + TPAKLL
        VOL_YR  = VOL_YR  + TVOL_SC
        VOLH_YR = VOLH_YR + HVOL_SC
        VOLK_YR = VOLK_YR + VOLK_SC
        BA_SP=BA_SP+TREE(ISIZ,1)*MSBA(ISIZ)*SPCLT(ISIZ)
        SPCL_TPA = SPCL_TPA + SPCLT(ISIZ)*TREE(ISIZ,1)
   20 CONTINUE
      DO 300 J=1,MXDWHC
         DO 320 K=1,MXDWHZ
            IPS_SLSH = IPS_SLSH + DWPHOS(J,K)
  320    CONTINUE
  300 CONTINUE
      WRITE(*,'(A,F16.7)') 'PreDispBKP ', OLDBKP
      WRITE(*,'(A,F16.7)') 'PostDispBKP', BKPA
      WRITE(*,'(A,F16.7)') 'StandRV    ', GRFSTD
      WRITE(*,'(A,F16.7)') 'StandBA    ', BASTD
      WRITE(*,'(A,F16.7)') 'BAH        ', BAH(3)+BAH(5)
      WRITE(*,'(A,F16.7)') 'BA_BtlKld  ', BAK_YR
      WRITE(*,'(A,F16.7)') 'TPA        ', TPA_YR
      WRITE(*,'(A,F16.7)') 'TPAH       ', TREE(3,1)+TREE(5,1)
      WRITE(*,'(A,F16.7)') 'TPA_BtlKld ', TPAK_YR
      WRITE(*,'(A,F16.7)') 'StandVol   ', VOL_YR
      WRITE(*,'(A,F16.7)') 'VolHost    ', VOLH_YR
      WRITE(*,'(A,F16.7)') 'VolBtlKld  ', VOLK_YR
      WRITE(*,'(A,F16.7)') 'BA_Special ', BA_SP
      WRITE(*,'(A,F16.7)') 'SpclTPA    ', SPCL_TPA
      WRITE(*,'(A,F16.7)') 'Ips_Slash  ', IPS_SLSH
C     per-size-class arrays for FVS_BM_Tree / FVS_BM_Vol golden (classes 3 and 5)
      DO 400 ISIZ=1,NSCL
        TPAKLL = PBKILL(ISIZ) + ALLKLL(ISIZ)
        TPA_SC = TREE(ISIZ,1) + TREE(ISIZ,2)
        TVOL_SC = (TVOL(ISIZ,1)*TREE(ISIZ,1) + TVOL(ISIZ,2)*TREE(ISIZ,2))
        HVOL_SC = TVOL(ISIZ,1)*TREE(ISIZ,1)
        VOLK_SC = TPAKLL * TVOL(ISIZ,1)
        IF (ISIZ.EQ.3 .OR. ISIZ.EQ.5) THEN
          WRITE(*,'(A,I2,5F14.6)') 'SC',ISIZ,TPA_SC,TREE(ISIZ,1),TPAKLL,
     >       SPCLT(ISIZ)*TREE(ISIZ,1),TVOL_SC
          WRITE(*,'(A,I2,2F14.6)') 'VL',ISIZ,HVOL_SC,VOLK_SC
        ENDIF
  400 CONTINUE
      END
