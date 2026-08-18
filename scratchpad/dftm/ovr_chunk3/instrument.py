#!/usr/bin/env python3
# Instrument pristine tmcoup.f with Float32-hex TRANSFER dumps at 4 anchors.
src = "/workspace/ForestVegetationSimulator/dftm/tmcoup.f"
out = "/tmp/claude-1000/-workspace/b4e1b3b1-495b-403e-810b-5db3604b56cc/scratchpad/dftmwork/ovr/tmcoup.f"
lines = open(src).read().split("\n")

DUMP_A = r"""C=== DFTM DUMP-A classification ===
      WRITE(0,'(A,6I7)') 'DBGA_HDR ',ITRN,IFIN,NACLAS(1),
     &   NACLAS(2),IDFCOD,IGFCOD
      WRITE(0,'(A,4I7)') 'DBGA_ISCT',ISCT(IDFCOD,1),
     &   ISCT(IDFCOD,2),ISCT(IGFCOD,1),ISCT(IGFCOD,2)
      DO 9011 I=1,ITRN
        WRITE(0,'(A,I4,I8)') 'DBGA_IPT ',I,IPT(I)
 9011 CONTINUE
      DO 9012 I=1,ITRN
        WRITE(0,'(A,I4,3(1X,Z8.8))') 'DBGA_REC ',I,
     &     TRANSFER(PROB(I),1),TRANSFER(PCNEWF(I),1),
     &     TRANSFER(FBIOMS(I),1)
 9012 CONTINUE
      DO 9013 I=1,IFIN
        WRITE(0,'(A,I4,2I8)') 'DBGA_ISC ',I,ISC(I,1),ISC(I,2)
 9013 CONTINUE"""

DUMP_B = r"""C=== DFTM DUMP-B ICON/eggs ===
      WRITE(0,'(A,I3,4(1X,Z8.8))') 'DBGB_HDR ',IEGTYP,
     &   TRANSFER(DFMEAN,1),TRANSFER(DFREGG(2),1),
     &   TRANSFER(GFMEAN,1),TRANSFER(GFREGG(2),1)
      DO 9021 I=1,IFIN
        WRITE(0,'(A,I4,I3,I4,7(1X,Z8.8))') 'DBGB_ICN ',
     &     I,IZ6(I),JCLAS2(I),TRANSFER(Z4(I),1),
     &     TRANSFER(Z5(I),1),TRANSFER(Z2(I),1),
     &     TRANSFER(Z3(I),1),TRANSFER(X5(I),1),
     &     TRANSFER(X6(I),1),TRANSFER(X7(I),1)
 9021 CONTINUE"""

DUMP_TK = r"""      WRITE(0,'(A,I6,1X,Z8.8)') 'DBGTK_PR ',IP,
     &   TRANSFER(PRTOPK,1)"""

DUMP_C = r"""C=== DFTM DUMP-C final treelist ===
      DO 9031 II=1,ITRN
        I=IPT(II)
        WRITE(0,'(A,I5,3I4,I8,5(1X,Z8.8))') 'DBGC_TRE ',
     &     I,IMC(I),ICR(I),NORMHT(I),ITRUNC(I),
     &     TRANSFER(WK2(I),1),TRANSFER(HT(I),1),
     &     TRANSFER(HTG(I),1),TRANSFER(DBH(I),1),
     &     TRANSFER(DG(I),1)
 9031 CONTINUE"""

# Insert bottom-up so line numbers stay valid.
# anchors (1-based line content):
#  line 311 '      ENDIF' (end GF block)  -> insert DUMP_A after
#  line 407 '  170 CONTINUE'             -> insert DUMP_B after
#  line 682 '          CALL TMRANN(PRTOPK)' -> insert DUMP_TK after
#  line 910 '  430 CONTINUE'             -> insert DUMP_C after
def check(idx, needle):
    assert needle in lines[idx], f"line {idx+1}={lines[idx]!r} lacks {needle!r}"
check(310, "ENDIF")
check(406, "170 CONTINUE")
check(681, "CALL TMRANN(PRTOPK)")
check(909, "430 CONTINUE")

lines[909:910] = [lines[909], DUMP_C]
lines[681:682] = [lines[681], DUMP_TK]
lines[406:407] = [lines[406], DUMP_B]
lines[310:311] = [lines[310], DUMP_A]

open(out, "w").write("\n".join(lines))
print("wrote", out, "lines", len(lines))
