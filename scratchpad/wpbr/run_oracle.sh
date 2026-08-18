#!/bin/bash
# Run stock FVSie + relinked FVSie_wpbr with a SEPARATE tree-data file (unit 02).
# stdin: line1 = keyword file, line2 = tree data file.
set -u
cd "$(dirname "$0")" || exit 2
STOCK=/workspace/ForestVegetationSimulator/bin/FVSie
WPBR=./FVSie_wpbr

cat > ie.tre <<'EOF'
   1      248112       0101   011WP 12014   0654   00111     0  0
   2      248112       0101   011WP 08010   0504   00111     0  0
   3      248112       0102   011DF 11014   0634   00111     0  0
   4      248112       0102   011WP 15019   0734   00111     0  0
   5      248112       0103   011DF 09011   0524   00111     0  0
EOF

mkhead () {  # $1 = tag/title -> keyfile body up to TREEDATA
cat <<EOF
SCREEN
NOAUTOES
NOTRIPLE
STATS
STDIDENT
S248112  IE WPBR $1
DESIGN                                        11.0       1.0
STDINFO     11406001     570.0      60.0     315.0      30.0      34.0
INVYEAR       1990.0
NUMCYCLE         5.0
TREEFMT
(T24,I4,T1,I4,T31,F2.0,I1,A3,F3.1,F2.1,T45,F3.0,T63,F3.0,T60,F3.1,T48,I1,
T52,I2,T66,5I1,T54,7I1,T75,F3.0)
TREEDATA
EOF
}

mkhead OFF > ie_off2.key
printf 'ECHOSUM\nPROCESS\nSTOP\n' >> ie_off2.key

mkhead ON > ie_on2.key
printf 'BRUST\nRUSTINDX       0.05\nEND\nECHOSUM\nPROCESS\nSTOP\n' >> ie_on2.key

printf 'ie_off2.key\nie.tre\n' | "$STOCK" > s2.log 2>&1
[ -f ie_off2.sum ] && mv -f ie_off2.sum stock_off2.sum
printf 'ie_off2.key\nie.tre\n' | "$WPBR" > w2off.log 2>&1
[ -f ie_off2.sum ] && mv -f ie_off2.sum relink_off2.sum
printf 'ie_on2.key\nie.tre\n' | "$WPBR" > w2on.log 2>&1
[ -f ie_on2.sum ] && mv -f ie_on2.sum relink_on2.sum

echo "=== TPA loaded? (stock off, row2) ==="; sed -n '2p' stock_off2.sum
echo "=== OFF: stock vs relink (WPBR-off) ==="
if diff stock_off2.sum relink_off2.sum >/dev/null; then echo IDENTICAL; else diff stock_off2.sum relink_off2.sum; fi
echo "=== ON vs OFF (drop timestamp line1) ==="
if diff <(tail -n +2 relink_off2.sum) <(tail -n +2 relink_on2.sum) >/dev/null; then echo "IDENTICAL (WPBR did NOT change projection)"; else echo "DIFFERS (WPBR engaged):"; diff <(tail -n +2 relink_off2.sum) <(tail -n +2 relink_on2.sum); fi
