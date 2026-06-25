#!/usr/bin/env python3
# key_to_structured_yaml.py — render an FVS .key as a STRUCTURED yaml: every
# keyword becomes a block with NAMED, typed parameters (numbers stay numbers,
# strings are quoted). Schema for the supported keywords is the SCHEMA dict below.
#   usage:  python3 key_to_structured_yaml.py <stand.key>  > stand.yaml
import sys, re

# keyword -> ordered (field_index_1based, name, type) ; STDIDENT/TREEFMT take following raw line
SCHEMA = {
  'DESIGN':  [(1,'basal_area_factor','f'),(2,'fixed_plot_area_inverse','f'),(3,'break_dbh','f'),
              (4,'number_of_plots','i'),(5,'nonstockable_code','i'),(6,'sample_weight','f'),
              (7,'stockable_proportion','f')],
  'STDINFO': [(1,'forest_code','i'),(2,'habitat','s'),(3,'stand_age','f'),(4,'aspect','f'),
              (5,'slope','f'),(6,'elevation','f'),(9,'stand_origin','i')],
  'SITECODE':[(1,'site_species','i'),(2,'site_index','f')],
  'INVYEAR': [(1,'year','i')],
  'NUMCYCLE':[(1,'cycles','i')],
  'TIMEINT': [(1,'cycle','i'),(2,'length','i')],
  'CYCLEAT': [(1,'year','i')],
  'MANAGED': [(1,'year','i')],
  'TFIXAREA':[(1,'area','f')],
  'NUMTRIP': [(1,'count','i')],
  'THINBBA': [(1,'year','i'),(2,'residual_basal_area','f'),(3,'cut_efficiency','f'),(4,'dbh_min','f'),(5,'dbh_max','f'),(6,'species','i'),(7,'plot','i')],
  'THINABA': [(1,'year','i'),(2,'residual_basal_area','f'),(3,'cut_efficiency','f'),(4,'dbh_min','f'),(5,'dbh_max','f'),(6,'species','i'),(7,'plot','i')],
  'THINBTA': [(1,'year','i'),(2,'residual_tpa','f'),(3,'cut_efficiency','f'),(4,'dbh_min','f'),(5,'dbh_max','f'),(6,'species','i'),(7,'plot','i')],
  'THINATA': [(1,'year','i'),(2,'residual_tpa','f'),(3,'cut_efficiency','f'),(4,'dbh_min','f'),(5,'dbh_max','f'),(6,'species','i'),(7,'plot','i')],
  'THINSDI': [(1,'year','i'),(2,'residual_sdi','f'),(3,'cut_efficiency','f'),(4,'dbh_min','f'),(5,'dbh_max','f'),(6,'species','i'),(7,'plot','i')],
  'THINCC':  [(1,'year','i'),(2,'residual_ccf','f'),(3,'cut_efficiency','f'),(4,'dbh_min','f'),(5,'dbh_max','f'),(6,'species','i'),(7,'plot','i')],
  'THINDBH': [(1,'year','i'),(2,'dbh_min','f'),(3,'dbh_max','f'),(4,'cut_efficiency','f'),(6,'residual_tpa','f'),(7,'species','i')],
  'THINHT':  [(1,'year','i'),(2,'ht_min','f'),(3,'ht_max','f'),(4,'cut_efficiency','f'),(6,'residual_tpa','f'),(7,'species','i')],
  'THINRDEN':[(1,'year','i'),(2,'residual_relsdi','f'),(3,'cut_efficiency','f'),(4,'dbh_min','f'),(5,'dbh_max','f'),(6,'species','i'),(7,'plot','i')],
  'THINAUTO':[(1,'year','i'),(2,'cut_efficiency','f')],
  'SERLCORR':[(1,'phi','f'),(2,'theta','f')],
  'RESETAGE':[(1,'year','i'),(2,'age','f')],
  'VOLUME':  [(1,'year','i'),(2,'species','i'),(3,'dbh_min','f'),(4,'top_diam','f'),(5,'stump','f'),(6,'form_class','f'),(7,'method','i'),(8,'scf_min_dbh','f'),(9,'scf_top_dib','f'),(10,'scf_stump','f')],
  'BFVOLUME':[(1,'year','i'),(2,'species','i'),(3,'bf_min_dbh','f'),(4,'bf_top_dib','f'),(5,'bf_stump','f')],
  'MCDEFECT':[(1,'year','i'),(2,'species','i'),(3,'defect_5','f'),(4,'defect_10','f'),(5,'defect_15','f'),(6,'defect_20','f'),(7,'defect_25','f')],
  'BFDEFECT':[(1,'year','i'),(2,'species','i'),(3,'defect_5','f'),(4,'defect_10','f'),(5,'defect_15','f'),(6,'defect_20','f'),(7,'defect_25','f')],
  'BAIMULT': [(1,'year','i'),(2,'species','i'),(3,'multiplier','f'),(4,'dbh_min','f'),(5,'dbh_max','f')],
  'HTGMULT': [(1,'year','i'),(2,'species','i'),(3,'multiplier','f'),(4,'dbh_min','f'),(5,'dbh_max','f')],
  'CRNMULT': [(1,'year','i'),(2,'species','i'),(3,'multiplier','f')],
  'REGDMULT':[(1,'year','i'),(2,'species','i'),(3,'multiplier','f')],
  'REGHMULT':[(1,'year','i'),(2,'species','i'),(3,'multiplier','f')],
  'MORTMULT':[(1,'year','i'),(2,'species','i'),(3,'multiplier','f'),(4,'dbh_min','f'),(5,'dbh_max','f')],
  'FIXMORT': [(1,'year','i'),(2,'species','i'),(3,'rate','f'),(4,'dbh_min','f'),(5,'dbh_max','f'),(6,'option','i')],
  'FIXDG':   [(1,'year','i'),(2,'species','i'),(3,'value','f'),(4,'dbh_min','f'),(5,'dbh_max','f')],
  'FIXHTG':  [(1,'year','i'),(2,'species','i'),(3,'value','f'),(4,'dbh_min','f'),(5,'dbh_max','f')],
  'BAMAX':   [(1,'bamax','f')],
  'SDIMAX':  [(1,'species','i'),(2,'sdimax','f'),(5,'pct_lo','f'),(6,'pct_hi','f')],
  'SDICALC': [(1,'method','i')],
  'FERTILIZ':[(1,'year','i'),(2,'nitrogen','f')],
  'SETSITE': [(1,'year','i'),(2,'habitat','s'),(3,'bamax','f'),(4,'species','i'),(5,'site_index','f'),(6,'si_flag','i'),(7,'sdimax','f')],
  'RANNSEED':[(1,'seed','i')],
  'NOCALIB': [(1,'species','i')],
  'DGSTDEV': [(1,'value','f')],
  'PLANT':   [(1,'year','i'),(2,'species','i'),(3,'tpa','f'),(4,'survival_pct','f'),(5,'age','i'),(6,'height','f'),(7,'shade','f')],
  'NATURAL': [(1,'year','i'),(2,'species','i'),(3,'tpa','f'),(4,'survival_pct','f'),(5,'age','i'),(6,'height','f'),(7,'shade','f')],
  'ESTAB':   [(1,'disturbance_date','i')],
  'VOLEQNUM':[(1,'species','i'),(2,'equation','i')],
}
NOPARAM = {'TREEDATA','PROCESS','STOP'}

def field(rec, i):                       # 10-col field i (1-based), cols 10*i+1..10*i+10
    s = rec[10*i:10*i+10] if len(rec) > 10*i else ''
    return s.strip()

def emit_num(v, typ):
    f = float(v)
    if typ == 'i' or f == int(f):
        return str(int(f))
    return ('%g' % f)

def main(path):
    lines = open(path).read().replace('\r\n','\n').replace('\r','\n').split('\n')
    out = ['# FVS keywords — structured form. Order is significant (list).',
           '# Numbers are numbers; strings are quoted. See README for the schema.',
           'keywords:']
    i = 0
    while i < len(lines):
        ln = lines[i]; name = ln[:8].strip().upper()
        if not name:
            i += 1; continue
        stripped = ln.strip()
        if stripped and not stripped[0].isalpha():   # free-form line (e.g. an IF condition)
            out.append('  - raw: "%s"' % stripped.replace('"','\\"')); i += 1; continue
        if name == 'STDIDENT':
            ident = lines[i+1].strip() if i+1 < len(lines) else ''
            out.append('  - STDIDENT:'); out.append('      id: "%s"' % ident); i += 2; continue
        if name == 'TREEFMT':
            fmt=''; j=i+1
            while j < len(lines):
                fmt += lines[j].strip()
                if fmt.endswith(')'): j+=1; break
                j+=1
            out.append('  - TREEFMT:'); out.append('      format: "%s"' % fmt); i=j; continue
        if name in NOPARAM:
            out.append('  - %s: {}' % name); i += 1; continue
        if name in SCHEMA:
            try:
                params = []
                for idx, pname, typ in SCHEMA[name]:
                    v = field(ln, idx)
                    if v == '': continue
                    if typ == 's': params.append('      %s: "%s"' % (pname, v))
                    else:          params.append('      %s: %s' % (pname, emit_num(v, typ)))
                out.append('  - %s:' % name)
                out += params if params else ['      {}']
                i += 1; continue
            except ValueError:
                pass   # misaligned/alpha field → fall through to a lossless raw passthrough
        # Unknown keyword OR a schema parse failure: emit the whole card verbatim as a raw
        # line so it ALWAYS round-trips (Task 8 adds named schemas / multi-record support).
        out.append('  - raw: "%s"' % ln.rstrip().replace('"','\\"')); i += 1
    return '\n'.join(out) + '\n'

print(main(sys.argv[1]), end='')
