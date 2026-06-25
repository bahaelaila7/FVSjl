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
  'THINBBA': [(1,'year','i'),(2,'residual_basal_area','f'),(3,'cut_efficiency','f'),
              (4,'dbh_min','f'),(5,'dbh_max','f'),(6,'species','i'),(7,'plot','i')],
  'THINSDI': [(1,'year','i'),(2,'residual_sdi','f'),(3,'cut_efficiency','f'),
              (4,'dbh_min','f'),(5,'dbh_max','f'),(6,'species','i'),(7,'plot','i')],
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
            params = []
            for idx, pname, typ in SCHEMA[name]:
                v = field(ln, idx)
                if v == '': continue
                if typ == 's': params.append('      %s: "%s"' % (pname, v))
                else:          params.append('      %s: %s' % (pname, emit_num(v, typ)))
            out.append('  - %s:' % name)
            out += params if params else ['      {}']
            i += 1; continue
        out.append('  - %s: {}  # (no schema)' % name); i += 1
    return '\n'.join(out) + '\n'

print(main(sys.argv[1]), end='')
