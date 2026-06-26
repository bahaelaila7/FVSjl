#!/usr/bin/env python3
# key_to_structured_yaml.py — render an FVS .key as an ORDER-AWARE, HIERARCHICAL
# yaml. The keyword stream stays an ordered sequence (FVS order is significant);
# the hierarchy only GROUPS it into named sections for readability, and the grouping
# is order-PRESERVING (consecutive same-section keywords form a block; blocks appear
# in the ORIGINAL order, so flattening top-to-bottom reproduces the exact sequence).
# Each keyword becomes a block with NAMED, typed params (numbers stay numbers,
# strings are quoted). Mirrors src/io/yaml_keywords.jl.
#   usage:  python3 key_to_structured_yaml.py <stand.key> [--flat] > stand.yaml
import sys

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
  'LEAVESP': [(1,'species','i')],
  'SPLEAVE': [(1,'species','i')],
  'SPECPREF':[(1,'species','i'),(2,'preference','i')],
  'CUTEFF':  [(1,'proportion','f')],
  'MINHARV': [(1,'min_volume','f')],
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
  'SPGROUP': [(1,'group','s')],
}
NOPARAM = {'TREEDATA','PROCESS','STOP'}

# keyword -> section (order-preserving grouping; mirrors _KW_SECTION in yaml_keywords.jl)
SECTION = {}
def _sec(names, s):
    for n in names.split(): SECTION[n] = s
_sec('STDIDENT STDINFO DESIGN SITECODE INVYEAR NUMCYCLE TIMEINT CYCLEAT MANAGED '
     'RESETAGE TFIXAREA TREEFMT TREEDATA NOTREES NOAUTOES PROCESS STOP', 'setup')
_sec('SDIMAX SDICALC BAMAX', 'density')
_sec('GROWTH NOCALIB READCORD REUSCORD READCORH REUSCORH READCORR REUSCORR BAIMULT '
     'HTGMULT CRNMULT REGDMULT REGHMULT DGSTDEV SERLCORR RANNSEED FIXDG FIXHTG FIXMORT '
     'MORTMULT TREESZCP SIZCAP', 'growth')
_sec('THINBBA THINABA THINBTA THINATA THINSDI THINCC THINHT THINQFA THINRDEN THINDBH '
     'THINPT SETPTHIN THINPRSC THINAUTO SPECPREF LEAVESP SPLEAVE CUTEFF MINHARV SALVAGE '
     'YARDLOSS', 'treatments')
_sec('ESTAB PLANT NATURAL SPROUT NOSPROUT', 'regeneration')
_sec('VOLUME BFVOLUME VOLEQNUM MCDEFECT BFDEFECT MCFDLN BFFDLN', 'volume')
_sec('SPGROUP', 'species_groups')
_sec('SETSITE FERTILIZ', 'site')
_sec('ECHOSUM DATABASE DSNOUT SUMMARY TREELIDB CUTLIST STRCLASS COMPUTDB', 'output')
_sec('COMPUTE IF', 'event_monitor')
_sec('COMPRESS NOTRIPLE NUMTRIP', 'control')

def section_of(name):
    return SECTION.get(name, 'other')

def field(rec, i):                       # 10-col field i (1-based), cols 10*i+1..10*i+10
    s = rec[10*i:10*i+10] if len(rec) > 10*i else ''
    return s.strip()

def emit_scalar(v, typ):
    # a non-numeric field (e.g. an alpha species code "LP" in PLANT) is quoted as a
    # string even where the schema expects a number, so the form stays lossless.
    if typ == 's':
        return '"%s"' % v
    try:
        f = float(v)
    except ValueError:
        return '"%s"' % v
    if typ == 'i' or f == int(f):
        return str(int(f))
    return ('%g' % f)

def is_plain(name):                      # name is a bare keyword (letters/digits)?
    return name != '' and name[0].isalpha() and name.replace(' ', '') == name and name.isalnum()

def parse_entries(path):
    """Parse the .key into a list of (kind, payload) entries in order.
       kind: 'raw' -> verbatim line ; 'kw' -> (NAME, [(pname,val_yaml)]) ; 'plain' -> NAME"""
    lines = open(path).read().replace('\r\n', '\n').replace('\r', '\n').split('\n')
    ents = []
    i = 0
    while i < len(lines):
        ln = lines[i]; name = ln[:8].strip().upper()
        if not name:
            i += 1; continue
        stripped = ln.strip()
        if stripped and not stripped[0].isalpha():   # free-form line (IF cond, etc.)
            ents.append(('raw', stripped)); i += 1; continue
        # a "name" with embedded spaces (e.g. the SPGROUP member list "SM HI") is not a
        # bare keyword — carry it verbatim, matching _is_plain_keyword in yaml_keywords.jl
        if not is_plain(name) and name not in SCHEMA and name not in NOPARAM:
            ents.append(('raw', stripped)); i += 1; continue
        if name == 'STDIDENT':
            ident = lines[i+1].strip() if i+1 < len(lines) else ''
            ents.append(('trail', ('STDIDENT', 'id', ident))); i += 2; continue
        if name == 'TREEFMT':
            fmt = ''; j = i + 1
            while j < len(lines):
                fmt += lines[j].strip()
                if fmt.endswith(')'): j += 1; break
                j += 1
            ents.append(('trail', ('TREEFMT', 'format', fmt))); i = j; continue
        if name in NOPARAM:
            ents.append(('plain', name)); i += 1; continue
        if name in SCHEMA:
            mapped = {idx for idx, _, _ in SCHEMA[name]}
            present = [k for k in range(1, 13) if field(ln, k) != '']
            # named form only if EVERY present field maps to a slot — else fall back to
            # the positional `params:` list so no field is silently dropped (lossless),
            # mirroring `named_ok` in src/io/yaml_keywords.jl.
            if all(k in mapped for k in present):
                params = []
                for idx, pname, typ in SCHEMA[name]:
                    v = field(ln, idx)
                    if v == '': continue
                    params.append((pname, emit_scalar(v, typ)))
                ents.append(('kw', (name, params))); i += 1; continue
            last = present[-1] if present else 0
            ents.append(('pos', (name, [field(ln, k) for k in range(1, last + 1)])))
            i += 1; continue
        # unknown keyword: keep any present fields positionally (lossless)
        present = [k for k in range(1, 13) if field(ln, k) != '']
        if present:
            ents.append(('pos', (name, [field(ln, k) for k in range(1, present[-1] + 1)])))
        else:
            ents.append(('plain', name))
        i += 1
    return ents

def render_entry(kind, payload, ind):
    out = []
    if kind == 'raw':
        out.append('%s- raw: "%s"' % (ind, payload.replace('"', '\\"')))
    elif kind == 'plain':
        out.append('%s- %s: {}' % (ind, payload))
    elif kind == 'trail':
        name, key, val = payload
        out.append('%s- %s:' % (ind, name))
        out.append('%s    %s: "%s"' % (ind, key, val.replace('"', '\\"')))
    elif kind == 'kw':
        name, params = payload
        out.append('%s- %s:' % (ind, name))
        if params:
            for pname, val in params:
                out.append('%s    %s: %s' % (ind, pname, val))
        else:
            out.append('%s    {}' % ind)
    elif kind == 'pos':                     # positional fallback (unmapped fields)
        name, vals = payload
        out.append('%s- keyword: "%s"' % (ind, name))
        out.append('%s  params: [%s]' % (ind, ', '.join('"%s"' % v for v in vals)))
    return out

def entry_section(kind, payload, cur):
    # named keyword -> its section; raw/continuation lines stay in the current block
    if kind in ('kw', 'trail', 'pos'):
        return section_of(payload[0])
    if kind == 'plain':
        return section_of(payload)
    return cur            # 'raw'

def main(path, flat):
    ents = parse_entries(path)
    if flat:
        out = ['# FVS keywords — order is significant. Flat ordered list.', 'keywords:']
        for kind, payload in ents:
            out += render_entry(kind, payload, '  ')
        return '\n'.join(out) + '\n'
    out = ['# FVS keywords — ORDER-AWARE hierarchical form. The keyword stream is',
           '# an ordered sequence; sections only GROUP it (flattening top-to-bottom',
           '# reproduces the exact order). `treatments`, `species_groups` and',
           '# `event_monitor` are order-significant (define-before-use, same-cycle',
           '# activities run in input order). See docs/KEYWORDS.md.',
           'stand:']
    cur = None
    for kind, payload in ents:
        sec = entry_section(kind, payload, cur)
        if sec != cur:
            out.append('  - %s:' % sec); cur = sec
        out += render_entry(kind, payload, '      ')
    return '\n'.join(out) + '\n'

if __name__ == '__main__':
    args = [a for a in sys.argv[1:] if not a.startswith('-')]
    flat = '--flat' in sys.argv
    print(main(args[0], flat), end='')
