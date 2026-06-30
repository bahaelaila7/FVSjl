import re
def logical_statements(path):
    """Join Fortran fixed-form lines into logical statements (handles col-6 continuation)."""
    out=[]; cur=""
    for raw in open(path, errors='replace'):
        line=raw.rstrip('\n')
        if not line.strip(): continue
        c1=line[0] if line else ' '
        if c1 in 'Cc*!': continue                      # full-line comment
        cont = len(line)>5 and line[5] not in (' ','0')  # col-6 continuation marker
        body = line[6:] if len(line)>6 else ''          # statement is cols 7+
        # also keep cols 1-5 region for the DATA keyword on a starting line (label area usually blank)
        start = (line[:5] if len(line)>=5 else line)
        if cont:
            cur += body
        else:
            if cur: out.append(cur)
            cur = (start + ' ' + body) if start.strip() else body
            # if the keyword sits in cols 1-6 region (e.g. "      DATA"), include from col 1
            if 'DATA' not in cur.upper() and 'DATA' in line.upper():
                cur = line[6:] if len(line)>6 else line
            if re.search(r'\bDATA\b', line, re.I):
                cur = line  # keep the whole physical line for DATA starts
    if cur: out.append(cur)
    return out

def grab(path,name):
    stmts=logical_statements(path)
    for s in stmts:
        m=re.search(rf'\bDATA\s+{name}\s*/(.*?)/', s, re.I|re.S)
        if m:
            body=m.group(1)
            vals=[]
            for t in re.split(r'[,\s]+', body.strip()):
                if not t: continue
                mm=re.match(r"(\d+)\*(.+)", t)
                if mm: vals+=[float(mm.group(2))]*int(mm.group(1))
                else:
                    try: vals.append(float(t))
                    except: pass
            return vals
    return None
