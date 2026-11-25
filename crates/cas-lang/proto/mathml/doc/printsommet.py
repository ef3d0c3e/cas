"""printsommet.py <giac src dir> -- list giac functions registered with a printsommet / texprint (src grep).

Output is printsommet.txt (giac 2.0.0.21, source fetched per C0009); mathml_proto.py reads that file."""
import glob, re, sys
SRC = sys.argv[1]
pat = re.compile(r'define_unary_function_(eval[0-9a-z_]*)\s*\(([^;]*)\);')
out = {}
for f in glob.glob(f'{SRC}/*.cc') + glob.glob(f'{SRC}/more/*.cc'):
    for m in pat.finditer(open(f, encoding='latin-1').read()):
        kind, args = m.group(1), [a.strip() for a in m.group(2).split(',')]
        if 'index' in kind:
            args = args[1:]
        name = re.sub(r'^(unary)?__', '', args[0])
        ps = tx = None
        if kind.startswith('eval2'):
            ps = args[3]
        elif kind.startswith('eval4'):
            ps, tx = args[3], args[4] if len(args) > 4 else None
        elif kind.startswith('eval5'):
            ps, tx = args[4], args[5]
        elif kind.startswith('eval_taylor2'):
            ps, tx = args[-2], args[-1]
        if (ps and ps != '0') or (tx and tx != '0'):
            out[name] = (ps, tx, f.split('/')[-1])
for k in sorted(out):
    print(k, *out[k])
