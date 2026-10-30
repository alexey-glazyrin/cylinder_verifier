"""Select exact rational arithmetic; no floating-point backend is allowed."""
import os

BACKEND = os.environ.get('CYLINDER_EXACT_BACKEND', 'fractions')
if BACKEND == 'fractions':
    from fractions import Fraction as F
    VERSION = 'Python standard library'
elif BACKEND == 'gmpy2':
    import gmpy2
    F = gmpy2.mpq
    VERSION = gmpy2.version()
else:
    raise RuntimeError('Unknown exact rational backend: ' + BACKEND)
