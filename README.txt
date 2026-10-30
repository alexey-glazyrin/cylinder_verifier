START HERE

1. Read the manuscript "Upper bounds for packings of cylinders around the
   unit sphere" (A. Glazyrin), including its verification appendix.
2. Read documentation.pdf for the mathematical explanation and
   the precise connection between the manuscript and the code.
3. Use Python 3.10 or later, standard library only, from this directory:

       python check_integrity.py
       python run_verification.py 7

   To verify every case N=4,...,24:

       python run_verification.py --jobs 4

   The jobs option runs several independent cases at once; omit it to run
   one case at a time. Long-running overlap checks may take substantial time;
   see MEASURED RUN TIMES below. The runner creates reports/ and logs/.
   Do NOT enable Python optimization (-O, -OO, or PYTHONOPTIMIZE).

WHAT IS PROVED

For each exact rational radius r and quadrilateral in Appendix A:
  * no triple intersection of the associated closed spherical regions;
  * connectedness of every directed difference Q_i minus Q_j;
  * pairwise intersection area at most the exact rational delta;
  * positive area margin N*A-(3*N-6)*delta-4*pi.

The mathematical theorem in the manuscript turns these verified hypotheses
into R_N < r. This is an exact-arithmetic computational certificate, not
a formal proof-assistant development. The proof decisions use integers and
rational numbers. Decimal diagnostic fields and elapsed times are not inputs
to any proof decision.

CONTENTS

bounds/             Checkers and data for N=4,...,12.
bounds_13_24/       Checkers and data for N=13,...,24.
parameters.csv      Exact r,u,v,rho,delta and lower angular endpoint for all N.
run_verification.py Runs complete case verifications and writes reports/ and logs/.
check_integrity.py  Checks the SHA-256 hashes listed in MANIFEST.sha256.json.
documentation.pdf   Mathematical documentation of the certificates.
LICENSE             Creative Commons Attribution 4.0 International (CC BY 4.0).

LICENSE

Copyright (c) 2026 Alexey Glazyrin. This package (code, data, and
documentation) is released under the Creative Commons Attribution 4.0
International License; see LICENSE. Please cite the manuscript when
using these certificates.

Each case JSON file names one overlap partition and one connectivity
partition. The historical name "cutoff" means connectivity in this release:
all four pieces must be marked. No small sliver is discarded or charged an
area penalty. The archive contains no historical unit-radius or optimization
experiments.

READING SUCCESS AND FAILURE

A successful full run prints "All requested cases passed COMPLETE
verification." A case report must have shard=null, checked_boxes equal to
overlap_boxes, cutoff_bound="0", and an empty failures list. The runner
requires these conditions. A successful partial shard is NOT a full proof.

Any failed assertion, nonzero exit code, missing report, or incomplete case
must be treated as failure. Inspect logs/nN.log. Never suppress an assertion
or change a numerical tolerance to obtain success. The checker is deliberately
conservative: subdivision may be required if parameters are changed.

MEASURED RUN TIMES

Complete run of all 21 cases on 30 October 2026: Python 3.13.5, Windows 11,
standard-library fractions backend, eight cases at a time; about 30 minutes
of wall-clock time in total. Seconds per case:

   N=4    288    N=11  1155    N=18   878
   N=5    580    N=12   132    N=19   644
   N=6    281    N=13   530    N=20   108
   N=7   1815    N=14   195    N=21   104
   N=8    867    N=15   208    N=22   226
   N=9    273    N=16   656    N=23   107
   N=10   478    N=17   207    N=24   388

Every case passed with checked_boxes equal to overlap_boxes, an empty
failures list, and all connectivity boxes fully marked.