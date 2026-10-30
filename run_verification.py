"""Run complete certificates, retaining subprocess logs and fresh reports."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
import argparse
import datetime
import hashlib
import json
import os
import platform
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parent

def run_case(n, backend):
    folder = ROOT/('bounds' if n<=12 else 'bounds_13_24')
    report = ROOT/'reports'/f'n{n}_verification.json'
    if report.exists():
        report.unlink()
    started = time.time()
    command = [sys.executable, '-u', str(folder/'verify_bounds.py'), str(n)]
    env=os.environ.copy()
    env['CYLINDER_EXACT_BACKEND']=backend
    with (ROOT/'logs'/f'n{n}.log').open('w', encoding='utf-8') as log:
        completed = subprocess.run(command, cwd=ROOT, env=env, stdout=log, stderr=subprocess.STDOUT)
    result={'n':n, 'returncode':completed.returncode, 'seconds':time.time()-started, 'passed':False}
    if completed.returncode==0 and report.exists():
        data=json.loads(report.read_text())
        result['passed']=(data['shard'] is None and data['checked_boxes']==data['overlap_boxes']
                          and not data['failures'] and data['cutoff_bound']=='0')
        result['overlap_boxes']=data['overlap_boxes']
        result['connectivity_boxes']=data['cutoff_boxes']
    print(('PASS' if result['passed'] else 'FAIL')+f' N={n}, {result["seconds"]:.1f}s', flush=True)
    return result

def main():
    if not __debug__:
        raise SystemExit('Run without -O or -OO.')
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('n', type=int, nargs='*', help='Cases 4 through 24; default: all')
    parser.add_argument('--jobs', type=int, default=1, help='Simultaneous complete case verifications')
    parser.add_argument('--backend', choices=['fractions','gmpy2'], default='fractions',
                        help='Exact rational arithmetic; gmpy2 is an optional installed accelerator')
    args=parser.parse_args()
    cases=sorted(set(args.n or range(4,25)))
    if not all(4<=n<=24 for n in cases) or args.jobs<1:
        parser.error('Use 4<=N<=24 and jobs>=1.')
    (ROOT/'reports').mkdir(exist_ok=True)
    (ROOT/'logs').mkdir(exist_ok=True)
    started=datetime.datetime.now(datetime.timezone.utc).isoformat()
    result={'started_utc':started, 'python':sys.version, 'platform':platform.platform(),
            'cases_requested':cases, 'jobs':args.jobs, 'arithmetic_backend':args.backend,
            'complete':False, 'results':[]}
    output=ROOT/'reports'/('run_all.json' if cases==list(range(4,25)) else 'run_selected.json')
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futures=[pool.submit(run_case,n,args.backend) for n in cases]
        for future in as_completed(futures):
            result['results'].append(future.result())
            output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    result['results'].sort(key=lambda x:x['n'])
    result['complete']=len(result['results'])==len(cases) and all(x['passed'] for x in result['results'])
    result['finished_utc']=datetime.datetime.now(datetime.timezone.utc).isoformat()
    result['checker_sha256']={str(p.relative_to(ROOT)).replace('\\','/'):hashlib.sha256(p.read_bytes()).hexdigest()
                             for group in ('bounds','bounds_13_24') for p in sorted((ROOT/group).glob('*.py'))}
    output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    if not result['complete']:
        raise SystemExit('FAIL: at least one complete case verification failed; see logs.')
    print('All requested cases passed COMPLETE verification.',flush=True)

if __name__=='__main__':
    main()
