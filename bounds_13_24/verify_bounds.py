"""Exact n=13,...,24 radius certificates. Geometry uses rational arithmetic.

Uses the arithmetic routines in the two earlier checkers, with the
parameters set by radius_case.json. Negative common-normal parameters
are included: they become possible when the cylinder radius is below 1.
"""
from pathlib import Path
from rational import F, BACKEND, VERSION
from functools import lru_cache
from math import isqrt
import json,gzip,time,sys
if not __debug__:
    raise SystemExit('Run without -O or -OO: assertions implement proof checks.')
import exact_geometry as core
import complement_geometry as cuts

HERE=Path(__file__).resolve().parent/'data'
N=int(sys.argv[1])
CASE=json.loads((HERE/f'n{N}_case.json').read_text())
assert N==CASE['n'] and 13<=N<=24
if not __debug__:
    raise SystemExit('Run without python -O; assertions are part of the proof checker.')
RADIUS=F(CASE['r']);U=F(CASE['u']);V=F(CASE['v']);R=F(CASE['R'])
DELTA=F(CASE['delta']);LOWER=F(CASE['angle_lower'])
H=1+RADIUS;K=2*RADIUS/H
D=core.D;AG=core.AG;P=core.ROOT_PI
assert RADIUS>0 and R>0 and DELTA>=0 and U>0 and V>0 and -F(1)<=LOWER<1
core.U=cuts.U=U;core.V=cuts.V=V
core.UI=int(U*D);core.VI=int(V*D)
assert F(core.UI,D)==U and F(core.VI,D)==V
fac=isqrt(D*D+core.UI*core.UI);core.FAC=(fac,fac+1)
positive_trig=core.trig_point

@lru_cache(None)
def signed_trig(x):
    if x<0:
        s,c=positive_trig(-x)
        return core.ni(s),c
    return positive_trig(x)

def trig_range(l,r):
    assert -F(1)<=l<=r<4
    sl,cl=signed_trig(l);sr,cr=signed_trig(r)
    s=(min(sl[0],sr[0]),max(sl[1],sr[1]))
    c=(min(cl[0],cr[0]),max(cl[1],cr[1]))
    if l<=core.PI_HI/2 and r>=core.PI_LO/2:s=(s[0],max(s[1],D))
    if l<=0<=r:c=(c[0],max(c[1],D))
    if l<=core.PI_HI and r>=core.PI_LO:c=(min(c[0],-D),c[1])
    return s,c

core.trig_point=signed_trig
core.trig_range=cuts.trig_range=trig_range


def norm_intervals(box):
    """Enclose every coefficient, including extrema at negative angles.

    A*sin(phi)+B*cos(phi) has critical points pi apart. On an interval
    of length less than pi, endpoint derivative signs locate its only
    possible interior extremum. The critical value is +/-sqrt(1+u^2).
    Ambiguous interval signs only enlarge the enclosure.
    """
    th0,th1,ph0,ph1,w0,w1=box
    assert 0<=ph1-ph0<core.PI_LO
    s,c=trig_range(th0,th1)
    t,d=trig_range(ph0,ph1)
    sw,cw=trig_range(w0,w1)
    t0,d0=signed_trig(ph0)
    t1,d1=signed_trig(ph1)

    def sinusoid(a,b):
        f0=core.ai(core.scale(t0,a),core.scale(d0,b))
        f1=core.ai(core.scale(t1,a),core.scale(d1,b))
        g0=core.si(core.scale(d0,a),core.scale(t0,b))
        g1=core.si(core.scale(d1,a),core.scale(t1,b))
        low=min(f0[0],f1[0])
        high=max(f0[1],f1[1])
        if g0[1]>=0 and g1[0]<=0:
            high=max(high,core.FAC[1])
        if g0[0]<=0 and g1[1]>=0:
            low=min(low,-core.FAC[1])
        return low,high

    ap=sinusoid(core.UI,D)
    bp=sinusoid(D,-core.UI)
    am=sinusoid(core.UI,-D)
    bm=sinusoid(-D,-core.UI)
    ans=[]
    for aa,bb in [(ap,bp),(am,bm)]:
        ans.append([
            core.ai(core.mi(aa,c),core.mi(core.mi(bb,s),cw)),
            core.si(core.mi(aa,s),core.mi(core.mi(bb,c),cw)),
            core.ni(core.mi(bb,sw)),
        ])
    for sign in (1,-1):
        mix=core.ai(core.scale(core.mi(d,cw),core.VI),
                    sw if sign==1 else core.ni(sw))
        ans.append([
            core.si(core.scale(core.mi(t,c),core.VI),core.mi(s,mix)),
            core.ai(core.scale(core.mi(t,s),core.VI),core.mi(c,mix)),
            core.si(core.scale(core.mi(d,sw),core.VI),
                    cw if sign==1 else core.ni(cw)),
        ])
    return ans

core.norm_intervals=cuts.norm_intervals=norm_intervals

def area_bound(poly):
    if len(poly)<3:return 0
    center=(sum(p[0] for p in poly)/len(poly),sum(p[1] for p in poly)/len(poly))
    lc=core.length_lower(center);total=0
    for p,q in zip(poly,poly[1:]+poly[:1]):
        lp=core.length_lower(p);lq=core.length_lower(q)
        den=(lc*lp*lq)//(D*D)+(core.dot_lower(center,p)*lq)//D+(core.dot_lower(p,q)*lc)//D+(core.dot_lower(q,center)*lp)//D
        assert den>0
        det=abs((p[0]-center[0])*(q[1]-center[1])-(p[1]-center[1])*(q[0]-center[0]))
        t=core.ceildiv(det.numerator*D*AG,det.denominator*den)
        # atan(x) <= x-x^3/3+x^5/5; the polynomial is increasing.
        if t<AG:
            total+=2*t-(2*t**3)//(3*AG**2)+core.ceildiv(2*t**5,5*AG**4)
        else:total+=2*t
    return total

core.area_bound=cuts.area_bound=area_bound
cuts.POINTS=[(x,z) for x in (-U,F(0),U) for z in (-V,F(0),V)]
cuts.CORNERS=[(x,z) for x in (-U,U) for z in (-V,V)]
cuts.QUADS=[]
for sx in (-1,1):
    for sz in (-1,1):
        x0,x1=sorted((F(0),sx*U));z0,z1=sorted((F(0),sz*V))
        cuts.QUADS.append((sx,sz,[(x0,z0),(x1,z0),(x1,z1),(x0,z1)]))

def radians(path):
    box=core.normalized_box(path)
    return tuple(LOWER+x*(P/2-LOWER) if j<4 else x*P for j,x in enumerate(box))

def scalar_checks():
    a,b=core.atan_bounds(F(1,5),80);c,d=core.atan_bounds(F(1,239),20)
    assert core.PI_LO<16*a-4*d and 16*b-4*c<core.PI_HI<P
    # The theta/phi lower endpoint covers arcsin(K-1).
    assert F(signed_trig(LOWER)[0][1],D)<K-1
    # Each sinusoidal coefficient has at most one critical point
    # on this interval. The new normal enclosure includes both kinds
    # of critical value, so negative phi needs no monotonicity assumption.
    assert P/2-LOWER<core.PI_LO
    S=1+U*U;T=S+V*V
    E=R*R+H*H-F(4,3)*RADIUS**2
    C=R*R*S/T+H*H-F(4,3)*RADIUS**2
    assert E>0 and C>0
    assert E*E*S<4*H*H*R*R
    assert C*C*T<4*H*H*R*R
    # sqrt(4/3)<7/6, so every tangency point is private.
    assert abs(R-H)+F(7,6)*RADIUS<2*RADIUS
    ratio_squared=U*U*V*V/T
    low=F(isqrt(ratio_squared.numerator*D*D//ratio_squared.denominator),D)
    if low<F(1,2):
        raw_area_low=4*core.atan_bounds(low,60)[0]
    else:
        # atan(x)=pi/4+atan((x-1)/(x+1)).
        z=(low-1)/(low+1)
        if z>0:angle_lower=core.atan_bounds(z,60)[0]
        elif z<0:angle_lower=-core.atan_bounds(-z,60)[1]
        else:angle_lower=F(0)
        raw_area_low=core.PI_LO+4*angle_lower
    scaled=raw_area_low*10**15
    area_low=F(scaled.numerator//scaled.denominator,10**15)
    margin=N*area_low-(3*N-6)*DELTA-4*core.PI_HI
    assert margin>0
    return {'area_lower':str(area_low),'area_lower_decimal':float(area_low),
            'margin_lower':str(margin),'margin_lower_decimal':float(margin)}

def read_partition(name, kind):
    with gzip.open(HERE/name,'rt') as f:rows=json.load(f)
    allowed = {'separation', 'cutoff'} if kind=='connectivity' else {'separation', 'order', 'overlap'}
    assert isinstance(rows,list) and rows, 'Empty or invalid partition'
    assert all(isinstance(row,list) and len(row)==2 and isinstance(row[0],str)
               and set(row[0])<=set('01') and row[1] in allowed for row in rows), 'Invalid partition leaf'
    rows.sort(key=lambda z:z[0]);paths=[p for p,_ in rows]
    assert all(not b.startswith(a) for a,b in zip(paths[:-1],paths[1:]))
    depth=max(map(len,paths))
    assert sum(1<<(depth-len(p)) for p in paths)==1<<depth
    return rows

def check_chunk(rows):
    types={};maximum=0;failures=[]
    for path,reason in rows:
        box=radians(path);types[reason]=types.get(reason,0)+1
        if reason=='order':assert box[0]>box[3]
        elif reason=='separation':
            upper=trig_range(box[0],box[1])[0][1]+trig_range(box[2],box[3])[0][1]
            assert upper*K.denominator<K.numerator*D
        elif reason=='overlap':
            value=core.bound(box);maximum=max(maximum,value)
            if F(value,AG)>DELTA:failures.append([path,str(F(value,AG))])
        elif reason=='cutoff':
            value,reached=cuts.bound(box)
            assert reached==4 and value==0,(path,value,reached)
        else:raise AssertionError(reason)
    core.normalized_box.cache_clear()
    return {'count':len(rows),'types':types,'max':maximum,'failures':failures}

def verify(shard=None):
    start=time.time();scalars=scalar_checks()
    print('Exact rational backend:', BACKEND, VERSION, flush=True)
    print('Scalar checks passed; area lower',scalars['area_lower_decimal'],'margin',scalars['margin_lower_decimal'],flush=True)
    cutoff=read_partition(CASE['cutoff_partition'], 'connectivity')
    cut_result=check_chunk(cutoff)
    print('Connected complements verified:',len(cutoff),'boxes.',flush=True)
    rows=read_partition(CASE['overlap_partition'], 'overlap')
    count=0;maximum=0;types={};failures=[]
    chunks=[rows[i:i+10000] for i in range(0,len(rows),10000)]
    if shard is not None:
        assert 0<=shard<4
        chunks=chunks[shard::4]
    expected=sum(map(len,chunks))
    for ch in chunks:
        got=check_chunk(ch);count+=got['count'];maximum=max(maximum,got['max']);failures+=got['failures']
        for k,v in got['types'].items():types[k]=types.get(k,0)+v
        if count%100000<10000 or count==expected:
            print('Verified',count,'of',expected,'max',str(F(maximum,AG)),'failures',len(failures),'seconds',round(time.time()-start,1),flush=True)
    assert not failures,failures[:10]
    result={'arithmetic_backend':BACKEND,'arithmetic_version':VERSION,'case':CASE,'scalar_bounds':scalars,'cutoff_boxes':len(cutoff),'cutoff_types':cut_result['types'],'cutoff_bound':'0',
            'overlap_boxes':len(rows),'checked_boxes':count,'shard':shard,'types':types,'max_overlap_bound':str(F(maximum,AG)),
            'failures':failures,'seconds':time.time()-start}
    name=f'n{N}_verification.json' if shard is None else f'n{N}_verification_part_{shard}.json'
    with open(Path(__file__).resolve().parent.parent/'reports'/name,'w') as f:json.dump(result,f,indent=2)
    print('Full radius certificate checks passed.' if shard is None else 'SHARD ONLY passed; this is not a complete verification.',flush=True)

if __name__=='__main__':verify(int(sys.argv[3]) if len(sys.argv)>3 and sys.argv[2]=='--shard' else None)
