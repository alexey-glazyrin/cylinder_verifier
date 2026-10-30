"""Check an overlap subdivision with integer and rational arithmetic.

No floating-point arithmetic is used in any geometric bound. Sine and
cosine are enclosed by Taylor polynomials with outward rounding, and a
spherical triangle is bounded using atan(t) <= t for t >= 0.
"""
from rational import F, BACKEND, VERSION
from functools import lru_cache
from math import isqrt,factorial
import json,time

D=10**40
G=10**15
AG=10**20
U=F(59,100);V=F(48,25)
UI=59*D//100;VI=48*D//25
ROOT_PI=F(3141592653589794,10**15)
PI_LO=F('3.14159265358979323846264338327950288419716939937510')
PI_HI=PI_LO+F(1,10**50)
assert PI_HI<ROOT_PI
assert 4**81*D<factorial(81) and 4**80*D<factorial(80)

def atan_bounds(x,terms):
    assert 0<x<1 and terms%2==0
    lower=sum((-1)**k*x**(2*k+1)/F(2*k+1) for k in range(terms))
    return lower,lower+x**(2*terms+1)/F(2*terms+1)

def verify_scalar_claims():
    # Machin's identity: pi = 16*atan(1/5)-4*atan(1/239).
    a,b=atan_bounds(F(1,5),80);c,d=atan_bounds(F(1,239),20)
    assert PI_LO<16*a-4*d and 16*b-4*c<PI_HI
    R=F(227,100);S=1+U*U;T=S+V*V
    # Exact rational forms of the two no-triple inequalities after squaring.
    assert (R*R+F(8,3))**2*S<16*R*R
    assert (R*R*S/T+F(8,3))**2*T<16*R*R
    qlower=F(25243,50000)
    assert qlower*qlower*T<U*U*V*V
    alower,_=atan_bounds(qlower,20)
    assert 4*alower>F(187,100)
    assert 7*F(187,100)-15*F(346,10000)>4*PI_HI
    print('Exact scalar checks passed: no triples, area > 1.87, and conditional Euler margin.',flush=True)

def ceildiv(a,b):return -((-a)//b)
def ai(a,b):return a[0]+b[0],a[1]+b[1]
def ni(a):return -a[1],-a[0]
def si(a,b):return ai(a,ni(b))
def mi(a,b):
    z=[a[0]*b[0],a[0]*b[1],a[1]*b[0],a[1]*b[1]]
    return min(z)//D,ceildiv(max(z),D)
def di(a,n):return a[0]//n,ceildiv(a[1],n)
def scale(a,x):return mi(a,(x,x))
def fi(x):return (x.numerator*D//x.denominator,ceildiv(x.numerator*D,x.denominator))

@lru_cache(None)
def trig_point(x):
    assert 0<=x<4
    xx=fi(x);x2=mi(xx,xx)
    st=xx;ct=(D,D);sn=st;cs=ct
    for k in range(1,40):
        st=di(mi(st,x2),2*k*(2*k+1))
        ct=di(mi(ct,x2),(2*k-1)*(2*k))
        sn=si(sn,st) if k%2 else ai(sn,st)
        cs=si(cs,ct) if k%2 else ai(cs,ct)
    return (sn[0]-1,sn[1]+1),(cs[0]-1,cs[1]+1)

def trig_range(l,r):
    sl,cl=trig_point(l);sr,cr=trig_point(r)
    s=(min(sl[0],sr[0]),max(sl[1],sr[1]))
    c=(min(cl[0],cr[0]),max(cl[1],cr[1]))
    if l<=PI_HI/2 and r>=PI_LO/2:s=(s[0],max(s[1],D))
    if l==0:c=(c[0],max(c[1],D))
    if l<=PI_HI and r>=PI_LO:c=(min(c[0],-D),c[1])
    return s,c

root=isqrt(D*D+UI*UI)
FAC=(root,root+1)
def norm_intervals(box):
    th0,th1,ph0,ph1,w0,w1=box
    s,c=trig_range(th0,th1);t,d=trig_range(ph0,ph1);sw,cw=trig_range(w0,w1)
    t0,d0=trig_point(ph0);t1,d1=trig_point(ph1)
    def endpoints(a,b):return min(a[0],b[0]),max(a[1],b[1])
    ap=endpoints(ai(d0,scale(t0,UI)),ai(d1,scale(t1,UI)))
    dp0=si(scale(d0,UI),t0);dp1=si(scale(d1,UI),t1)
    if dp0[1]>=0 and dp1[0]<=0:ap=(ap[0],max(ap[1],FAC[1]))
    bp=endpoints(si(t0,scale(d0,UI)),si(t1,scale(d1,UI)))
    am=endpoints(si(scale(t0,UI),d0),si(scale(t1,UI),d1))
    bm=endpoints(ni(ai(t0,scale(d0,UI))),ni(ai(t1,scale(d1,UI))))
    dm0=si(scale(t0,UI),d0);dm1=si(scale(t1,UI),d1)
    if dm0[0]<=0 and dm1[1]>=0:bm=(min(bm[0],-FAC[1]),bm[1])
    ans=[]
    for aa,bb in [(ap,bp),(am,bm)]:
        ans.append([ai(mi(aa,c),mi(mi(bb,s),cw)),si(mi(aa,s),mi(mi(bb,c),cw)),ni(mi(bb,sw))])
    for sign in (1,-1):
        mix=ai(scale(mi(d,cw),VI),sw if sign==1 else ni(sw))
        ans.append([si(scale(mi(t,c),VI),mi(s,mix)),ai(scale(mi(t,s),VI),mi(c,mix)),
                    si(scale(mi(d,sw),VI),cw if sign==1 else ni(cw))])
    return ans

def clip(poly,a,b,c):
    if not poly:return []
    out=[];p=poly[-1];fp=a*p[0]+b*p[1]-c
    for q in poly:
        fq=a*q[0]+b*q[1]-c
        if (fq>=0)!=(fp>=0):
            ratio=fp/(fp-fq)
            out.append((p[0]+ratio*(q[0]-p[0]),p[1]+ratio*(q[1]-p[1])))
        if fq>=0:out.append(q)
        p=q;fp=fq
    return out

def fl(x):return x.numerator*D//x.denominator
def length_lower(p):
    x=1+p[0]*p[0]+p[1]*p[1]
    return isqrt(x.numerator*D*D//x.denominator)
def dot_lower(p,q):return fl(1+p[0]*q[0]+p[1]*q[1])

def area_bound(poly):
    if len(poly)<3:return 0
    center=(sum(p[0] for p in poly)/len(poly),sum(p[1] for p in poly)/len(poly))
    lc=length_lower(center);total=0
    for p,q in zip(poly,poly[1:]+poly[:1]):
        lp=length_lower(p);lq=length_lower(q)
        den=(lc*lp*lq)//(D*D)+(dot_lower(center,p)*lq)//D+(dot_lower(p,q)*lc)//D+(dot_lower(q,center)*lp)//D
        assert den>0
        det=abs((p[0]-center[0])*(q[1]-center[1])-(p[1]-center[1])*(q[0]-center[0]))
        total+=ceildiv(2*det.numerator*D*AG,det.denominator*den)
    return total

def bound(box):
    normals=norm_intervals(box);total=0;unit=D//G
    for sx in (-1,1):
        for sz in (-1,1):
            x0,x1=sorted((F(0),sx*U));z0,z1=sorted((F(0),sz*V))
            p=[(x0,z0),(x1,z0),(x1,z1),(x0,z1)]
            for n in normals:
                a=ceildiv(n[0][1],unit) if sx>0 else n[0][0]//unit
                b=ceildiv(n[2][1],unit) if sz>0 else n[2][0]//unit
                c=n[1][0]//unit
                p=clip(p,a,b,c)
                if not p:break
            total+=area_bound(p)
    return total

@lru_cache(None)
def normalized_box(path):
    if not path:return (F(0),F(1),F(0),F(1),F(0),F(1))
    old=normalized_box(path[:-1]);box=list(old)
    widths=[box[1]-box[0],box[3]-box[2],2*(box[5]-box[4])]
    k=max(range(3),key=lambda j:widths[j]);mid=(box[2*k]+box[2*k+1])/2
    box[2*k+(1 if path[-1]=='0' else 0)]=mid
    return tuple(box)
