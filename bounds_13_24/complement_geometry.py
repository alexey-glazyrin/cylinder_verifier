r"""Exact bound on components of Q1\Q2 not containing Q1's center.

Each complement is the union of four affine caps in a gnomonic rectangle.
Caps proved connected to the center whenever nonempty are marked.
The release verifier accepts ONLY when all four caps are marked.
The historical area fallback is not used to accept any case.
"""
from exact_geometry import (F,D,AG,U,V,ROOT_PI,PI_HI,ceildiv,
    norm_intervals,trig_range,clip,area_bound,normalized_box)
import json,time

POINTS=[(x,z) for x in (-U,F(0),U) for z in (-V,F(0),V)]
CORNERS=[(x,z) for x in (-U,U) for z in (-V,V)]
QUADS=[]
for sx in (-1,1):
    for sz in (-1,1):
        x0,x1=sorted((F(0),sx*U));z0,z1=sorted((F(0),sz*V))
        QUADS.append((sx,sz,[(x0,z0),(x1,z0),(x1,z1),(x0,z1)]))

def upper(n,p):
    x,z=p
    return n[0][1 if x>=0 else 0]*x+n[2][1 if z>=0 else 0]*z-n[1][0]

def lower(n,p):
    x,z=p
    return n[0][0 if x>=0 else 1]*x+n[2][0 if z>=0 else 1]*z-n[1][1]

def common_witness(a,b):
    for sx,sz,quad in QUADS:
        p=quad
        for n in (a,b):
            p=clip(p,-n[0][1 if sx>0 else 0],-n[2][1 if sz>0 else 0],-n[1][0])
        if len(p)>=3:
            q=(sum(v[0] for v in p)/len(p),sum(v[1] for v in p)/len(p))
            if upper(a,q)<0 and upper(b,q)<0:return True
    return False

def cap_bound(n):
    total=0
    for sx,sz,quad in QUADS:
        p=clip(quad,-n[0][0 if sx>0 else 1],-n[2][0 if sz>0 else 1],-n[1][1])
        total+=area_bound(p)
    return total

def bound(box):
    ns=norm_intervals(box)
    reached={i for i,n in enumerate(ns) if n[1][0]>0}
    assert reached, 'No uniformly private center certified'
    adj=[set() for _ in range(4)]
    for p in POINTS:
        neg=[i for i,n in enumerate(ns) if upper(n,p)<0]
        for i in neg:adj[i].update(neg)
    while True:
        new=reached|set().union(*(adj[i] for i in reached))
        for i in set(range(4))-reached:
            # Every nonempty affine cap has a strictly negative vertex.
            # Each potentially negative vertex must already be in a
            # center-connected cap. This also handles an empty cap.
            possible=[p for p in CORNERS if lower(ns[i],p)<=0]
            if all(any(upper(ns[j],p)<0 for j in reached) for p in possible):
                new.add(i)
        if new==reached:break
        reached=new
    if len(reached)<4:
        for _ in range(3):
            new={i for i in set(range(4))-reached
                 if any(common_witness(ns[i],ns[j]) for j in reached)}
            if not new:break
            reached.update(new)
    return sum(cap_bound(ns[i]) for i in set(range(4))-reached),len(reached)
