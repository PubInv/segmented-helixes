#!/usr/bin/env python3
"""Independent numerical verification; Python 3, no third-party packages.

Run beside catalog_data.py and ExactCatalog.json. Uses Rodrigues matrices,
projected three-point circumcircles, explicit face vertices, and a transcript
of the 2022 catalog. It does not use quaternions to obtain its reference axes.
"""
import ast
from collections import Counter
import json
import math
from pathlib import Path
from catalog_data import ROWS, FACES, SOURCE_REVISION, vertices, source_edge, inradius

ROOT = Path(__file__).resolve().parent
W = math.sqrt(5)
G = (1+W)/2
BETA = math.acos((3*W-1)/8)

def dot(a,b): return sum(x*y for x,y in zip(a,b))
def cross(a,b):
    return (a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0])
def add(a,b): return tuple(x+y for x,y in zip(a,b))
def sub(a,b): return tuple(x-y for x,y in zip(a,b))
def scale(a,k): return tuple(x*k for x in a)
def norm(a): return math.sqrt(dot(a,a))
def unit(a): return scale(a,1/norm(a))
def clamp(a): return max(-1.0,min(1.0,a))
def mv(M,a): return tuple(dot(row,a) for row in M)
def mm(A,B): return tuple(tuple(dot(row,col) for col in zip(*B)) for row in A)

def rotation(axis,angle):
    x,y,z=u=unit(axis); c=math.cos(angle);s=math.sin(angle)
    K=((0,-z,y),(z,0,-x),(-y,x,0))
    return tuple(tuple(c*(i==j)+(1-c)*u[i]*u[j]+s*K[i][j]
                 for j in range(3)) for i in range(3))

def shortest(a,b):
    a,b=unit(a),unit(b); v=cross(a,b); sn=norm(v);co=clamp(dot(a,b))
    if sn<1e-14:
        if co>0:return ((1,0,0),(0,1,0),(0,0,1))
        k=min(range(3),key=lambda i:abs(a[i]))
        return rotation(cross(a,tuple(float(i==k) for i in range(3))),math.pi)
    return rotation(v,math.atan2(sn,co))

def baseline(nb,nc,e):
    return mm(shortest(e,nc),shortest(nb,scale(e,-1)))

def face_geometry(solid):
    vs=vertices(solid);fs=[]
    for ids in FACES[solid]:
        ps=[vs[i] for i in ids]
        c=tuple(sum(p[i] for p in ps)/len(ps) for i in range(3))
        fs.append({'center':c,'normal':unit(c),'corners':ps})
    return fs

GEOMETRY={s:face_geometry(s) for s in FACES}

def direct(solid,face,tau):
    """Axis from matrix skew; radius from a projected circumcircle."""
    f0,fc=GEOMETRY[solid][0],GEOMETRY[solid][face]
    b,c=f0['center'],fc['center'];nb,nc=f0['normal'],fc['normal']
    delta=sub(c,b);L=norm(delta);e=unit(delta)
    R=mm(rotation(nc,tau),baseline(nb,nc,e));t=sub(c,mv(R,b))
    T=lambda x:add(mv(R,x),t)
    D=T(c);step2=sub(D,c)
    ct=clamp((sum(R[i][i] for i in range(3))-1)/2)
    theta=math.acos(ct)
    # Straight joint points do not determine a directed rotational axis.
    if norm(cross(delta,step2))<1e-12:
        assert norm(sub(step2,delta))<1e-11
        r=0.0;d=L;cp=1.0
    else:
        h=unit((R[2][1]-R[1][2],R[0][2]-R[2][0],R[1][0]-R[0][1]))
        d=dot(delta,h);cp=d/L
        P=sub(delta,scale(h,d));Q=sub(step2,scale(h,d))
        # Circumradius = abc/(4 area); cross(P,Q) has magnitude 2 area.
        r=norm(P)*norm(Q)*norm(add(P,Q))/(2*norm(cross(P,Q)))
    return {'r':r,'d':d,'cos_phi':cp,'cos_theta':ct,'L':L,
            'theta':math.degrees(theta),'phi':math.degrees(math.acos(clamp(cp))),
            'R':R,'T':T,'cos_alpha':-dot(nb,nc)}

def formula(solid,face,tau):
    fs=GEOMETRY[solid];b,c=fs[0]['center'],fs[face]['center']
    L=norm(sub(c,b));a=-dot(fs[0]['normal'],fs[face]['normal'])
    Q=1-(1+a)/2*math.cos(tau/2)**2
    ct=1-2*Q
    if abs(a-1)<1e-12 or abs(abs(tau)-math.pi)<1e-12:
        return dict(r=0,d=L,cos_phi=1,cos_theta=ct)
    cp=math.sin(tau/2)/math.sqrt(Q)
    return dict(r=L*math.sqrt((1-a)/2)*abs(math.cos(tau/2))/(2*Q),
                d=L*cp,cos_phi=cp,cos_theta=ct)

def expression(text):
    """Evaluate only the arithmetic grammar used in the exact data file."""
    tree=ast.parse(text,mode='eval')
    allowed=(ast.Expression,ast.BinOp,ast.UnaryOp,ast.Call,ast.Name,ast.Load,
             ast.Constant,ast.Add,ast.Sub,ast.Mult,ast.Div,ast.Pow,ast.USub,ast.UAdd)
    for node in ast.walk(tree):
        if not isinstance(node,allowed):raise ValueError('Unexpected expression')
        if isinstance(node,ast.Call) and not (isinstance(node.func,ast.Name)
             and node.func.id=='sqrt' and len(node.args)==1 and not node.keywords):
            raise ValueError('Unexpected function')
    vals={'sqrt':math.sqrt,'w':W,'g':G,'p':math.sqrt(6)+math.sqrt(2),
          'q':math.sqrt(6)-math.sqrt(2),'K':math.sqrt(5+2*W)/(2*W)}
    return eval(compile(tree,'<exact radical>','eval'),{'__builtins__':{}},vals)

def source_base(solid,face):
    if solid=='Oct' and face==5:return 60.0
    if solid=='Dod' and face==8:return 36.0
    if solid=='Ico':
        if face==13:return 60.0
        return 45.0 if face in (8,9,10,11,17,19) else 120.0
    return 0.0

def corrected_base(solid,face):
    if solid=='Ico' and face in (8,11,19):return math.degrees(BETA)
    if solid=='Ico' and face in (9,10,17):return -math.degrees(BETA)
    return source_base(solid,face)

def matching_error(solid,face,tau):
    z=direct(solid,face,tau)
    return max(min(norm(sub(z['T'](p),q)) for q in GEOMETRY[solid][face]['corners'])
               for p in GEOMETRY[solid][0]['corners'])

def check_geometry():
    for solid,fs in GEOMETRY.items():
        vs=vertices(solid)
        edges=sorted(norm(sub(a,b)) for i,a in enumerate(vs) for b in vs[i+1:])
        assert abs(edges[0]-1)<1e-12
        for f in fs:
            assert abs(norm(f['center'])-inradius(solid))<1e-12
            assert max(abs(dot(sub(p,f['center']),f['normal'])) for p in f['corners'])<1e-12
        counts=Counter(round(-dot(fs[0]['normal'],f['normal']),9) for f in fs[1:])
        print(f'{solid}: unit edges, equal face inradii, cos(alpha) counts {dict(counts)}')

def check_catalog(data):
    worst={'r':0.0,'d':0.0,'theta':0.0,'phi':0.0};exact_error=0.0
    print('\n28 catalog checks; lengths converted back to the original constructor scale:')
    print(f"{'Name':23s} {'r (edge 1)':>12s} {'d (edge 1)':>12s} {'theta_R':>11s} {'phi':>11s}")
    for row,record in zip(ROWS,data['original_catalog']):
        name,solid,_,face,twist,r0,th0,d0,ph0=row
        z=direct(solid,face,math.radians(twist))
        for field in ('r','d','cos_phi','cos_theta'):
            val=expression(record['expressions'][field])
            err=abs(val-z[field]);exact_error=max(exact_error,err)
            assert err<2e-11,(name,field,err)
        comparisons={'r':(z['r']*source_edge(solid),r0),
            'd':(z['d']*source_edge(solid),d0),
            'theta':(0.0 if record['straight'] else z['theta'],th0),
            'phi':(z['phi'],ph0)}
        for field,(got,want) in comparisons.items():
            worst[field]=max(worst[field],abs(got-want))
            assert abs(got-want)<.0005,(name,field,got,want)
        print(f"{name:23s} {z['r']:12.8f} {z['d']:12.8f} {z['theta']:11.6f} {z['phi']:11.6f}")
    print('PASS: 28/28 match every printed parameter after scale/convention conversion.')
    print('Maximum printed-value discrepancies:',worst)
    print('Maximum exact-radical versus independent matrix/circumcircle discrepancy:',exact_error)

def check_all_rules(data):
    n=0;bad=[];worst_match=0.0;worst_formula=0.0
    # Handedness is paired. Straight joint orbits do not distinguish twists
    # of the full rigid body; this reproduces the catalog's multiplicities.
    representatives={r['name']:r for r in data['original_catalog']}
    representatives.update({r['name']:r for r in data['corrected_matching_rows']})
    multiplicities=Counter()
    for solid,fs in GEOMETRY.items():
        count={'Tet':3,'Cub':4,'Oct':3,'Dod':5,'Ico':3}[solid]
        for face in range(1,len(fs)):
            for j in range(count):
                old=source_base(solid,face)+(j-count//2)*360/count
                new=corrected_base(solid,face)+(j-count//2)*360/count
                old=((old+180)%360)-180;new=((new+180)%360)-180
                err=matching_error(solid,face,math.radians(old))
                if err>1e-10:bad.append((solid,face,old,err))
                err=matching_error(solid,face,math.radians(new))
                worst_match=max(worst_match,err)
                assert err<1e-10,(solid,face,new,err)
                for deg in (old,new):
                    z=direct(solid,face,math.radians(deg));f=formula(solid,face,math.radians(deg))
                    for field in f:worst_formula=max(worst_formula,abs(z[field]-f[field]))
                candidates=[]
                for name,record in representatives.items():
                    if record['solid']!=solid:continue
                    v=record['values']
                    err=max(abs(z['r']-v['r']),abs(abs(z['d'])-abs(v['d'])))
                    if not record['straight']:
                        err=max(err,abs(z['cos_theta']-v['cos_theta']))
                    if err<1e-10:candidates.append(name)
                assert len(candidates)==1,(solid,face,new,candidates)
                multiplicities[candidates[0]]+=1
                n+=1
    assert n==162 and len(bad)==18
    assert {s for s,*_ in bad}=={'Ico'}
    assert {f for _,f,*_ in bad}=={8,9,10,11,17,19}
    assert worst_formula<2e-11
    print('\nPASS: all 162 corrected distinct-face rules match their face vertices.')
    print('Original rules failing exact matching:',len(bad),'(the six class-T icosahedral faces).')
    print('Maximum corrected vertex error:',worst_match)
    print('Maximum matrix/circumcircle versus generic-formula discrepancy:',worst_formula)
    for name,record in representatives.items():
        assert multiplicities[name]==record['analogs'],(name,multiplicities[name])
    print('PASS: all 28 representative multiplicities, pairing handedness and grouping straight orbits.')

def check_continuous():
    n=0;worst=0.0
    for solid,fs in GEOMETRY.items():
        for bi,b in enumerate(fs):
            for ci,c in enumerate(fs):
                a=-dot(b['normal'],c['normal'])
                if bi==ci or abs(abs(a)-1)<1e-12:continue
                # Re-index the same solid temporarily so arbitrary incoming faces
                # test the actual matrices; the formulas do not assume face 0.
                saved=GEOMETRY[solid]
                GEOMETRY[solid]=[b,c]
                try:
                    for deg in (-165,-120,-72,0,45,120,165):
                        z=direct(solid,1,math.radians(deg));f=formula(solid,1,math.radians(deg))
                        for field in f:worst=max(worst,abs(z[field]-f[field]))
                        n+=1
                finally:GEOMETRY[solid]=saved
    assert n==3948 and worst<2e-11
    print('\nPASS:',n,'ordered-face continuous-twist cases; maximum discrepancy',worst)

def check_corrections(data):
    print('\nExact icosahedral matching correction: beta =',math.degrees(BETA),'degrees')
    for rec in data['corrected_matching_rows']:
        tau=BETA + (2*math.pi/3 if rec['twist_label']=='beta+120'
                    else -2*math.pi/3 if rec['twist_label']=='beta-120' else 0)
        z=direct('Ico',8,tau)
        for field in ('r','d','cos_phi','cos_theta'):
            assert abs(expression(rec['expressions'][field])-z[field])<2e-11
        assert matching_error('Ico',8,tau)<1e-10
        print(f"{rec['name']:23s} tau={math.degrees(tau):12.8f} r={z['r']:.9f} d={z['d']:.9f} theta={z['theta']:.8f} phi={z['phi']:.8f}")
    # beta's independently derived sine/cosine are a unit-circle point.
    co=(3*W-1)/8;si=math.sqrt(3)*(1+W)/8
    assert abs(co*co+si*si-1)<1e-14
    print('PASS: corrected radicals, matching vertices, and exact beta sine/cosine.')

def main():
    data=json.loads((ROOT/'ExactCatalog.json').read_text())
    assert data['source_revision']==SOURCE_REVISION
    print('Closed-form Platonic helix verification; standard-library Python.')
    check_geometry();check_catalog(data);check_all_rules(data);check_continuous();check_corrections(data)
    print('\nAll checks passed. Numerical checks support, and do not replace, the manuscript proofs.')

if __name__=='__main__':main()
