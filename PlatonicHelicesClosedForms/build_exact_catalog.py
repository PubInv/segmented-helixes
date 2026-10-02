#!/usr/bin/env python3
"""Regenerate exact data and LaTeX tables. Requires SymPy 1.14 or later.

The manuscript is self-contained: running this file is optional. The independent
verifier uses only Python's standard library. Every displayed radical is checked
symbolically against the theorem before it is written.
"""
import json
from pathlib import Path
import sympy as sp
from catalog_data import ROWS, SOURCE_REVISION

ROOT = Path(__file__).resolve().parent
w = sp.sqrt(5)
g = (1+w)/2
p = sp.sqrt(6)+sp.sqrt(2)
q = sp.sqrt(6)-sp.sqrt(2)
K = sp.sqrt(5+2*w)/(2*w)
CONSTANTS = {'w': w, 'g': g, 'p': p, 'q': q, 'K': K}
SYMBOLS = {name: sp.Symbol(name, positive=True) for name in CONSTANTS}

RHO = {'Tet': sp.sqrt(6)/12, 'Cub': sp.Rational(1,2),
       'Oct': sp.sqrt(6)/6, 'Dod': sp.sqrt(250+110*w)/20,
       'Ico': sp.sqrt(3)*(3+w)/12}
ALPHA = {'Tet': {1: sp.Rational(1,3)},
         'Cub': {1: sp.Integer(1), 2: sp.Integer(0)},
         'Oct': {1: -sp.Rational(1,3), 2: sp.Rational(1,3), 5: sp.Integer(1)},
         'Dod': {1: -w/5, 2: w/5, 8: sp.Integer(1)},
         'Ico': {1: -w/3, 2: -sp.Rational(1,3), 8: sp.Rational(1,3),
                 12: w/3, 13: sp.Integer(1)}}

# r, d, cos(phi). Constants are defined in the manuscript, not free inputs.
DISPLAY = {
 'Tetrahelix': ('sqrt(3)/30', '-1/sqrt(10)', '-3/sqrt(10)'),
 'Tetratorus': ('sqrt(3)/6', '0', '0'),
 'Boxbeam': ('0', '1', '1'),
 'Staircase': ('0', '1/sqrt(2)', '1'),
 'Blockhelix': ('sqrt(2)/6', '-1/sqrt(3)', '-sqrt(2/3)'),
 'Cubatorus': ('1/2', '0', '0'),
 'Octabeam': ('0', 'sqrt(6)/3', '1'),
 'Octaspiky': ('2*sqrt(3)/33', '-sqrt(22)/11', '-3/sqrt(11)'),
 'Octamedium': ('sqrt(3)/15', '-sqrt(10)/5', '-3/sqrt(10)'),
 'Octagear': ('sqrt(3)/6', '0', '0'),
 'Treestar': ('sqrt(3)/3', '0', '0'),
 'Dodecabeam': ('0', 'sqrt(250+110*w)/10', '1'),
 'Dodecadoubler': ('2*w*K/(15+2*w)', '-g**2/w*sqrt(10/(13-w))', '-sqrt(10/(13-w))'),
 'The Alternater': ('4*w*K/(35+w)', '-g**3/w*sqrt(10/(17-3*w))', '-sqrt(10/(17-3*w))'),
 'Dodecashaft': ('4*w*g**2*K/(35-w)', '-g**2/w*sqrt(10/(17+3*w))', '-sqrt(10/(17+3*w))'),
 'Dodecagear': ('K', '0', '0'),
 'Dodecacorkscrew': ('2*w*g**2*K/(15-2*w)', '-g**3/w*sqrt(10/(13+w))', '-sqrt(10/(13+w))'),
 'Dodecadonut': ('g**2*K', '0', '0'),
 'Pearlshaft': ('0', 'g**2/sqrt(3)', '1'),
 'Quasi-planar': ('g**2*sqrt(4-p)/(sqrt(3)*(8+p))', 'g**2*sqrt((4+p)/(3*(8+p)))', 'sqrt(3*(4+p)/(2*(8+p)))'),
 'Two Strands': ('2*g**2/(sqrt(3)*(21+w))', 'g*sqrt(2/(21+w))', 'sqrt(18/(21+w))'),
 'Slow Twist': ('2*g**2/(sqrt(3)*(21-w))', 'g**3*sqrt(2/(21-w))', 'sqrt(18/(21-w))'),
 'Rock Candy': ('g**2*sqrt(6)/33', 'g**2/sqrt(11)', '3/sqrt(11)'),
 'Icosa Tree Star': ('sqrt(3)/6', '0', '0'),
 'Icosacorkscrew': ('g**2*sqrt(4+q)/(sqrt(3)*(8-q))', '-g**2*sqrt((4-q)/(3*(8-q)))', '-sqrt(3*(4-q)/(2*(8-q)))'),
 'Planar point cluster': ('g**2/(2*sqrt(6))', '0', '0'),
 'Big Icosacorkscrew': ('g**2*sqrt(26+17*sqrt(2))/(14*sqrt(3))', 'g**2*sqrt((3-sqrt(2))/21)', 'sqrt((9-3*sqrt(2))/14)'),
 'The Wheel': ('g**4/(2*sqrt(3))', '0', '0'),
}
CORRECTED = {
 'Quasi-planar': ('beta+120', -(3*w+1)/8,
   '2/(sqrt(3)*(17+3*w))', 'g**3*sqrt(2/(17+3*w))', '3*g/sqrt(17+3*w)'),
 'Icosacorkscrew': ('beta-120', sp.Rational(1,4),
   'g**2*sqrt(15)/21', '-g**2/sqrt(7)', '-3/sqrt(14)'),
 'Big Icosacorkscrew': ('beta', (3*w-1)/8,
   '2*g**4/(sqrt(3)*(17-3*w))', 'g*sqrt(2/(17-3*w))', '3/(g*sqrt(17-3*w))'),
}

def exact(text):
    return sp.sympify(text, locals=CONSTANTS)

def latex_display(text):
    # Preserve a compact ratio under a radical rather than expanding it into
    # reciprocal roots. The evaluated expression is checked independently.
    result=sp.latex(sp.sympify(text, locals=SYMBOLS, evaluate=False))
    return (result.replace(r'1 \cdot ', '')
                  .replace(r'\left(-1\right) ', '- ')
                  .replace(r'1 \frac{1}', r'\frac{1}'))

def length_display(name, field, text):
    # Keep the dodecahedral product outside the radical rather than putting
    # a fraction under a radical inside the numerator of another fraction.
    powers={'Dodecadoubler':(2,'13-w'), 'The Alternater':(3,'17-3*w'),
            'Dodecashaft':(2,'17+3*w'), 'Dodecacorkscrew':(3,'13+w')}
    if field=='d' and name in powers:
        power,den=powers[name]
        return rf'-\frac{{g^{power}}}{{w}}\sqrt{{\frac{{10}}{{{latex_display(den)}}}}}'
    return latex_display(text)

def zero(expr):
    # After squaring, every expression belongs to Q(sqrt(2),sqrt(3),sqrt(5)).
    # Clearing denominators and expanding gives an exact coefficient check,
    # without expensive general-purpose radical denesting.
    return sp.expand(sp.together(expr).as_numer_denom()[0]) == 0

def parameter_record(row, replacement=None):
    name, solid, analogs, face, tau, old_r, old_theta, old_d, old_phi = row
    a = ALPHA[solid][face]
    x = sp.simplify(sp.cos(sp.pi*sp.Rational(tau,180)))
    label = str(tau)
    texts = DISPLAY[name]
    signed = -1 if tau < 0 else 1
    if replacement:
        label, x, *texts = replacement
        signed = -1 if 'minus' in label or label == 'beta-120' else 1
    L2 = sp.expand(2*RHO[solid]**2*(1+a))
    D = sp.expand(3-a-(1+a)*x)
    ctheta = sp.expand((1+a)*(1+x)/2-1)
    theta_text = str(sp.factor(ctheta)).replace('sqrt(5)', 'w')
    if not replacement and name == 'Quasi-planar':
        theta_text = '-(2+p)/6'
    elif not replacement and name == 'Icosacorkscrew':
        theta_text = '(q-2)/6'
    if not zero(exact(theta_text)-ctheta):
        raise AssertionError((name, 'cos(theta) display'))
    r2 = L2*(1-a)*(1+x)/D**2
    straight = bool(zero(r2))
    cp2 = sp.Integer(1) if straight else 2*(1-x)/D
    d2 = L2 if straight else L2*cp2
    r, d, cp = map(exact, texts)
    for field, expr, squared in [('r',r,r2),('d',d,d2),('cos_phi',cp,cp2)]:
        if not zero(expr**2-squared):
            raise AssertionError((name, field, sp.simplify(expr**2-squared)))
    if float(r) < 0 or (not straight and float(cp)*signed < 0):
        raise AssertionError((name, 'sign'))
    # The squared identities above and equal signs establish d=L*cos(phi).
    # Avoid asking a heuristic simplifier to denest unrelated square roots.
    if abs(float(d)-float(sp.sqrt(L2)*cp)) > 1e-12:
        raise AssertionError((name,'d=L*cos(phi)'))
    values = {'cos_theta': float(ctheta), 'cos_phi': float(cp), 'r': float(r), 'd':float(d)}
    return {'name': name, 'solid':solid, 'analogs':analogs, 'face':face,
       'twist_label':label, 'twist_degrees': tau if not replacement else None,
       'cos_alpha': str(a), 'cos_tau': str(x), 'L_squared': str(L2),
       'straight':straight,
       'expressions':dict(zip(('r','d','cos_phi'),texts), cos_theta=theta_text),
       'latex':dict(zip(('r','d','cos_phi'),
                    (length_display(name,field,text) for field,text in zip(('r','d','cos_phi'),texts))),
                    cos_theta=latex_display(theta_text)),
       'values':values,
       'source_printed':{'r':old_r,'theta':old_theta,'d':old_d,'phi':old_phi}}

def tables(records):
    out=[]
    cell=lambda text: r'$\displaystyle '+text+'$'
    for solid in ('Tet','Cub','Oct','Dod','Ico'):
        rs=[r for r in records if r['solid']==solid]
        if solid in ('Tet','Cub','Oct'):
            out += [r'\begin{center}',r'\small',r'\renewcommand{\arraystretch}{1.7}',
                r'\begin{tabular}{@{}lrrcccc@{}}',r'\toprule',
                r'Name & $C$-face & $\tau$ (deg.) & $\cos\theta$ & $\cos\phi$ & $r$ & $d$\\',r'\midrule']
            for r in rs:
                cp=r['latex']['cos_phi']+(r'^{\dagger}' if r['straight'] else '')
                fields=(r['latex']['cos_theta'],cp,r['latex']['r'],r['latex']['d'])
                out.append(f"{r['name']} & {r['face']} & {r['twist_label']} & "+' & '.join(map(cell,fields))+r'\\')
            out += [r'\bottomrule',r'\end{tabular}',r'\end{center}']
            (ROOT/f'tables_{solid}.tex').write_text('\n'.join(out)+'\n')
            out=[]
            continue
        out += [r'\begin{center}',r'\small',r'\renewcommand{\arraystretch}{1.55}',
                r'\begin{tabular}{@{}lrrcc@{}}',r'\toprule',
                r'Name & $C$-face & $\tau$ (deg.) & $\cos\theta$ & $\cos\phi$\\',r'\midrule']
        for r in rs:
            cp=r['latex']['cos_phi']; ct=r['latex']['cos_theta']
            if r['straight']:cp+=r'^{\dagger}'
            out.append(f"{r['name']} & {r['face']} & {r['twist_label']} & "+cell(ct)+' & '+cell(cp)+r'\\')
        out += [r'\bottomrule',r'\end{tabular}',r'\end{center}',
                r'\begin{center}',r'\small',r'\renewcommand{\arraystretch}{1.6}',
                r'\begin{tabular}{@{}lcc@{}}',r'\toprule',r'Name & $r$ & $d$\\',r'\midrule']
        for r in rs:
            out.append(r['name']+' & '+cell(r['latex']['r'])+' & '+cell(r['latex']['d'])+r'\\')
        out += [r'\bottomrule',r'\end{tabular}',r'\end{center}']
        (ROOT/f'tables_{solid}.tex').write_text('\n'.join(out)+'\n')
        out=[]

def check_icosahedral_matching_phase():
    """Check the face-0/face-8 vertex alignment in exact coordinates."""
    nb=sp.Matrix([-1,1,1])/sp.sqrt(3)
    nc=sp.Matrix([-1/g,0,-g])/sp.sqrt(3)
    B=g**2*sp.Matrix([-1,1,1])/6
    C=sp.Matrix([-g,0,-g**3])/6
    PB=sp.Matrix([-g/2,0,sp.Rational(1,2)])
    PC=sp.Matrix([-g/2,0,-sp.Rational(1,2)])
    v=nb.cross(-nc)
    V=sp.Matrix([[0,-v[2],v[1]],[v[2],0,-v[0]],[-v[1],v[0],0]])
    A=sp.eye(3)+V+sp.Rational(3,4)*V**2
    a=A*(PB-B);b=PC-C
    claims=(nb.dot(nc)+sp.Rational(1,3),
            a.dot(a)-sp.Rational(1,3),b.dot(b)-sp.Rational(1,3),
            3*a.dot(b)-(3*w-1)/8,
            3*nc.dot(a.cross(b))-sp.sqrt(3)*(1+w)/8)
    assert all(zero(expr) for expr in claims)
    print('Verified exact vertex derivation of cos(beta) and sin(beta).')

def main():
    check_icosahedral_matching_phase()
    original=[]
    for row in ROWS:
        original.append(parameter_record(row))
        print('Verified exact expressions:',row[0],flush=True)
    corrected=[parameter_record(row,CORRECTED[row[0]]) for row in ROWS if row[0] in CORRECTED]
    obj={'source_revision':SOURCE_REVISION,'normalization':'unit edge length',
         'constants':{n:str(v) for n,v in CONSTANTS.items()},
         'angle_convention':'principal rigid-motion theta; forward axis in straight rows',
         'original_catalog':original,'corrected_matching_rows':corrected}
    (ROOT/'ExactCatalog.json').write_text(json.dumps(obj,indent=2)+'\n')
    tables(original)
    print('Symbolic checks passed for 28 original catalog rows and 3 matching corrections.')
    for r in original+corrected:
        print(r['name'],r['twist_label'],r['expressions'])

if __name__ == '__main__':
    main()
