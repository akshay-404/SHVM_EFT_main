"""
Exact finite-group check of the printed SHVM charge assignment.
Verifies perturbative congruences and instanton charge sums. It is not a
modern full discrete anomaly classification or a UV charge completion.
"""

from pathlib import Path
from itertools import combinations_with_replacement
from collections import deque
import json
import math


G = (2, 4, 14, 12, 4)
'''
Complete group, G = Z2 x Z4 x Z14 x Z12 x Z4
'''
'''
Q   Left-handed quark doublets
u	Right-handed up-type quarks
d	Right-handed down-type quarks
L	Left-handed lepton doublets
e	Right-handed charged leptons
n	Right-handed neutrinos
chi	Seven scalar flavon fields
'''
base = {
    'Q': [(0, 0, 11), (0, 0, 9), (1, 0, 8)],
    'u': [(1, 1, 9), (0, 0, 4), (0, 0, 0)],
    'd': [(1, 1, 9), (0, 0, 12), (1, 0, 2)],
    'L': [(0, 3, 12), (0, 3, 10), (0, 3, 6)],
    'e': [(1, 0, 10), (0, 3, 13), (0, 3, 1)],
    'n': [(0, 0, 8), (1, 1, 3), (1, 1, 3)],
    'chi': [(1, 3, 2), (0, 0, 5), (1, 0, 8), (0, 0, 13), (0, 0, 11), (0, 0, 6), (1, 2, 8)]
}
'''Original first three charge assignment'''

f12 = {
    'Q': [10, 0, 4],
    'u': [1, 5, 7],
    'd': [1, 9, 5],
    'L': [8, 10, 0],
    'e': [11, 7, 5],
    'n': [11, 3, 3],
    'chi': [9, 7, 9, 1, 3, 11, 6]
}
'''Fourth charge assignment'''

f4 = {
    'Q': [1, 3, 0],
    'u': [0, 2, 0],
    'd': [0, 2, 3],
    'L': [0, 2, 0],
    'e': [3, 1, 3],
    'n': [1, 3, 3],
    'chi': [1, 1, 0, 3, 1, 1, 0]
}
'''Fifth charge assignment'''

q = {
    key: [
        tuple(a)+(b, c)
        for a, b, c in zip(base[key], f12[key], f4[key])
    ]
    for key in base
}
'''Complete five component charge vector'''


def add(a, b, g=G): return tuple((x+y) % n for x, y, n in zip(a, b, g))
def mul(k, a, g=G): return tuple((k*x) % n for x, n in zip(a, g))
def charge(label, g=G): return mul(1 if label > 0 else -1, q['chi'][abs(label)-1], g)


def subgroup(generators, g=G):
    '''
    Generates subgroup from generators
    '''
    zero = (0,) * len(g)
    found = {zero}
    todo = deque([zero])
    while todo:
        u = todo.popleft()
        for a in generators:
            z = add(u, a, g)
            if z not in found:
                found.add(z)
                todo.append(z)
    return found


def charges(g):
    '''
    Charges for the specified group g. 
    Arrays truncate/retarget the fourth modulus for the successive groups.
    '''
    return {
        k: [
            tuple(x % n for x, n in zip(a, g))
            for a in arr
        ]
        for k, arr in q.items()
    }


expected = {
    'u': {(1, 1): 1, (1, 3): 5, (2, 2): 2, (2, 3): -2, (3, 3): 3},
    'd': {(1, 1): 1, (1, 2): 4, (2, 2): 5, (3, 3): 6},
    'e': {(1, 1): 1, (1, 2): 4, (1, 3): 5, (2, 2): 5, (2, 3): -2, (3, 3): 2}
}
'''
Expected leading-order Yukawa operator assignments.
The keys are: (i,j) where i denotes the left-handed fermion generation and 
j denotes the right-handed fermion generation.
The value identifies the scalar field inserted into the operator.
'''

light = [1, 2, 4, 5, 6]
'''Scalar fields considered in the Yukawa and bilinear scans'''

counts = []
for g in [(2, 4, 14, 4), (2, 4, 14, 12), G]:
    qq = charges(g)         # Charges reduced to the current group
    zero = (0,)*len(g)      # Identity element of the current group
    lead = {}               # Leading-order Yukawa operators found
    classes = []            # Higher-order Yukawa operators found
    bil = set()             # Scalar bilinear operator classes

    for sector in ['u', 'd', 'e']:
        '''Yukawa Sector for up-type quarks, down-type quark and charged lepton'''
        ll = 'L' if sector == 'e' else 'Q'
        lead[sector] = {}

        for i in range(3):
            for j in range(3):
                '''Calculate bare Yukawa Charge'''
                bare = add(mul(-1, qq[ll][i], g), qq[sector][j], g)

                '''
                Searching for the leading flavon insertion. All r for which
                -qL + qR + qXr = 0G
                '''
                labels = [
                    r for r in [
                        *range(1, 7), 
                        *range(-1, -7, -1)
                    ]
                    if add(bare, charge(r, g), g) == zero
                ]

                '''Recording the leading operator'''
                if labels:
                    lead[sector][(i+1, j+1)] = labels[0]
                    assert len(labels) == 1

                '''Searching for higher-order Yukawa operators'''
                for r in light:
                    '''Flavon and its conjugate'''
                    for sg in (1, -1):
                        '''Any number k of additional X3 insertions'''
                        for k in range(28):
                            total = add(
                                add(bare, charge(sg*r, g), g), 
                                mul(k, qq['chi'][2], g),
                                g
                            )
                            if total == zero and not (
                                k == 0 and expected[sector].get((i+1, j+1)) == sg*r
                            ):
                                '''
                                Searching for only higher order operators. All +/-r, k for which
                                -qL + qR +/- qXr + kqX3 = 0G
                                '''
                                classes.append((sector, i+1, j+1, sg*r, k))

    assert lead == expected, (g, lead)

    '''
    Search for bilinear operators.
    Scalar classes are identified with their Hermitian conjugates.
    '''
    signed = sorted([*light, *[-r for r in light]])
    for pair in combinations_with_replacement(signed, 2):
        if pair[0] == -pair[1]:
            continue
        for k in range(28):
            total = add(
                add(charge(pair[0], g), charge(pair[1], g), g),
                mul(k, qq['chi'][2], g), 
                g
            )
            if total == zero:
                conj = tuple(sorted([-r for r in pair]))
                bil.add(min((pair, k), (conj, (-k) % 28)))

    counts.append(
        dict(
            group=list(g), # Current group moduli
            single_channel_nonleading=len(classes), # Number of additional Yukawa operator solutions found
            bilinear_hc_classes=len(bil), # Number of distinct scalar bilinear classes after identifying Hermitian conjugates
            surviving_classes=classes if g == G else None # Complete list of higher-order Yukawa solutions is stored only for the full group \(G\)
        )
    )

'''Checking expected operator count'''
assert [a['single_channel_nonleading'] for a in counts] == [67, 14, 2]
assert [a['bilinear_hc_classes'] for a in counts] == [10, 2, 0]

D37 = subgroup([q['chi'][2], q['chi'][6]]) # Subgroup generated by X3 and X7
Dall = subgroup(q['chi']) # Subgroup generated by all 7 scalar field
assert len(D37) == 56 and len(Dall) == math.prod(G) == 5376

'''Soft scalar bilinear and phase-invariant interaction'''
Bcharge = mul(-2, q['chi'][2]) # Soft bilinear charge
assert Bcharge == (0, 0, 12, 6, 0)
assert subgroup([q['chi'][2], q['chi'][6], Bcharge]) == D37 # Checking whether the soft bilinear changes the generated subgroup

'''
Neutrino mass operator verification
Degree <=2 neutrino vertices, without assuming a chi7 insertion.
'''
labels = sorted([*range(1, 8), *range(-1, -8, -1)])
monomials = [()]+[(a,) for a in labels] + list(combinations_with_replacement(labels, 2))
'''All possible scalar monomials with zero, one or two scalar insertions'''

nu = {}
for i in range(3):
    for j in range(3):
        bare = add(mul(-1, q['L'][i]), q['n'][j])
        good = []
        for mono in monomials:
            total = bare
            for r in mono:
                total = add(total, charge(r))
            if total == (0,)*5:
                good.append(mono)
        if good:
            nu[(i+1, j+1)] = good

nuex = {
    (1, 1): [(-7, -1)], 
    (1, 2): [(-4, 7)], 
    (1, 3): [(-4, 7)], 
    (2, 2): [(4, 7)], 
    (2, 3): [(4, 7)], 
    (3, 2): [(-7, 5)], 
    (3, 3): [(-7, 5)]
}
assert nu == nuex, (nu, nuex)

'''
Instanton and gravitational charge sums. Calculates three charge sums:
1. Colour-sector charge sum.
2. Weak-isospin-sector charge sum.
3. Gravitational charge sum.
'''
def chsum(terms):
    '''
    terms ~ (field sector, multiplicity)
    returns = sum(mult_i * sum(charge(sector)))
    '''
    total = (0,)*5
    for key, mult in terms:
        for a in q[key]:
            total = add(total, mul(mult, a))
    return total


anomalies = dict(
    color=chsum([('Q', 2), ('u', -1), ('d', -1)]), # Color-charge sum
    weak=chsum([('Q', 3), ('L', 1)]), # Weak-isospin charge sum
    gravity_charge_sum=chsum([
        ('Q', 6),   # +6 = 3 colours × 2 weak components
        ('u', -3),  # -3 = 3 colours, conjugate-charge convention
        ('d', -3),  # -3 = 3 colours, conjugate-charge convention
        ('L', 2),   # 2 = 2 weak components
        ('e', -1),  # -1 = One charged-lepton singlet
        ('n', -1)   # -1 = One right-handed neutrino
    ]) # Gravitational charge sum
)
assert anomalies == dict(
    color=(1, 2, 6, 0, 1),
    weak=(1, 1, 0, 0, 2),
    gravity_charge_sum=(0, 0, 8, 8, 1)
)
assert add(
    mul(-4, q['chi'][2]),
    mul(4, q['chi'][6])
) == (0,)*5

out = dict(
    status='PASS', 
    group_order=math.prod(G), 
    D37_order=len(D37), 
    quotient_order=math.prod(G)//len(D37), 
    all_VEV_subgroup_order=len(Dall),
    sequential_operator_counts=counts, 
    neutrino_vertices={str(k): v for k, v in nu.items()}, 
    instanton_and_gravity_sums=anomalies,
    phase_invariant='(chi3_dagger)^4 chi7^4 allowed',
    soft_B_chi_charge=Bcharge, soft_B_preserves_D37_quotient=True,
    scope='Exact perturbative charge congruences and stated measure sums; not a complete gauged discrete or constituent completion'
)
'''
status: All programmed assertions passed
group_order: Order of the complete finite group
D37_order: Order of subgroup generated by X3, X7
quotient_order: Order of G/D37
all_VEV_subgroup_order: Order of subgroup generated by all seven scalar charges
sequential_operator_counts: Higher-order Yukawa and scalar bilinear counts
neutrino_vertices: Allowed neutrino operators up to degree two
instanton_and_gravity_sums: Three calculated charge sums
phase_invariant: Allowed phase-sensitive scalar interaction
soft_B_chi_charge: Charge of the specified soft bilinear
soft_B_preserves_D37_quotient: Whether including the bilinear changes the generated subgroup
scope: Explicit statement of the audit's limitations
'''

Path(__file__).with_name('flavour_charge_audit_v9.json').write_text(json.dumps(out, indent=2)+'\n')
print(json.dumps(out, indent=2))
