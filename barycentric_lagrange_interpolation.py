import random

def fp_inv(a: int, p: int) -> int:
    # Extended euclidean algorithm to find modular inverses for integers
    a %= p
    if a == 0:
        return 0
    lm, hm = 1, 0
    low, high = a % p, p
    while low > 1:
        r = high // low
        nm, new = hm - lm * r, high - low * r
        lm, low, hm, high = nm, new, lm, low
    return lm % p

def fp_div(x, y, p:int):
    return x * fp_inv(y, p) % p

def fp_div_polys(a, b, p: int):
    """
    Long polynomial difivion for two polynomials in coefficient form
    """
    a = [x for x in a]
    o = []
    apos = len(a) - 1
    bpos = len(b) - 1
    diff = apos - bpos
    while diff >= 0:
        quot = fp_div(a[apos], b[bpos], p)
        o.insert(0, quot)
        for i in range(bpos, -1, -1):
            a[diff + i] -= b[i] * quot
        apos -= 1
        diff -= 1
    return [x % p for x in o]

'''
def lagrange_polynomial(mdata: list, p: int):
    x_coords = [mdata[i][0] for i in range(len(mdata))]
    y_coords = [mdata[i][1] for i in range(len(mdata))]
    base_poly = [0]
    for i in range(len(x_coords)):
        basis = [1]
        for j in range(len(x_coords)):
            if j == i:
                continue
            basis = poly_mul_modp(basis, [-x_coords[j], 1], p)
            basis = poly_rmul_modp(basis, fp_inv(x_coords[i] - x_coords[j], p), p)
        base_poly = poly_add(base_poly, poly_rmul_modp(basis, y_coords[i], p))
    return [coef %p for coef in base_poly]
'''

def padding_zero(lst: list, itr: int):
    return [lst[i] if i in range(len(lst)) else 0  for i in range(len(lst) + max(itr, 0))]

def poly_add(poly_a: list, poly_b: list):
    difc = max(len(poly_a), len(poly_b)) - min(len(poly_a), len(poly_b))
    return [i + j for i, j in zip(padding_zero(poly_a, difc), padding_zero(poly_b, difc))]

def poly_add_modp(poly_a: list, poly_b: list, p: int):
    difc = max(len(poly_a), len(poly_b)) - min(len(poly_a), len(poly_b))
    return [i + j %p for i, j in zip(padding_zero(poly_a, difc), padding_zero(poly_b, difc))]

def poly_eval_modp(poly: list, r: int, p: int):
    r_pow = [(r**i)%p for i in range(len(poly))]
    return sum([(poly[i]*r_pow[i])%p for i in range(len(poly))])%p

def poly_rmul_modp(poly: list, a: int, p: int):
    return [(item*a)%p for item in poly]

def poly_mul_modp(poly_a: list, poly_b: list, p: int):
    prod = [0 for i in range(len(poly_a) + len(poly_b) - 1)]
    for i in range(len(poly_a)):
        for j in range(len(poly_b)):
            prod[i + j] += (poly_a[i]*poly_b[j]) %p
    return [coef%p for coef in prod]

def z_poly_modp(ind_list: list, p: int):
    zp = [1]
    for i in ind_list:
        zp = poly_mul_modp(zp, [-i, 1], p)
    assert len(zp) == len(ind_list) + 1
    return list(zp)

'''
def z_poly_modp_np(ind_list: list, p: int):
    zp = np.poly1d([1])
    for i in ind_list:
        zp = np.polymul(zp, np.poly1d([-i, 1])) #too large anyway
        zp = np.poly1d([zp[i]%p for i in range(len(zp))])
    #assert len(zp) == len(ind_list)+1
    return list(zp)
'''

'''
def fp_div_polys_np(a, b, p): # slower than above
    qp = np.polydiv(a[::-1], b[::-1])[0]
    return [int(qp[i])%p for i in range(len(qp))][::-1]
'''

def z_inv_modp(ind_list: list, t: int, p: int):
    zp = [1]
    for i in range(t):
        if i in ind_list:
            continue
        zp = poly_mul_modp(zp, [-i, 1], p)
    return list(zp)

def gen_ranindlst(lngth: int, rnge: int):
    temp = set()
    while len(temp) < lngth:
        temp.add(random.randint(0, rnge-1))
    temlst = list(temp)
    temlst.sort()
    return temlst

'''
def diff_modp(poly: list, p: int):
    return [i*poly[i+1] %p for i in range(len(poly)-1)]
'''

def poly_sum_modp(poly_list: list, p: int):
    ps = [0]
    for poly in poly_list:
        ps = poly_add_modp(ps, poly, p)
    return ps

def prod_modp(ind_list: list, p: int):
    prod = 1
    for i in ind_list:
        prod = prod* i %p
    return prod


import time
# 2**7 ~ 0.2s, 2**8 ~ 0.7s, 2**9 ~ 2.6s, 2**10 ~ 11.5s, 2**11 ~ 50s
p = 65000549695646603732796438742359905742570406053903786389881062969044166799969
n = 9
l = 2**n
print('length: ', l)

### Barycentric Lagrange Interpolation, SIAM REVIEW
t1=time.time()

indice = gen_ranindlst(l, p-1)
corr_y = [random.randint(1, p-1) for i in range(l)]

A_poly = z_poly_modp(indice, p)
indice_list = [indice[:i] + indice[i+1:] for i in range(l)]
A_poly_list = [fp_div_polys(A_poly, [-indice[i], 1], p) for i in range(len(indice))]

t2=time.time()
print('setup: ', t2-t1)

w_list = [fp_inv(prod_modp(poly_add_modp(indice_list[i], [-indice[i] for _ in range(l-1)] , p), p), p) for i in range(len(indice))]
d2_list = [poly_rmul_modp(A_poly_list[i], (-1)*w_list[i]*corr_y[i], p) for i in range(len(indice))]
P2 = poly_sum_modp(d2_list, p)

t3=time.time()
print('interpolation: ', t3-t2)

'''
# test part
for i in range(len(indice)//4):
    if poly_eval_modp(P2, indice[i], p) != corr_y[i]:
        print('bugs', i)
'''

'''
t1=time.time()
indice = gen_ranindlst(l, p-1)
corr_y = [random.randint(1, p-1) for i in range(l)]

A_poly = z_poly_modp(indice, p)
Ap_poly = diff_modp(A_poly, p)
indice_list = [indice[:i] + indice[i+1:] for i in range(l)] #len correct
A_poly_list = [fp_div_polys(A_poly, [-indice[i]%p, 1], p) for i in range(len(indice))]
t2=time.time()

#t1=time.time()
#A_poly_list = [z_poly_modp(indice_list[i], p) for i in range(l)]
#A_poly_list = [fp_div_polys(A_poly, [-indice[i], 1], p) for i in range(len(indice))]
#A_poly_list_np = [fp_div_polys(A_poly, [-indice[i]%p, 1], p) for i in range(len(indice))]

c_list = [corr_y[i]* fp_inv(poly_eval_modp(Ap_poly, indice[i], p),p) %p for i in range(len(indice))]
d_list = [poly_rmul_modp(A_poly_list[i], c_list[i], p) for i in range(len(indice))]
P = poly_sum_modp(d_list, p)
t3=time.time()

for i in range(len(indice)):
    if poly_eval_modp(P, indice[i], 1033) != corr_y[i]:
        print('bugs')
        break
print('done')
'''
