import random, time

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
    return lm %p

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

def padding_zero(lst: list, itr: int):
    return [ lst[i] if i in range(len(lst)) else 0  for i in range(len(lst) + max(itr, 0)) ]

def poly_add_modp(poly_a: list, poly_b: list, p: int):
    difc = max(len(poly_a), len(poly_b)) - min(len(poly_a), len(poly_b))
    return [ i + j %p  for i, j in zip(padding_zero(poly_a, difc), padding_zero(poly_b, difc)) ]

def poly_rmul_modp(poly: list, a: int, p: int):
    return [(item*a)%p for item in poly]

def poly_mul_modp(poly_a: list, poly_b: list, p: int):
    prod = [0 for i in range(len(poly_a) + len(poly_b) - 1)]
    for i in range(len(poly_a)):
        for j in range(len(poly_b)):
            prod[i + j] += (poly_a[i]*poly_b[j]) %p
    return [item %p for item in prod]

def poly_eval_modp(poly: list, r: int, p: int):
    r_pw = [r**i %p for i in range(len(poly))]
    return sum([(poly[i]*r_pw[i])%p for i in range(len(poly))])%p


def z_poly_modp(ind_list: list, p: int):
    zp = [1]
    for i in ind_list:
        zp = poly_mul_modp(zp, [-i, 1], p)
    return list(zp)

def z_root_eval_modp(root: list, r: int, p: int):
    value = 1
    for i in range(len(root)):
        value *= (r - root[i]) %p
    return value %p

def gen_ranindlst(lngth: int, rnge: int):
    temp = set()
    while len(temp) < lngth:
        temp.add(random.randint(0, rnge-1))
    temlst = list(temp)
    temlst.sort()
    return temlst

def poly_sum_modp(poly_list: list, p: int):
    ps = [0]
    for poly_i in poly_list:
        ps = poly_add_modp(ps, poly_i, p)
    return ps


### Newton Interpolation
def fdd_ep_modp(m: int, dpts: list, p: int):
    xpnts = [dpts[i][0] for i in range(len(dpts))][:m+1]
    val_sum = 0
    for j in range(len(xpnts)):
        xpnts_j = [xpnts[i] for i in list(range(j))+list(range(j+1, len(xpnts)))]
        val_sum += (dpts[j][1] *fp_inv(z_root_eval_modp(xpnts_j, xpnts[j], p), p))%p
    return val_sum %p

def newton_interpolation(mdata: list, p: int):
    xpnts = [mdata[i][0] for i in range(len(mdata))]
    enum_znx = [z_poly_modp(xpnts[:i], p) for i in range(len(mdata))]
    coef_znx = [poly_rmul_modp(enum_znx[i], fdd_ep_modp(i, mdata, p), p) for i in range(len(mdata))]
    return poly_sum_modp(coef_znx, p)

def newton_interpolation_append(nip: list, mdata: list, apndata: list, p: int):
    ttldata = mdata + apndata
    xpnts = [ttldata[i][0] for i in range(len(ttldata))]
    enum_znx = [z_poly_modp(xpnts[:i], p) for i in range(len(ttldata))]
    coef_znx_apnd = [poly_rmul_modp(enum_znx[i], fdd_ep_modp(i, ttldata, p), p) for i in range(len(mdata), len(ttldata))]
    nx_apnd = poly_sum_modp(coef_znx_apnd, p)
    return poly_add_modp(nip, nx_apnd, p)


from bn256 import order
dq = order
n = 4
l = 2**n

#xpts = gen_ranindlst(l, dq-1)
xpts = list(range(l))
datapts = [(xpts[i], random.randint(0, dq-1)) for i in range(l)]

xpts_apd = list(range(l, 2*l))
datapts_apd = [(xpts_apd[i], random.randint(0, dq-1)) for i in range(l)]
ttldatapts = datapts + datapts_apd

xpts_apd2 = list(range(2*l, 4*l))
datapts_apd2 = [(xpts_apd2[i], random.randint(0, dq-1)) for i in range(2*l)]
ttldatapts2 = datapts + datapts_apd + datapts_apd2

ts=time.time()
Nx = newton_interpolation(datapts, dq)
tn=time.time()
print('Newton Poly: ', tn-ts)
for i in range(l):
    if not poly_eval_modp(Nx, datapts[i][0], dq) == datapts[i][1]:
        print('error', i)

tsapn=time.time()
Nx_apnd = newton_interpolation_append(Nx, datapts, datapts_apd, dq)
tnapn=time.time()
print('Appended Newton Poly: ', tnapn-tsapn)
print('total time: ', tnapn-ts)
for i in range(2*l):
    if not poly_eval_modp(Nx_apnd, ttldatapts[i][0], dq) == ttldatapts[i][1]:
        print('error', i)

tsapn2=time.time()
Nx_apnd2 = newton_interpolation_append(Nx_apnd, ttldatapts, datapts_apd2, dq)
tnapn2=time.time()
print('2nd-Appended Newton Poly: ', tnapn2-tsapn2)
print('total time: ', tnapn2-ts)
for i in range(4*l):
    if not poly_eval_modp(Nx_apnd2, ttldatapts2[i][0], dq) == ttldatapts2[i][1]:
        print('error', i)
