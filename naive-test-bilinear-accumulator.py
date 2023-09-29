import bn256, random

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

def prod_modp(lst: list):
    pd = 1
    for pt in lst:
        pd = pd*pt %p
    return pd

def gf(m:int, Xlst: list):
    Xs = [xi+s for xi in Xlst]
    return P.scalar_mul((prod_modp(Xs)*m)%p)

def accu_add(Vf, xp):
    assert xp not in X
    X.append(xp)
    return Vf.scalar_mul(xp+s)

def accu_del(Vf, xp):
    assert xp in X
    return Vf.scalar_mul(fp_inv(xp+s, p))


p = bn256.order
s = random.randint(1, p-1)
h = random.randint(1, p-1)

P = bn256.curve_G
Q = bn256.twist_G
Qpub = s*Q
H = h*P

u = random.randint(1, p-1)
X = [1,3,5]
V = [P]

# user 1, 2, 3

# steps 1~4
user_cred = [(random.randint(1, p-1), random.randint(1, p-1), random.randint(1, p-1), random.randint(1, p-1)) for i in range(3)]
user_op = [(user_cred[i][2]*user_cred[i][0]+user_cred[i][3], user_cred[i][2]*user_cred[i][1]) for i in range(3)]
lhs_c = user_cred[1][3] + user_cred[1][2]*(user_cred[1][0]+user_cred[1][1]*h)
rhs_c = user_op[1][0] + user_op[1][1]*h
print('steps 1~4 are verified: ', lhs_c == rhs_c)

# steps 5~7
a = [random.randint(1, p-1), random.randint(1, p-1), random.randint(1, p-1)]
for i in range(len(a)):
    V.append(accu_add(V[-1], a[i]))

S = [fp_inv(a[i]+s, p)*(user_op[1][0]*P + P) for i in range(3)]

lhs_s = bn256.optimal_ate(a[1]*Q+Qpub, S[1])
rhs_s = bn256.optimal_ate(Q, user_op[1][0]*P+P)
print('steps 5~7 are verified: ', lhs_s==rhs_s)

# steps 8~10, can be taken into (dynamic) pairing-based accumulator solely
lhs_v = bn256.optimal_ate(a[1]*Q+Qpub, V[1])
rhs_v = bn256.optimal_ate(Q, V[2])
print('steps 8~10 are verified: ', lhs_v==rhs_v)
