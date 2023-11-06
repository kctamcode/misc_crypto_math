import math, random

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

def gcdExtended(a, b):
    # Base Case
    if a == 0 :
        return b,0,1

    gcd,x1,y1 = gcdExtended(b%a, a)

    # Update x and y using results of recursive
    # call
    x = y1 - (b//a) * x1
    y = x1

    return gcd, x, y


p = 47
q = 53
n = p*q
phi_n = (p-1)*(q-1)
g = 3

Accu = g
sc = [] # shortcut, required to have gdc(m, phi(n)) == 1


# add sc into argument
def accu_AddElem(x, At): #At equiv u/sc
    if x not in sc:
        assert gcdExtended(x, n)[0] == 1
        sc.append(x)
        return At**x %n
    else:
        return At

def accu_DelElem(x, At): #At equiv u/sc
    if x not in sc:
        return At
    else:
        sc.remove(x)
        return At**fp_inv(x, phi_n) %n

def accu_MemWitCreate(x, At): #At equiv u/sc
    if x not in sc:
        return At
    else:
        return At**fp_inv(x, phi_n) %n

def accu_VerMem(x, wx, At):
    return wx**x %n == At

def accu_NonMemWitCreate(x, At): #At equiv u/sc
    u = 1 # this should be stored in the system
    for i in sc:
        u *= i

    _, a, b = gcdExtended(u, x)
    while a < 0:
        a = a + x
        b = b - u #negative
    assert a*u + b*x == 1
    
    #a_inv_exp = (1+ (-a*u %phi_n))* fp_inv(x, (p-1)*(q-1)) %phi_n
    b0 = -b %phi_n

    t = g**u %n
    d = g**b0 %n

    assert t**a %n == d**x * g %n
    return a, d, x

def accu_VerNonMem(x, ux, At):
    return At**ux[0] %n == ux[1]**x * g %n


X = [5, 7, 11]
for i in X:
    assert gcdExtended(i, phi_n)[0] == 1
    Accu = accu_AddElem(i, Accu)
Accu = accu_DelElem(11, Accu)
Accu = accu_AddElem(11, Accu)
print('Add and Del are: ', Accu == 3**(5*7*11)%n)

w5 = accu_MemWitCreate(5, Accu)
v5 = accu_VerMem(5, w5, Accu)
print('VerMem Part: ', v5)

nm17 = accu_NonMemWitCreate(17, Accu)
nvm17 = accu_VerNonMem(17, nm17, Accu)
print('NonVerMem Part: ', nvm17)


def accu_update_add_memwit(w, x):
    return w**x %n

def accu_update_del_nonmemwit(x, u): # for small x's
    a = u[0]*x
    d = u[1]
    return a, d

def accu_update_del_memwit(w, x, accun): ### could generate bugs
    _, a, b = gcdExtended(x, sc[0])
    while a < 0:
        a = a + sc[0]
        b = b - x #negative
    assert a*x + b*sc[0] == 1

    #a_inv_exp = (1+ (-a*x %phi_n))* fp_inv(sc[0], phi_n) %phi_n
    b0 = b %phi_n

    s = g**x %n
    t = g**b0 %n

    assert s**a * t**sc[0] %n == g
    return w**a * accun**b0 %n

def accu_update_add_nonmemwit(x, u, accuo): ##true in theory
    _, a0, r0 = gcdExtended(x, u[2])
    while a0 < 0:
        a0 = a0 + u[2]
        r0 = r0 - x #negative
    assert a0*x + r0*u[2] == 1

    u0 = u[0]*a0 %u[2]
    r = (u0*x-u[0])//u[2] %phi_n
    u1 = u[1]*accuo**r %n
    return u0, u1, u[2]


'''
nm17_7 = accu_update_del_nonmemwit(7, nm17)
print(nm17, nm17_7)
accu_7 = accu_DelElem(7, Accu)
print(Accu, accu_7)
print(accu_VerNonMem(17, nm17_7, accu_7))
'''
'''
accu17=accu_AddElem(17, Accu)
w5p17=accu_update_add_memwit(w5, 17)
print('update memwit with adding: ', accu_VerMem(5, w5p17, accu17))
'''

'''
accu7 = accu_DelElem(7, Accu)
w5_7 = accu_update_del_memwit(w5, 7, accu7)
print('update memwit with deleting: ', accu_VerMem(5, w5_7, accu7))
'''
'''
nm17p3 = accu_update_add_nonmemwit(3, nm17, Accu)
Accu = accu_AddElem(3, Accu)
print('update nonmemwit with adding: ', accu_VerNonMem(17, nm17p3, Accu))
'''
