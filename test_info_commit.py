from sympy.matrices import Matrix
from secp256k1 import G, N, P, add, multiply

n = 8 # padding to n=4 at least
sz = int(math.sqrt(2**n))
print("size of matrix: ", sz) # expect to have maximal 16, in practise n=2

poly_coef = [random.randint(1, N) for _ in range(2**n)] # A
Poly_Matrix = Matrix(sz, sz, poly_coef) # A
RD_Matrix = Matrix(sz, sz, [random.randint(1, N) for _ in range(2**n)]) # B
while RD_Matrix.det() == 0: # compulsory
    RD_Matrix = Matrix(sz, sz, [random.randint(1, N) for _ in range(2**n)])

c = random.randint(1, sz)
print('value of c: ', c)
K_v = [random.randint(1, N) for _ in range(2*c)]

# default
lv = Matrix(1, sz, [x**(sz*i) for i in range(sz)])
rv = Matrix(sz, 1, [x**i for i in range(sz)])

def gnr_vmd(kvlst, cinp, szinp):
    lambda_kvl_vamd = [[kvlst[j]**(szinp*i) for i in range(szinp)] for j in range(cinp)]
    theta_kvl_vdmd = [[kvlst[cinp+j]**i for i in range(szinp)] for j in range(cinp)]
    return Matrix(lambda_kvl_vamd), Matrix(theta_kvl_vdmd), lambda_kvl_vamd, theta_kvl_vdmd

# gnr rdmatr for btoken, kvlst to be the commited poly; grant v, u to btoken; x:=usercred
def gnr_infocom_key(kvlst, cinp, szinp, polymatr, rdmatr): # define carrying data; cinp,szinp
    fLambda_Matrix, fTheta_Matrix, _, _ = gnr_vmd(kvlst, cinp, szinp)
    fGamma_Matrix = fLambda_Matrix * (polymatr + rdmatr) # prover's work
    fOmega_Matrix = rdmatr * fTheta_Matrix.transpose() # prover's work
    # evaluate x, or not
    fv = (polymatr + rdmatr) * rv # prover's work
    fu = lv * rdmatr # prover's work
    return fLambda_Matrix, fTheta_Matrix, fGamma_Matrix, fOmega_Matrix, fv, fu

def infocom_verify(gmam, lamm, omgm, tham, vv, vu): # using secp256k1
    fscrt = random.randint(1, N)
    fgrv = gmam * rv
    flbv = lamm * vv
    fr1 = fgrv[0].as_poly().eval(x, fscrt)%N
    fl1 = flbv[0].as_poly().eval(x, fscrt)%N
    v1 = (multiply(G, fr1)==multiply(G, fl1))

    flvom = lv*omgm
    futm = vu*tham.transpose()
    fr2 = flvom[0].as_poly().eval(x, fscrt)%N
    fl2 = futm[0].as_poly().eval(x, fscrt)%N
    v2 = (multiply(G, fr2)==multiply(G, fl2))
    return v1 & v2

def recover_poly(v_vec, u_vec):
    hx = lv * v_vec
    gx = u_vec * rv
    trp = sympy.expand((hx-gx)[0]).as_poly()
    return trp.all_coeffs()[::-1]


Lambda_Matrix, Theta_Matrix, Gamma_Matrix, Omega_Matrix, v, u = gnr_infocom_key(K_v, c, sz, Poly_Matrix, RD_Matrix)

# print("v: ", v)
print("verification 1: ", Gamma_Matrix * rv == Lambda_Matrix * v) # in ecc format, another curve
# print("u: ", u)
print("verification 2: ", lv*Omega_Matrix == u*Theta_Matrix.transpose()) # in ecc format, another curve
### test recovering
print("poly coef recoverable: ", recover_poly(v, u) == poly_coef) #done
# test final verifying
vf = infocom_verify(Gamma_Matrix, Lambda_Matrix, Omega_Matrix, Theta_Matrix, v, u)
print("infocommit verification: ", vf)
