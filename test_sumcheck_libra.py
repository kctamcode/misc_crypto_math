import sympy, time, random, math

x, y, z = sympy.symbols('x y z')
polyt = sympy.poly(x**2+y**3+x*z**2+z*y) # make up one from f(x), naive homogenuous

def naive_homogns(polycoef: list):
    deg = len(polycoef) - 1
    assert deg >= 4
    nhpoly = sympy.poly(0, x)
    for i in range(deg):
        if i % 2 == 0:
            nhpoly = nhpoly + sympy.poly(polycoef[i] * x**i * y**(i+1) * z**(3*deg-2*i-1))
        else:
            nhpoly = nhpoly + sympy.poly(polycoef[i] * x**i * y**(deg+i) * z**(3*deg-2*i-deg))
    return nhpoly


ranlist = [5, 17]
binlist = [0, 1]

def sumcheck_libra(spolyt, rlst):
    binlist = [0, 1]
    # round 0
    sH0 = 0
    for b0 in binlist:
        st0Hx = spolyt.eval(x, b0)
        for b1 in binlist:
            st1Hxy = st0Hx.eval(y, b1)
            for b2 in binlist:
                st2Hxyz = st1Hxy.eval(z, b2)
                sH0 += int(st2Hxyz)
    stemp_f1 = sympy.poly(0, x)
    # round 1
    for b1 in binlist:
        stemp_Hxz = spolyt.eval(y, b1)
        for b2 in binlist:
            stemp_Hx = stemp_Hxz.eval(z, b2)
            #print("temp eval", temp_Hx)
            stemp_f1 = stemp_f1 + stemp_Hx
    sH1 = stemp_f1.eval(0) + stemp_f1.eval(1)
    # round 2
    stemp_f2 = sympy.poly(0, y)
    stemp_2Hx = spolyt.eval(x, rlst[0])
    for b2 in binlist:
        stemp_2Hxz = stemp_2Hx.eval(z, b2)
        stemp_f2 = stemp_f2 + stemp_2Hxz
    sH2 = stemp_f2.eval(y, 0) + stemp_f2.eval(y, 1)
    # round 3
    stemp_3Hx = spolyt.eval(x, rlst[0])
    stemp_f3 = stemp_3Hx.eval(y, rlst[1])
    sH3 = stemp_f3.eval(z,0) + stemp_f3.eval(z,1)

    return [sH0, sH1, sH2, sH3]

SCTHL = [12, 12, 113, 9898]
SCTNHL = [7, 7, 3906, 107114502385142082]
print("test sumcheck: ", sumcheck_libra(polyt, ranlist)==SCTHL)
print("test sumcheck and naive hmgns: ", sumcheck_libra(naive_homogns([1,1,1,1,1,1,1]),ranlist)==[7,7,3906,107114502385142082])
