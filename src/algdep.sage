def Qp2split(a):
    # return x,y in Qp where a = x + ky for k a generator of K = a.parent()
    K = a.parent()
    p = K.prime()
    m = a.precision_absolute()
    F = Qp(p, m)
    x = F(0)
    y = F(0)
    e = a.expansion()
    for n in range(len(e)):
        # for some ridiculous reason the length of a[n] depends on
        # whether or not the corresponding coefficient is non-zero,
        # instead of simply having entries equal to zero
        if len(e[n]) > 0:
            x += e[n][0] * p ^ n
        if len(e[n]) == 2:
            y += e[n][1] * p ^ n

    return (x, y)

def is_power_of_p(a,p):
    assert p.is_prime()
    assert a != 0, "will cause philosophical uproar"
    if a == 1:
        return True
    if a in Integers():
        return ((len(factor(a)) == 1) and (factor(a)[0][0] == p))
    if 1/a in Integers():
        return ((len(factor(ZZ(1/a))) == 1) and (factor(ZZ(1/a))[0][0] == p))
    return False
    
    

def algdepQp(a, d):
    # algebraic recognition of algebraic number/Q of degree d from
    # approximation in Qp

    P = 0
    Qp = a.parent()
    p = Qp.prime()
    m = a.precision_absolute()
    N = p ^ ((5 * m / 7).floor())
    M = copy(zero_matrix(RR, d + 2, d + 2))

    for j in [0..d]:
        M[j, j] = 1
        M[j, d + 1] = QQ(a ^ j)
    M[d + 1, d + 1] = N

    short_vec = M.LLL()[0][:]
    return sum(x ^ i * short_vec[i] for i in [0..d])


def algdepQp2(a, d):
    # algebraic recognition of algebraic number/Q of degree d from
    # approximation in Qp2 (deg 2 unramified ext of Qp)

    P = 0
    K = a.parent()
    p = K.prime()
    m = a.precision_absolute()
    N = p ^ ((5 * m / 7).floor())
    M = copy(zero_matrix(RR, d + 3, d + 3))

    for j in [0..d]:
        M[j, j] = 1
        c1, c2 = Qp2split(a ^ j)

        M[j, d + 1] = QQ(c1)
        M[j, d + 2] = QQ(c2)

    M[d + 1, d + 1] = N
    M[d + 2, d + 2] = N
    print(M)
    short_vec = M.LLL()[0][:]
    return sum(x ^ i * short_vec[i] for i in [0..d])


def GS_algdep(a, deg, Lvals):
    P = 0
    Fp = a.parent()
    p = Fp.prime()
    ZZx.<x> = PolynomialRing(ZZ)
    print(f"p equals {p}")
    m = a.precision_absolute()
    # N = p ^ ((5 * m / 7).floor())
    N = p^m
    Lpos = sorted([l for l in Lvals if l >= 0],
                  reverse=True)  # eg [1,3,8] from [1, 3, -3, 8, -1]

    partial_sums = [sum(Lpos[-i:]) for i in [1..len(Lpos)]]  # eg [12,11,8]

    d = ZZ(deg / 2)
    assert len(partial_sums) == d
    # primitive p^2-1 th root of 1 in Qp^2
    zeta = cyclotomic_polynomial(p ^ 2 - 1).roots(Fp)[0][0]
    print(f"zeta equals {zeta}")
    for k in [0..p ^ 2 - 1]:
        print(f"k = {k}")
        b = a * zeta ^ k

        M = copy(zero_matrix(ZZ, d + 3, d + 3))
        for j in [0..d]:
            # print(f"j = {j}")
            # print(f"partial sums = {partial_sums}, \n {partial_sums[d-j-1]}")
            M[j, j] = 1
            if j == d:
                cj = Fp(b ^ d)
            else:
                cj = Fp(p ^ partial_sums[d - j - 1] * (b ^ j + b ^ deg - j))
            # print("precision of cj =", cj.precision_absolute(), "m = ", m)
            X, Y = Qp2split(cj)
            M[j, d + 1] = ZZ(X)  # get 1st component

            M[j, d + 2] = ZZ(Y)  # and 2nd component

        M[d + 1, d + 1] = N
        M[d + 2, d + 2] = N
        # print("M is given by:")
        # print(M)
        # print("\n")
        Mnew = M.LLL()  # LLL outputs new basis which hopefully is good
        # print("\nLLL ran successfully")
        # print(Mnew)
        

        P1 = x ^ d * Mnew[0][d]  # start with the middle coeff
        
        for j in [0..d - 1]:
            c = Mnew[0][j]
            P1 += QQ(p ^ partial_sums[d - j - 1]) * c * (x ^ j + x ^ (deg - j))

        P1 = P1 / gcd(list(P1))
        if P1[0] != 0:
        
            print(f"P1 = {P1} has first coeff {P1[0].factor()}\n")
            if P1[0] < 0:
                P1 = -P1
            if is_power_of_p(P1[0], p) and P1[d] != 0:
                P = P1
    # print(algdep(a,deg))
    return P
