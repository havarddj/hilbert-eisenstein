def Qp2split(a):
    # return x,y in Qp where a = x + ky for k a generator of K = a.parent()
    k = a.parent().gen()
    x = a.trace()/2

    y = (a -x)/k
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

def GS_val_vec (D):
    # returns a list of decreasing integers which bound the valuations of coefficients of the minimal polynomial of Gross-Stark units
    Lvals = []
    zero_flag = true                   # to ensure we adjoin a value of zero only once per pair  
    e = genus_field_roots_of_1(D)
    for Q in BinaryQF_reduced_representatives(D):
        val = Meyer(Q)*e
        if val > 0:
            Lvals.append(val)
        else:
            if val == 0:
                if zero_flag:   # if we didn't count a zero already
                    Lvals.append(0)   # then we can add 0
                    zero_flag = false   # but prevent counting the next one by setting the flag
                else:                   # if we did count a zero, then we don't add zero but reset the flag 
                    zero_flag = true
    # next sort things
    Lsort = sorted([v / gcd(Lvals) for v in Lvals])
    print(f"Lsort = {Lsort}")
    # next add up 

    return list(reversed([sum(Lsort[:i]) for i in [1..len(Lsort)]]))  # eg [12,11,8]

def algdepQp(a, d):
    # algebraic recognition of algebraic number/Q of degree d from
    # approximation in Qp

    P = 0
    Qp = a.parent()
    p = Qp.prime()
    m = a.precision_absolute()
    N = p ^ ((5 * m / 7).floor())
    M = copy(zero_matrix(QQ, d + 2, d + 2))

    for j in [0..d]:
        M[j, j] = 1
        M[j, d + 1] = QQ(a ^ j)
    M[d + 1, d + 1] = N

    short_vec = M.LLL()[0]
    return sum(x ^ i * short_vec[i] for i in [0..d])


def algdepQp2(a, d):
    # algebraic recognition of algebraic number/Q of degree d from
    # approximation in Qp2 (deg 2 unramified ext of Qp)

    P = 0
    K = a.parent()
    p = K.prime()
    m = a.precision_absolute()
    N = p ^ ((5 * m / 7).floor())
    M = copy(zero_matrix(QQ, d + 3, d + 3))

    for j in [0..d]:
        M[j, j] = 1
        c1, c2 = Qp2split(a ^ j)

        M[j, d + 1] = QQ(c1)
        M[j, d + 2] = QQ(c2)

    M[d + 1, d + 1] = N
    M[d + 2, d + 2] = N

    short_vec = M.LLL()[0]
    return sum(x ^ i * short_vec[i] for i in [0..d])


def GS_algdep(a, deg, val_list):
    P = 0

    Fp = a.parent()
    p = Fp.prime()
    ZZx.<x> = PolynomialRing(ZZ)

    m = a.precision_absolute()
    N = p ^ ((5 * m / 7).floor())
    # N = p^m

    d = ZZ(deg / 2)
        
    assert len(val_list) == d #if it has zeros this might throw errors
    # note that the extra zeros beyond the first don't contribute to any terms.
    
    # primitive p^2-1 th root of 1 in Qp^2
    zeta = cyclotomic_polynomial(p ^ 2 - 1).roots(Fp)[0][0]
    # print(f"zeta equals {zeta}")
    print(f"partial sums = {val_list}")
    
    for k in [0..ZZ((p ^ 2 - 1)/2)]:
        print(f"k = {k}")
        b = a * zeta ^ k

        M = copy(zero_matrix(ZZ, d + 3, d + 3))
        for j in [0..d]:
            # print(f"j = {j}")
            M[j, j] = 1
            if j == d:
                cj = Fp(b ^ d)
            else:
                cj = Fp(p ^ val_list[j] * (b ^ j + b ^ deg - j))

            X, Y = Qp2split(cj)
            # M[j, d + 1] = ZZ(X/p^valuation(X,p))*p^valuation(X,p)  # get 1st component
            M[j, d + 1] = QQ(X).floor()
            # M[j, d + 2] = ZZ(Y/p^valuation(Y,p))*p^valuation(Y,p)  # and 2nd component
            M[j, d + 2] = QQ(Y).floor()

            
        M[d + 1, d + 1] = N
        M[d + 2, d + 2] = N
        # print("M is given by:")
        # print(M)
        # print("\n")
        Mnew = M.BKZ(block_size=d, float_type=qd1)  # LLL outputs new basis which hopefully is good
        # print("\nLLL ran successfully")
        # print(Mnew)
        

        P1 = x ^ d * Mnew[0][d]  # start with the middle coeff
        
        for j in [0..d - 1]:
            c = Mnew[0][j]
            P1 += QQ(p ^ val_list[j]) * c * (x ^ j + x ^ (deg - j))

        P1 = P1 / gcd(list(P1))
        if P1[0] != 0:
        
            print(f"P1 = {P1} has first coeff {P1[0].factor()}\n")
            if P1[0] < 0:
                P1 = -P1
            if is_power_of_p(P1[0], p) and P1[d] != 0:
                P = P1
    # print(algdep(a,deg))
    return P
