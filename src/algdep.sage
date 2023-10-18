def Qp2split(a):
    # return x,y in Qp where a = x + ky for k a generator of K = a.parent()
    k = a.parent().gen()
    x = a.trace() / 2

    y = (a - x) / k
    return (x, y)


def is_power_of_p(a, p):
    assert p.is_prime()
    assert a != 0, "will cause philosophical uproar"
    if a == 1:
        return True
    if a in Integers():
        return ((len(factor(a)) == 1) and (factor(a)[0][0] == p))
    if 1 / a in Integers():
        return ((len(factor(ZZ(1 / a))) == 1)
                and (factor(ZZ(1 / a))[0][0] == p))
    return False


def GS_val_vec(D):
    # returns a list of decreasing integers which bound the valuations of coefficients of the minimal polynomial of Gross-Stark units
    Lvals = []
    zero_flag = true  # to ensure we adjoin a value of zero only once per pair
    e = genus_field_roots_of_1(D)
    for Q in BinaryQF_reduced_representatives(D):
        val = Meyer(Q) * e
        if val > 0:
            Lvals.append(val)
        else:
            if val == 0:
                if zero_flag:  # if we didn't count a zero already
                    Lvals.append(0)  # then we can add 0
                    zero_flag = false  # but prevent counting the next one by setting the flag
                else:  # if we did count a zero, then we don't add zero but reset the flag
                    zero_flag = true
    # Next sort things
    # Lsort = sorted([v / gcd(Lvals) for v in Lvals])
    Lsort = sorted([v for v in Lvals])
    print(f"Lsort = {Lsort}")
    # next add up

    return list(reversed([sum(Lsort[:i])
                          for i in [1..len(Lsort)]]))  # eg [12,11,8]


def algdepQp(a, d):
    """Basic algebraic recognition over Qp:

    Given element a in Qp with finite precision, return integer
    polynomial of degree d which has roots approximately equal to a.

    There's no proof of correctness, because a is given to finite precision.
    """

    Qp = a.parent()
    p = Qp.prime()
    m = a.precision_absolute()
    N = p ^ ((5 * m / 7).floor())

    M = copy(zero_matrix(ZZ, d + 2, d + 1))
    for j in [0..d]:
        M[j, j] = 1
        M[j, d + 1] = ZZ(a ^ j)
    M[d + 1, d] = N
    print(f"M = {M}")

    short_vec = M.LLL()[0]
    return sum(x ^ i * short_vec[i] for i in [0..d])


def algdepQp2(a, d):
    # algebraic recognition of algebraic number/Q of degree d from
    # approximation in deg n unramified ext of Qp
    K = a.parent()
    n = K.degree()
    p = K.prime()
    m = a.precision_absolute()
    N = p ^ ((5 * m / 7).floor())
    M = copy(zero_matrix(ZZ, n, d + 2))
    for j in [0..d]:

        c1, c2 = Qp2split(a ^ j)
        M[0, j] = ZZ(c1)
        M[1, j] = ZZ(c2)
    M[0, d + 1] = N
    M[1, d + 1] = N
    Mdash = M.right_kernel().basis()
    # print(Mdash)
    Mred = Matrix(Mdash).LLL()
    # print(Mred)
    return QQ['x'](list(Mred[0])[:d + 1])
    # M = copy(zero_matrix(ZZ,d+n+1,d+n+1))

    # for j in [0..d]:
    #     M[j, j] = 1
    #     for i,c in enumerate(Qp2split(a^j)):
    #         M[j, d+1 + i] = ZZ(c)
    #         M[d+1+i,d+1+i] = N

    # Mnew = M.BKZ()
    # print(Mnew)
    # short_vec = Mnew[1]
    # print(f"short vector: {short_vec}")
    # return QQ['x'](short_vec[:d+1].list())


def p_adic_eltseq(a):
    """Given p-adic number a in field K, return list of
    Qp-coefficients obtained by writing a as a linear combination of
    basis vectors with respect to the power basis of Kp
    """
    Kp = a.parent()
    _, _, psi = Kp.vector_space()
    l = psi(a).list()
    g = Kp.gens()[0]
    # assert Kp(a) == sum(l[i]*g^i for i in range(len(l))), f"a = {a} != {sum(l[i]*g^i for i in range(len(l)))}"
    # note: this assertion fails for inexplicable reasons!
    return l


def algdep_p_adic(a, s):
    """Algebraic recognition based on algorithm in [Gaudry et al,
    §4.2]

    Recognise an element `a' of a p-adic extension as a root of
    a degree s integer polynomial.
    """

    if a.valuation() < 0:
        return ZZ["x"](list(algdep_p_adic(1 / a, s))[::-1])
    Kp = a.parent()
    p = Kp.prime()
    d = Kp.degree()

    N = floor(Kp.precision_cap())
    # Define the matrix A
    A = copy(zero_matrix(ZZ, s + 1, d))
    for i in range(s + 1):
        A[i] = p_adic_eltseq(a ^ i)

    M = block_matrix([[A], [p ^ N * identity_matrix(ZZ, d)]])
    K = M.left_kernel()

    return ZZ["x"](K.matrix().LLL()[0].list()[:s + 1])


def GS_algdep(a, s, D):

    if a.valuation() < 0:
        return GS_algdep(1 / a, s, D)
    ZZx.<x> = PolynomialRing(ZZ)
    
    Kp = a.parent()
    d = Kp.degree()
    p = Kp.prime()

    m = a.precision_absolute()
    N = floor(Kp.precision_cap())
    # N = p^m

    val_list = GS_val_vec(D)
    assert len(val_list) == s / 2  # if it has zeros this might throw errors
    # note that the extra zeros beyond the first don't contribute to any terms.

    # primitive p^2-1 th root of 1 in Qp^2
    zeta = cyclotomic_polynomial(p ^ 2 - 1).roots(Kp)[0][0]

    print(f"Lower bound on p-valuation of coefficients  = {val_list}")

    P = 0

    for k in [0..ZZ(ceil((p ^ 2 - 1)/2))]:
        print(f"Running algdep on exp(log(u))*zeta^{k}")
        b = a * zeta ^ k
        A = copy(zero_matrix(ZZ, s / 2 + 1, d))
        for i in range(s / 2 + 1):

            if i == s / 2:
                A[i] = p_adic_eltseq(b ^ i)
            else:
                A[i] = p_adic_eltseq((b ^ i + b ^ (s - i)) * p ^ val_list[i])

        M = block_matrix([[A], [p ^ N * identity_matrix(ZZ, d)]])
        K = M.left_kernel()
        short_vec = K.matrix().LLL()[0].list()[:s / 2 + 1]
        # print(f"short vector: {short_vec}")
        P1 = x ^ (s / 2) * short_vec[s / 2]  # start with the middle coeff

        for j in range(s / 2):
            c = short_vec[j]
            P1 += QQ(p ^ val_list[j]) * c * (x ^ j + x ^ (s - j))

        P1 = P1 / gcd(list(P1))
        if P1[0] != 0:
            # print(f"P1 = {P1} has first coeff {P1[0].factor()}\n")
            if P1[0] < 0:
                P1 = -P1
            if P1[d - 1] != 0 and is_power_of_p(P1[0], p) and is_HCF(D,P1,p):
                return P1
    # print(algdep(a,deg))
    return P
