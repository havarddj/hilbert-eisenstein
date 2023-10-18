def p_adic_eltseq(a):
    """Given p-adic number a in field K, return list of
    Qp-coefficients obtained by writing a as a linear combination of
    basis vectors with respect to the power basis of Kp
    """
    Kp = a.parent()
    _, _, psi = Kp.vector_space()
    l = psi(a).list()
    g = Kp.gens()[0]
    
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


def is_HCF(D,P,p):
    """Given polynomial P, check

    - if extension of F = Q(sqrt D) generated by P has discriminant a
    power of D,

    - if extension is Galois,

    - if roots of P has correct valuations in extension with respect to L-values above the fixed prime p
    """

    F = QuadraticField(D)
    R.<x> = PolynomialRing(F)
    Pirr = R(P).factor()[0][0]

    H.<h> = F.extension(Pirr)

    if not log(H.absolute_discriminant(),D) in ZZ:
        print(f"Discriminant of H not power of {D}")
        return False

    if not H.is_galois_relative():
        print("H is not Galois over Q(sqrt D)")
        return False
    
    
    frakp = H.prime_factors(p)[0]

    root_valn_list = [r.ord(frakp) for r,_ in P.roots(H)]
    L_valn_list = [Q.Zagier_L_value()*genus_field_roots_of_1(D) for Q in BinaryQF_reduced_representatives(D)]
    if not (sorted(root_valn_list) == sorted(L_valn_list)):
        print(f"valuation lists are different: {root_valn_list} != {L_valn_list}")
        return False
    print("found candidate narrow Hilbert class field")
    return True    

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

