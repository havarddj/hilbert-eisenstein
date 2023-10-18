# attach("quadforms.sage")
# load("modforms.sage")
# from sage.modular.overconvergent.hecke_series import katz_expansions
# from sage.modular.overconvergent.hecke_series import higher_level_katz_exp

def divides_exactly(d,n):
    """return True if d || n (d divides n exactly), else False"""
    return True if (n % d == 0 and n/d % d != 0) else False

def diagonal_restriction(F,k,m, pStab=1):
    """Compute the diagonal restriction of the Hilbert Eisenstein
    series E_(k,k)(1,psi) when psi corresponds to an 'odd indicator
    function' on the class group, determined by the quadratic form F;
    or the p-stabilisation E^(p)(z) = E(z) - E(pz)
    """
    
    Fs, Forms = Delta_data(F, m)

    D = F.discriminant()
    N = F.conductor()
    
    R.<q> = PowerSeriesRing(ZZ, m)
    q = R.gen()
    Diag_F = 0
        
    for n in [1..m-1]:
        Fp, Fm = Fs[n-1]
        assert len(Fp) == len(Fm)
        coeff_n = 0
        for i in range(len(Fp)):
            Ap = Fp[i][0]
            Am = Fm[i][0]
            Bp = Fp[i][1]
            Bm = Fm[i][1]
            if gcd(Ap, pStab) == 1:
                coeff_n +=  Ap**(k-1)

            if gcd(Am, pStab) == 1:
                coeff_n += (-1)**k*(-Am)**(k-1)
        Diag_F += 4*coeff_n*q^n
    # for Q in BinaryQF_reduced_representatives(D/N^2):
    #     M2Z = MatrixSpace(ZZ,2,2)
    #     for H in Hecke_matrices(N):
    #         if N > 1 and F.is_equivalent(Q.matrix_action_right(M2Z(H))):
    #             print(f"F = {F} is a lift of {Q}")
    #             Diag_F -= 2*diagonal_restriction(Q,k,m,pStab=pStab)(q^N)
    return Diag_F

    M = ModularForms(Gamma0(N*pStab),2*k)
    if Diag_F == 0:
        return M(0)
    else:    
        ct = find_const_term(M.q_expansion_basis(m), Diag_F + O(q^m))    
        return M((ct+Diag_F).add_bigoh(m))
   

def diagonal_restriction_derivative(Q,p,m, pprec=None):
    """Compute the p-adic overconvergent modular form coming from the
    first order deformation of the Hilbert Eisenstein series E_(1,1)(1,psi) when psi
    corresponds to an 'odd indicator function' on the class group
    determined by the quadratic form Q. This is described in [DPV2].
    """
    if pprec == None:
        pprec = m
    D = Q.discriminant()    
    f = Q.conductor()
    D0 = ZZ(D/f^2)

    # for now, assume f = 0
    # assert f == 1, "Only implemented for fundamental discriminants"
    assert kronecker_symbol(D,p) == -1, "p must be inert in Q(sqrt(D))"
    F = QuadraticField(D)
    assert F.narrow_class_group().order() != 1
    
    ZZx = PolynomialRing(ZZ,"x")
    QQp = pAdicField(p, pprec, type = 'capped-rel', print_mode="val-unit")
    if p == 2:
        FFp.<unr2gen> = QQp.extension(x^2+x+1)
        sqrtD = sqrt(FFp(D))
    else:
        FFp.<sqrtD> = QQp.extension(x^2-D)
    assert FFp.degree() == 2
    

    Fs, Forms = Delta_data(Q, m)

    Rp.<q> = PowerSeriesRing(FFp, m)
    Diag_F = 0
        
    for n in [1..m-1]:
        Fp, Fm = Fs[n-1]
        prod_n = 1
        for i in range(len(Fp)):
            a,b,c = Fp[i]
            if gcd(a, p) == 1:
                prod_n *= (-b+n*sqrtD)/(2*a)
            a,b,c = Fm[i]
            if gcd(a, p) == 1:
                prod_n /= (-b+n*sqrtD)/(2*a)
        # This is p-adic because sqrtD is in the field Fp
        Diag_F += log(prod_n/p^prod_n.valuation(p))*q^n
    # print(Diag_F)

    if Diag_F == 0:
        return 0
    bnd = ModularForms(weight=2+(p-1)*floor(pprec*(p+1)/p)).dimension()-1
    # print("number of terms needed for overconvergent basis:", bnd)  # 
    # if m < bnd:
    #     print(f"Need at least {bnd} terms to compute katz basis")
    #     return None
    M = overconvergent_basis(2,p, m*p, m, base_ring=FFp)

    ct = is_overconvergent(M, Diag_F)
    return (ct+Diag_F).add_bigoh(m)


def is_trace_0(f,N,p):
    """
    Test if the trace of a modular form f from level Np to level p vanishes. What we actually do is check whether <f,g> for each g coming from level p, which we can do explicitly by writing them both in terms of an eigenbasis {f_i} for M_k(Gamma0(Np)) and computing <f_i,f_i>^2 using built-in magma functionality.
    """

    k = f.weight()
    m = f.prec()
    M = ModularForms(Gamma0(N*p),k,prec=m)
    try:
        f = M(f)
    except:
        print("f not found in M_k(Gamma0(Np))")
    pNewforms = [M(g) for g in
                 ModularForms(Gamma0(N),k).q_expansion_basis(m)] + [M(h) for h in
                                                                    M.newforms()]
    pNewspace = M.span_of_basis(pNewforms)

    return M.eisenstein_subspace()


def trace_of_diagonal_restriction(F,k,m):
    """
    Given quadratic form F of discriminant D and conductor f, return
    the first m coefficients of the modular form Tr^f_1 E_{k,F}
    """
    D = F.discriminant()
    Ek = diagonal_restriction(F,k,m)
    
    f = F.conductor()
    D0 = ZZ(D/f^2)
 
    F0 = F.reduced_form()
    fibre = []
    for Q0, l in ring_class_fibres(D0,f):
        if F0 in l:
            fibre = l
            break


    if len(fibre) == 0:         # sum is empty, so return Ek
        return Ek
    print("Number of terms in sum", len(fibre))
    series = [ diagonal_restriction(Q,k,m,flip=True) for Q in fibre]
    print("series considered:\n",  [Ek] + series)
    return Ek + (-1)^k * sum(series)
    # compute phi^A for each A equivalent to F
        
    
