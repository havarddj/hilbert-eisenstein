def diagonal_restriction(F,k, pStab=1):
    """Compute the diagonal restriction of the Hilbert Eisenstein
    series E_(k,k)(1,psi) when psi corresponds to an 'odd indicator
    function' on the class group, determined by the quadratic form F;
    If pStab != 1, returns instead the p-stabilisation, where p is the
    value of pStab.

    """

    D = F.discriminant()
    N = F.conductor()

    M = ModularForms(Gamma0(N*pStab),2*k)
    m = M.sturm_bound()

    Fs, Forms = Delta_data(F, m)

    # assert N == 1, "nontrivial conductor not currently supported"
    
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
            if gcd(Ap, N*pStab) == 1:
                coeff_n +=  Ap**(k-1)

            if gcd(Am, N*pStab) == 1:
                coeff_n += (-1)**k*(-Am)**(k-1)
        Diag_F += 4*coeff_n*q^n

    if Diag_F == 0:
        return M(0)
    else:    
        ct = find_const_term(M.q_expansion_basis(m), Diag_F + O(q^m))    
        return M((ct+Diag_F).add_bigoh(m))
   

def diagonal_restriction_derivative(Q,p,m, pprec=None):
    """Compute the first m nonconstant terms of the p-adic
    overconvergent modular form coming from the first order
    deformation of the Hilbert Eisenstein series E_(1,1)(1,psi) when
    psi corresponds to an 'odd indicator function' on the class group
    determined by the quadratic form Q.

    """
    if pprec == None:
        pprec = m
    D = Q.discriminant()    
    f = Q.conductor()

    assert f == 1, "Only implemented for fundamental discriminants"
    assert kronecker_symbol(D,p) == -1, "p must be inert in Q(sqrt(D))"

    if has_negative_fundamental_unit(D):
        return 0

    F = QuadraticField(D)
    
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
    M = overconvergent_basis(2,p, m*p, m, base_ring=FFp)

    ct = is_overconvergent(M, Diag_F)
    return (ct+Diag_F).add_bigoh(m)    
