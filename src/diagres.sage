# attach("quadforms.sage")
# load("modforms.sage")
from sage.modular.overconvergent.hecke_series import katz_expansions
from sage.modular.overconvergent.hecke_series import higher_level_katz_exp

def DiagonalRestriction(F,k,m, pStab=1):
    """Compute the diagonal restriction of the Hilbert Eisenstein
    series E_(k,k)(1,psi) when psi corresponds to an 'odd indicator
    function' on the class group, determined by the quadratic form F;
    or the p-stabilisation E^(p)(z) = E(z) - E(pz)
    """
    
    Fs, Forms = Delta_data(F, m)

    D = F.discriminant()
    D0 = fundamental_discriminant(D)
    f = F.conductor()
    
    R.<q> = PowerSeriesRing(ZZ, m)
    Diag_F = 0
        
    for n in [1..m-1]:
        Fp, Fm = Fs[n-1]
        coeff_n = 0
        for i in range(len(Fp)):
            Ap = Fp[i][0]
            Am = Fm[i][0]
            if gcd(Ap, f*pStab) == 1:
                coeff_n +=  Ap**(k-1)
            if gcd(Am, f*pStab) == 1:
                coeff_n +=  (-1)**k*(-Am)**(k-1)
        Diag_F += 4*coeff_n*q^n

    M = ModularForms(Gamma0(f*pStab),2*k)
    if Diag_F == 0:
        return M(0)
    else:    
        ct = find_const_term(M.q_expansion_basis(m), Diag_F)    
        return M((ct+Diag_F).add_bigoh(m))
   

def DiagonalRestrictionDerivative(Q,p,m, pprec=None):
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

    ZZx.<x> = PolynomialRing(ZZ)
    QQp = pAdicField(p,pprec, print_mode="val-unit")
    FFp.<sqrtD> = QQp.extension(x^2-D0)
    assert FFp.degree() == 2
    

    Fs, Forms = Delta_data(Q, m)

    Rp.<q> = PowerSeriesRing(FFp, m)
    Diag_F = 0
        
    for n in [1..m-1]:
        Fp, Fm = Fs[n-1]
        prod_n = 1
        for i in range(len(Fp)):
            a,b,c = Fp[i]
            if gcd(a, f*p) == 1:
                prod_n *= (-b+n*sqrtD)/(2*a)
            a,b,c = Fm[i]
            if gcd(a, f*p) == 1:
                prod_n /= (-b+n*sqrtD)/(2*a)
        # This is p-adic because sqrtD is in the field Fp
        Diag_F += log(prod_n/p^prod_n.valuation(p))*q^n
    # print(Diag_F)

    if Diag_F == 0:
        return 0
    bnd = ModularForms(weight=2+(p-1)*floor(pprec*(p+1)/p)).dimension()-1
    print("number of terms needed for overconvergent basis:", bnd)
    if m < bnd:
        print(f"Need at least {bnd} terms to compute katz basis")
        return None
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

