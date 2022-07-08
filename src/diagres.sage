# attach("quadforms.sage")
# load("modforms.sage")
from sage.modular.overconvergent.hecke_series import katz_expansions

def DiagonalRestriction(F,k,m, pStab=1):
    """Compute the diagonal restriction of the Hilbert Eisenstein
    series E_(k,k)(1,psi) when psi corresponds to an 'odd indicator function' on the class group, determined by the quadratic form F; or the p-stabilisation E^(p)(z) = E(z) - E(pz)
    """
    
    Fs, Forms = Delta_data(F, m)

    D = F.discriminant()
    D0 = fundamental_discriminant(D)
    f = sqrt(D/D0)
    
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
   

def DiagonalRestrictionDerivative(Q,p,m):
    """Compute the p-adic overconvergent modular form coming from the
    first order deformation of the Hilbert Eisenstein series E_(1,1)(1,psi) when psi
    corresponds to an 'odd indicator function' on the class group
    determined by the quadratic form Q. This is described in [DPV2];
    NB: Given that we are not taking p-adic logarithms, this is probably DPV1,
    i.e. parallel weight direction (PWD)
    """

    D = Q.discriminant()
    D0 = fundamental_discriminant(D)
    f = sqrt(D/D0)
    # for now, assume f = 0
    assert f == 1, "Only implemented for fundamental discriminants"
    assert kronecker_symbol(D,p) == -1, "p must be inert in Q(sqrt(D))"
    assert p >= 5, "for now overconvergent forms in sage only work when p > 3"
    F = QuadraticField(D)
    assert F.ideal(p).is_prime()

    ZZx.<x> = PolynomialRing(ZZ)
    QQp = pAdicField(p,m)
    Fp.<sqrtD> = QQp.extension(x^2-D0)
    assert Fp.degree() == 2
    

    Fs, Forms = Delta_data(Q, m)

    Rp.<q> = PowerSeriesRing(Fp, m)
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
        Diag_F += (2*q^n)*log(prod_n/p^prod_n.valuation())
    # print(Diag_F)
    # return Diag_F
    if Diag_F == 0:
        return 0
    elif p >= 5:
        M = katz_expansions(2,p,m-1,m,m)[0]  # katz_expansions returns ([q_exp_basis],E_p-1)
        # print(M)
        ct = is_overconvergent(M, Diag_F)
        return (ct+Diag_F).add_bigoh(m)


    

# def Diagonal_Restriction_DerivativeAPWD(F,k,m, pStab=1):
#     """Compute the p-adic overconvergent modular form coming from a
#     first order deformation of the Hilbert Eisenstein series
#     E_(k,k)(1,psi) when psi corresponds to an 'odd indicator function'
#     on the class group determined by the quadratic form F. This is
#     described in [DPV2]; this is the *anti-parallel weight deformation*!

#     """
        
#     As, Forms = Delta_data(F, m)

#     D = F.discriminant()
#     D0 = fundamental_discriminant(D)
#     f = sqrt(D/D0)
    
#     R.<q> = PowerSeriesRing(ZZ, m)
#     Diag_F = 0
        
#     for n in [1..m-1]:
#         Ap, Am = As[n-1]
#         coeff_n = 0
#         for i in range(len(Ap)):
#             if gcd(Ap[i], f*pStab) == 1:
#                 coeff_n +=  Ap[i]**(k-1)
#             if gcd(Am[i], f*pStab) == 1:
#                 coeff_n +=  (-1)**k*(-Am[i])**(k-1)
#         Diag_F += 4*coeff_n*q^n

#     M = ModularForms(Gamma0(f*pStab),2*k)
#     if Diag_F == 0:
#         return M(0)
#     else:    
#         ct = find_const_term(M.q_expansion_basis(m), Diag_F)    
#         return M((ct+Diag_F).add_bigoh(m))


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

