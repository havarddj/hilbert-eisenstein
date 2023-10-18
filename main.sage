"""
Compute Gross-Stark units using p-adic modular forms

AUTHORS:

- Håvard Damm-Johnsen (2023): Initial version

"""

attach("./src/algdep.sage")
attach("./src/oc.sage")
attach("./src/diagres.sage")
attach("./src/modforms.sage")
attach("./src/quadforms.sage")
attach("./src/test.sage")
attach("./src/stark_heegner.sage")
# attach("./src/lvals.sage")


def GS_unit(D, p, nterms=0, pprec=0):
    """
    Compute Brumer--Stark unit in the narrow Hilbert class field of
    $\mathbb Q(\sqrtD)$, which is a `p`-unit.

    """
    # pick quadratic form with minimal special value <=> GS-unit has
    # minimal p-valuation
    assert kronecker_symbol(D, p) == -1, f"{p} is not inert in Q(sqrt({D}))"
    assert is_fundamental_discriminant(
        D), "{D} is not a fundamental discriminant"
    _, Q = sorted([(Q.Zagier_L_value().abs(), Q)
                   for Q in BinaryQF_reduced_representatives(D)])[0]
    return GS_unit_BQF(Q, p, nterms, pprec)


def GS_unit_BQF(F, p, nterms=0, pprec=0):

    D = F.discriminant()
    assert kronecker_symbol(D, p) == -1, f"{p} is not inert in Q(sqrt({D}))"

    if nterms == 0:
        # number of terms necessary to compute diagonal restriction derivative
        # to desired precision
        # heuristic, based on examples; should be a bit bigger for p small,
        # but scale about linearly in p
        nterms = 6 * p + ceil(10 * 1 / p)

    if pprec == 0:
        target_prec = max(10, 3 * GS_val_vec(F.discriminant())[0])
        bound = floor(nterms * (p + 1) / (hasse_power(p) * p))
        loss_factor = floor(bound * hasse_power(p) * p / (p + 1))
        pprec = target_prec + loss_factor
        print(f"p-adic precision = {pprec}")

        bnd = ModularForms(weight=2 +
                           (p - 1) * floor(pprec *
                                           (p + 1) / p)).dimension() - 1
        nterms = max(bnd, nterms)
        print(f"Number of modular form coefficients to be computed = {nterms}")

    drd = diagonal_restriction_derivative(F, p, nterms, pprec=pprec)
    ct = drd[0]

    h = len(BinaryQF_reduced_representatives(D))
    e = genus_field_roots_of_1(D)

    print(f"number of roots of unity in HCF={e}")
    # To compute the valuation of the constant term, our choice of
    # normalisation is precisely so that the constant term is p^R where
    # R is the sum of the positive slopes of the vector of L-values
    Qs = BinaryQF_reduced_representatives(D)
    Lvals = [ZZ(Q.Zagier_L_value() * e) for Q in Qs]
    Lvals = [v for v in Lvals]
    for i in range(len(Qs)):

        print(f"{Qs[i]} has L-value equal to", Lvals[i])
    u = exp(e * ct) * p ^ (-F.Zagier_L_value() * e)
    print(f"Attempting to find algebraic relations for {u}:\n")
    # return algdep_p_adic(u, 2 * h)
    if F.conductor() > 1:
        P = algdep_p_adic(u, h)
        K = QuadraticField(D)
        PF = PolynomialRing(K, "x")(P)
        print("Is irreducible?", PF.factor())
        print("Discriminant factors:", K.discriminant().factor())
        # for f, _ in PolynomialRing(K, "x")(P).factor():
        #     H.<h> = K.extension(f)
        #     print("Absolute disc:", H.absolute_discriminant().factor())

        return P
    return GS_algdep(u, h, D)


# def trace_test(D, p, bd=20, m=30):
#     # assert kronecker_symbol(
#     #     D, p) == -1, f"p = {p} should be split in Q(\sqrt {D})"
#     # trace_test(69,17, bd=3) shows that the trace is not zero in general
#     assert is_discriminant(
#         D), f"D = {D} should be a (fundamental) discriminant"
#     if kronecker_symbol(D, p) == 1:
#         print(f"p = {p} split in Q(sqrt {D}), coherent case")
#     elif kronecker_symbol(D, p) == -1:
#         print(f"p = {p} inert in Q(sqrt {D}), incoherent case")
#     else:
#         print(f"p = {p} ramified in Q(sqrt {D}), wacky case")
#     if is_fundamental_discriminant(D):
#         print("D is a fundamental discriminant")
#     for N in range(2, bd):
#         print(f"Testing discriminants of conductor {N}")
#         if N != p and D % N != 0 and is_prime(N):
#             for Q in BinaryQF_reduced_representatives(N ^ 2 * D):
#                 f = diagonal_restriction(Q, 1, m, pStab=p)
#                 print("f =", f)
#                 print("trace of f =", modform_trace(f, p), "\n")


#     return 0
def find_difference(f, g):
    # print("Normalised Eis:", f[1] / g[1] * g)
    return f - f[1] / g[1] * g


# def test_hyp(D,p,m):
#     R.<q> = PowerSeriesRing(q,m)
#     E = ModularForms(5).basis()[0].q_expansion(50)
#     for Q in BinaryQF_reduced_representatives(D):
