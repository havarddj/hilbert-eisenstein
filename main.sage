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


def GS_unit(D, p, nterms=0, pprec=30):
    """
    Compute Brumer--Stark unit in the narrow Hilbert class field of
    $\mathbb Q(\sqrtD)$, which is a `p`-unit.

    """
    # pick quadratic form with minimal special value <=> GS-unit has
    # minimal p-valuation
    _, Q = sorted([(Q.Zagier_L_value().abs(), Q)
                   for Q in BinaryQF_reduced_representatives(D)])[0]
    return GS_unit_BQF(Q, p, nterms, pprec)


def GS_unit_BQF(F, p, nterms=0, pprec=30):
    if pprec == 0:
        pprec = 2 * ceil(log(F.discriminant()) * p)
        print(f"p-adic precision = {pprec}")
    if nterms == 0:
        # number of terms necessary to compute diagonal restriction derivative
        # to desired precision
        nterms = ceil(log(p) * pprec)
        print(f"Number of modular form coefficients computed = {nterms}")

    D = F.discriminant()
    assert kronecker_symbol(D, p) == -1, f"{p} is not inert in Q(sqrt({D}))"
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
    print(f"running algdep on {u}")
    # return algdep_p_adic(u, 2 * h)
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
