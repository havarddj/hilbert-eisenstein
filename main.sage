attach("./src/algdep.sage")
attach("./src/diagres.sage")
attach("./src/modforms.sage")
attach("./src/quadforms.sage")


def GS_unit(F, p, m=50):
    D = F.discriminant()
    assert kronecker_symbol(D, p) == -1, f"{p} is not inert in Q(sqrt({D}))"
    drd = DiagonalRestrictionDerivative(F, p, m)
    ct = drd[0]

    h = len(BinaryQF_reduced_representatives(D))
    e = genus_field_roots_of_1(D)
    print(f"number of roots of unity in HCF={e}")
    # To compute the valuation of the constant term, our choice of
    # normalisation is precisely so that the constant term is p^R where
    # R is the sum of the positive slopes of the vector of L-values
    Qs = BinaryQF_reduced_representatives(D)
    Lvals = [ZZ(Meyer(Q) * e) for Q in Qs]
    Lvals = [v / gcd(Lvals) for v in Lvals]
    for i in range(len(Qs)):
        print(f"{Qs[i]} has Meyer special value equal to", Lvals[i])
    # print(f"ct = {ct}")

    return GS_algdep(
        exp(ct) * QQ(p ^ (Meyer(F) * e / gcd(Lvals))), h, GS_val_vec(D))


def trace_test(D, p, bd=20, m=30):
    # assert kronecker_symbol(
    #     D, p) == -1, f"p = {p} should be split in Q(\sqrt {D})"
    # trace_test(69,17, bd=3) shows that the trace is not zero in general
    assert is_discriminant(
        D), f"D = {D} should be a (fundamental) discriminant"
    if kronecker_symbol(D, p) == 1:
        print(f"p = {p} split in Q(sqrt {D}), coherent case")
    elif kronecker_symbol(D, p) == -1:
        print(f"p = {p} inert in Q(sqrt {D}), incoherent case")
    else:
        print(f"p = {p} ramified in Q(sqrt {D}), wacky case")
    if is_fundamental_discriminant(D):
        print("D is a fundamental discriminant")
    for N in range(2, bd):
        print(f"Testing discriminants of conductor {N}")
        if N != p and D % N != 0 and is_prime(N):
            for Q in BinaryQF_reduced_representatives(N ^ 2 * D):
                f = DiagonalRestriction(Q, 1, m, pStab=p)
                print("f =", f)
                print("trace of f =", modform_trace(f, p), "\n")

    return 0


def p_new_test(D, p, bd=20, m=30):
    assert kronecker_symbol(
        D, p) == -1, f"p = {p} should be split in Q(\sqrt {D})"
    for N in range(2, bd):
        print(f"Testing discriminants of conductor {N}")
        if N != p:
            for Q in BinaryQF_reduced_representatives(N ^ 2 * D):
                f = DiagonalRestriction(Q, 1, m, pStab=p)
                # does not work
                # M = f.parent()
                # S = ModularSymbols(f.level()).cuspidal_submodule()
                # if S.dimension() > 0:
                #     S_new = S.new_subspace(N).q_expansion_basis()
                # else:

                # E_new = [
                #     M.eisenstein_subspace().old_submodule(
                #         ell).q_expansion_basis()
                #     for ell in divisors(f.level() / N)
                # ]

                # Mnew = M.span(Ebasis_new + Sbasis_new)

    return None
