attach("src/algdep.sage")
attach("src/diagres.sage")
attach("src/modforms.sage")
attach("src/quadforms.sage")


def GS_unit(F, p, m=50):
    D = F.discriminant()
    assert kronecker_symbol(D, p) == -1
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
    for i in range(len(Qs)):
        print(f"{Qs[i]} has Meyer special value equal to", Lvals[i])
    print(f"ct = {ct}")

    return GS_algdep(exp(ct) * QQ(p ^ (Meyer(F) * e)), h, Lvals)
