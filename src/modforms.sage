def find_const_term(M, f, K=QQ):
    """given a q-expansion basis M, and f a q-expansion in M with
    missing constant term, find constant term of f.
    """
    if M == []:
        print("Basis is empty")
        return(0)
    R.<Z> = PowerSeriesRing(K)
    A = [R(M[j]).padded_list()[1:] for j in [0..len(M)-1]]
    Mat = Matrix(K, A)
    try:
        # print("f has coeffs", f.padded_list())
        # print("basis matrix given by", Mat)
        soln = Mat.transpose() \ vector(K, f.padded_list()[1:])
    except ValueError as e:
        # print(e)
        if str(e) != "matrix equation has no solutions":
            raise
        else:
            print("not in space of forms")
            return(0)
    CT = (sum(soln[n]*M[n] for n in [0..len(M)-1]))[0]
    return(CT)

def is_overconvergent(M,f):
    """given a q-expansion basis M, and f a q-expansion in space of
    katz expansions M without constant term, find constant term
    of f if f is in M, or complain that f is not in M

    We are essentially solving the linear system Ax = f, where A is
    the matrix constructed by removing the constant terms of the forms
    in M.

    """
    Kp = f[0].parent()
    Rp = Kp.integer_ring()
    m = min([len(f.padded_list())]
            + [len(M[i].padded_list()) for i in range(len(M))])
    # print(Kp)
    if M == []:
        print("Basis is empty")
        return(0)
    R.<Z> = PowerSeriesRing(Rp)
    A = []
    for g in M: 
        # make row vector of lifts of coefficients to Kp, and append to A
        # A.append([Rp(c.lift()) for c in row[1:]])
        A.append([Rp(Kp(c)) for c in g.padded_list()[1:m]])
    print("Solving linear system")

    try:
        # print("f has coeffs", f.padded_list())
        # print("basis matrix given by", Mat)
        # the \ operator solves Ax = B for x, same as solve_right
        soln = Matrix(Rp,A).transpose() \ vector(Rp, f.padded_list()[1:m])

        # soln = Matrix(Rp, A).transpose().solve_right(fvec)
        # soln = Matrix(Kp,A).kernel().basis()[0]
        # print(f"kernel is given by {Y.basis()[0]} etc")

    except ValueError as e:
        # print(e)
        if str(e) != "matrix equation has no solutions":
            raise
        else:
            print("not in space of forms")
            return(0)
    print("found constant term!")

    # for i, x in enumerate(soln):
    #     print(f"valuation coordinate {i} = {valuation(x,Rp.prime())}")
    # print(f"length of soln vec = {len(soln)} and length of M = {len(M)}")
    CT = sum(soln[n]*Rp(M[n][0]) for n in range(len(soln)))

    return(CT)

def find_in_space(f, M, K = QQ):
    """given a set of q-expansions M, and f a q-expansion, write f as
    linear combination of forms in M.
    """
    if K == QQ:
        K = f[0].parent()
    if "Power Series" not in str(f.parent()):
        print("Warning: f should be power series. Currently f is in ", f.parent())
    if M == []:
        return 0 
    
    R.<Z> = PowerSeriesRing(K)
    m = min([len(f.padded_list())]
            + [len(M[i].padded_list()) for i in range(len(M))])
    
    A = [R(M[i]).padded_list()[:m] for i in range(len(M))]

    Mat = Matrix(K, A)

    try:
        v = vector(K, f.padded_list()[:m])
        soln = Mat.transpose() \ v

    except ValueError as e:
        if str(e) != "matrix equation has no solutions":
            raise
        else:
            print("f is not in K-span of M for K = ",K)
            return 0
    return(soln)


def modform_trace(f,p, m=30):
    """compute trace of f from level Np to level p
    RMK: I'm not sure if this is correct if N is composite!
    """
    N = f.level()/p
    assert N in ZZ, f"level of f = {f.level()} must be divisible by p = {p}"
    N = ZZ(N)


    M = f.parent()

    # problem with newforms:
    if "Newform" in str(type(f)):
        f = M(f.q_expansion(m))
    # confusing point: if M has level Np, then Mold_submodule(N)
    # returns the images of the degeneracy maps from level p, i.e. is
    # precisely the p-new space in my notation.
    Mold = M.old_submodule(N) # everything coming from level p
    try:
        Mnew = Mold.complement() # everything else, i.e. coming from
                                 # level N = Q.conductor()
    except NotImplementedError as e:
        print("failed to compute new/old decomposition")
        return None
                             
    B = Mnew.q_expansion_basis(m) + Mold.q_expansion_basis(m)

    combo = find_in_space(f.q_expansion(m),B)

    dimold = Mold.dimension()
    dimnew = Mnew.dimension()
    dim = M.dimension()
    assert dimold + dimnew == dim
    fnew = sum(Mnew.basis()[i]*combo[i] for i in
                              range(dimnew))
    fold = f - fnew
    print(f"fnew = {fnew}")
    print(f"fold = {fold}")

    # now to compute the composite of the degeneracy map and the
    # trace, we need to find a preimage of fold under the product map
    # (f,g) ~> (f(q)+g(q^N)). "ll" means "lower level"
    Mll = ModularForms(Gamma0(p))
    ll_dim = Mll.dimension()
    R.<q> = PowerSeriesRing(QQ)
    im_basis = [g for g in Mll.q_expansion_basis(m)] + [g(q^N) for g in Mll.q_expansion_basis(m)]

    ll_combo = find_in_space(fold.q_expansion(m), im_basis)

    f1 = sum([Mll.basis()[i]*ll_combo[i] for i in range(ll_dim)])
    f2 = sum([Mll.basis()[i]*ll_combo[i+ll_dim] for i in range(ll_dim)])
    print(f"fold equals f_1(q) + f_2(q^N), where f1 = {f1} and f2 = {f2}")
    TN = Mll.hecke_operator(N)
    tr = Gamma0(N).index()*f1 + TN(f2)  # (N+1)f_1 + T(N)f_2
    return Mll(tr)



