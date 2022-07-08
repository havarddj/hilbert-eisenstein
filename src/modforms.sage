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
    # print(Kp)
    if M == []:
        print("Basis is empty")
        return(0)
    R.<Z> = PowerSeriesRing(Rp)
    A = []
    for j in [0..len(M)-1]:   # for each row of M do:
        # make row vector of lifts of coefficients to Kp, and append to A
        A.append([Rp(c.lift()) for c in M[j][1:]])
    # assert len(A[0]) == len(f.padded_list()[1:])  ## should be m-1
    print("trying to solve linear system")
    # fvec = [Rp(Kp(c)) for c in f.padded_list()[1:]]
    try:
        # print("f has coeffs", f.padded_list())
        # print("basis matrix given by", Mat)
        # the \ operator solves Ax = B for x, same as solve_right
        soln = Matrix(Rp,A).transpose() \ vector(Rp, f.padded_list()[1:])

        # soln = Matrix(Rp, A).transpose().solve_right(fvec)

        # print(f"kernel is given by {Y}")

    except ValueError as e:
        # print(e)
        if str(e) != "matrix equation has no solutions":
            raise
        else:
            print("not in space of forms")
            return(0)
    print("found constant term!")

    for l in range(len(soln)):
        print(f"valuation coordinate {l} = {valuation(soln[l],Kp.prime())}")
    # print(len(soln),len(A))
    assert len(soln) == len(A)
    CT = (sum(soln[n]*Rp(M[n][0].lift()) for n in range(len(A))))

    return(CT)
