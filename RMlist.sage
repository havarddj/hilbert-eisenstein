mfrak = 17
QFlist = BinaryQF_reduced_representatives(mfrak)
print(QFlist)
for Q in QFlist:
    RMlist = []
    for n in [1..S-1]:
        # print(n)
        RMset = get_RM_set(n, Q)
        # print(RMset)
        if len(RMset) !=0:
            RMlist.append(RMset)
    print("Computed RM sets \n")
    for j in [0..deltam]:
        Delta = 0
        for n in [1..S-1]:
            for l in [0..len(RMlist)-1]:
                RM = RMlist[l]
                print(RM[0])
                if RM != []:
                    Q2 = RM[0]
                    tau = stable_root(Q2[0])
                    print(tau)
                    Delta += q^n *2^d*psi(tau)*Q2[0][0]^(kj[j]-1)
    DeltaList.append(Delta)
print("Computed Delta Q = ", Q)

def get_RM_set(n, Q):
    """
    Compute the set of RM-points RM(n, tau)_f using algorithm from
    [lauder-vonk?]_

    Input:
    - $n$ an integer
    - a reduced binary quadratic form with stable root tau
    Output:
    - the set of augmented RM points of discriminant n^2 D, as a list
      of lists [F,gamma_n] where F is a quadratic form and gamma_n a
      determinant n matrix with integer coeffs

    Note that the norm of the ideal corresponding to F is precisely
    the first coefficient, a.
    """
    D = Q.discriminant()
    D0 = D.squarefree_part()
    f = sqrt(D/D0)
    # First, we do a stupid computation to find possible elements of
    # stabilisers of RMpts
    Stab = []
    for a in [1..-1,step=-1]:
        for b in [1..-1,step=-1]:
            for c in [1..-1,step=-1]:
                for d in [1..-1,step=-1]:
                    if a*d-b*c == 1:
                        s = Matrix(ZZ,[[a,b],[c,d]])
                        if s^(12) == Matrix(ZZ,[[1,0],[0,1]]) and -s not in Stab:
                            Stab.append(s)

    Mn = []
    for d in divisors(n):
        if gcd(d, n/d) == 1:
            for j in [0..d-1]:
                Mn.append(Matrix(ZZ,[[n/d,j], [0,d]]))

    def is_in_SL2Z(Q, gamman, gnp, Stab):
        # """Returns True if $\gamma_n'\mathrm{Stab(\tau)}\gamma_n^-1
        # \subset \Sl_2(\Z), and false otherwise"""    
        flag = True
        for s in Stab:
            if Q.matrix_action_left(s) == Q:
                Mat = gnp*s*(gamman^-1)
                for i in [0..1]:
                    for j in [0..1]:
                        if not Mat[i][j].is_integer():
                            flag = False
        return(flag)

    RMn = []
    for k in [0..len(Mn)-1]:
        gamman = Mn[k]
        flag = False                # if True, then there exists some gamma_n'
                                    # for which is_in_SL2Z is true
        for gnp in Mn[:k-1]:
            if is_in_SL2Z(Q,gamman, gnp, Stab):
                flag = True
        if not flag:
            Q2 = Q.matrix_action_left(gamman)
            if gcd(f, Q2[0]) == 1:
                RMn.append([Q2.reduced_form(), gamman])

    return(RMn)

def stable_root(Q):
    """ Returns stable root of indefinite quadratic form Q
    """
    if Q[0] == 0:
        return(0)
    return((-Q[1]+sqrt(Q.discriminant()))/(2*Q[0]))
