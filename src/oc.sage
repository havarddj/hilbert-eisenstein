"""
Overconvergent modular forms of tame level 1

AUTHORS:

- Håvard Damm-Johnsen (2023): Initial version

Ported from Vonk's magma code.

"""

def hasse_power(p):
    """
    Returns the power of the Hasse invariant that has been lifted.
    """
    if p == 2:
        return 4
    elif p == 3:
        return 3
    else:
        return 1


def hasse_Eis_weight(p):
    """Returns the weight of the Eisenstein series to be used for
    lifting the Hasse invariant.
    """
    return hasse_power(p) * (p - 1)


def complementary_spaces(k, p, delta, delta_j, E4, E6):
    """Returns basis for a subspace of weight k complementary to
    image of multiplication by the Eisenstein series.
    """
    a = k % 3
    b = (k // 2) % 2

    d = dimension_modular_forms(1, k) -1
    e = dimension_modular_forms(1, k - hasse_Eis_weight(p)) - 1
    compl_space = []
    for j in [e + 1 .. d]:
        aj = delta_j * E6 ^ (2 * (d - j) + b) * E4 ^ a
        delta_j *= delta
        compl_space.append(aj)
    return compl_space, delta_j

def normalised_eis(wt,nterms):
    """Returns `nterms` terms of the q-expansion of wt `wt`
    Eisenstein series, normalised to have constant term 1

    """
    f = eisenstein_series_qexp(wt, prec=nterms)
    return f/f.coefficients()[0]

def overconvergent_basis(wt,p, mp, nterms, base_ring=False):
    """Computes Katz basis of overconvergent modular forms using
    Lauder's algorithms.
    """
    if not base_ring:
        base_ring = pAdicField(p, print_mode="val-unit")

    R.<q> = PowerSeriesRing(base_ring.integer_ring())

    Ep1 = R(normalised_eis(hasse_Eis_weight(p), mp))
    E4 = R(normalised_eis(4, mp))
    E6 = R(normalised_eis(6, mp))
    Delta = R(CuspForms(1,12).0.q_expansion(mp))
    delta_j = R(1)

    basis = []

    bound = floor(nterms * (p + 1) / (hasse_power(p) * p))
    for i in [0..bound]:
        compl_space, delta_j = complementary_spaces(
            wt + i * hasse_Eis_weight(p),
            p,
            Delta,
            delta_j,
            E4,E6)
        B = [p^(floor(i*hasse_power(p)*p/(p+1)))*Ep1^(-i)*f for f in compl_space]
        basis += [R(b) for b in B if b != 0]
        # print(f"i = {i}, len(B) = {len(basis)}")
    return basis


def ordinary_projection(G):
    """Compute the ordinary projection of G, in the space of
    overconvergent modular forms of weight 2 and tame level 1
    """
    Rp = G.coefficients()[0].parent()
    p = Rp.prime()

    Rpq = G.parent()

    prec = Rp.precision_cap()
    m = prec

    # Rp = ZqFM(p^2,m)              # make fixed precision p-adic ring
    
    
    bnd = ModularForms(weight=2+(p-1)*floor(prec*(p+1)/p)).dimension()-1
    print("number of terms needed for overconvergent basis:", bnd)
    
    B = overconvergent_basis(2,p,p*m,m,base_ring = Rp)
    print(f"number of terms of B[0] is {len(B[0].coefficients())}")
    # this is currently a bit wasteful when trying to compute SH points:
    # we computed this basis already to find constant term of DRD
    # TODO: cache this in a clever way

    ell = len(B)
    print(f"ell = {ell}")
    print(f"ell*p = {p*ell}" )
    Up = copy(zero_matrix(Rp, ell+1, ell+1))

    # image of basis under Up:

    for j in range(ell):
        print(p*(j-1))
    T = copy(zero_matrix(Rp,ell,ell))
    for i in range(ell):
        for j in range(ell):
            T[i,j] = B[i].coefficients()[p*j]
        
    print("found basis under Up")
    for i in range(ell):
        Ti = T[i]
        print(f"i = {i}, Ti = {Ti}")
        for j in range(ell):
            Bj = Ti.parent()(B[j].coefficients()[:ell])
            lj = ZZ(Bj[j])
            u = ZZ(Ti[j])/lj
            assert valuation(u,p) >= 0, f"encountered p in the denominator! for i = {i}, j = {j}";
            if lj == 0:
                print(j, Bj)
            Up[i,j] = Rp(u)
            Ti  -=  Up[i,j]*Bj;

    print("Computed Up matrix")
    print("random coeff:", Up[1,3])    
    Up_power = Up**(2*m)

    print("Computed power of Up matrix")

    print("first coeff:", Up_power[0,0])
    comb = find_in_space(G,B, K=Rp)
    
    comb_ord = comb*Up_power
    print("computed image of G under Up^2m")
    print(f"equals: {vector(comb_ord)}") 
    # Gord = Rpq(0)
        # for i in range(ell):
    #     Gord -= Rp(comb_ord[i]/comb[ell])*B[i]
    return sum(comb_ord[i]/comb[ell-1]*B[i] for i in range(ell))
