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
    return f/f.padded_list()[0]

def overconvergent_basis(wt,p, mp, nterms, base_ring=False):
    """Computes Katz basis of overconvergent modular forms using
    Lauder's algorithms.
    """
    if not base_ring:
        base_ring = pAdicField(p, print_mode="val-unit")

    R.<q> = PowerSeriesRing(base_ring.integer_ring())

    Ep1 = R(normalised_eis(hasse_power(p)*(p-1), mp))
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
        # B = [p^(floor(i*hasse_power(p)*p/(p+1)))*Ep1^(-i)*f for f in compl_space]
        B = [p^(floor(i*hasse_power(p)*p/(p+1)))*Ep1^(-i)*f for f in compl_space]  # cf alan's code
        basis += [R(b) for b in B if b != 0]
        # print(f"i = {i}, len(B) = {len(basis)}")
    return basis


def ordinary_projection(G):
    """Compute the ordinary projection of G, in the space of
    overconvergent modular forms of weight 2 and tame level 1
    """
    Rp = G[0].parent()
    p = Rp.prime()

    Rpq = G.parent()

    prec = Rp.precision_cap()
    m = bnd

    
    bnd = ModularForms(weight=2+(p-1)*floor(prec*(p+1)/p)).dimension()-1
    
    assert m >= bnd, f"number of terms needed for overconvergent basis: {bnd}"
    
    B = overconvergent_basis(2,p,p*m,m,base_ring = Rp)

    # # this is currently a bit wasteful when trying to compute SH points:
    # # we computed this basis already to find constant term of DRD
    # # TODO: cache this in a clever way

    ell = len(B)

    Up = copy(zero_matrix(Rp, ell, ell))

    # image of basis under Up:
 
    T = matrix(Rp, [B[i].padded_list()[0:p*ell:p] for i in range(ell)])
    assert T.dimensions() == (ell,ell)
    
    for i in range(ell):
        Ti = T[i]
        for j in range(ell):
            Bj = T[i].parent()(B[j].padded_list()[:ell])
            lj = Bj[j]
            u = ZZ(Ti[j])/ZZ(lj)
            assert valuation(u,p) >= 0, f"encountered p in the denominator! for i = {i}, j = {j}, u = {u}, Tij = {T[i,j]}"
            
            Up[i,j] = Rp(u)
            Ti  -=  Up[i,j]*Bj;

    ord_proj = Up**(2*m)

    # write G in terms of basis 
    G_vector = find_in_space(G,B, K=Rp)

    # apply ordinary projection
    comb_ord = G_vector*ord_proj

    # rewrite in terms of basis vectors
    return sum(comb_ord[i]*B[i] for i in range(ell))
