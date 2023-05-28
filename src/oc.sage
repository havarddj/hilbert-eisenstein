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
        
    R.<q> = PowerSeriesRing(base_ring, default_prec = nterms)

    Ep1 = R(normalised_eis(hasse_Eis_weight(p), nterms))
    E4 = R(normalised_eis(4, nterms))
    E6 = R(normalised_eis(6, nterms))
    Delta = R(CuspForms(1,12).0.q_expansion(nterms))
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
        basis += [R(b).add_bigoh(nterms) for b in B if b != 0]
    return basis
