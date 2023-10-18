# Compute Stark-Heegner points using ordinary projection


def SH_point(D, p, pprec=30):
    """
    Compute Stark--Heegner point on the elliptic curve of conductor p, when
    $X_0(p)$ has genus 1.
    
    """
    assert CuspForms(
        Gamma0(p)).dimension() == 1, f"p = {p} is not a genus 1 prime."

    assert kronecker_symbol(D, p) == -1, f"p = {p} should be inert in Q(sqrtD)"
    E = EllipticCurve(f"{p}a")
    # note: curve N.a in the Cremona label is always
    # Gamma0(N)-optimal, aka. "best quotient", so this is the one we
    # want.
    MF_space = ModularForms(p)
    assert MF_space.dimension() == 2, "genus 1 primes only, for now!"

    Eis = MF_space.eisenstein_series()[0].qexp(pprec)  # arbitrary prec
    f = E.modular_form().qexp(pprec)
    ETate = E.tate_curve(p)
    Lalg = E.sha().an() * prod(E.tamagawa_numbers()) / E.torsion_order() ^ 2

    for Q in BinaryQF_reduced_representatives(D):
        G = ordinary_projection(diagonal_restriction_derivative(Q, p, pprec))
        Rp = G.parent()
        f_vector = find_in_space(G, [Rp(f), Rp(Eis)])
        print("found f in space of modular forms!!")
        print("f_vector:", f_vector)

        u = exp(2 * f_vector[0] / Lalg)
        u_vec = p_adic_eltseq(u)
        assert u_vec[
            1] == 0, f"Sage only supports Tate curve over Qp, and lambda_f is not defined over Qp"
        P = ETate.parametrisation_onto_original_curve(u_vec[0])
        print([algdep(P[i], 2) for i in range(2)])

    return P

    # Todo:
    #  - [x] implement ordinary projection
    #  - [x] find OP in MF space
    #  - [x] compute spectral expansion
    #  - [x] compute L_alg for E
    #  - detect ratio as element Qp (or Fp?) and map to Tate curve
    #  - try to lift point to Q or use algdep to solve for coordinates and plug
    #    in roots over C or Qbar to check if OK
