def SH_point(D, p, pprec=50):
    """
    Compute Stark--Heegner point on the elliptic curve of conductor p, when
    $X_0(p)$ has genus 1.
    
    """
    assert CuspForms(
        Gamma0(p)).dimension() == 1, f"p = {p} is not a genus 1 prime."
    E = EllipticCurve(f"{p}a")
    # note: curve N.a in the Cremona label is always
    # Gamma0(N)-optimal, aka. "best quotient", so this is the one we
    # want.
    ETate = E.tate_curve(p)
    # Todo:
    #  - implement ordinary projection
    #  - find OP in MF space
    #  - compute spectral expansion
    #  - compute L_alg for E
    #  - detect ratio as element Qp (or Fp?) and map to Tate curve
    #  - try to lift point to Q or use algdep to solve for coordinates and plug
    #    in roots over C or Qbar to check if OK
