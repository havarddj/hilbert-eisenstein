def has_negative_fundamental_unit(D):
    """ Test if Q(sqrt D) has negative fundamental unit
    """
    F = QuadraticField(D)
    return F(F.unit_group().gens()[-1]).norm() == -1


def is_GS_disc(D, p):
    return (is_fundamental_discriminant(D) and kronecker_symbol(D, p) == -1
            and not has_negative_fundamental_unit(D))


# useful function for looking for q-expansions modulo constant terms:
def find_difference(f, g):
    # print("Normalised Eis:", f[1] / g[1] * g)
    return f - f[1] / g[1] * g


def BQFs(D):
    """Shorthand for BinaryQF_reduced_representatives, which is a mouthfull"""
    return BinaryQF_reduced_representatives(D, primitive_only=True)
