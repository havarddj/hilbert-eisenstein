attach("src/quadforms.spyx")


def Meyer(F):
    # based on Rademacher-Grosswald and DIT18, computes Psi(F)/6, which by Meyer's theorem equals a difference of partial L-values attached to associated narrow ideal classes over Q(\sqrt{D}) where D is the discriminant of F.
    M = automorph(F)
    a, b = M[0]
    c, d = M[1]
    if c == 0:
        Phi = b / d
    else:
        Phi = (a + d) / c - 12 * sign(c) * my_dedekind_sum(
            a, c
        )  # apparently dedekind sums exist in sage, but they don't agree with ours! (probably means ours are wrong)
    Psi = Phi - 3 * sign(c * (a + d))
    return Psi / 6  # uses normalisation L(0,1_A-A*) = 2zeta_-(0,A) in DIT18


def my_dedekind_sum(a, c):
    assert c != 0 and gcd(a, c) == 1
    return sum(B1(a * k / c) * B1(k / c) for k in [1..c.abs()])


def B1(a):
    # Bernoulli sawtooth function
    if a == a.floor():
        return 0
    return a - a.floor() - 1 / 2


def automorph(F):
    # first, find fundamental unit in order associated to F
    D = F.discriminant()
    f = F.conductor()
    eps = fundamental_unit_in_order(D, f)
    if eps.norm() == 1:
        t = eps.trace().abs()
    else:
        t = (eps ^ 2).trace().abs()

    u = ZZ(sqrt((t ^ 2 - 4) / D))
    assert t ^ 2 - D * u ^ 2 == 4
    a, b, c = list(F)
    M = Matrix([[(t + b * u) / 2, c * u], [-a * u, (t - b * u) / 2]])
    assert M.determinant() == 1
    return M


def genus_field_roots_of_1(D):
    m = D.squarefree_part()
    if m % 4 == 1:
        disc = m
    else:
        disc = 4 * m
    e = 2

    for fac in factor(disc):
        p = fac[0]
        if p == 2 and m % 4 == 3:
            e *= 2
        if p == 3:
            e *= 3

    return e


def BQFconductor(Q):
    D = Q.discriminant()
    return ZZ(sqrt(D / fundamental_discriminant(D)))


def fundamental_unit_in_order(D, f):
    """ Compute fundamental unit in the order O of Q(\sqrtD) of conductor f;
    naively, by looking for power of fundamental unit which lands in O.
    """
    F = QuadraticField(D)
    eps = F.unit_group().fundamental_units()[0]
    O = F.order(F.gen())
    while eps not in O:
        eps *= eps
    return eps


def is_discriminant(D):
    return sqrt(D) not in ZZ and D == QuadraticField(D).discriminant()


sage.quadratic_forms.binary_qf.BinaryQF.automorph = automorph
sage.quadratic_forms.binary_qf.BinaryQF.conductor = BQFconductor
