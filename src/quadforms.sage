from sage.all import *

attach("src/quadforms.spyx")

## Meyer stuff
# computes special values based on Meyer's formula.
# Below you will find an alternative method based on Zagier reduction

def Meyer(F):
    # based on Rademacher-Grosswald and DIT18, computes Psi(F)/6, which by Meyer's theorem equals a difference of partial L-valueqs attached to associated narrow ideal classes over Q(\sqrt{D}) where D is the discriminant of F.
    M = F.automorph()

    a, b = M[0]
    c, d = M[1]
    if c == 0:
        Phi = b / d
    else:
        Phi = (a + d) / c - 12 * sign(c) * my_dedekind_sum(
            a, c
        )  # apparently dedekind sums exist in sage, but they don't agree with ours! (probably means ours are wrong)
    Psi = Phi - 3 * sign(c * (a + d))
    return Psi / 12  # as in magma code


def my_dedekind_sum(a, c):
    assert c != 0 and GCD(a, c) == 1
    if c < 0:
        return my_dedekind_sum(a, -c)

    if a < 0 or a > c:
        return my_dedekind_sum(a % c, c)

    if c == 1:
        return 0  # then k/c is an integer

    # now we're in the situation of Apostol Mod forms ex. 3.10
    # run euclidean algo:
    x = a
    y = c
    R = []
    while y != 0:
        t = y  # r_k-1
        y = x % y  # r_k+1
        x = t  # r_k
        R.append(t)

    return sum(
        (-1) ^ (j + 1) * (R[j] ^ 2 + R[j - 1] ^ 2 + 1) / (R[j] * R[j - 1])
        for j in [1..len(R) - 1]) / 12 - ((-1) ^ len(R) + 1) / 8


def B1(a):
    # Bernoulli sawtooth function
    if a == a.floor():
        return 0
    return a - a.floor() - 1 / 2


def automorph_via_unit(F):
    """Compute automorph by computing fundamental unit in order"""
    if not Q.is_indefinite():
        return False
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

def BQF_automorph(Q):
    """ Compute automorph by Zagier reduction
    """
    D = Q.discriminant()
    reduction_matrix = identity_matrix(ZZ,2)
    # matrix sending Q to a reduced form Q'

    Q1 = Q
    while not Q.is_Zagier_reduced():
         a, b, c = list(Q)
         n = ceil((b + sqrt(D.n())) / (2 * a))
         Q = Q.parent()([a * n ^ 2 - b * n + c, 2 * a * n - b, a])
         reduction_matrix = reduction_matrix*Matrix(ZZ,[[n,1],[-1,0]])

    assert Q1.matrix_action_right(reduction_matrix) == Q
    Q2 = Q
    

    stab_matrix = identity_matrix(ZZ,2)
    # matrix stabilising the reduced form Q'
    while True:
        a, b, c = list(Q2)
        n = ceil((b + sqrt(D.n())) / (2 * a))
        Q2 = Q2.parent()([a * n ^ 2 - b * n + c, 2 * a * n - b, a])
        stab_matrix = stab_matrix*Matrix(ZZ,[[n,1],[-1,0]])
        if Q2 == Q:
            break
    gamma = reduction_matrix* stab_matrix*reduction_matrix^-1
    assert Q1.matrix_action_right(gamma) == Q1
    return gamma





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


def BQF_conductor(Q):
    D = Q.discriminant()
    return ZZ(sqrt(D / fundamental_discriminant(D)))



# this is provided by I.quadratic_form()
def ideal_to_BQF(I):
    F = I.parent().ring()
    sqrtD = F.gen()

    a,b = I.basis()

    aconj = list(a)[0] - list(a)[1]*sqrtD
    bconj = list(b)[0] - list(b)[1]*sqrtD
    if  ( aconj*b-bconj*a )/sqrtD > 0:
        a,b = b,a               # this is legal in python!
        aconj,bconj = bconj, aconj
    N = I.norm()
    return BinaryQF(a*aconj/N, ( a*bconj +b*aconj)/N, b*bconj/N)


def BQF_to_ideal(Q):
    """
    return ideal corresponding to Q of discriminant D in Q(\sqrt D)

    """
    D = Q.discriminant()
    a,b,c = list(Q)
    F.<sqrtD> = QuadraticField(D)
    
    if a > 0:
        return F.ideal([1, (b + sqrtD)/(2*a)])
    else:
        return F.ideal([sqrtD, sqrtD*(b + sqrtD)/(2*a)])
    




def fundamental_unit_in_order(D, f):
    """ Compute fundamental unit in the order O of Q(\sqrtD) of conductor f;
    naively, by looking for power of fundamental unit which lands in O.
    """
    F = QuadraticField(D)
    OF = F.maximal_order()
    eps = F.unit_group().fundamental_units()[0]

    O = F.order(f * OF.gen(0), OF.gen(1))
    assert O.index_in(OF) == f
    while eps not in O:
        eps *= eps
    return eps


def is_discriminant(D):
    return sqrt(D) not in ZZ and D == QuadraticField(D).discriminant()

def stable_root(Q):
    F = QuadraticField(Q.discriminant())
    return (-Q[1] + F.0)/(2*Q[0])



## Zagier reduction
# This was suggested by a very helpful anonymous reviewer. Thank you very much!


def is_zreduced(Q):
    """ Test if indefinite BQF Q  is Zagier-reduced, as described in [Zag, §13]
    """
    assert Q.is_indefinite(),  f"Q = {Q} is not an indefinite binary quadratic form"
    a, b, c = list(Q)
    return a > 0 and c > 0 and b > a + c



def zred_apply_T(Q, n_out=False):
    """ Apply operator T = Sn from eq. 13.1 in [Zag] to Q
    Optionally, output the number n = ceil(r_Q), where r_Q is the 'first root' of Q
    """
    
    a, b, c = list(Q)

    D = Q.discriminant()

    n = ceil((b + sqrt(D.n())) / (2 * a+.0001))  # nasty fix for bug for D=225
    if n_out:
        return Q.parent()([a * n ^ 2 - b * n + c, 2 * a * n - b, a]), n
    else:
        return Q.parent()([a * n ^ 2 - b * n + c, 2 * a * n - b, a])


def zreduce(Q):
    assert Q.is_indefinite(), f"Q = {Q} is not an indefinite binary quadratic form"
    while not Q.is_Zagier_reduced():
        Q = zred_apply_T(Q)
    return Q


def zred_detect_cycle(Q):
    """Compute reduced cycle of Q, and output the cycle as a list of
    quadratic forms Q_i as well as the list of transition numbers n_i.
    """
    
    assert Q.is_indefinite(), f"Q = {Q} is not an indefinite binary quadratic form"
    Q = Q.Zagier_reduce()

    Q0 = Q
    Q_list = []
    n_list = []
    while True:
        Q_list.append(Q)
        Q, n = zred_apply_T(Q, n_out=True)
        n_list.append(n)
        if Q == Q0:
            break

    return Q_list, n_list


def zred_L_value(Q):
    """Compute L-value of partial zeta function with narrow ideal
    class corresponding to Q, using Satz 2 in §14 of [Zag]
    """
    assert Q.is_indefinite(), f"Q = {Q} is not an indefinite binary quadratic form"
    return sum(ni - 3 for ni in Q.Zagier_reduced_cycle()[1]) / 12



# For some benchmarking, run
# > timeit("[Meyer(Q) for Q in BinaryQF_reduced_representatives(400024)]")
# > timeit("[Q.Zagier_L_value() for Q in BinaryQF_reduced_representatives(400024)]")
# for various values of 400024. They're roughly equally fast!

# define some useful methods: (is this bad practise?)

sage.quadratic_forms.binary_qf.BinaryQF.conductor = BQF_conductor
sage.quadratic_forms.binary_qf.BinaryQF.ideal = BQF_to_ideal
sage.quadratic_forms.binary_qf.BinaryQF.stable_root = stable_root
sage.quadratic_forms.binary_qf.BinaryQF.is_Zagier_reduced = is_zreduced
sage.quadratic_forms.binary_qf.BinaryQF.Zagier_reduce = zreduce
sage.quadratic_forms.binary_qf.BinaryQF.Zagier_reduced_cycle = zred_detect_cycle
sage.quadratic_forms.binary_qf.BinaryQF.Zagier_L_value = zred_L_value
sage.quadratic_forms.binary_qf.BinaryQF.automorph = BQF_automorph
