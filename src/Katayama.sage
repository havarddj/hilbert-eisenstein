# Here we implement Theorems 1 and 2 of [Katayama 1976], computing special values of Hecke L-functions attached to ray class characters of real quadratic fields.
# NB: this relies on the Hecke character branch of sage, which stopped working for me:<


def Katayama_Lvalue(chi, k):

    K = chi.parent().number_field()
    assert chi.is_primitive(), "chi must be primitive"
    assert K.degree() == 2, f"K = {K} must be real quadratic"
    k_neg = False
    l = k
    if ZZ(k) <= 0:
        k_neg = True
        l = 1 - k

    if chi.signature() == range(K.signature()[0]):
        if (l % 2 == 0):
            if k_neg:
                return fnl_eqn_factors(chi, k) * Katayama_Lvalue_pos_even(
                    chi, l)
            else:
                return Katayama_Lvalue_pos_even(chi, l)
        else:
            if k_neg:
                return 0
            else:
                raise ValueError(
                    "Not implemented, try Dokchitser for approximation")
    else:
        if (l % 2 == 1):
            if k_neg:
                return fnl_eqn_factors(chi, k) * Katayama_Lvalue_pos_odd(
                    chi, l)
            else:
                return Katayama_Lvalue_pos_odd(chi, l)
        else:
            if k_neg:
                return 0
            else:
                raise ValueError(
                    "Not implemented, try Dokchitser for approximation")


def Katayama_Lvalue_pos_even(chi, l):
    return None


def Katayama_Lvalue_pos_odd(chi, l):
    return None


def Gauss_sum():
    return None


def sigma_B():
    return None


def epsilon_f():
    return None


def Katayama_Dedekind_sum():
    return None
