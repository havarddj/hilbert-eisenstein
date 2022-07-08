def Barner_Lval(F, k, reflect=False):
    # for k in ks:
    #     print(f"argument s = {k}")
    assert k % 2 == 1, "k should be odd"
    D = F.discriminant()
    D0 = fundamental_discriminant(D)
    f = sqrt(D/D0)
    
    K.<a> = QuadraticField(D0)
    sqrtD = sqrt(K(D))
    if F[0] <= 0:
        F = BinaryQF(F[2],-F[1],F[0])
    rho = (-F[1]+sqrtD)/(2*F[0])
    print(rho)
    # find fundamental unit
    r,s = solve_pell(D)
    u = r + s*sqrtD
    assert u > 0
    c = s/(rho[1])
    d = r-c*(rho[0])
    print(u, c*rho + d, f"norm of u is {u.norm()}")
    assert u == c*rho + d, "wrong fundamental unit"
    
    # Lvals = [0 for k in ks]
    Lval = 0
    # for i in range(len(ks)):
    # k = ks[i]
    # print(f"k = {k}")
    for m in [0..2*k]:
        bin1 = binomial(2*k,m)
        gDs = gen_Dedekind_sum(2*k,m,d,c)
        Lval += (-1)^m *bin1*gDs * sum(binomial(m-1,mu)*binomial(2*k-1-m,k-1 -mu)*K(u^(2*mu+1-m)).trace() for mu in [0..k-1])
    
    Lval*= 2^(2*k-2)*pi^(2*k)/(factorial(2*k)*sqrt(K(D)^(2*k-1))*K.ideal(1,rho).norm()^(k-1))
    if reflect == True:
        C = (sqrt(K(D)*K.ideal(f).norm())/pi)^(2*k-1)
        # TODO: insert root number sign W(chi).
        Lval*=C*gamma((1+k)/2)*gamma(k/2)/gamma(1-k/2)/gamma(1/2-k/2)
    return Lval

# def reflect_Lval(Lval,s):
#     return Lval*



def gen_Dedekind_sum(n, m, d, c):
    # see eq (3.60) in [Barner 69]
    assert c > 0, "c should be positive"
    assert gcd(c, d) == 1, "c and d should be coprime"
    assert m <= 2 * n, "bad choice of m and n"
    gDs = 0
    for mu in range(0, c):
        gDs += bernoulli_polynomial(d * mu / c, m) * bernoulli_polynomial(
            mu / c, 2 * n - m)
    return gDs

def solve_pell(N, numTry = 100):
    cf = continued_fraction(sqrt(N))
    for i in range(numTry):
        denom = cf.denominator(i)
        numer = cf.numerator(i)
        if numer ^2 - N * denom ^2 == 1:
            return numer , denom
    return None , None
